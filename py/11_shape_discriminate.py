# -*- coding: utf-8 -*-
"""
11_shape_discriminate.py

为 Fig 5（方向一致性热图）与 Fig 7（分布形状判别）准备数据。

口径（与 课题故事线与框架_v1.md 一致）：
  - 配对 7 例：01/09/12/13/16/18/23
  - 稳健时点：每例「最后一个 <= T4 的活检时点」（T5 为手术标本，排除）
  - 恶性细胞：CNV class == tumor
  - 高表达阈值：该 MP 在全部 7 例基线恶性细胞上的 Q75，固定不变

形状判别（位移以基线 IQR 标准化）：
  z25 / z50 / z75 = ΔQ25 / ΔQ50 / ΔQ75 ÷ IQR_base
  asym = z75 - z25
  disp = IQR_post / IQR_base

  按顺序判定（先形态、后位置）：
    1. disp > 1.6 且 |z50| < 0.5  -> dispersion up   分布展宽（细胞被推向两端）
    2. disp < 0.63 且 |z50| < 0.5 -> dispersion down 分布收缩（趋于同质）
    3. z25,z50,z75 同向且 > 0.3    -> whole shift up   整体右移（诱导）
    4. z25,z50,z75 同向且 < -0.3   -> whole shift down 整体左移
    5. z75 > 0.5 且 |z25|,|z50| 小  -> right-tail expansion 右尾扩张（亚群扩增）
    6. z25 > 0.5 且 |z75| 小        -> left truncation    左端截断（克隆内状态选择）
    7. 其余                        -> no change

说明：z 值以 IQR 为尺度，IQR 小的程序会被放大，故另附 dq50 绝对值列供核对。
"""
import csv
import os
from collections import defaultdict

BASE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                    "data", "raw", "nc2025")
TAB = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "results", "tables")

PAIRED = ["01", "09", "12", "13", "16", "18", "23"]
LAST = {"01": "T3", "09": "T3", "12": "T3", "13": "T3",
        "16": "T4", "18": "T4", "23": "T2"}
RCB = {"01": "favorable", "09": "unfavorable", "12": "unfavorable",
       "13": "unfavorable", "16": "unfavorable", "18": "unfavorable",
       "23": "unfavorable"}
ARM = {"01": "Chemo", "09": "ICI", "12": "ICI", "13": "Chemo",
       "16": "Chemo", "18": "ICI", "23": "ICI"}
MPS = ["MP2", "MP3", "MP4", "MP5", "MP6", "MP7", "MP8", "MP9", "MP11"]


def quantile(v, p):
    """线性插值分位数（与 R type=7 一致）"""
    n = len(v)
    if n == 0:
        return None
    v = sorted(v)
    h = (n - 1) * p
    lo = int(h)
    hi = min(lo + 1, n - 1)
    return v[lo] + (h - lo) * (v[hi] - v[lo])


def main():
    # ---------- 1. CNV 表：细胞的 class 与亚克隆 ----------
    kv = {}
    with open(os.path.join(BASE, "gex_cnv_SCEVAN_karyotype.csv")) as f:
        r = csv.reader(f)
        next(r)
        for row in r:
            if len(row) < 6:
                continue
            kv[f"{row[0]}_{row[4]}_{row[5]}"] = (row[1], row[3])
    print(f"CNV 表细胞 {len(kv)}")

    # ---------- 2. 评分表：收集恶性细胞评分 ----------
    # cells[(patient, tp)] = {mp: [scores]}
    cells = defaultdict(lambda: defaultdict(list))
    n_row = n_match = n_tumor = 0
    with open(os.path.join(BASE, "signature_scores.csv")) as f:
        r = csv.reader(f)
        header = next(r)
        mp_order = header[1:]
        for row in r:
            if not row or not row[0]:
                continue
            n_row += 1
            bc = row[0]
            if bc not in kv:
                continue
            n_match += 1
            cls, _ = kv[bc]
            if cls != "tumor":
                continue
            n_tumor += 1
            parts = bc.split("_")
            if len(parts) < 4 or not parts[1].lower().startswith("ccg"):
                continue
            pid, tp = parts[2].zfill(2), parts[3]
            if pid not in PAIRED:
                continue
            if tp != "T1" and tp != LAST[pid]:
                continue
            for i, mp in enumerate(mp_order, start=1):
                if mp not in MPS:
                    continue
                cells[(pid, tp)][mp].append(float(row[i]))
    print(f"评分表 {n_row} 行；匹配 CNV {n_match}；tumor {n_tumor}")

    # ---------- 3. 每个 MP 的固定高表达阈值（基线 Q75，跨 7 例合并） ----------
    thresh = {}
    for mp in MPS:
        allb = []
        for p in PAIRED:
            allb.extend(cells[(p, "T1")][mp])
        thresh[mp] = quantile(allb, 0.75)
        print(f"  {mp:<5} 基线合并 n={len(allb):>6}  阈值 Q75={thresh[mp]:.3f}")

    # ---------- 4. 分位数位移 + 形状判定 ----------
    shape_rows = []
    heat_rows = []
    for p in PAIRED:
        tp2 = LAST[p]
        for mp in MPS:
            b = cells[(p, "T1")][mp]
            a = cells[(p, tp2)][mp]
            if len(b) < 30 or len(a) < 30:
                continue
            qb = {k: quantile(b, v) for k, v in
                  (("10", .10), ("25", .25), ("50", .50), ("75", .75), ("90", .90))}
            qa = {k: quantile(a, v) for k, v in
                  (("10", .10), ("25", .25), ("50", .50), ("75", .75), ("90", .90))}
            iqr = qb["75"] - qb["25"]
            if iqr <= 0:
                continue
            z = {k: (qa[k] - qb[k]) / iqr for k in qb}
            asym = z["75"] - z["25"]
            disp = (qa["75"] - qa["25"]) / iqr
            big = max(abs(z["25"]), abs(z["50"]), abs(z["75"]))

            if disp > 1.6 and abs(z["50"]) < 0.5:
                call = "dispersion up"
            elif disp < 0.63 and abs(z["50"]) < 0.5:
                call = "dispersion down"
            elif z["25"] > 0.3 and z["50"] > 0.3 and z["75"] > 0.3:
                call = "whole shift up"
            elif z["25"] < -0.3 and z["50"] < -0.3 and z["75"] < -0.3:
                call = "whole shift down"
            elif z["75"] > 0.5 and abs(z["25"]) < 0.35 and abs(z["50"]) < 0.35:
                call = "right-tail expansion"
            elif z["25"] > 0.5 and abs(z["75"]) < 0.35:
                call = "left truncation"
            else:
                call = "no change"
            shape_rows.append(dict(
                patient=p, tp=tp2, mp=mp, n_base=len(b), n_post=len(a),
                q25_base=f"{qb['25']:.3f}", q50_base=f"{qb['50']:.3f}",
                q75_base=f"{qb['75']:.3f}",
                q25_post=f"{qa['25']:.3f}", q50_post=f"{qa['50']:.3f}",
                q75_post=f"{qa['75']:.3f}",
                dq25=f"{qa['25'] - qb['25']:.3f}",
                dq50=f"{qa['50'] - qb['50']:.3f}",
                dq75=f"{qa['75'] - qb['75']:.3f}",
                iqr_base=f"{iqr:.3f}", disp=f"{disp:.2f}",
                z25=f"{z['25']:.3f}", z50=f"{z['50']:.3f}", z75=f"{z['75']:.3f}",
                asym=f"{asym:.3f}", call=call,
                rcb=RCB[p], arm=ARM[p]))

            # 热图用：Δ高表达占比 与 Δ中位
            hb = sum(1 for v in b if v >= thresh[mp]) / len(b) * 100
            ha = sum(1 for v in a if v >= thresh[mp]) / len(a) * 100
            heat_rows.append(dict(
                patient=p, tp=tp2, mp=mp, rcb=RCB[p], arm=ARM[p],
                hi_base=f"{hb:.1f}", hi_post=f"{ha:.1f}",
                d_hi=f"{ha - hb:.1f}",
                d_med=f"{qa['50'] - qb['50']:.3f}",
                z50=f"{z['50']:.3f}",
                shape_call=call))

    # ---------- 5. 密度曲线用的细胞级长表 ----------
    # 只保留「有实质位移」的 患者×程序 组合，外加 MP7/MP11 全部（作对照）
    keep = {(r["patient"], r["mp"]) for r in shape_rows if abs(float(r["z50"])) > 0.45}
    keep |= {(p, "MP7") for p in PAIRED} | {(p, "MP11") for p in PAIRED}
    dense_rows = []
    for p in PAIRED:
        for tp in ("T1", LAST[p]):
            for mp in MPS:
                if (p, mp) not in keep:
                    continue
                for v in cells[(p, tp)][mp]:
                    dense_rows.append((p, tp, mp, f"{v:.4f}",
                                       RCB[p], ARM[p],
                                       "baseline" if tp == "T1" else "on-treatment"))
    print(f"密度长表保留 {len(keep)} 个 患者×程序 组合")

    # ---------- 6. 写出 ----------
    os.makedirs(TAB, exist_ok=True)

    def dump(name, rows, cols=None):
        if not rows:
            print(f"  !! {name} 为空")
            return
        cols = cols or list(rows[0].keys())
        with open(os.path.join(TAB, name), "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=cols)
            w.writeheader()
            for r in rows:
                w.writerow({c: r[c] for c in cols})
        print(f"  -> {name}  ({len(rows)} 行)")

    dump("nc2025_v3_shape.csv", shape_rows)
    dump("nc2025_v3_heatmap.csv", heat_rows)

    # 每例的「主导响应程序」：|z50| 最大者
    dom = []
    for p in PAIRED:
        rs = [r for r in shape_rows if r["patient"] == p]
        if not rs:
            continue
        top = max(rs, key=lambda r: abs(float(r["z50"])))
        ranked = sorted(rs, key=lambda r: -abs(float(r["z50"])))
        dom.append(dict(patient=p, tp=LAST[p], rcb=RCB[p], arm=ARM[p],
                        top_mp=top["mp"], top_z50=top["z50"], top_dq50=top["dq50"],
                        top_call=top["call"],
                        second_mp=ranked[1]["mp"] if len(ranked) > 1 else "",
                        second_z50=ranked[1]["z50"] if len(ranked) > 1 else "",
                        n_mp_abs_z50_gt05=sum(1 for r in rs if abs(float(r["z50"])) > 0.5)))
    dump("nc2025_v3_dominant.csv", dom)

    with open(os.path.join(TAB, "nc2025_v3_density.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["patient", "tp", "mp", "score", "rcb", "arm", "phase"])
        w.writerows(dense_rows)
    print(f"  -> nc2025_v3_density.csv  ({len(dense_rows)} 行)")

    # ---------- 7. 汇总打印 ----------
    print("\n=== Fig 7 形状判别（7 例 × 9 程序）===")
    print("患者 时点 MP     z25    z50    z75   disp  判定")
    for r in shape_rows:
        print(f"{r['patient']}   {r['tp']}   {r['mp']:<5} "
              f"{float(r['z25']):+6.2f} {float(r['z50']):+6.2f} {float(r['z75']):+6.2f} "
              f"{float(r['disp']):5.2f}   {r['call']}")

    print("\n=== 判定分布（按程序）===")
    for mp in MPS:
        rs = [r for r in shape_rows if r["mp"] == mp]
        cnt = defaultdict(int)
        for r in rs:
            cnt[r["call"]] += 1
        print(f"  {mp:<5} " + "  ".join(f"{k}={v}" for k, v in sorted(cnt.items())))

    print("\n=== 每例主导响应程序 ===")
    print("患者 时点  分组         主导程序   z50     次位        z50     实质位移程序数")
    for d in dom:
        print(f"P{d['patient']}  {d['tp']}  {d['rcb'][:3]}/{d['arm']:<6} "
              f"{d['top_mp']:<8} {float(d['top_z50']):+6.2f}  "
              f"{d['second_mp']:<8} {float(d['second_z50']):+6.2f}   {d['n_mp_abs_z50_gt05']}")


if __name__ == "__main__":
    main()
