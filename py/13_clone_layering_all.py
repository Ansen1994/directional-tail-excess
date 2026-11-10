# -*- coding: utf-8 -*-
"""
13_clone_layering_all.py

把「克隆内协同转变」从 P16 / P18 两例推广检验到全部配对病例。

背景
----
第三轮（11_shape_discriminate.py）确认 7 例全部是 whole shift up，即整体平移，
而非右尾扩张 / 左端截断。整体平移在细胞层面有两种互斥解释：

  H_pre  预存耐药克隆：治疗前已存在一个高表达克隆，治疗后它占比升高。
         可检验推论：(a) 基线各克隆在主导程序上应有显著差异；
                     (b) 位移应随克隆占比变化而出现，克隆内分布不动。
  H_ind  诱导性状态转变：治疗把细胞从低状态推到高状态，各克隆一起动。
         可检验推论：(a) 基线各克隆无显著差异（无预存）；
                     (b) 同一克隆治疗前后自身位移显著且方向一致；
                     (c) 位移幅度应超过基线克隆间极差（预存解释无法覆盖）。

本脚本对每例、每个 MP 计算：
  base_spread   基线各克隆中位数的极差 ÷ 基线全体 IQR（标准化克隆间异质性）
  separation    随访「最低」留存克隆中位数 − 基线「最高」克隆中位数
                > 0 表示治疗后最差的克隆都已高于治疗前最好的克隆，即位移
                完全越出基线克隆异质性的范围 —— 预存假说无法解释（主判据）
  effect_ratio  |Δ中位| ÷ (base_spread × IQR_base / 2) > 2 亦表示位移超越克隆间差异
  icc_post      随访时点克隆间方差 ÷ 总方差（克隆能解释多少细胞间差异）
  perm_p        打乱数值-标签配对重算极差的置换检验 p（仅作参考，见下）

关于 perm_p 的方法学限制
------------------------
置换检验把每个细胞当作独立观测，但同一克隆内的细胞高度相关（共享核型与
克隆内状态），属于伪重复，会系统性低估 p 值。本脚本仍输出该值以便与其他
分析口径衔接，但**不以 p 值作判定依据**，判定以 separation / effect_ratio /
符号一致性 / ICC 四个效应量指标为准。

克隆身份可靠性
--------------
SCEVAN 亚克隆编号跨时点全局一致（P16 在 T1 编 1-4，T4 编 1-2，编号同源；
P18 在 T2/T3 编 1-5，T1 也编 1-5），因此同编号跨时点视为同一克隆。
局限：不排除治疗后同一编号内部发生亚克隆演化，故本分析只作「克隆层面
是否整体同步」的判据，不作克隆系统发生推断。

不可分层病例
------------
  P12  T1 有 1281 个肿瘤细胞但 SCEVAN 未给出任何亚克隆 → 只能做整体层面
  P01  T3 每克隆仅 8-86 个细胞 → 低于 min_cells，如实标注 low_n 不强行判定
"""
import csv
import os
import random
from collections import defaultdict

BASE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                    "data", "raw", "nc2025")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(ROOT, "results", "tables")

PAIRED = ["01", "09", "12", "13", "16", "18", "23"]
LAST = {"01": "T3", "09": "T3", "12": "T3", "13": "T3",
        "16": "T4", "18": "T4", "23": "T2"}
RCB = {"01": "favorable", "09": "unfavorable", "12": "unfavorable",
       "13": "unfavorable", "16": "unfavorable", "18": "unfavorable",
       "23": "unfavorable"}
ARM = {"01": "Chemo", "09": "ICI", "12": "ICI", "13": "Chemo",
       "16": "Chemo", "18": "ICI", "23": "ICI"}
MPS = ["MP2", "MP3", "MP4", "MP5", "MP6", "MP7", "MP8", "MP9", "MP11"]

MIN_CELLS = 30      # 单个 克隆×时点 的最低细胞数
N_PERM = 2000       # 置换次数
SEED = 20261108


def quantile(v, p):
    n = len(v)
    if n == 0:
        return None
    v = sorted(v)
    h = (n - 1) * p
    lo = int(h)
    hi = min(lo + 1, n - 1)
    return v[lo] + (h - lo) * (v[hi] - v[lo])


def median(v):
    return quantile(v, 0.5)


def spread_of(groups):
    """组中位数极差；组数 <2 返回 0"""
    meds = [median(g) for g in groups if g]
    if len(meds) < 2:
        return 0.0
    return max(meds) - min(meds)


def perm_p(labels, values, obs_stat, rng, n_perm=N_PERM):
    """打乱数值与克隆标签的配对，重算组中位数极差，返回单侧 p（越大越异质）"""
    n = len(values)
    if obs_stat <= 0 or n < 20 or len(set(labels)) < 2:
        return 1.0
    v = list(values)
    hit = 0
    for _ in range(n_perm):
        rng.shuffle(v)                      # 只打乱数值，标签位置固定
        buckets = defaultdict(list)
        for lab_i, val in zip(labels, v):
            buckets[lab_i].append(val)
        if spread_of(list(buckets.values())) >= obs_stat:
            hit += 1
    return (hit + 1) / (n_perm + 1)


def boot_ci(a, b, rng, n_boot=300):
    """中位数差 (b - a) 的 bootstrap 95% CI"""
    if len(a) < 10 or len(b) < 10:
        return float("nan"), float("nan")
    ds = []
    na, nb = len(a), len(b)
    for _ in range(n_boot):
        sa = [a[rng.randrange(na)] for _ in range(na)]
        sb = [b[rng.randrange(nb)] for _ in range(nb)]
        ds.append(median(sb) - median(sa))
    ds.sort()
    return ds[max(0, int(0.025 * n_boot))], ds[min(n_boot - 1, int(0.975 * n_boot))]


def main():
    rng = random.Random(SEED)

    # ---------- 1. CNV：细胞 -> (patient, tp, clone) ----------
    # 显式要求 Epithelial：SCEVAN 的 tumor class 里混有少量被误判的基质/免疫
    # 细胞（P16 T4 的 1181 个 tumor 中有 744 个实为 Stromal/Immune），
    # cNMF 评分表只覆盖上皮细胞，两层过滤后留下的才是真肿瘤上皮。
    comp = {}
    with open(os.path.join(BASE, "GEX_OBS_Lineage.csv")) as f:
        r = csv.reader(f)
        next(r)
        for row in r:
            if len(row) >= 5:
                comp[row[0]] = row[3]

    cell2clone = {}
    cell2meta = {}
    nosub = defaultdict(int)
    n_cnv = 0
    n_nonepi = 0
    with open(os.path.join(BASE, "gex_cnv_SCEVAN_karyotype.csv")) as f:
        r = csv.reader(f)
        next(r)
        for row in r:
            if len(row) < 6:
                continue
            n_cnv += 1
            if row[1] != "tumor":
                continue
            parts = row[4].split("_")
            if len(parts) < 4:
                continue
            pid, tp = parts[1].zfill(2), parts[2]
            if pid not in PAIRED or tp not in ("T1", LAST[pid]):
                continue
            key = f"{row[0]}_{row[4]}_{row[5]}"
            if comp.get(key, "") != "Epithelial":
                n_nonepi += 1
                continue
            nosub[(pid, tp)] += 1        # 上皮层内的肿瘤细胞总数（含未分克隆者）
            if not row[3]:
                continue        # tumor 但无亚克隆（P12 T1 的 1281 个）
            cell2clone[key] = row[3]
            cell2meta[key] = (pid, tp)
    print(f"CNV 表 {n_cnv} 细胞；非上皮误判 {n_nonepi}；"
          f"上皮层肿瘤细胞 {len(cell2clone)}")

    # ---------- 2. 评分表 ----------
    # vals[(pid, tp)][clone][mp] = [scores]
    vals = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    n_hit = 0
    with open(os.path.join(BASE, "signature_scores.csv")) as f:
        r = csv.reader(f)
        header = next(r)
        mp_order = header[1:]
        cmp_idx = {mp: i + 1 for i, mp in enumerate(mp_order) if mp in MPS}
        for row in r:
            if not row or not row[0]:
                continue
            key = row[0]
            if key not in cell2clone:
                continue
            n_hit += 1
            pid, tp = cell2meta[key]
            cl = cell2clone[key]
            for mp, i in cmp_idx.items():
                vals[(pid, tp)][cl][mp].append(float(row[i]))
    print(f"评分表命中 {n_hit} 细胞（应为 {len(cell2clone)}）")

    # ---------- 3. 分层明细 + 汇总量 ----------
    layer_rows = []
    sum_rows = []
    story = []
    dom_cells = []

    for p in PAIRED:
        tp0, tp1 = "T1", LAST[p]
        cl0 = sorted(vals[(p, tp0)].keys(), key=lambda x: float(x))
        cl1 = sorted(vals[(p, tp1)].keys(), key=lambda x: float(x))
        n_base_tot = sum(len(vals[(p, tp0)][c][MPS[0]]) for c in cl0)
        n_post_tot = sum(len(vals[(p, tp1)][c][MPS[0]]) for c in cl1)
        if len(cl0) < 2:
            story.append(dict(patient=p, verdict="no_subclone",
                              note=f"基线 {nosub[(p, tp0)]} 个肿瘤细胞中 "
                                   f"{n_base_tot} 个被分配亚克隆（{len(cl0)} 个克隆），"
                                   f"SCEVAN 未能解析出可比较的亚克隆结构，仅做整体层面"))
            print(f"\n[P{p}] 无亚克隆可比（基线克隆 {len(cl0)}，已分配细胞 {n_base_tot}"
                  f"/{nosub[(p, tp0)]}）→ 跳过分层")
            continue

        # ---------- 每个 MP 的克隆分层 ----------
        mps_use = []
        for mp in MPS:
            g0 = [vals[(p, tp0)][c][mp] for c in cl0 if len(vals[(p, tp0)][c][mp]) >= MIN_CELLS]
            if len(g0) < 2:
                continue
            all0 = [v for g in g0 for v in g]
            iqr0 = quantile(all0, .75) - quantile(all0, .25)
            if iqr0 <= 0:
                continue
            med0 = [median(g) for g in g0]
            obs = max(med0) - min(med0)
            bs = obs / iqr0

            # 置换检验：只对「有实质位移」的程序做（无位移的程序不需要判预存）
            all1_pre = [v for c in cl1 for v in vals[(p, tp1)][c][mp]
                        if len(vals[(p, tp1)][c][mp]) >= MIN_CELLS]
            d_pre = (median(all1_pre) - median(all0)) if len(all1_pre) >= MIN_CELLS else 0.0
            if len(all0) <= 8000 and abs(d_pre) >= 0.15:
                labels = [c for c in cl0 if len(vals[(p, tp0)][c][mp]) >= MIN_CELLS
                          for _ in range(len(vals[(p, tp0)][c][mp]))]
                pp = perm_p(labels, all0, obs, rng, n_perm=1000)
            else:
                pp = float("nan")

            # 留存克隆的自身位移
            ret, up, down = [], 0, 0
            detail = []
            for i, c in enumerate(cl0):
                if len(vals[(p, tp0)][c][mp]) < MIN_CELLS:
                    continue
                if len(vals[(p, tp1)][c][mp]) < MIN_CELLS:
                    layer_rows.append(dict(
                        patient=p, tp=tp1, mp=mp, clone=c,
                        n_base=len(vals[(p, tp0)][c][mp]),
                        med_base=f"{median(vals[(p, tp0)][c][mp]):.4f}",
                        n_post=len(vals[(p, tp1)][c][mp]),
                        med_post="", d_med="", ci_lo="", ci_hi="", status="dropped"))
                    continue
                mb = median(vals[(p, tp0)][c][mp])
                ma = median(vals[(p, tp1)][c][mp])
                d = ma - mb
                ret.append(d)
                detail.append((c, vals[(p, tp0)][c][mp], vals[(p, tp1)][c][mp], d))
                if d > 0:
                    up += 1
                elif d < 0:
                    down += 1
                layer_rows.append(dict(
                    patient=p, tp=tp1, mp=mp, clone=c,
                    n_base=len(vals[(p, tp0)][c][mp]), med_base=f"{mb:.4f}",
                    n_post=len(vals[(p, tp1)][c][mp]), med_post=f"{ma:.4f}",
                    d_med=f"{d:.4f}", ci_lo="", ci_hi="", status="retained"))

            # 新发克隆
            for c in cl1:
                if c not in cl0 and len(vals[(p, tp1)][c][mp]) >= MIN_CELLS:
                    layer_rows.append(dict(
                        patient=p, tp=tp1, mp=mp, clone=c,
                        n_base=0, med_base="",
                        n_post=len(vals[(p, tp1)][c][mp]),
                        med_post=f"{median(vals[(p, tp1)][c][mp]):.4f}",
                        d_med="", ci_lo="", ci_hi="", status="new"))

            # 整体位移（全体细胞中位数差）
            all1 = all1_pre
            if len(all1) < MIN_CELLS:
                continue
            d_overall = median(all1) - median(all0)

            er = abs(d_overall) / (bs * iqr0 / 2) if bs > 0 and iqr0 > 0 else float("inf")

            # post 时点克隆内/克隆间方差分解（ICC）
            post_g = [vals[(p, tp1)][c][mp] for c in cl1
                      if len(vals[(p, tp1)][c][mp]) >= MIN_CELLS]
            icc = float("nan")
            if len(post_g) >= 2:
                tot = [v for g in post_g for v in g]
                gm = sum(tot) / len(tot)
                ss_tot = sum((v - gm) ** 2 for v in tot)
                ss_b = sum(len(g) * (sum(g) / len(g) - gm) ** 2 for g in post_g)
                icc = ss_b / ss_tot if ss_tot > 0 else float("nan")

            # 分离度：随访最低留存克隆 vs 基线最高克隆（含已消失克隆，取保守口径）
            base_meds_all = [(c, median(vals[(p, tp0)][c][mp])) for c in cl0
                             if len(vals[(p, tp0)][c][mp]) >= MIN_CELLS]
            post_meds_ret = []
            for c, bl, pl, d in detail:
                post_meds_ret.append((c, median(pl)))
            sep = float("nan")
            if base_meds_all and post_meds_ret:
                sep = (min(m for _, m in post_meds_ret)
                       - max(m for _, m in base_meds_all))

            # 消失克隆是否系统性偏低（克隆选择的另一种形态）
            ret_base = [median(bl) for _, bl, _, _ in detail]
            drop_base = [median(vals[(p, tp0)][c][mp]) for c in cl0
                         if len(vals[(p, tp0)][c][mp]) >= MIN_CELLS
                         and len(vals[(p, tp1)][c][mp]) < MIN_CELLS]
            drop_bias = float("nan")
            if ret_base and drop_base:
                drop_bias = (sum(drop_base) / len(drop_base)
                             - sum(ret_base) / len(ret_base))

            mps_use.append(dict(mp=mp, bs=bs, pp=pp, n_ret=len(ret),
                                up=up, down=down, d_overall=d_overall,
                                z50=d_overall / iqr0, sep=sep,
                                n_drop=len(drop_base), drop_bias=drop_bias,
                                er=er, icc=icc, iqr0=iqr0, n_base=len(all0),
                                n_post=len(all1), g0=len(g0), g1=len(post_g),
                                detail=detail))

        if not mps_use:
            story.append(dict(patient=p, verdict="low_n",
                              note="各克隆细胞数均低于门槛，无法分层"))
            print(f"\n[P{p}] 全部分层单元低于 {MIN_CELLS} 细胞 → 不下结论")
            continue

        # ---------- 主导程序（按 |z50| 最大，与 v3 dominant 同口径） ----------
        top = max(mps_use, key=lambda x: abs(x["z50"]))

        # 对主导程序的留存克隆补 bootstrap 95% CI（克隆内自身位移是否稳健）
        for c, bl, pl, d in top["detail"]:
            lo, hi = boot_ci(bl, pl, rng)
            for row in layer_rows:
                if (row["patient"] == p and row["mp"] == top["mp"]
                        and row["clone"] == c and row["status"] == "retained"):
                    row["ci_lo"] = f"{lo:.4f}"
                    row["ci_hi"] = f"{hi:.4f}"
                    break
        top["n_sig"] = sum(1 for row in layer_rows
                           if row["patient"] == p and row["mp"] == top["mp"]
                           and row["status"] == "retained"
                           and row["ci_lo"]
                           and (float(row["ci_lo"]) > 0 or float(row["ci_hi"]) < 0))

        # 主导程序的细胞级长表（供作图，只含有足够细胞的克隆）
        for c in cl0:
            n0 = len(vals[(p, tp0)][c][top["mp"]])
            n1 = len(vals[(p, tp1)][c][top["mp"]])
            if n0 < MIN_CELLS and n1 < MIN_CELLS:
                continue
            for v in vals[(p, tp0)][c][top["mp"]]:
                dom_cells.append((p, c, "baseline", f"{v:.4f}"))
            for v in vals[(p, tp1)][c][top["mp"]]:
                dom_cells.append((p, c, "on-treatment", f"{v:.4f}"))

        for m in mps_use:
            sum_rows.append(dict(
                patient=p, tp=tp1, rcb=RCB[p], arm=ARM[p], mp=m["mp"],
                n_base=m["n_base"], n_post=m["n_post"],
                n_clone_base=m["g0"], n_clone_post=m["g1"],
                base_spread=f"{m['bs']:.2f}", perm_p=f"{m['pp']:.4f}",
                n_retained=m["n_ret"], n_up=m["up"], n_down=m["down"],
                n_dropped=m["n_drop"],
                drop_bias=f"{m['drop_bias']:.4f}" if m["drop_bias"] == m["drop_bias"] else "",
                d_overall=f"{m['d_overall']:.4f}", z50=f"{m['z50']:.3f}",
                separation=f"{m['sep']:.4f}" if m["sep"] == m["sep"] else "",
                sep_z=f"{m['sep'] / m['iqr0']:.3f}" if m["sep"] == m["sep"] else "",
                effect_ratio=f"{m['er']:.2f}", icc_post=f"{m['icc']:.3f}",
                is_dominant="yes" if m is top else ""))

        # ---------- 该例判定（效应量口径，不用 p 值） ----------
        verdicts = []
        sep_z = top["sep"] / top["iqr0"] if top["sep"] == top["sep"] else float("nan")
        if top["bs"] >= 0.5:
            verdicts.append("baseline_heterogeneous")
        else:
            verdicts.append("baseline_homogeneous")
        if sep_z == sep_z:
            if sep_z > 0:
                verdicts.append("separated_beyond_baseline")
            elif sep_z > -0.5:
                verdicts.append("not_separated")
            else:
                verdicts.append("clearly_overlapping_baseline")
        if top["n_ret"] >= 2:
            if top["up"] == top["n_ret"]:
                verdicts.append("all_retained_same_direction")
            elif max(top["up"], top["down"]) / top["n_ret"] >= 0.75:
                verdicts.append("mostly_same_direction")
            else:
                verdicts.append("mixed_direction")
        if top["er"] > 2:
            verdicts.append("shift_exceeds_clonal_spread")

        story.append(dict(patient=p, verdict="|".join(verdicts),
                          note=f"主导 {top['mp']}（z50={top['z50']:+.2f}），"
                               f"基线 {top['g0']} 克隆 / 随访 {top['g1']} 克隆，"
                               f"留存 {top['n_ret']} 个（升 {top['up']} 降 {top['down']}），"
                               f"分离度 {top['sep']:+.3f}（{sep_z:+.2f} IQR），"
                               f"克隆异质性 {top['bs']:.2f} IQR，ICC={top['icc']:.3f}"))
        print(f"\n[P{p}] 主导 {top['mp']}  z50={top['z50']:+.2f}  "
              f"留存 {top['n_ret']} 升{top['up']}降{top['down']}  "
              f"分离度={top['sep']:+.3f} ({sep_z:+.2f} IQR)  "
              f"克隆异质性={top['bs']:.2f} IQR  ICC={top['icc']:.3f}  "
              f"置换p={top['pp']:.3f}(伪重复,仅参考)")
        print(f"       判定: {' | '.join(verdicts)}")

    # ---------- 4. 写出 ----------
    os.makedirs(TAB, exist_ok=True)

    def dump(name, rows):
        if not rows:
            print(f"  !! {name} 为空")
            return
        cols = list(rows[0].keys())
        with open(os.path.join(TAB, name), "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=cols)
            w.writeheader()
            for r in rows:
                w.writerow({c: r[c] for c in cols})
        print(f"  -> {name}  ({len(rows)} 行)")

    dump("nc2025_v4_clone_layer.csv", layer_rows)
    dump("nc2025_v4_clone_summary.csv", sum_rows)
    dump("nc2025_v4_clone_verdict.csv", story)

    with open(os.path.join(TAB, "nc2025_v4_clone_cells.csv"), "w",
              newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["patient", "clone", "phase", "score"])
        w.writerows(dom_cells)
    print(f"  -> nc2025_v4_clone_cells.csv  ({len(dom_cells)} 行)")

    # ---------- 5. 汇总 ----------
    print("\n=== 分层判定总览 ===")
    for s in story:
        print(f"P{s['patient']}  {s['verdict']:<60} {s['note']}")

    ok = [s for s in story if "separated_beyond_baseline" in s["verdict"]]
    nos = [s for s in story if "not_separated" in s["verdict"]
           or "clearly_overlapping_baseline" in s["verdict"]]
    same = [s for s in story if "all_retained_same_direction" in s["verdict"]]
    print(f"\n位移完全越出基线克隆异质性: {len(ok)}/{len(story)} 例 "
          f"({', '.join('P' + s['patient'] for s in ok)})")
    print(f"未分离（无法排除克隆选择）: {len(nos)}/{len(story)} 例 "
          f"({', '.join('P' + s['patient'] for s in nos)})")
    print(f"留存克隆全部同向: {len(same)}/{len(story)} 例 "
          f"({', '.join('P' + s['patient'] for s in same)})")


if __name__ == "__main__":
    main()
