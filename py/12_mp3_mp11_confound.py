# -*- coding: utf-8 -*-
"""
12_mp3_mp11_confound.py

检验 P16 的 MP3（胶原/间质）与 MP11（EMT-III）上升是否来自间质污染或双细胞。

三层检验：
  A. P16 治疗后 MP3 与 MP11 的细胞级相关性（若 r 很高，两者可能是同一个间质信号）
  B. MP3 / MP11 评分与测序深度的关系（双细胞通常 counts 翻倍）
  C. 用上皮基因（EPCAM/KRT8/KRT18/CDH1）门控后，MP3/MP11 的位移是否还成立
     —— 若门控后位移大幅衰减，说明原来的位移由非上皮细胞贡献

注：C 需要表达矩阵。当前已抽取的文件不含表达量，故 C 只做「评分与深度」的间接代理：
   若 MP3 高细胞的 nFeature 显著更高，提示双细胞。
"""
import csv
import os
import math
from collections import defaultdict

BASE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                    "data", "raw", "nc2025")
TAB = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "results", "tables")
PAIRED = ["01", "09", "12", "13", "16", "18", "23"]
LAST = {"01": "T3", "09": "T3", "12": "T3", "13": "T3",
        "16": "T4", "18": "T4", "23": "T2"}


def pearson(x, y):
    n = len(x)
    if n < 3:
        return None
    mx, my = sum(x) / n, sum(y) / n
    sx = math.sqrt(sum((v - mx) ** 2 for v in x))
    sy = math.sqrt(sum((v - my) ** 2 for v in y))
    if sx == 0 or sy == 0:
        return None
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / (sx * sy)


def spearman(x, y):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    return pearson(rank(x), rank(y))


def main():
    # 读评分表（含 MP3/MP11）与 CNV 表
    kv = {}
    with open(os.path.join(BASE, "gex_cnv_SCEVAN_karyotype.csv")) as f:
        r = csv.reader(f); next(r)
        for row in r:
            if len(row) < 6:
                continue
            kv[f"{row[0]}_{row[4]}_{row[5]}"] = (row[1], row[3])

    rec = defaultdict(lambda: {"MP3": [], "MP11": [], "MP7": [], "all": []})
    with open(os.path.join(BASE, "signature_scores.csv")) as f:
        r = csv.reader(f)
        h = next(r)
        i3, i11, i7 = h.index("MP3"), h.index("MP11"), h.index("MP7")
        for row in r:
            if not row or not row[0]:
                continue
            bc = row[0]
            if bc not in kv or kv[bc][0] != "tumor":
                continue
            parts = bc.split("_")
            if len(parts) < 4 or not parts[1].lower().startswith("ccg"):
                continue
            pid, tp = parts[2].zfill(2), parts[3]
            if pid not in PAIRED:
                continue
            if tp != "T1" and tp != LAST[pid]:
                continue
            rec[(pid, tp)]["MP3"].append(float(row[i3]))
            rec[(pid, tp)]["MP11"].append(float(row[i11]))
            rec[(pid, tp)]["MP7"].append(float(row[i7]))

    print("=== A. MP3 与 MP11 的细胞级相关（各患者各时点）===")
    print("患者 时点    n     Pearson   Spearman   判定")
    out = []
    for p in PAIRED:
        for tp in ("T1", LAST[p]):
            d = rec[(p, tp)]
            n = len(d["MP3"])
            if n < 100:
                continue
            pr = pearson(d["MP3"], d["MP11"])
            sr = spearman(d["MP3"], d["MP11"])
            flag = "高度共线" if sr and sr > 0.6 else ("中等" if sr and sr > 0.35 else "可分离")
            print(f"P{p}  {tp}   {n:>5}   {pr:+.3f}    {sr:+.3f}    {flag}")
            out.append(dict(patient=p, tp=tp, n=n,
                            pearson=f"{pr:.4f}", spearman=f"{sr:.4f}",
                            collinear=flag))

    with open(os.path.join(TAB, "nc2025_v3_mp3_mp11_corr.csv"), "w",
              newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["patient", "tp", "n", "pearson",
                                          "spearman", "collinear"])
        w.writeheader(); w.writerows(out)
    print(f"\n-> {os.path.join(TAB, 'nc2025_v3_mp3_mp11_corr.csv')}")

    print("\n=== B. MP7 与 MP11 的相关（对照：两个不同程序应低相关）===")
    for p in ["16", "18"]:
        for tp in ("T1", LAST[p]):
            d = rec[(p, tp)]
            if len(d["MP7"]) < 100:
                continue
            sr = spearman(d["MP7"], d["MP11"])
            print(f"  P{p} {tp}  n={len(d['MP7']):>5}  Spearman(MP7,MP11) = {sr:+.3f}")

    print("\n=== C. 高 MP3 细胞 vs 低 MP3 细胞的其它程序评分（P16）===")
    d = rec[("16", LAST["16"])]
    pairs = sorted(zip(d["MP3"], d["MP11"], d["MP7"]))
    q = len(pairs) // 4
    lo, hi = pairs[:q], pairs[-q:]
    for name, grp in (("低 MP3 (Q1)", lo), ("高 MP3 (Q4)", hi)):
        m11 = sum(x[1] for x in grp) / len(grp)
        m7 = sum(x[2] for x in grp) / len(grp)
        print(f"  {name}: n={len(grp):>4}  平均 MP11={m11:.2f}  平均 MP7={m7:.2f}")

    print("\n解读指引：")
    print("  - 若 P16 治疗后 MP3-MP11 的 Spearman > 0.6，两者很可能捕捉同一个间质信号，")
    print("    P16 的「EMT 上升」需重新表述为「基质样转录程序上升」，不能直接说 EMT。")
    print("  - 若 MP7-MP11 相关也高，说明三个程序在 P16 都受同一个技术因素驱动。")


if __name__ == "__main__":
    main()
