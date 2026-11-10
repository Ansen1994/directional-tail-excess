#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
补算：全部 190 对患者的细胞级 KS 距离，与伪 bulk 余弦相似度的关系。

论点：如果 bulk 能反推细胞组成，余弦越高就该 KS 越低（强负相关）。
实测若相关性弱，说明「伪 bulk 相同」推不出「细胞组成相同」。

产出 results/tables/nc2025_fig2_sim_vs_ks.csv
列：p1, p2, cosine, ks_median, ks_max, ks_argmax
"""

import csv
import math
import os
import statistics
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(ROOT, "results", "tables")
CELLS = os.path.join(TAB, "nc2025_fig2_cells.csv")
SIM = os.path.join(TAB, "nc2025_fig2_pseudobulk_sim.csv")
OUT = os.path.join(TAB, "nc2025_fig2_sim_vs_ks.csv")

MPS = ["MP2", "MP8", "MP7", "MP3", "MP11", "MP6", "MP9", "MP4", "MP5"]


def ks(x, y):
    """两样本 KS 距离；输入需已排序"""
    i = j = 0
    nx, ny = len(x), len(y)
    d = 0.0
    while i < nx and j < ny:
        if x[i] <= y[j]:
            i += 1
        else:
            j += 1
        v = abs(i / nx - j / ny)
        if v > d:
            d = v
    return d


def main():
    # patient -> mp -> 排序后的得分列表
    data = defaultdict(lambda: defaultdict(list))
    with open(CELLS, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            data[r["patient"]][r["mp"]].append(float(r["score"]))
    for p in data:
        for m in data[p]:
            data[p][m].sort()
    pats = sorted(data)
    print(f"患者 {len(pats)} 例，每组细胞数中位 "
          f"{int(statistics.median([len(data[p]['MP2']) for p in pats]))}")

    sim = {(r["p1"], r["p2"]): float(r["cosine"])
           for r in csv.DictReader(open(SIM, encoding="utf-8"))}

    rows = []
    for i, a in enumerate(pats):
        for b in pats[i + 1:]:
            c = sim.get((a, b))
            if c is None:
                continue
            kv = {m: ks(data[a][m], data[b][m]) for m in MPS}
            mks = statistics.median(kv.values())
            amax = max(kv, key=lambda m: kv[m])
            rows.append(dict(p1=a, p2=b, cosine=f"{c:.5f}",
                             ks_median=f"{mks:.4f}",
                             ks_max=f"{kv[amax]:.4f}", ks_argmax=amax))
    rows.sort(key=lambda r: -float(r["cosine"]))

    with open(OUT, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"-> {os.path.basename(OUT)}  {len(rows)} 行")

    # 相关性
    x = [float(r["cosine"]) for r in rows]
    y = [float(r["ks_median"]) for r in rows]

    def pearson(a, b):
        n = len(a)
        ma, mb = sum(a) / n, sum(b) / n
        num = sum((p - ma) * (q - mb) for p, q in zip(a, b))
        da = math.sqrt(sum((p - ma) ** 2 for p in a))
        db = math.sqrt(sum((q - mb) ** 2 for q in b))
        return num / (da * db)

    print(f"\n余弦 vs KS 中位：Pearson r = {pearson(x, y):+.3f}")
    print(f"余弦范围 {min(x):.3f}-{max(x):.3f}，KS 中位范围 {min(y):.3f}-{max(y):.3f}")
    print("\n最相似的 5 对：")
    for r in rows[:5]:
        print(f"  {r['p1']} vs {r['p2']}  cos={r['cosine']}  "
              f"KS中位={r['ks_median']}  最大 KS 在 {r['ks_argmax']}={r['ks_max']}")
    print("\n最不相似的 3 对：")
    for r in rows[-3:]:
        print(f"  {r['p1']} vs {r['p2']}  cos={r['cosine']}  KS中位={r['ks_median']}")


if __name__ == "__main__":
    main()
