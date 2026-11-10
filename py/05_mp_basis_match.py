"""
在 nc2025 逐样本 NMF basis 中定位 MP7 / MP11 的代表程序，
并导出其完整基因权重排名（突破 gmt 只给 50 个基因的上限）。

basis 结构：dict[sample] -> DataFrame(7000 HVG x (k x n_programs))
，列名为 "<sample>.k<K>.p<I>"。
"""
import pickle, os, csv, warnings
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "raw", "nc2025")
TAB = os.path.join(ROOT, "results", "tables")

with open(os.path.join(RAW, "nmf_basis.pickle"), "rb") as f:
    basis = pickle.load(f)

gmt = {}
with open(os.path.join(RAW, "MP_Programs.gmt")) as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) >= 4:
            gmt[p[0]] = [g.upper() for g in p[2:] if g]

print(f"basis 样本数 {len(basis)}；每样本程序数示例 "
      f"{basis[list(basis)[0]].shape[1]}\n")

rows, best = [], {}
for mp, genes in gmt.items():
    hits = []
    for s, df in basis.items():
        for c in df.columns:
            top = {g.upper() for g in df[c].nlargest(50).index}
            k = len(top & set(genes))
            hits.append((k, len(top & set(genes)), s, c))
    hits.sort(reverse=True)
    k, n_avail, s, c = hits[0]
    best[mp] = (s, c, k, n_avail)
    print(f"{mp:<5} {c:<30} top50 重叠 {k}/50")

# 导出 MP7 / MP11 代表程序的完整排名
FOCUS = ["MP7", "MP11"]
KEY = ["AGR2", "AGR3", "SCGB2A2", "PIP", "TFF1", "TFF3", "ESR1", "GATA3",
       "AREG", "ERBB4", "INPP4B", "TBC1D9", "STC2", "ANKRD30A", "KRT15", "SOX9"]
for mp in FOCUS:
    s, c, k, n = best[mp]
    v = basis[s][c].sort_values(ascending=False)
    print(f"\n=== {mp} 代表程序 {c}（样本 {s}）top 25 ===")
    print(", ".join(v.index[:25]))
    print("--- 关键基因排名（0 = 无权重）---")
    vn = {g.upper(): g for g in v.index}
    for g in KEY:
        if g in vn:
            val = v[vn[g]]
            rank = int((v > val).sum()) + 1
            print(f"   {g:<10} rank {rank:>5}  weight {val:.4f}")
        else:
            print(f"   {g:<10} 不在 7000 HVG 内")
    for i, (gene, val) in enumerate(v.items()):
        if val <= 0:
            break
        rows.append(dict(mp=mp, program=c, sample=s, rank=i + 1,
                         gene=gene, weight=round(float(val), 6)))

with open(os.path.join(TAB, "nc2025_mp_full_genes.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=["mp", "program", "sample", "rank", "gene", "weight"])
    w.writeheader(); w.writerows(rows)
print(f"\n-> {TAB}/nc2025_mp_full_genes.csv ({len(rows)} 行)")
