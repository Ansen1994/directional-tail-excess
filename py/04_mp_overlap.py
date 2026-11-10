"""
两套 NMF 程序的基因集交叉比对
------------------------------------------------
目的：判断 nc2025 的 MP7 (ER-I) / MP11 (EMT-III) 与 GSE176078 的
      MP1（腔面-ER+）/ 逐样本 k6c1（腔面分泌）是否同构。
判据：超几何富集（背景 20000 编码基因）+ Jaccard + 共享基因清单。
输出：results/tables/mp_overlap_matrix.csv, mp_overlap_pairs.csv
"""
import csv, os
import math

def hyper_sf(k, bg, nb, na):
    """P(X >= k)，X ~ Hypergeom(bg, nb, na)；用精确整数组合数"""
    if k <= 0:
        return 1.0
    lim = min(na, nb)
    tail = sum(math.comb(nb, i) * math.comb(bg - nb, na - i) for i in range(k, lim + 1))
    return tail / math.comb(bg, na)

def bh_fdr(p):
    n = len(p); order = sorted(range(n), key=lambda i: p[i])
    q = [0.0] * n; prev = 1.0
    for rank, i in enumerate(reversed(order), start=1):
        idx = n - rank + 1
        prev = min(prev, p[i] * n / idx)
        q[i] = prev
    return q

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(ROOT, "results", "tables")
GMT = os.path.join(ROOT, "data", "raw", "nc2025", "MP_Programs.gmt")
BG = 20000

def read_ours():
    out = {}
    for fn, label in [("epi_mp_genes.csv", "G176078"), ("ps_nmf_genes.csv", "G176078ps")]:
        with open(os.path.join(TAB, fn)) as f:
            r = csv.reader(f); next(r)
            for row in r:
                if len(row) < 3:
                    continue
                prog, _, gene = row[0].strip('"'), row[1], row[2].strip('"')
                out.setdefault(f"{label}:{prog}", set()).add(gene.upper())
    return out

def read_nc():
    out = {}
    with open(GMT) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            out[f"nc2025:{p[0]}"] = {g.upper() for g in p[2:] if g}
    return out

ours, nc = read_ours(), read_nc()
print(f"我们的程序集 {len(ours)} 个；nc2025 的 MP {len(nc)} 个\n")

rows, pairs = [], []
for a, ga in sorted(ours.items()):
    for b, gb in sorted(nc.items()):
        inter = ga & gb
        k, na, nb = len(inter), len(ga), len(gb)
        p = hyper_sf(k, BG, nb, na)
        jac = k / len(ga | gb)
        rows.append(dict(ours=a, nc=b, n_ours=na, n_nc=nb, n_shared=k,
                         jaccard=round(jac, 4), p=float(f"{p:.3e}")))
        pairs.append((a, b, inter, p, jac))

pv = bh_fdr([r["p"] for r in rows])
for r, q in zip(rows, pv):
    r["fdr"] = float(f"{q:.3e}")

rows.sort(key=lambda r: (-r["n_shared"], r["p"]))
with open(os.path.join(TAB, "mp_overlap_matrix.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)

print(f"{'我们的程序':<22}{'nc MP':<10}{'共享':>4}{'Jaccard':>9}{'FDR':>11}   共享基因")
for r in rows[:18]:
    key = (r["ours"], r["nc"])
    genes = [g for (a, b, g, _, _) in pairs if a == key[0] and b == key[1]][0]
    print(f"{r['ours']:<22}{r['nc'].replace('nc2025:',''):<10}{r['n_shared']:>4}"
          f"{r['jaccard']:>9.3f}{r['fdr']:>11.1e}   {', '.join(sorted(genes))[:70]}")

print("\n--- 我们 4 个全局程序的身份对照 ---")
for key in ["G176078:MP1", "G176078:MP2", "G176078:MP3", "G176078:MP4"]:
    sub = [r for r in rows if r["ours"] == key][:3]
    print(key.replace("G176078:", ""), "→", ", ".join(
        f"{r['nc'].replace('nc2025:','')}({r['n_shared']},{r['fdr']:.0e})" for r in sub))

print("\n--- MP7 与 MP11 在全部 7 个我方程序中的最佳匹配 ---")
for mp in ["MP7", "MP11"]:
    sub = [r for r in rows if r["nc"] == f"nc2025:{mp}"]
    sub.sort(key=lambda r: -r["n_shared"])
    print(f"{mp}:", ", ".join(f"{r['ours']}({r['n_shared']})" for r in sub))
