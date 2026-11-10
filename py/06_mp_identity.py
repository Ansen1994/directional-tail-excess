"""
程序身份比对 v2：用 nc2025 的完整 basis（7000 HVG 全排名），
而不是 gmt 的 50 基因上限，判断 MP7 / MP11 与 GSE176078 各程序是否同构。

背景：gmt 只保留每个 MP 的 50 个基因，导致 MP7 里 rank 31 的 TFF1、
      rank 52 的 SCGB2A2 被切掉，造成"两边无重叠"的假象。
输出：results/tables/mp_identity_v2.csv
"""
import csv, os, math

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(ROOT, "results", "tables")
BG = 20000


def hyper_sf(k, bg, nb, na):
    if k <= 0:
        return 1.0
    lim = min(na, nb)
    tail = sum(math.comb(nb, i) * math.comb(bg - nb, na - i) for i in range(k, lim + 1))
    return tail / math.comb(bg, na)


def bh(p):
    n = len(p); order = sorted(range(n), key=lambda i: p[i])
    q = [0.0] * n; prev = 1.0
    for rank, i in enumerate(reversed(order), start=1):
        prev = min(prev, p[i] * n / (n - rank + 1))
        q[i] = prev
    return q


def read_ours():
    out = {}
    for fn, lab in [("epi_mp_genes.csv", "G:MP"), ("ps_nmf_genes.csv", "PS:")]:
        with open(os.path.join(TAB, fn)) as f:
            r = csv.reader(f); next(r)
            for row in r:
                if len(row) < 3:
                    continue
                prog, gene = row[0].strip('"'), row[2].strip('"').upper()
                out.setdefault(lab + prog, set()).add(gene)
    return out


def read_nc():
    out = {}
    with open(os.path.join(TAB, "nc2025_mp_full_genes.csv")) as f:
        for row in csv.DictReader(f):
            if float(row["weight"]) <= 0:
                continue
            out.setdefault(row["mp"], []).append((int(row["rank"]), row["gene"].upper()))
    return {k: [g for _, g in sorted(v)] for k, v in out.items()}


ours, nc = read_ours(), read_nc()

print("=== 完整排名下，腔面/ER 基因在 MP7 与 MP11 中的位次 ===")
KEYSET = ["TFF1", "SCGB2A2", "GATA3", "AREG", "TBC1D9", "AGR3", "AGR2", "TFF3",
          "ERBB4", "ESR1", "INPP4B", "KRT15", "AZGP1", "PIP", "CLDN4", "TSPAN13"]
print(f"{'gene':<10}{'MP7':>8}{'MP11':>8}")
for g in KEYSET:
    a = nc["MP7"].index(g) + 1 if g in nc["MP7"] else None
    b = nc["MP11"].index(g) + 1 if g in nc["MP11"] else None
    print(f"{g:<10}{(a if a else '-'):>8}{(b if b else '-'):>8}")

rows = []
for N in (30, 50, 100, 200):
    print(f"\n=== top-{N} vs 我方 top-30（nc 取前 {N}，我方取前 30）===")
    for mp in ["MP7", "MP11"]:
        top = set(nc[mp][:N])
        for ours_name, gs in sorted(ours.items()):
            ov = top & gs
            k, na, nb = len(ov), len(gs), min(N, len(nc[mp]))
            p = hyper_sf(k, BG, nb, na)
            rows.append(dict(nc_mp=mp, topN=N, ours=ours_name, n_shared=k,
                             jaccard=round(k / len(top | gs), 4), p=float(f"{p:.2e}"),
                             shared=";".join(sorted(ov))))
            if N == 100 and k >= 3:
                print(f"  {mp} vs {ours_name:<14} 共享 {k:>2}  Jaccard {k/len(top|gs):.3f}  "
                      f"p={p:.1e}  {', '.join(sorted(ov))}")

q = bh([r["p"] for r in rows])
for r, v in zip(rows, q):
    r["fdr"] = float(f"{v:.2e}")

with open(os.path.join(TAB, "mp_identity_v2.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)

print("\n=== 我方腔面分泌程序（PS:k6c1）在 MP7 前 200 中的共享基因 ===")
top = set(nc["MP7"][:200])
print(", ".join(sorted(top & ours["PS:k6c1"])))
print("\n=== 我方 MP1 在 MP7 前 200 中的共享基因 ===")
print(", ".join(sorted(top & ours["G:MPMP1"])))
