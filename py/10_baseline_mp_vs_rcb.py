"""
横断面检验：基线 MP 水平 vs RCB 应答
------------------------------------------------
背景：配对分析只有 7 例、其中仅 1 例 favorable，无法检验"MP7 升高与应答"的关系。
      但 12 例有基线样本且 RCB 已知，可做横断面检验（患者级中位数，非细胞级）。
注意：n=12（3 favorable / 9 unfavorable），功效极低，结果只作方向性参考。
输出：results/tables/nc2025_baseline_mp_rcb.csv
"""
import csv, os, statistics as st
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "raw", "nc2025")
TAB = os.path.join(ROOT, "results", "tables")
MPS = ["MP2", "MP3", "MP4", "MP5", "MP6", "MP7", "MP8", "MP9", "MP11"]

rcb = {r["patient"]: (r["response"], r["rcb"], r["arm"])
       for r in csv.DictReader(open(os.path.join(TAB, "nc2025_rcb_labels.csv")))}

kv = {}
with open(os.path.join(RAW, "gex_cnv_SCEVAN_karyotype.csv")) as f:
    r = csv.reader(f); next(r)
    for row in r:
        if len(row) >= 6:
            kv[f"{row[0]}_{row[4]}_{row[5]}"] = row[1]

cells = defaultdict(lambda: defaultdict(list))
with open(os.path.join(RAW, "signature_scores.csv")) as f:
    r = csv.reader(f); h = next(r)
    idx = {n: i for i, n in enumerate(h) if n}
    for row in r:
        if not row or not row[0]:
            continue
        k = row[0]
        if kv.get(k) != "tumor":
            continue
        p = k.split("_")
        if p[3] != "T1":
            continue
        pid = f"P{p[2]}"
        for m in MPS:
            cells[pid][m].append(float(row[idx[m]]))

MIN = 50
rows = []
for pid in sorted(cells, key=lambda x: int(x[1:])):
    if pid not in rcb or rcb[pid][0] == "unknown":
        continue
    n = len(cells[pid]["MP7"])
    if n < MIN:
        continue
    r = dict(patient=pid, response=rcb[pid][0], rcb=rcb[pid][1],
             arm=rcb[pid][2], n_cells=n)
    for m in MPS:
        r[m] = round(st.median(cells[pid][m]), 3)
    rows.append(r)

with open(os.path.join(TAB, "nc2025_baseline_mp_rcb.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)

print(f"{'患者':<6}{'应答':<13}{'细胞数':>7}" + "".join(f"{m:>7}" for m in MPS))
for r in rows:
    print(f"{r['patient']:<6}{r['response']:<13}{r['n_cells']:>7}"
          + "".join(f"{r[m]:>7.2f}" for m in MPS))


def mw(x, y):
    """Mann-Whitney U 精确近似（正态近似，含并列校正），返回 U 与双侧 p"""
    import math
    n1, n2 = len(x), len(y)
    allv = [(v, 0) for v in x] + [(v, 1) for v in y]
    allv.sort()
    ranks = {}
    i = 0
    while i < len(allv):
        j = i
        while j + 1 < len(allv) and allv[j+1][0] == allv[i][0]:
            j += 1
        avg = (i + j) / 2 + 1
        for k in range(i, j + 1):
            ranks[k] = avg
        i = j + 1
    r1 = sum(ranks[k] for k in range(len(allv)) if allv[k][1] == 0)
    u1 = r1 - n1 * (n1 + 1) / 2
    mu = n1 * n2 / 2
    sd = math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
    z = (u1 - mu) / sd if sd else 0
    p = 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))
    return u1, p


fav = [r for r in rows if r["response"] == "favorable"]
unf = [r for r in rows if r["response"] == "unfavorable"]
print(f"\n=== 基线 MP7：favorable(n={len(fav)}) vs unfavorable(n={len(unf)}) ===")
for m in ["MP7", "MP11"]:
    a = [r[m] for r in fav]
    b = [r[m] for r in unf]
    u, p = mw(a, b)
    print(f"{m}: favorable 中位 {st.median(a):.2f} {sorted(round(x,2) for x in a)} | "
          f"unfavorable 中位 {st.median(b):.2f} {sorted(round(x,2) for x in b)} | "
          f"U={u:.1f} p={p:.3f}")
print("\n注：n=12（3 vs 9），功效极低；p 值仅作方向性参考，不用于声明显著性。")
