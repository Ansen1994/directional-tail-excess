"""
检验 P26 是否与 P16 为同一患者的重复取样（P06 与 P08 同理）
--------------------------------------------------------------
理由：Fig1B 把 P06/P26 标在 P08/P16 下方（"同患者样本"标注），
      且补充数据 Baseline n=12 而数据里有基线样本的编号有 13 个，
      差 1 例正好可由 P26≡P16 解释。
检验：比较各样本对的 9 个 MP 评分分布（中位数绝对差 + KS 距离），
      以 P26 vs P16 为假设对，其余样本对为背景分布。
"""
import csv, os
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "raw", "nc2025")

vals = defaultdict(lambda: defaultdict(list))
with open(os.path.join(RAW, "signature_scores.csv")) as f:
    r = csv.reader(f); h = next(r)
    mps = h[1:]
    for row in r:
        bc = row[0]
        parts = bc.split("_")
        # barcode 形如 <BC>_CCG1112_16_T4_A1_CCG1112_MO_Batch1
        # parts[1]=批次前缀, parts[2]=患者号, parts[3]=时间点
        if len(parts) < 4 or not parts[1].lower().startswith("ccg"):
            continue
        pid, tp = parts[2], parts[3]
        key = f"P{int(pid):02d}.{tp}"
        for i, v in enumerate(row[1:]):
            if v:
                vals[mps[i]][key].append(float(v))

samples = sorted({k for m in vals for k in vals[m]})
print(f"样本 {len(samples)} 个；MP {len(mps)} 个\n")


def ks(a, b):
    a, b = sorted(a), sorted(b)
    n, m = len(a), len(b)
    i = j = 0
    d = 0.0
    while i < n and j < m:
        if a[i] <= b[j]:
            i += 1
        else:
            j += 1
        d = max(d, abs(i / n - j / m))
    return d


def pair_dist(s1, s2):
    md, ksd = [], []
    for m in mps:
        a, b = vals[m].get(s1, []), vals[m].get(s2, [])
        if len(a) < 30 or len(b) < 30:
            continue
        a_sorted, b_sorted = sorted(a), sorted(b)
        md.append(abs(a_sorted[len(a)//2] - b_sorted[len(b)//2]))
        ksd.append(ks(a, b))
    return (sum(md)/len(md) if md else None), (sum(ksd)/len(ksd) if ksd else None)


targets = [("P16.T1", "P26.T1"), ("P08.T1", "P06.T4"), ("P16.T1", "P16.T4")]
print("=== 假设的重复取样对 ===")
for a, b in targets:
    if a in samples and b in samples:
        md, kd = pair_dist(a, b)
        print(f"  {a} vs {b}: 中位差均值 {md:.3f}  KS 均值 {kd:.3f}  "
              f"(n={len(vals['MP7'].get(a,[]))}/{len(vals['MP7'].get(b,[]))})")

print("\n=== 背景：同患者不同时点 ===")
same = [("P01.T1", "P01.T2"), ("P09.T1", "P09.T3"), ("P18.T1", "P18.T3"),
        ("P13.T1", "P13.T3"), ("P23.T1", "P23.T2")]
for a, b in same:
    if a in samples and b in samples:
        md, kd = pair_dist(a, b)
        print(f"  {a} vs {b}: 中位差均值 {md:.3f}  KS 均值 {kd:.3f}")

print("\n=== 背景：不同患者 ===")
diff = [("P16.T1", "P09.T1"), ("P26.T1", "P09.T1"), ("P08.T1", "P16.T1"),
        ("P16.T1", "P18.T1"), ("P26.T1", "P01.T1")]
for a, b in diff:
    if a in samples and b in samples:
        md, kd = pair_dist(a, b)
        print(f"  {a} vs {b}: 中位差均值 {md:.3f}  KS 均值 {kd:.3f}")
