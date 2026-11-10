import csv, statistics as st
from collections import defaultdict

BASE = "./nc2025"

def parse(bc):
    p = bc.split("_")
    return p[2], p[3]          # patient, timepoint

# ---- MP scores ----
mps = ["MP2", "MP4", "MP7", "MP11"]
scores = {}                     # bc -> {mp: val}
with open(f"{BASE}/signature_scores.csv") as f:
    r = csv.reader(f); hdr = next(r)
    idx = {n: i for i, n in enumerate(hdr) if n}
    for row in r:
        if not row or not row[0]:
            continue
        scores[row[0]] = {m: float(row[idx[m]]) for m in mps}

# ---- cellstate ----
state = {}
with open(f"{BASE}/Annotation.csv") as f:
    r = csv.reader(f); next(r)
    for row in r:
        if row and row[0]:
            state[row[0]] = row[1]

print(f"cells with MP scores: {len(scores)}   with cellstate: {len(state)}\n")

# ---- per patient x timepoint ----
grp = defaultdict(list)
for bc, s in scores.items():
    p, t = parse(bc)
    grp[(p, t)].append((s, state.get(bc, "NA")))

pts = sorted({p for p, _ in grp})
tps = ["T1", "T2", "T3", "T4", "T5"]

print("=== 恶性细胞数（患者 × 时间点）===")
print("pat  " + "".join(f"{t:>9}" for t in tps))
for p in pts:
    line = f"{p:>3}  "
    for t in tps:
        n = len(grp.get((p, t), []))
        line += f"{n:>9}" if n else f"{'.':>9}"
    print(line)

n_pair = sum(1 for p in pts if ("T1" in [t for q, t in grp if q == p]) and
             any(t in [tt for q, tt in grp if q == p] for t in ["T2", "T3", "T4", "T5"]))
paired = [p for p in pts if any(q == p and t == "T1" for q, t in grp) and
          any(q == p and t in ("T2", "T3", "T4") for q, t in grp)]
print(f"\n有基线(T1)且有后续时间点的患者：{len(paired)} 例 -> {paired}")

# ---- 选择 vs 诱导 的判别指标 ----
print("\n=== MP7 (ER-I) 配对变化：分布位移 vs 高表达占比 ===")
print("pat   tp     n    mean   median  pct>Q75(global)")
allvals = [s["MP7"] for s in scores.values()]
q75 = st.quantiles(allvals, n=4)[2]
print(f"(global MP7 Q75 = {q75:.3f})")
for p in paired:
    for t in ["T1", "T2", "T3", "T4"]:
        rows = grp.get((p, t))
        if not rows or len(rows) < 30:
            continue
        v = [x[0]["MP7"] for x in rows]
        frac = sum(1 for x in v if x > q75) / len(v) * 100
        print(f"{p:>3}  {t}  {len(v):>5}  {st.mean(v):6.2f}  {st.median(v):6.2f}  {frac:5.1f}%")
    print()

print("=== MP11 (EMT-III) 配对变化 ===")
q75b = st.quantiles([s["MP11"] for s in scores.values()], n=4)[2]
print("pat   tp     n    mean   median  pct>Q75")
for p in paired:
    for t in ["T1", "T2", "T3", "T4"]:
        rows = grp.get((p, t))
        if not rows or len(rows) < 30:
            continue
        v = [x[0]["MP11"] for x in rows]
        frac = sum(1 for x in v if x > q75b) / len(v) * 100
        print(f"{p:>3}  {t}  {len(v):>5}  {st.mean(v):6.2f}  {st.median(v):6.2f}  {frac:5.1f}%")
    print()

# ---- cellstate 组成 ----
print("=== Tumor cellstate 组成（按时间点汇总，%）===")
comp = defaultdict(lambda: defaultdict(int))
for bc, s in scores.items():
    p, t = parse(bc)
    comp[t][state.get(bc, "NA")] += 1
allstates = sorted({k for d in comp.values() for k in d})
print("state".ljust(18) + "".join(f"{t:>8}" for t in tps))
for sname in allstates:
    line = sname.ljust(18)
    for t in tps:
        tot = sum(comp[t].values()) or 1
        line += f"{comp[t][sname]/tot*100:7.1f}%"
    print(line)
