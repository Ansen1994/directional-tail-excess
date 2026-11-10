"""
A5 主臂 v2：诱导 vs 选择 判别（修正版）

相对 v1 的修正：
  1. 患者 12 基线亚克隆全为空值 -> 仅做 L1，剔除 L2/L3
  2. 不再只取最后一个时间点，输出全部可用时间点的轨迹
  3. 把"状态位移"（中位数）与"尾部扩张"（高表达占比）拆成两个独立指标
     —— 选择的核心特征是尾部扩张，中位数往往不动，v1 会误判为"无变化"
  4. 对治疗后消失的克隆，回溯其基线 MP 水平，检验选择的方向是否与程序变化一致

输出：results/tables/nc2025_v2_*.csv
"""
import csv, os, random, statistics as st
from collections import defaultdict

RAW, OUT = "data/raw/nc2025", "results/tables"
os.makedirs(OUT, exist_ok=True)
MPS = ["MP2", "MP3", "MP4", "MP5", "MP6", "MP7", "MP8", "MP9", "MP11"]
TPS = ["T1", "T2", "T3", "T4", "T5"]
MIN_CELLS = 100      # 时间点进入分析的最低恶性细胞数
MIN_CLONE = 30       # 克隆内比较的最低细胞数

kv = {}
with open(f"{RAW}/gex_cnv_SCEVAN_karyotype.csv") as f:
    r = csv.reader(f); next(r)
    for row in r:
        if len(row) >= 6:
            kv[f"{row[0]}_{row[4]}_{row[5]}"] = (row[1], row[3])

cells = []
with open(f"{RAW}/signature_scores.csv") as f:
    r = csv.reader(f); hdr = next(r)
    idx = {n: i for i, n in enumerate(hdr) if n}
    for row in r:
        if not row or not row[0]:
            continue
        k = row[0]
        if k not in kv:
            continue
        cls, sc = kv[k]
        if cls != "tumor":
            continue
        p = k.split("_")
        cells.append((p[2], p[3], {m: float(row[idx[m]]) for m in MPS}, sc))

def q(v, frac):
    v = sorted(v); k = (len(v) - 1) * frac
    lo = int(k); hi = min(lo + 1, len(v) - 1)
    return v[lo] + (v[hi] - v[lo]) * (k - lo)

# 固定参照：全部基线细胞的 Q75
REF = {m: q([c[2][m] for c in cells if c[1] == "T1"], .75) for m in MPS}

by_pt = defaultdict(list)
for c in cells:
    by_pt[(c[0], c[1])].append(c)

pats = sorted({c[0] for c in cells})
paired = []
for p in pats:
    tps = [t for t in TPS if len(by_pt.get((p, t), [])) >= MIN_CELLS]
    if "T1" in tps and len(tps) > 1:
        paired.append((p, tps))

print(f"恶性细胞 {len(cells)}；配对患者 {len(paired)} 例")
print("基线参照 Q75: " + ", ".join(f"{m}={REF[m]:.3f}" for m in ["MP7", "MP11"]))

# ---------------- L1 轨迹 ----------------
l1 = []
for p, tps in paired:
    base = [c[2] for c in by_pt[(p, "T1")]]
    for t in tps:
        cs = by_pt[(p, t)]
        for m in MPS:
            v = [c[2][m] for c in cs]
            vb = [x[m] for x in base]
            iqr_b = q(vb, .75) - q(vb, .25)
            l1.append(dict(
                patient=p, tp=t, mp=m, n=len(v),
                med=round(st.median(v), 3), q25=round(q(v, .25), 3), q75=round(q(v, .75), 3),
                med_base=round(st.median(vb), 3), iqr_base=round(iqr_b, 3),
                d_med=round(st.median(v) - st.median(vb), 3),
                state_z=round((st.median(v) - st.median(vb)) / iqr_b, 3) if iqr_b else "",
                hi_frac=round(sum(1 for x in v if x > REF[m]) / len(v) * 100, 1),
                hi_frac_base=round(sum(1 for x in vb if x > REF[m]) / len(vb) * 100, 1),
                d_hi=round((sum(1 for x in v if x > REF[m]) / len(v)
                            - sum(1 for x in vb if x > REF[m]) / len(vb)) * 100, 1)))

# ---------------- L2 克隆组成 ----------------
l2 = []
for p, tps in paired:
    base_cs = by_pt[(p, "T1")]
    if any(c[3] == "" for c in base_cs):
        continue                                    # 基线无克隆标注，跳过
    b = defaultdict(int)
    for c in base_cs:
        b[c[3]] += 1
    tb = sum(b.values())
    for t in tps:
        a = defaultdict(int)
        for c in by_pt[(p, t)]:
            a[c[3]] += 1
        ta = sum(a.values())
        tvd = sum(abs(a[s]/ta*100 - b.get(s, 0)/tb*100) for s in set(a) | set(b)) / 2
        for s in sorted(set(a) | set(b)):
            l2.append(dict(patient=p, tp=t, subclone=s,
                           n_base=b.get(s, 0), n_tp=a.get(s, 0),
                           pct_base=round(b.get(s, 0)/tb*100, 1),
                           pct_tp=round(a.get(s, 0)/ta*100, 1),
                           d_pct=round(a.get(s, 0)/ta*100 - b.get(s, 0)/tb*100, 1),
                           tvd=round(tvd, 1),
                           lost=(b.get(s, 0) >= 50 and a.get(s, 0) == 0)))

# 消失克隆的基线 MP 水平（检验选择方向）
lost_mp = []
for p, tps in paired:
    base_cs = by_pt[(p, "T1")]
    for s in sorted({c[3] for c in base_cs}):
        vb = [c[2]["MP7"] for c in base_cs if c[3] == s]
        if len(vb) < 50:
            continue
        alive = [t for t in tps if t != "T1" and by_pt[(p, t)]]
        if not alive:
            continue
        last = alive[-1]
        na = sum(1 for c in by_pt[(p, last)] if c[3] == s)
        lost_mp.append(dict(patient=p, subclone=s, tp_last=last,
                            n_base=len(vb), n_last=na,
                            med_MP7_base=round(st.median(vb), 3),
                            fate=("消失" if na == 0 else ("扩张" if na/len(vb) > 1.3 else "维持"))))

# ---------------- L3 克隆内 ----------------
random.seed(42)
def perm_p(vb, va, n=1000):
    obs = abs(st.median(va) - st.median(vb))
    pool = vb + va; nb = len(vb); cnt = 0
    for _ in range(n):
        random.shuffle(pool)
        if abs(st.median(pool[nb:]) - st.median(pool[:nb])) >= obs:
            cnt += 1
    return (cnt + 1) / (n + 1)

l3 = []
for p, tps in paired:
    for m in ["MP7", "MP11"]:
        for t in tps[1:]:
            cs_b = by_pt[(p, "T1")]; cs_a = by_pt[(p, t)]
            vb = [c[2][m] for c in cs_b]; va = [c[2][m] for c in cs_a]
            d_all = st.median(va) - st.median(vb)
            grp = defaultdict(lambda: ([], []))
            for c in cs_b:
                grp[c[3]][0].append(c[2][m])
            for c in cs_a:
                grp[c[3]][1].append(c[2][m])
            cl = []
            for s, (x, y) in grp.items():
                if len(x) < MIN_CLONE or len(y) < MIN_CLONE:
                    continue
                pp = perm_p(x, y)
                cl.append(dict(patient=p, tp=t, mp=m, subclone=s,
                               n_base=len(x), n_tp=len(y),
                               med_base=round(st.median(x), 3),
                               med_tp=round(st.median(y), 3),
                               d_med=round(st.median(y) - st.median(x), 3),
                               p_perm=round(pp, 3)))
            l3.extend(cl)
            if cl:
                w = st.median([c["d_med"] for c in cl])
                l3.append(dict(patient=p, tp=t, mp=m, subclone="ALL",
                               n_base=len(vb), n_tp=len(va),
                               med_base=round(st.median(vb), 3),
                               med_tp=round(st.median(va), 3),
                               d_med=round(d_all, 3),
                               p_perm="",
                               within_med=round(w, 3),
                               within_share=round(w / d_all, 2) if d_all else "",
                               n_clone=len(cl),
                               n_sig=sum(1 for c in cl if c["p_perm"] < 0.05)))

# ---------------- 判定（最后一时间点） ----------------
call = []
for p, tps in paired:
    last = tps[-1]
    for m in ["MP7", "MP11"]:
        r1 = [x for x in l1 if x["patient"] == p and x["tp"] == last and x["mp"] == m][0]
        r3 = [x for x in l3 if x["patient"] == p and x["tp"] == last
              and x["mp"] == m and x["subclone"] == "ALL"]
        share = r3[0]["within_share"] if r3 else ""
        sz, dh = r1["state_z"], r1["d_hi"]
        dmed = r1["d_med"]
        arrow = "↑" if dmed > 0 else "↓"
        if abs(dmed) < 0.20 and abs(dh) < 10:
            lab = "无变化"
        elif share == "":
            lab = f"仅L1({arrow})"
        elif share >= 0.5:
            lab = f"诱导{arrow}"
        elif share < 0.3:
            lab = f"选择{arrow}"
        else:
            lab = f"混合{arrow}"
        call.append(dict(patient=p, tp=last, mp=m, d_med=dmed,
                         state_z=sz, d_hi=dh,
                         hi_base=r1["hi_frac_base"], hi_tp=r1["hi_frac"],
                         within_share=share, call=lab))

def dump(name, rows):
    if not rows:
        return
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(f"{OUT}/{name}", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    print(f"  -> {OUT}/{name}  ({len(rows)} 行)")

print("\n=== 判定（二元指标：state_z = 中位位移/IQR；d_hi = 高表达占比变化 pp）===")
print("患者  时点  MP     Δ中位   state_z  d_hi(pp)  高表达%        克隆内share  判定")
for c in call:
    print(f"{c['patient']:>3}  {c['tp']}   {c['mp']:<5} {c['d_med']:+6.2f}  "
          f"{c['state_z']:+6.2f}  {c['d_hi']:+6.1f}   "
          f"{c['hi_base']:5.1f}→{c['hi_tp']:5.1f}   "
          f"{str(c['within_share']):>6}      {c['call']}")

print("\n=== 消失克隆的基线 MP7（检验选择方向）===")
for r in sorted(lost_mp, key=lambda x: (x["patient"], x["subclone"])):
    if r["fate"] != "维持":
        print(f"  患者{r['patient']} 克隆{r['subclone']}  基线{r['n_base']:>4} → 末次{r['n_last']:>4}  "
              f"基线MP7中位={r['med_MP7_base']:.2f}  {r['fate']}")

dump("nc2025_v2_l1_trajectory.csv", l1)
dump("nc2025_v2_l2_clone.csv", l2)
dump("nc2025_v2_l3_within.csv", l3)
dump("nc2025_v2_call.csv", call)
dump("nc2025_v2_lostclone.csv", lost_mp)
print("\n完成。")
