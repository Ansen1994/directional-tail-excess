"""
时点敏感性分析：判定是否依赖"取最后一个时间点"
------------------------------------------------
对每个 (患者, MP, 所有治疗后时点) 重算判定，比较方向与机制标签的稳定性。
输入：nc2025_v2_l1_trajectory.csv / nc2025_v2_l3_within.csv
输出：results/tables/nc2025_v3_call_all_tp.csv
"""
import csv, os
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(ROOT, "results", "tables")

l1 = list(csv.DictReader(open(os.path.join(TAB, "nc2025_v2_l1_trajectory.csv"))))
l3 = list(csv.DictReader(open(os.path.join(TAB, "nc2025_v2_l3_within.csv"))))

share = {(r["patient"], r["tp"], r["mp"]): r["within_share"]
         for r in l3 if r["subclone"] == "ALL"}

post = defaultdict(list)
for r in l1:
    post[r["patient"]].append(r["tp"])
    post[r["patient"]] = sorted(set(post[r["patient"]]))

rows = []
for r in l1:
    if r["tp"] == "T1":
        continue
    p, tp, m = r["patient"], r["tp"], r["mp"]
    dmed, sz, dh = float(r["d_med"]), float(r["state_z"]), float(r["d_hi"])
    sh = share.get((p, tp, m), "")
    shv = float(sh) if sh not in ("", None) else None
    arrow = "up" if dmed > 0 else "down"
    # 机制判定需先有位移：|Δ中位| < 0.20 时比值（within_share）分母近零，不做机制归类
    if abs(dmed) < 0.20:
        lab = "no change (below effect threshold)"
    elif shv is None:
        lab = f"L1 only ({arrow})"
    elif shv >= 0.5:
        lab = f"induction ({arrow})"
    elif shv < 0.3:
        lab = f"selection ({arrow})"
    else:
        lab = f"mixed ({arrow})"
    rows.append(dict(patient=p, tp=tp, mp=m, d_med=dmed, state_z=sz, d_hi=dh,
                     within_share="" if shv is None else round(shv, 2),
                     direction=arrow, call=lab))

with open(os.path.join(TAB, "nc2025_v3_call_all_tp.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)

print("=== 全部治疗后时点的 Δ中位（正=上升）===")
for m in ["MP7", "MP11"]:
    print(f"\n### {m}")
    hdr = f"{'患者':<6}" + "".join(f"{t:>10}" for t in ["T2", "T3", "T4", "T5"])
    print(hdr)
    byp = defaultdict(dict)
    for r in rows:
        if r["mp"] == m:
            byp[r["patient"]][r["tp"]] = r["d_med"]
    for p in sorted(byp, key=int):
        line = f"{p:<6}"
        for t in ["T2", "T3", "T4", "T5"]:
            v = byp[p].get(t)
            line += f"{v:>+10.2f}" if v is not None else f"{'-':>10}"
        # 方向是否在时点间翻转
        vs = [v for v in byp[p].values() if v is not None]
        flip = (max(vs) > 0.20 and min(vs) < -0.20) if vs else False
        print(line + ("   <-- 方向翻转" if flip else ""))

print("\n=== 判定标签在各时点间的一致性 ===")
n_flip = 0
for m in ["MP7", "MP11"]:
    print(f"\n{m}")
    for p in sorted({r["patient"] for r in rows}, key=int):
        labs = [r["call"] for r in rows if r["mp"] == m and r["patient"] == p]
        if not labs:
            continue
        uniq = sorted(set(labs))
        stable = len(uniq) == 1
        if not stable:
            n_flip += 1
        print(f"  患者{p:<3} " + " | ".join(f"{r['tp']}:{r['call']}"
                                          for r in rows if r["mp"] == m and r["patient"] == p)
              + ("" if stable else "   <-- 不稳定"))
print(f"\n不稳定的 (患者, MP) 组合：{n_flip} / {sum(1 for m in ['MP7','MP11'] for p in {r['patient'] for r in rows})}")
