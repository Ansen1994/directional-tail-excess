"""
RCB 疗效标签的建立与验证
------------------------------------------------
来源：Fig 1B（逐样本临床注释矩阵）读图，非公开表格（Supplementary Data 4 只有汇总）。
验证：用补充数据 4 的"基线可分析 12 例 = 3 favorable + 9 unfavorable"做交叉核对。
输出：results/tables/nc2025_rcb_labels.csv
"""
import csv, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(ROOT, "results", "tables")

# 读自 Fig 1B：列按 "先 RCB 0-I 后 RCB II-III、组内按治疗臂" 排序
FAV = ["P01", "P07", "P15", "P08", "P27"]
UNF = ["P10", "P13", "P16", "P24", "P02", "P03", "P04", "P09", "P12", "P14", "P18", "P21", "P23"]
ARM = {p: "Chemo->Combo" for p in ["P01", "P07", "P15", "P08", "P27",
                                   "P10", "P13", "P16", "P24"]}
ARM.update({p: "ICI->Combo" for p in ["P02", "P03", "P04", "P09", "P12", "P14", "P18", "P21", "P23"]})

# 数据中实际存在的患者及其时点
tp = {}
with open(os.path.join(ROOT, "data", "raw", "nc2025", "GEX_OBS_Lineage.csv")) as f:
    r = csv.reader(f); next(r)
    for row in r:
        pid, t = row[1].split(".")
        tp.setdefault(pid, set()).add(t)

rows = []
for p in sorted(tp, key=lambda x: int(x[1:])):
    if p in FAV:
        rcb, grp = "RCB 0-I", "favorable"
    elif p in UNF:
        rcb, grp = "RCB II-III", "unfavorable"
    else:
        rcb, grp = "not shown", "unknown"
    rows.append(dict(patient=p, rcb=rcb, response=grp, arm=ARM.get(p, ""),
                     timepoints=";".join(sorted(tp[p])), has_T1="yes" if "T1" in tp[p] else "no"))

with open(os.path.join(TAB, "nc2025_rcb_labels.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)

print(f"{'患者':<6}{'RCB':<12}{'应答':<13}{'治疗臂':<15}{'时点':<20}T1基线")
for r in rows:
    print(f"{r['patient']:<6}{r['rcb']:<12}{r['response']:<13}{r['arm']:<15}"
          f"{r['timepoints']:<20}{r['has_T1']}")

# ---- 交叉验证：有基线且 RCB 已知者应为 12 例 = 3 favorable + 9 unfavorable ----
base_t1 = [r for r in rows if r["has_T1"] == "yes"]
base = [r for r in base_t1 if r["response"] != "unknown"]
nf = sum(1 for r in base if r["response"] == "favorable")
nu = sum(1 for r in base if r["response"] == "unfavorable")
print(f"\n=== 交叉验证（对照 Supplementary Data 4） ===")
print(f"有基线样本的编号 {len(base_t1)} 个，其中 RCB 已知 {len(base)} 个："
      f"favorable {nf} / unfavorable {nu}")
print("Supplementary Data 4 报告：Baseline n=12，3 favorable (RCB 0-I) / 9 unfavorable (RCB II-III)")
ok = (len(base) == 12 and nf == 3 and nu == 9)
print("判定：" + ("一致（差出的 1 例为图中未标注 RCB 的 P26）" if ok else "不一致，需复核"))

print("\n=== 配对分析所用 7 例的应答归属 ===")
paired = ["P01", "P09", "P12", "P13", "P16", "P18", "P23"]
for p in paired:
    r = [x for x in rows if x["patient"] == p][0]
    print(f"  {p}  {r['rcb']:<12}{r['response']:<13}{r['arm']}")
