#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fig 1 队列与设计：生成所需汇总表。

产出 results/tables/nc2025_fig1_cohort.csv
列：patient, tp, n_total, n_epi, n_tumor, rcb, response, arm, is_paired

口径：
  tp   取自 Sample 字段（形如 CCG1112_01_T1_A）
  细胞类型取自 Lineage 列（Tumor / Epithelial / 其余）
  RCB 与治疗臂取自 nc2025_rcb_labels.csv（已双向核验）
"""

import csv
import os
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "raw", "nc2025")
TAB = os.path.join(ROOT, "results", "tables")

TUM = os.path.join(RAW, "GEX_OBS_Lineage.csv")
RCB = os.path.join(TAB, "nc2025_rcb_labels.csv")
OUT = os.path.join(TAB, "nc2025_fig1_cohort.csv")

TP_ORDER = ["T1", "T2", "T3", "T4", "T5"]

# 已确认可做克隆分层的 7 例（见 A5第四轮文档）
PAIRED = ["01", "09", "12", "13", "16", "18", "23"]

# 每例的「末次可分析时点」：最后一个 <=T4 且肿瘤细胞足够的活检时点。
# 与 py/13_clone_layering_all.py 的 LAST 保持完全一致——Fig 1 不能画成统一 T4。
LAST = {"01": "T3", "09": "T3", "12": "T3", "13": "T3",
        "16": "T4", "18": "T4", "23": "T2"}


def main():
    rcb, arm, resp = {}, {}, {}
    with open(RCB, encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            p = r["patient"].lstrip("P").zfill(2)
            rcb[p] = r["rcb"]
            arm[p] = r["arm"]
            resp[p] = r["response"]

    cnt = defaultdict(lambda: [0, 0, 0])   # (pat,tp) -> [total, epi, tumor]
    with open(TUM, encoding="utf-8") as f:
        rd = csv.reader(f)
        h = next(rd)
        i_sh, i_lin = h.index("Sample"), h.index("Lineage")
        for row in rd:
            if len(row) <= max(i_sh, i_lin):
                continue
            parts = row[i_sh].split("_")
            if len(parts) < 3:
                continue
            pat, tp = parts[1].zfill(2), parts[2]
            lin = row[i_lin]
            c = cnt[(pat, tp)]
            c[0] += 1
            if lin == "Epithelial":
                c[1] += 1
            if lin == "Tumor":
                c[2] += 1

    rows = []
    pats = sorted({k[0] for k in cnt})
    for p in pats:
        for tp in TP_ORDER:
            c = cnt.get((p, tp))
            base = dict(patient="P" + p, tp=tp,
                        rcb=rcb.get(p, ""), response=resp.get(p, ""),
                        arm=arm.get(p, "") or "not reported",
                        is_paired="yes" if p in PAIRED else "",
                        is_last="yes" if LAST.get(p) == tp else "")
            if not c or c[0] == 0:
                rows.append(dict(base, n_total=0, n_normal_epi=0, n_tumor=0))
                continue
            rows.append(dict(base, n_total=c[0], n_normal_epi=c[1], n_tumor=c[2]))

    os.makedirs(TAB, exist_ok=True)
    with open(OUT, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    # ------- 终端核对 -------
    tot = sum(r["n_total"] for r in rows)
    print(f"患者 {len(pats)} 例，时点 {TP_ORDER}")
    print(f"合计细胞 {tot:,}（应为 249,379）")
    print(f"Tumor 合计 {sum(r['n_tumor'] for r in rows):,}（应为 65,206）")
    print(f"正常上皮合计 {sum(r['n_normal_epi'] for r in rows):,}（应为 28,415）")
    print(f"上皮层合计 {sum(r['n_tumor'] + r['n_normal_epi'] for r in rows):,}"
          f"（应为 93,621）")
    print()
    print(f"{'patient':8s} " + " ".join(f"{tp:>7s}" for tp in TP_ORDER)
          + "   rcb           arm")
    for p in pats:
        line = [r for r in rows if r["patient"] == "P" + p]
        vals = " ".join(f"{r['n_tumor']:7d}" if r["n_total"] else "      ."
                        for r in line)
        print(f"P{p:6s} {vals}   {line[0]['rcb']:<13s} {line[0]['arm']}")

    # 配对可分析性
    print()
    print("配对 7 例（基线 T1 -> 各自末次可分析时点）：")
    for p in PAIRED:
        t1 = next(r["n_tumor"] for r in rows
                  if r["patient"] == "P" + p and r["tp"] == "T1")
        tl = next(r["n_tumor"] for r in rows
                  if r["patient"] == "P" + p and r["tp"] == LAST[p])
        print(f"  P{p}: T1={t1:6d}  {LAST[p]}={tl:6d}")

    print(f"\n-> {os.path.basename(OUT)}  {len(rows)} 行")


if __name__ == "__main__":
    main()
