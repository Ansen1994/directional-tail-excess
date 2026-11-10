#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fig 2 C/D/E：批次效应的量级与配对设计是否免疫。

背景：nc2025 的 6 个批次与患者高度共线——
  Batch1: P03/P26/P08/P13 (+P16 T4)
  Batch2: P16 (T1, T5)
  Batch3: P04/P01/P02
  Batch4: P21/P12/P14/P07/P27
  Batch5: P18/P15/P06
  Batch6: P23/P09/P10/P24
7 例配对中，6 例的全部时点落在同一批次（批次效应在配对差值中被抵消），
只有 P16 跨 Batch2(T1) -> Batch1(T4)。

本脚本量化三件事：
  1. 每个 MP 的批次间 SD，与配对 Δ 的量级比较（batch_sd vs |delta|）
  2. 批次中心化后，7 例 x 9 MP 的 Δ 方向是否改变
  3. P16 单例的批次归因（Δ 中有多少可由批均值差解释）

产出：
  nc2025_fig2_batch_effect.csv   每 MP：批次间 SD、患者间 SD、|Δ| 中位、比值
  nc2025_fig2_batch_delta.csv    7 例 x 9 MP：原始 Δ vs 批次校正 Δ
  nc2025_fig2_batch_map.csv      批次-患者-时点 细胞数（画共线性图用）
"""

import csv
import math
import os
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(ROOT, "results", "tables")
CELLS = os.path.join(TAB, "nc2025_fig2_cells.csv")

MPS = ["MP2", "MP8", "MP7", "MP3", "MP11", "MP6", "MP9", "MP4", "MP5"]
PAIRED = ["01", "09", "12", "13", "16", "18", "23"]
LAST = {"01": "T3", "09": "T3", "12": "T3", "13": "T3",
        "16": "T4", "18": "T4", "23": "T2"}
MIN_CELLS = 30


def dump(name, rows):
    path = os.path.join(TAB, name)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"-> {name}  {len(rows)} 行")


def median(v):
    s = sorted(v)
    n = len(s)
    return s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2


def sd(v):
    if len(v) < 2:
        return float("nan")
    m = sum(v) / len(v)
    return math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1))


def main():
    # (patient, tp, batch, mp) -> [scores]
    bucket = defaultdict(list)
    # (patient, tp, batch) -> 细胞数
    celln = defaultdict(int)

    with open(CELLS, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            bucket[(r["patient"], r["tp"], r["batch"], r["mp"])].append(
                float(r["score"]))
            celln[(r["patient"], r["tp"], r["batch"])] += 1

    # ---------- 批次-患者映射 ----------
    map_rows = [dict(batch=k[2], patient=k[0], tp=k[1], n_cells=v)
                for k, v in celln.items()]
    map_rows.sort(key=lambda r: (r["batch"], r["patient"], r["tp"]))
    dump("nc2025_fig2_batch_map.csv", map_rows)

    # ---------- 批次均值（按细胞数加权）----------
    # batch -> mp -> 加权均值
    bsum = defaultdict(lambda: defaultdict(float))
    bwt = defaultdict(lambda: defaultdict(int))
    for (p, tp, b, m), v in bucket.items():
        bsum[b][m] += sum(v)
        bwt[b][m] += len(v)
    bmean = {b: {m: bsum[b][m] / bwt[b][m] for m in MPS if bwt[b][m] > 0}
             for b in bsum}

    # ---------- 每例每时点的 MP 中位数 ----------
    med = defaultdict(dict)
    for (p, tp, b, m), v in bucket.items():
        if len(v) >= MIN_CELLS:
            med[(p, tp)][m] = median(v)

    # ---------- 1. batch_sd vs |delta| ----------
    eff_rows = []
    for m in MPS:
        bvals = [bmean[b][m] for b in bmean if m in bmean[b]]
        b_sd = sd(bvals)
        if not b_sd:
            b_sd = 0.0
        # 患者间（基线时点，跨患者）
        base = [med[(p, "T1")][m] for p in {k[0] for k in med}
                if (p, "T1") in med and m in med[(p, "T1")]]
        p_sd = sd(base) if len(base) >= 3 else float("nan")
        # 配对 Δ
        deltas = []
        for p in PAIRED:
            pp = "P" + p
            a = med.get((pp, "T1"), {}).get(m)
            bb = med.get((pp, LAST[p]), {}).get(m)
            if a is not None and bb is not None:
                deltas.append(bb - a)
        d_abs = median([abs(x) for x in deltas]) if deltas else float("nan")
        d_max = max([abs(x) for x in deltas]) if deltas else float("nan")
        eff_rows.append(dict(
            mp=m,
            batch_sd=f"{b_sd:.4f}",
            patient_sd=f"{p_sd:.4f}" if p_sd == p_sd else "",
            delta_median_abs=f"{d_abs:.4f}" if d_abs == d_abs else "",
            delta_max_abs=f"{d_max:.4f}" if d_max == d_max else "",
            ratio_delta_to_batch=f"{d_max / b_sd:.2f}" if b_sd > 0 and d_max == d_max else ""))
    dump("nc2025_fig2_batch_effect.csv", eff_rows)

    # ---------- 2. 批次中心化后的 Δ ----------
    d_rows = []
    for p in PAIRED:
        pp = "P" + p
        for m in MPS:
            b0 = _batch_of(celln, pp, "T1")
            b1 = _batch_of(celln, pp, LAST[p])
            v0 = bucket.get((pp, "T1", b0, m))
            v1 = bucket.get((pp, LAST[p], b1, m))
            if not v0 or not v1 or len(v0) < MIN_CELLS or len(v1) < MIN_CELLS:
                continue
            raw = median(v1) - median(v0)
            corr = (median(v1) - bmean[b1][m]) - (median(v0) - bmean[b0][m])
            d_rows.append(dict(
                patient="P" + p, tp=LAST[p], mp=m,
                batch_base=b0, batch_last=b1, cross_batch="yes" if b0 != b1 else "",
                delta_raw=f"{raw:+.4f}", delta_batchcorr=f"{corr:+.4f}",
                flip="yes" if raw * corr < 0 else ""))
    dump("nc2025_fig2_batch_delta.csv", d_rows)

    # ---------- 3. P16 的批次归因 ----------
    print("\n== 批次效应量级（每 MP）==")
    print(f"{'MP':<5s} {'批次SD':>8s} {'患者SD':>8s} {'|Δ|中位':>8s} {'|Δ|最大':>8s} {'最大Δ/批次SD':>12s}")
    for r in eff_rows:
        print(f"{r['mp']:<5s} {r['batch_sd']:>8s} {r['patient_sd']:>8s} "
              f"{r['delta_median_abs']:>8s} {r['delta_max_abs']:>8s} "
              f"{r['ratio_delta_to_batch']:>12s}")

    n_flip = sum(1 for r in d_rows if r["flip"] == "yes")
    print(f"\n批次中心化后 Δ 方向翻转：{n_flip}/{len(d_rows)}")
    print("P16 的 Δ（唯一跨批次病例）：")
    for r in d_rows:
        if r["patient"] == "P16":
            print(f"  {r['mp']:<5s} raw={r['delta_raw']:>9s}  "
                  f"batch-corr={r['delta_batchcorr']:>9s}  {r['batch_base']}->{r['batch_last']}")

    # P16 MP3 归因
    b0, b1 = _batch_of(celln, "P16", "T1"), _batch_of(celln, "P16", "T4")
    shift = bmean[b1]["MP3"] - bmean[b0]["MP3"]
    raw = med[("P16", "T4")]["MP3"] - med[("P16", "T1")]["MP3"]
    print(f"\nP16 MP3：Δ={raw:+.3f}，批次均值差({b0}->{b1})={shift:+.3f}，"
          f"批次可解释 {100 * shift / raw:.1f}%")


def _batch_of(celln, p, tp):
    hits = [(v, k[2]) for k, v in celln.items() if k[0] == p and k[1] == tp]
    return max(hits)[1] if hits else ""


def _find_batch(bucket, p, tp):
    for k in bucket:
        if k[0] == p and k[1] == tp:
            return k[2]
    return ""


if __name__ == "__main__":
    main()
