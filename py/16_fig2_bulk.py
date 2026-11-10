#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fig 2 测量效度：为什么 bulk 反推不了瘤内机制。

四个输出（都写到 results/tables/）：

nc2025_fig2_icc.csv          每个元程序的患者间/患者内方差分解与 ICC
nc2025_fig2_pseudobulk.csv   每个患者的 9 维伪 bulk 向量（细胞均值）
nc2025_fig2_mpcorr.csv       9 x 9 细胞级 Spearman 相关矩阵
nc2025_fig2_cells.csv        细胞级长表（patient, tp, batch, mp, score），供分布图使用

口径：
  - 只用 Lineage=="Tumor" 且带得分的细胞（65,206）
  - pseudo-bulk = 该患者所有肿瘤细胞在 9 个程序上的平均得分。
    NMF 得分是基因表达的线性投影，因此这一均值与"把同样细胞混成 bulk 再打分"等价。
  - batch 取自 barcode 末段（..._Batch1）
"""

import csv
import math
import os
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "raw", "nc2025")
TAB = os.path.join(ROOT, "results", "tables")
SC = os.path.join(RAW, "signature_scores.csv")

MPS = ["MP2", "MP8", "MP7", "MP3", "MP11", "MP6", "MP9", "MP4", "MP5"]


def dump(name, rows):
    if not rows:
        raise SystemExit(f"{name}: 无数据")
    path = os.path.join(TAB, name)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"-> {name}  {len(rows)} 行")


def ranks(v):
    """平均秩，处理并列"""
    idx = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(idx):
        j = i
        while j + 1 < len(idx) and v[idx[j + 1]] == v[idx[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1
        for k in range(i, j + 1):
            r[idx[k]] = avg
        i = j + 1
    return r


def spearman(x, y):
    n = len(x)
    rx, ry = ranks(x), ranks(y)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    dy = math.sqrt(sum((b - my) ** 2 for b in ry))
    return num / (dx * dy) if dx > 0 and dy > 0 else float("nan")


def main():
    # ---------- 读得分 ----------
    meta = {}          # barcode -> (patient, tp, batch)
    vals = defaultdict(list)     # (patient, mp) -> [score]
    allv = defaultdict(list)     # mp -> [(patient, score)]
    pats = set()

    with open(SC, encoding="utf-8") as f:
        rd = csv.reader(f)
        h = next(rd)
        col = {m: h.index(m) for m in MPS if m in h}
        if len(col) != 9:
            raise SystemExit(f"列不匹配：{h}")
        for row in rd:
            if len(row) < len(h):
                continue
            bc = row[0]
            parts = bc.split("_")
            # CTCTAGCTCCGTGACA-1_CCG1112_16_T4_A1_CCG1112_MO_Batch1
            if len(parts) < 7:
                continue
            pat = parts[2].zfill(2)
            tp = parts[3]
            batch = parts[-1]
            pats.add(pat)
            for m, i in col.items():
                v = float(row[i])
                vals[(pat, m)].append(v)
                allv[m].append((pat, v))

    print(f"患者 {len(pats)} 例；总得分条目 {sum(len(v) for v in allv['MP2']):,}")

    # ---------- 1. 方差分解 / ICC ----------
    icc_rows = []
    for m in MPS:
        groups = defaultdict(list)
        for p, v in allv[m]:
            groups[p].append(v)
        groups = {k: g for k, g in groups.items() if len(g) >= 30}
        k = len(groups)
        N = sum(len(g) for g in groups.values())
        grand = sum(sum(g) for g in groups.values()) / N
        ssb = sum(len(g) * (sum(g) / len(g) - grand) ** 2 for g in groups.values())
        ssw = sum(sum((x - sum(g) / len(g)) ** 2 for x in g) for g in groups.values())
        dfb, dfw = k - 1, N - k
        msb, msw = ssb / dfb, ssw / dfw
        n0 = (N - sum(len(g) ** 2 for g in groups.values()) / N) / (k - 1)
        s2b = max((msb - msw) / n0, 0.0)
        s2w = msw
        icc = s2b / (s2b + s2w) if (s2b + s2w) > 0 else float("nan")
        allv_m = [v for _, v in allv[m]]
        mu = sum(allv_m) / len(allv_m)
        sd = math.sqrt(sum((x - mu) ** 2 for x in allv_m) / (len(allv_m) - 1))
        icc_rows.append(dict(
            mp=m, n_cells=N, n_patients=k,
            var_between=f"{s2b:.4f}", var_within=f"{s2w:.4f}",
            icc=f"{icc:.4f}",
            pct_between=f"{100 * icc:.1f}", pct_within=f"{100 * (1 - icc):.1f}",
            cell_sd=f"{sd:.3f}"))
    dump("nc2025_fig2_icc.csv", icc_rows)

    # ---------- 2. 伪 bulk ----------
    pb = {}
    for p in sorted(pats):
        vec = {}
        for m in MPS:
            v = vals.get((p, m), [])
            vec[m] = sum(v) / len(v) if v else float("nan")
        pb[p] = vec

    pb_rows = []
    for p in sorted(pats):
        d = dict(patient="P" + p)
        d.update({m: f"{pb[p][m]:.4f}" for m in MPS})
        d["n_cells"] = len(vals.get((p, "MP2"), []))
        pb_rows.append(d)
    dump("nc2025_fig2_pseudobulk.csv", pb_rows)

    # 伪 bulk 两两余弦相似度
    def cos(a, b):
        num = sum(a[m] * b[m] for m in MPS)
        da = math.sqrt(sum(a[m] ** 2 for m in MPS))
        db = math.sqrt(sum(b[m] ** 2 for m in MPS))
        return num / (da * db)

    plist = sorted(pats)
    best = None
    sim_rows = []
    for i, a in enumerate(plist):
        for b in plist[i + 1:]:
            s = cos(pb[a], pb[b])
            sim_rows.append(dict(p1="P" + a, p2="P" + b, cosine=f"{s:.5f}"))
            if best is None or s > best[2]:
                best = (a, b, s)
    sim_rows.sort(key=lambda r: -float(r["cosine"]))
    dump("nc2025_fig2_pseudobulk_sim.csv", sim_rows)

    # 细胞级差异：同两例，各程序的 KS 距离均值
    def ks(x, y):
        x, y = sorted(x), sorted(y)
        i = j = 0
        d = 0.0
        while i < len(x) and j < len(y):
            if x[i] <= y[j]:
                i += 1
            else:
                j += 1
            d = max(d, abs(i / len(x) - j / len(y)))
        return d

    ks_rows = []
    for a, b in [(best[0], best[1])]:
        for m in MPS:
            ks_rows.append(dict(patient_pair=f"P{a}_vs_P{b}", mp=m,
                                ks=f"{ks(vals[(a, m)], vals[(b, m)]):.4f}",
                                mean_a=f"{pb[a][m]:.3f}", mean_b=f"{pb[b][m]:.3f}"))
    dump("nc2025_fig2_pair_ks.csv", ks_rows)
    print(f"伪 bulk 最相似的一对：P{best[0]} vs P{best[1]}  余弦={best[2]:.4f}")

    # ---------- 3. 程序间相关 ----------
    # 抽样控制计算量
    idx = list(range(len(allv["MP2"])))
    step = max(1, len(idx) // 12000)
    idx = idx[::step]
    mat = {m: [allv[m][i][1] for i in idx] for m in MPS}
    cor_rows = []
    for a in MPS:
        for b in MPS:
            cor_rows.append(dict(mp1=a, mp2=b,
                                 rho=f"{spearman(mat[a], mat[b]):.4f}"))
    dump("nc2025_fig2_mpcorr.csv", cor_rows)
    print(f"相关矩阵用 {len(idx)} 个细胞抽样计算")

    # ---------- 4. 细胞级长表（含 batch） ----------
    cell_rows = []
    with open(SC, encoding="utf-8") as f:
        rd = csv.reader(f)
        h = next(rd)
        col = {m: h.index(m) for m in MPS}
        for row in rd:
            if len(row) < len(h):
                continue
            parts = row[0].split("_")
            if len(parts) < 7:
                continue
            p, tp, batch = parts[2].zfill(2), parts[3], parts[-1]
            for m in MPS:
                cell_rows.append(dict(patient="P" + p, tp=tp, batch=batch,
                                      mp=m, score=f"{float(row[col[m]]):.4f}"))
    dump("nc2025_fig2_cells.csv", cell_rows)

    # 批次 × MP 的均值（查技术信号）
    bm = defaultdict(list)
    for r in cell_rows:
        bm[(r["batch"], r["mp"])].append(float(r["score"]))
    bat_rows = []
    for (b, m), v in sorted(bm.items()):
        bat_rows.append(dict(batch=b, mp=m, n=len(v),
                             mean=f"{sum(v)/len(v):.4f}"))
    dump("nc2025_fig2_batch_mp.csv", bat_rows)

    # ---------- 终端核对 ----------
    print("\nICC 一览：")
    for r in icc_rows:
        print(f"  {r['mp']:<5s} 患者间 {r['pct_between']:>5s}%  "
              f"细胞间 {r['pct_within']:>5s}%  |  细胞 SD = {r['cell_sd']}")
    med = sorted(float(r["icc"]) for r in icc_rows)[len(icc_rows) // 2]
    print(f"  ICC 中位 = {med:.3f}")


if __name__ == "__main__":
    main()
