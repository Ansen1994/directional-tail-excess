# -*- coding: utf-8 -*-
"""
14_counterfactual_composition.py

把「克隆组成变化」对总体位移的贡献定量分解出来。

问题
----
13 号脚本用分离度判断位移是否越出基线克隆范围，P09 是唯一不分离的病例。
但「不分离」只说明不能排除克隆选择，没有回答「克隆选择能解释多少」。

反事实
------
构造一个「只有克隆组成变化、克隆内部完全不转变」的虚拟随访样本：
  对每个克隆 c，从该克隆的基线分布中有放回抽 n_post(c) 个细胞，
  n_post(c) 是随访时点该克隆的实际细胞数。
合并后算中位，即为反事实中位。
  反事实位移 = 反事实中位 - 基线全体中位
  实测位移   = 随访全体中位 - 基线全体中位
  组成解释率 = 反事实位移 / 实测位移

读法
----
  解释率 ≈ 0     位移几乎全部来自克隆内状态转变
  解释率 0.3-0.7 两种机制共存
  解释率 ≈ 1     位移可完全由克隆组成变化解释（即克隆选择）
  解释率 < 0     克隆组成变化的方向与实测位移相反，克隆选择不可能是原因

注意：反事实保留了随访时点的克隆占比，所以它同时吸收了取样偏差
（某克隆没被取到会被当成占比下降）。对本队列 20 例规模的活检，这个
偏差不可忽略，故解释率只作上限意义下的参考，不作精确分解。
"""
import csv
import os
import random
from collections import defaultdict

BASE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                    "data", "raw", "nc2025")
TAB = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "results", "tables")

PAIRED = ["01", "09", "12", "13", "16", "18", "23"]
LAST = {"01": "T3", "09": "T3", "12": "T3", "13": "T3",
        "16": "T4", "18": "T4", "23": "T2"}
RCB = {"01": "favorable", "09": "unfavorable", "12": "unfavorable",
       "13": "unfavorable", "16": "unfavorable", "18": "unfavorable",
       "23": "unfavorable"}
ARM = {"01": "Chemo", "09": "ICI", "12": "ICI", "13": "Chemo",
       "16": "Chemo", "18": "ICI", "23": "ICI"}
MIN_CELLS = 30
N_CF = 400
SEED = 20261108


def quantile(v, p):
    n = len(v)
    if n == 0:
        return None
    v = sorted(v)
    h = (n - 1) * p
    lo = int(h)
    hi = min(lo + 1, n - 1)
    return v[lo] + (h - lo) * (v[hi] - v[lo])


def median(v):
    return quantile(v, 0.5)


def main():
    rng = random.Random(SEED)

    # 主导程序：直接读 13 号脚本的产物
    dom = {}
    with open(os.path.join(TAB, "nc2025_v4_clone_summary.csv")) as f:
        for r in csv.DictReader(f):
            if r["is_dominant"] == "yes":
                dom[r["patient"]] = r["mp"]
    print("主导程序:", dom)

    # ---- 上皮细胞 -> 克隆 ----
    comp = {}
    with open(os.path.join(BASE, "GEX_OBS_Lineage.csv")) as f:
        r = csv.reader(f)
        next(r)
        for row in r:
            if len(row) >= 5:
                comp[row[0]] = row[3]

    cell2 = {}
    with open(os.path.join(BASE, "gex_cnv_SCEVAN_karyotype.csv")) as f:
        r = csv.reader(f)
        next(r)
        for row in r:
            if len(row) < 6 or row[1] != "tumor" or not row[3]:
                continue
            p = row[4].split("_")
            if len(p) < 4:
                continue
            pid, tp = p[1].zfill(2), p[2]
            if pid not in PAIRED or tp not in ("T1", LAST[pid]):
                continue
            key = f"{row[0]}_{row[4]}_{row[5]}"
            if comp.get(key, "") != "Epithelial":
                continue
            cell2[key] = (pid, tp, row[3])

    # ---- 评分 ----
    cells = defaultdict(list)   # (pid, tp, clone) -> [score]
    mp_target = {(p, dom[p]) for p in dom}
    with open(os.path.join(BASE, "signature_scores.csv")) as f:
        r = csv.reader(f)
        header = next(r)
        mp_order = header[1:]
        idx = {mp: i + 1 for i, mp in enumerate(mp_order)}
        for row in r:
            if not row or not row[0]:
                continue
            if row[0] not in cell2:
                continue
            pid, tp, cl = cell2[row[0]]
            if pid not in dom:
                continue        # P12 无亚克隆、无主导程序
            cells[(pid, tp, cl)].append(float(row[idx[dom[pid]]]))

    out = []
    for p in PAIRED:
        if p not in dom:
            out.append(dict(patient=p, mp="", verdict="not_testable",
                            obs_delta="", cf_delta="", cf_share="", note="不可分层"))
            continue
        mp = dom[p]
        tp0, tp1 = "T1", LAST[p]
        base = {cl: v for (q, t, cl), v in cells.items() if q == p and t == tp0}
        post = {cl: v for (q, t, cl), v in cells.items() if q == p and t == tp1}
        base = {c: v for c, v in base.items() if len(v) >= MIN_CELLS}

        base_all = [x for c in sorted(base) for x in base[c]]
        post_all = [x for c in sorted(post) for x in post[c] if len(post[c]) >= MIN_CELLS]
        if len(base_all) < MIN_CELLS or len(post_all) < MIN_CELLS:
            out.append(dict(patient=p, mp=mp, verdict="not_testable",
                            obs_delta="", cf_delta="", cf_share="", note="细胞不足"))
            continue

        obs_delta = median(post_all) - median(base_all)

        # 反事实：按随访克隆占比，从各克隆基线分布重采样
        plan = [(c, len(post[c])) for c in sorted(base) if len(post.get(c, [])) >= MIN_CELLS]
        new_cl = [c for c in post if c not in base and len(post[c]) >= MIN_CELLS]
        n_new = sum(len(post[c]) for c in new_cl)
        cf_meds = []
        for _ in range(N_CF):
            samp = []
            for c, n in plan:
                pool = base[c]
                samp.extend(pool[rng.randrange(len(pool))] for _ in range(n))
            if n_new:   # 新发克隆没有基线分布，用基线全体代替
                samp.extend(base_all[rng.randrange(len(base_all))] for _ in range(n_new))
            cf_meds.append(median(samp))
        cf_delta = median(cf_meds) - median(base_all)

        share = cf_delta / obs_delta if abs(obs_delta) > 1e-9 else float("nan")

        if share != share:
            verdict = "not_testable"
        elif share < 0.15:
            verdict = "composition_minor"       # 克隆组成几乎不贡献
        elif share < 0.5:
            verdict = "composition_partial"     # 两种机制共存，转变为主
        elif share < 0.85:
            verdict = "composition_major"       # 共存，选择为主
        else:
            verdict = "composition_dominant"    # 可完全由克隆组成解释

        out.append(dict(patient=p, mp=mp, rcb=RCB[p], arm=ARM[p], verdict=verdict,
                        n_base=len(base_all), n_post=len(post_all),
                        n_clone_base=len(base), n_clone_post=len(plan),
                        med_base=f"{median(base_all):.4f}",
                        med_post=f"{median(post_all):.4f}",
                        obs_delta=f"{obs_delta:.4f}",
                        cf_delta=f"{cf_delta:.4f}",
                        cf_share=f"{share:.3f}",
                        cf_lo=f"{quantile(cf_meds, .025) - median(base_all):.4f}",
                        cf_hi=f"{quantile(cf_meds, .975) - median(base_all):.4f}",
                        note=""))

    os.makedirs(TAB, exist_ok=True)
    cols = []
    for r in out:
        for k in r:
            if k not in cols:
                cols.append(k)
    with open(os.path.join(TAB, "nc2025_v4_composition_counterfactual.csv"),
              "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in out:
            w.writerow({c: r.get(c, "") for c in cols})
    print("-> nc2025_v4_composition_counterfactual.csv")

    print("\n患者 程序  实测位移  反事实位移  组成解释率  判定")
    for r in out:
        if r["obs_delta"] == "":
            print(f"P{r['patient']}  {r['mp']:<5} {r['note']}")
            continue
        print(f"P{r['patient']}  {r['mp']:<5} {float(r['obs_delta']):+9.3f}  "
              f"{float(r['cf_delta']):+10.3f}  {float(r['cf_share']):>9.2f}  {r['verdict']}")


if __name__ == "__main__":
    main()
