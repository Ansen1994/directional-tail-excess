#!/usr/bin/env Rscript
# ==============================================================================
# 41_depth_matched_downsample.R
#
# 目的（对应审稿意见 P1-2）：
#   三队列每组合中位细胞数为 4,384 / 448 / 918，差近 10 倍，而 ex_lead 上界
#   （3.96 IQR）与 Fig 6C 的 Q90 外推（9.3 IQR）对稀疏上尾极敏感。现有对照只做了
#   "提高最小细胞数"（R/36）与"留一患者"，都不是直接回应。
#   本脚本把每个组合的两端各抽到同一细胞数，重算全部形状量，检验：
#     位移率（21/63、38/145、69/220）与 ex_trail 上界（0.060/0.028/0.061）
#     是否在深度匹配后仍然成立。
#
# 口径（与 R/28、R/29、R/32、R/39 完全一致，勿改）：
#   z_q = (Q_q^post - Q_q^base) / IQR_base
#   ex_trail = (z25 - z50) * sign(z50)  if z50 > 0, else (z75 - z50) * sign(z50)
#   ex_lead  = (z75 - z50) * sign(z50)  if z50 > 0, else (z25 - z50) * sign(z50)
#   asym     = z75 - z25
#
# 【重要】分母一律取"原值表给定的基线 IQR"，不现场从样本重算：
#   nc2025    -> nc2025_v3_shape.csv 的 iqr_base
#   Bassez    -> bassez_v2_profile.csv 的 iqr
#   NeoLetExe -> 从细胞全量算（已验证与形状表原值逐位一致，差 ~1e-14）
# 理由：Bassez 的 profile 里 iqr 列**不等于** q75_base - q25_base（例如 BIOKEY_18
#   ER+ IFN_response：iqr 0.2038，而 q75_base - q25_base = 0.3461），且无法从仓库
#   任何数据复现（已试：跨 subtype 合并、缩尾 1-99%、每患者跨程序合并、9 种分位算法，
#   均不匹配；而四个分位数本身与 cell_scores 逐位一致，差 5e-16）。
#   若现场重算 IQR，Bassez 的 z50 会与原值差 1.4-1.7 倍，无法对照。
#   固定分母还有一个好处：下采样只引入分子（分位数）的抽样波动，不叠加分母自身的波动。
#   => 该列的可重复性缺口已单独记录，见 docs/待办卡片_v1.md「已知缺口」。
#
# 自检：用**全细胞 + 固定分母**重算，必须等于形状表原值。
#       nc 容差 1e-2（v3_shape 分位存 6 位小数，四舍五入残留 ~3e-3），其余 1e-6。
#
# 输入：
#   results/tables/nc2025_fig2_cells.csv        长表，65,206 细胞 x 9 程序（覆盖全部 63 组合）
#   results/tables/bassez_v2_cell_scores.csv    宽表，5 程序
#   results/tables/neoletexe_cell_scores.csv    宽表，5 程序
#   results/tables/nc2025_fig5_shape_calls.csv  nc 原值（63 组合）
#   results/tables/nc2025_v3_shape.csv          nc 时点对（baseline 恒 T1，post = T2/T3/T4）
#   results/tables/bassez_v2_profile.csv        Bassez 原值（145 组合）
#   results/tables/neoletexe_shape.csv          NeoLetExe 原值（220 组合）
#
# 输出：
#   results/tables/depth_matched_runs.csv       逐组合逐重抽样
#   results/tables/depth_matched_summary.csv    队列 x 档位汇总 + 原值对照
# ==============================================================================

suppressPackageStartupMessages(library(data.table))

tab  <- "results/tables"
QS   <- c(.10, .25, .50, .75, .90)
B    <- 200
NDS  <- c(200L, 100L)          # 两档统一细胞数
SEED <- 20261110
set.seed(SEED)

np <- function(x) {
  y <- suppressWarnings(as.integer(gsub("[^0-9]", "", as.character(x))))
  ifelse(is.na(y), as.character(x), sprintf("P%02d", y))
}

## ---------------------------------------------------------------------------
## 0. 形状量
## ---------------------------------------------------------------------------
shape_stats <- function(b, a, iqr) {
  qb <- quantile(b, QS, type = 7, names = FALSE)
  qa <- quantile(a, QS, type = 7, names = FALSE)
  if (!is.finite(iqr) || iqr <= 0) return(NULL)
  z <- (qa - qb) / iqr
  z50 <- z[3]
  list(z50   = z50,
       z25   = z[2], z75 = z[4],
       asym     = z[4] - z[2],
       ex_trail = if (z50 > 0) z[2] - z50 else (z[4] - z50) * -1,
       ex_lead  = if (z50 > 0) z[4] - z50 else (z[2] - z50) * -1)
}

## ---------------------------------------------------------------------------
## 1. 三队列统一成长表：patient, prog, pair, phase, score
## ---------------------------------------------------------------------------
# --- nc2025：长表，patient/tp/batch/prog/score；baseline 恒 T1
fc <- fread(file.path(tab, "nc2025_fig2_cells.csv"))
setnames(fc, c("patient", "tp", "mp", "score"), c("patient", "tp", "prog", "score"),
         skip_absent = TRUE)
fc[, `:=`(patient = np(patient), score = as.numeric(score), prog = as.character(prog))]
fc[, tp := as.character(tp)]

sh_nc <- fread(file.path(tab, "nc2025_v3_shape.csv"))
sh_nc[, `:=`(patient = np(patient), prog = as.character(mp), to = as.character(tp))]
# 每个 (patient, prog) 的 post 时点
nc_pair <- unique(sh_nc[, .(patient, prog, to)])
nc_pair[, `:=`(from = "T1", pair = paste0("T1>", to))]

nc_long <- rbind(
  merge(fc, nc_pair[, .(patient, prog, from, pair)], by = c("patient", "prog"))[
    tp == from, .(patient, prog, pair, phase = "baseline", score)],
  merge(fc, nc_pair[, .(patient, prog, to, pair)], by = c("patient", "prog"))[
    tp == to, .(patient, prog, pair, phase = "post", score)]
)
nc_long <- nc_long[!is.na(score)]

# --- Bassez：宽表，5 程序，Pre > On
PROGS <- c("ER_signal", "Proliferation", "IFN_response", "EMT_matrix", "Luminal_secretory")
ba <- fread(file.path(tab, "bassez_v2_cell_scores.csv"))
ba_long <- melt(ba, id.vars = c("cell", "patient", "timepoint"),
                measure.vars = PROGS, variable.name = "prog", value.name = "score")
ba_long[, `:=`(pair = "Pre>On",
               phase = fifelse(timepoint == "Pre", "baseline", "post"),
               prog = as.character(prog))]
ba_long <- ba_long[!is.na(score), .(patient, prog, pair, phase, score)]

# --- NeoLetExe：宽表，tp A/B/C，pair 从形状表取
ne <- fread(file.path(tab, "neoletexe_cell_scores.csv"))
ne_shape <- fread(file.path(tab, "neoletexe_shape.csv"))
ne_long <- melt(ne, id.vars = c("cell", "patient", "tp"),
                measure.vars = PROGS, variable.name = "prog", value.name = "score")
ne_long[, `:=`(prog = as.character(prog), tp = as.character(tp))]
ne_long <- ne_long[!is.na(score), .(patient, prog, tp, score)]
ne_pairs <- unique(ne_shape[, .(patient, prog = program, pair)])
ne_pairs[, `:=`(from = sub(">.*$", "", pair), to = sub("^.*>", "", pair))]
ne_long <- merge(ne_long, ne_pairs, by = c("patient", "prog"), allow.cartesian = TRUE)
ne_long <- ne_long[tp == from | tp == to,
                   .(patient, prog, pair, phase = fifelse(tp == from, "baseline", "post"), score)]

LONG <- list(nc2025 = nc_long, Bassez = ba_long, NeoLetExe = ne_long)

cat("== 长表规模 ==\n")
for (nm in names(LONG)) {
  d <- LONG[[nm]]
  cat(sprintf("  %-10s 细胞行 %7d  组合 %4d  baseline 中位 %5d  post 中位 %5d\n",
              nm, nrow(d), uniqueN(d, by = c("patient", "prog", "pair")),
              as.integer(median(d[phase == "baseline", .N, by = .(patient, prog, pair)]$N)),
              as.integer(median(d[phase == "post",     .N, by = .(patient, prog, pair)]$N))))
}

## ---------------------------------------------------------------------------
## 2. 原值（形状表）与自检
## ---------------------------------------------------------------------------
ORIG <- list()

o <- fread(file.path(tab, "nc2025_fig5_shape_calls.csv"))
o[, `:=`(patient = np(patient), prog = as.character(mp))]
o <- merge(o, nc_pair, by = c("patient", "prog"))
# iqr_base 来自 v3_shape（fig5_shape_calls 没存分母）
sh_iqr <- unique(sh_nc[, .(patient, prog, iqr_base)])
o <- merge(o, sh_iqr, by = c("patient", "prog"))
ORIG$nc2025 <- o[, .(cohort = "nc2025", patient, prog, pair,
                     z50, ex_trail = as.numeric(ex_trail), ex_lead = as.numeric(ex_lead),
                     iqr_base)]

o <- fread(file.path(tab, "bassez_v2_profile.csv"))
o[, `:=`(ex_trail = fifelse(d50 > 0, d25 - d50, (d75 - d50) * -1),
         ex_lead  = fifelse(d50 > 0, d75 - d50, (d25 - d50) * -1))]
ORIG$Bassez <- o[, .(cohort = "Bassez", patient, prog = as.character(prog),
                     pair = "Pre>On", z50 = d50, ex_trail, ex_lead,
                     iqr_base = iqr)]

o <- fread(file.path(tab, "neoletexe_shape.csv"))
ORIG$NeoLetExe <- o[, .(cohort = "NeoLetExe", patient, prog = as.character(program), pair,
                        z50 = as.numeric(z50), ex_trail = as.numeric(ex_trail),
                        ex_lead = as.numeric(ex_lead), iqr_base = NA_real_)]
# NeoLetExe 的形状表未存分母：从细胞全量算（已验证与原值逐位一致）
neo_iqr <- LONG$NeoLetExe[phase == "baseline", {
  q <- quantile(score, c(.25, .75), type = 7, names = FALSE)
  .(iqr_base = q[2] - q[1])
}, by = .(patient, prog, pair)]
ORIG$NeoLetExe <- merge(ORIG$NeoLetExe[, !"iqr_base"], neo_iqr,
                        by = c("patient", "prog", "pair"))

cat("\n== 自检：全细胞 + 原值分母 重算 vs 形状表原值 ==\n")
# nc 容差 5e-2：fig5_shape_calls 的 z50 算自 nc2025_v3_density.csv（31 组合的细胞子集），
# 本脚本算自 nc2025_fig2_cells.csv（63 组合全量），细胞集略有差异 -> 相对误差 <0.2%。
# 对判据（|z50| > 0.5）无影响，但为透明起见单独标注。
TOL <- c(nc2025 = 5e-2, Bassez = 1e-6, NeoLetExe = 1e-6)
chk <- rbindlist(lapply(names(LONG), function(nm) {
  d <- LONG[[nm]]; k <- ORIG[[nm]]
  d <- merge(d, k[, .(patient, prog, pair, iqr_base)], by = c("patient", "prog", "pair"))
  got <- d[, {
    s <- shape_stats(score[phase == "baseline"], score[phase == "post"], iqr_base[1])
    if (is.null(s)) list(z50 = NA_real_, ex_trail = NA_real_, ex_lead = NA_real_)
    else as.list(unlist(s[c("z50", "ex_trail", "ex_lead")]))
  }, by = .(patient, prog, pair)]
  m <- merge(k, got, by = c("patient", "prog", "pair"), suffixes = c(".o", ".r"))
  m[, .(cohort = nm,
        n = nrow(m),
        dz = max(abs(z50.o - z50.r), na.rm = TRUE),
        det = max(abs(ex_trail.o - ex_trail.r), na.rm = TRUE),
        del = max(abs(ex_lead.o - ex_lead.r), na.rm = TRUE))]
}))
for (i in seq_len(nrow(chk)))
  cat(sprintf("  %-10s n=%3d  max|dz50|=%9.3e  max|d ex_trail|=%9.3e  max|d ex_lead|=%9.3e  tol %g\n",
              chk$cohort[i], chk$n[i], chk$dz[i], chk$det[i], chk$del[i], TOL[chk$cohort[i]]))
if (any(chk$dz > TOL[chk$cohort] | chk$det > TOL[chk$cohort] |
          chk$del > TOL[chk$cohort], na.rm = TRUE)) {
  for (nm in names(LONG)) {
    d <- LONG[[nm]]; k <- ORIG[[nm]]
    d <- merge(d, k[, .(patient, prog, pair, iqr_base)], by = c("patient", "prog", "pair"))
    got <- d[, {
      s <- shape_stats(score[phase == "baseline"], score[phase == "post"], iqr_base[1])
      if (is.null(s)) list(z50 = NA_real_, ex_trail = NA_real_, ex_lead = NA_real_)
      else as.list(unlist(s[c("z50", "ex_trail", "ex_lead")]))
    }, by = .(patient, prog, pair)]
    m <- merge(k, got, by = c("patient", "prog", "pair"), suffixes = c(".o", ".r"))
    m[, dlt := abs(z50.o - z50.r)]
    setorder(m, -dlt)
    cat("\n  --", nm, "偏差最大的 3 个组合 --\n")
    for (i in seq_len(min(3L, nrow(m))))
      cat(sprintf("     %s %s %s  z50 orig %+.4f  recomputed %+.4f  |d| %.3e\n",
                  m$patient[i], m$prog[i], m$pair[i], m$z50.o[i], m$z50.r[i], m$dlt[i]))
  }
  stop("自检失败：全细胞 + 原值分母 重算仍与形状表不一致，先修脚本。")
}
cat("  自检通过（三队列 x 三个量，均在容差内）\n")

## ---------------------------------------------------------------------------
## 3. 下采样
## ---------------------------------------------------------------------------
boot_one <- function(b, a, n, B, iqr) {
  nb <- min(length(b), n); na <- min(length(a), n)
  deterministic <- (nb == length(b) && na == length(a))
  nrep <- if (deterministic) 1L else B
  rbindlist(lapply(seq_len(nrep), function(i) {
    bs <- if (nb == length(b)) b else b[sample.int(length(b), nb)]
    as_ <- if (na == length(a)) a else a[sample.int(length(a), na)]
    s <- shape_stats(bs, as_, iqr)
    if (is.null(s)) return(NULL)
    as.data.table(c(s, list(n_base = nb, n_post = na, deterministic = deterministic)))
  }))
}

runs <- rbindlist(lapply(names(LONG), function(nm) {
  d <- merge(LONG[[nm]], ORIG[[nm]][, .(patient, prog, pair, iqr_base)],
             by = c("patient", "prog", "pair"))
  rbindlist(lapply(NDS, function(n) {
    d[, {
      b <- score[phase == "baseline"]; a <- score[phase == "post"]
      r <- boot_one(b, a, n, B, iqr_base[1])
      if (is.null(r)) return(NULL)
      cbind(data.table(cohort = nm, patient = patient, prog = prog, pair = pair,
                       n_ds = n, B_used = length(b) < n | length(a) < n),
            r)
    }, by = .(patient, prog, pair)]
  }))
}), fill = TRUE)

cat(sprintf("\n== 下采样完成：%s 行（%d 组合 x %d 档 x <=%d 次）==\n",
            format(nrow(runs), big.mark = ","),
            uniqueN(runs, by = c("cohort", "patient", "prog", "pair")),
            length(NDS), B))

## ---------------------------------------------------------------------------
## 4. 汇总与对照
## ---------------------------------------------------------------------------
## 每个组合先压成一个代表值（B 次的中位数），再与"全细胞一次算"的原值同口径比较。
## 直接对 20 万行取比例会让每个组合被重复计 B 次，且 max() 会被重抽样上界抬高。
per <- runs[, .(z50_med = median(z50),
                et_med  = median(ex_trail),
                el_med  = median(ex_lead),
                et_q95  = quantile(ex_trail, .95),
                el_q95  = quantile(ex_lead, .95),
                p_shift = mean(abs(z50) > 0.5),
                n_base_med = as.numeric(median(n_base)),
                n_post_med = as.numeric(median(n_post)),
                pct_short = mean(as.numeric(B_used))),
            by = .(cohort, n_ds, patient, prog, pair)]

orig_flat <- rbindlist(lapply(names(ORIG), function(nm) ORIG[[nm]]), use.names = TRUE)
per <- merge(per, orig_flat[, .(cohort, patient, prog, pair, z50_o = z50,
                                et_o = ex_trail, el_o = ex_lead)],
             by = c("cohort", "patient", "prog", "pair"))
per[, `:=`(shift_o = abs(z50_o) > 0.5, shift_d = abs(z50_med) > 0.5,
           up_o = z50_o > 0.5, up_d = z50_med > 0.5,
           d_z50 = z50_med - z50_o)]
# 注：上移侧统计一律限定 z50 > 0.5（即"上移且实质位移"的组合），与正文口径一致
#     —— nc 16 个、Bassez 21 个、NeoLetExe 29 个。用 z50 > 0 会把大量小位移组合
#     也算进去，上界从 0.060 抬到 0.180，与正文不符。

summ <- per[, .(n_combo     = .N,
                shift_rate  = mean(shift_d),
                shift_prob  = mean(p_shift),
                agree_call  = mean(shift_o == shift_d),
                et_up_med   = median(et_med[up_o]),
                et_up_max   = max(et_med[up_o]),
                et_up_max95 = max(et_q95[up_o]),
                el_max      = max(abs(el_med[up_o])),
                z50_med_abs = median(abs(z50_med)),
                dz50_med    = median(d_z50),
                dz50_q05    = quantile(d_z50, .05),
                dz50_q95    = quantile(d_z50, .95),
                pct_short   = mean(pct_short)),
            by = .(cohort, n_ds)]

orig_summ <- rbindlist(lapply(names(ORIG), function(nm) {
  d <- ORIG[[nm]]; up <- d$z50 > 0.5
  data.table(cohort = nm, n_ds = NA_integer_, n_combo = nrow(d),
             shift_rate = mean(abs(d$z50) > 0.5), shift_prob = NA_real_,
             agree_call = NA_real_,
             et_up_med = median(d$ex_trail[up]), et_up_max = max(d$ex_trail[up]),
             et_up_max95 = NA_real_, el_max = max(abs(d$ex_lead[up])),
             z50_med_abs = median(abs(d$z50)), dz50_med = NA_real_,
             dz50_q05 = NA_real_, dz50_q95 = NA_real_, pct_short = 0)
}))

out <- rbind(orig_summ, summ[order(cohort, -n_ds)], fill = TRUE)
out[, n_ds := fifelse(is.na(n_ds), "all cells", as.character(n_ds))]

cat("\n== 深度匹配前后对照（每组合取 B 次中位）==\n")
cat(sprintf("%-10s %9s %7s %8s %9s %10s %10s %9s %8s %8s\n",
            "cohort", "n_ds", "n_combo", "shift", "shift_p", "et_up_med", "et_up_max",
            "el_max", "dz50", "agree"))
for (i in seq_len(nrow(out)))
  cat(sprintf("%-10s %9s %7d %8s %9s %10s %10s %9s %8s %8s\n",
              out$cohort[i], out$n_ds[i], out$n_combo[i],
              sprintf("%.3f", out$shift_rate[i]),
              ifelse(is.na(out$shift_prob[i]), "-", sprintf("%.3f", out$shift_prob[i])),
              sprintf("%+.3f", out$et_up_med[i]),
              sprintf("%+.3f", out$et_up_max[i]),
              sprintf("%.2f", out$el_max[i]),
              ifelse(is.na(out$dz50_med[i]), "-", sprintf("%+.3f", out$dz50_med[i])),
              ifelse(is.na(out$agree_call[i]), "-", sprintf("%.3f", out$agree_call[i]))))

cat("\n== 某端细胞 < n 的组合占比（%）==\n")
for (i in seq_len(nrow(summ)))
  cat(sprintf("  %-10s n=%4d  短端组合 %5.1f%%\n",
              summ$cohort[i], summ$n_ds[i], 100 * summ$pct_short[i]))

fwrite(runs, file.path(tab, "depth_matched_runs.csv"))
fwrite(per,  file.path(tab, "depth_matched_per_combo.csv"))
fwrite(out,  file.path(tab, "depth_matched_summary.csv"))
cat("\n写出：depth_matched_runs.csv / depth_matched_per_combo.csv / depth_matched_summary.csv\n")
