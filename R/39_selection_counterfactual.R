# 39 — 选择反事实的逐组合配对检验 + 后缘超额阈值扫描
#
# 目的：
#   ① 原 Fig 5D 里 "selection would land here" 的阴影矩形 x 起点 0.55 是按 Q10 直觉画的，
#      与 Methods / ex_trail 实际用的 Q25 口径不符。该矩形已按方案 3 删除；R/28 现直接
#      读取本脚本的 paired 表，画"16 个反事实的中位竖线 + 10th-90th percentile 带"。
#      本脚本同时对每个上移组合用自己的基线细胞重算 top-20% 反事实，做配对比较。
#   ② 扫描 z50 阈值 0 / 0.3 / 0.5 / 0.7 / 0.9，检验"上移侧后缘超额上界"
#      是否依赖 0.5 这个约定阈值（S2 敏感性原先只扫了位移率，没扫上界）。
#
# 口径（与 Methods / R/28 / R/29 / R/32 一致，勿改）：
#   z_q  = (Q_q^post - Q_q^pre) / IQR_pre
#   ex_trail = (z25 - z50) * sign(z50)  [上移] / (z75 - z50) * sign(z50) [下移]
#   Bassez 的 profile 表没有 ex_trail 列，用 d25/d50/d75 按上式现算（d* 已是 IQR 标准化位移）
#   选择反事实：只保留基线端 top 20% 细胞，重算 z25、z50，取 cf_ex_trail = z25 - z50
#
# 输出（results/tables/）：
#   trailing_excess_threshold_scan.csv     阈值 × 队列的上界汇总
#   selection_counterfactual_paired.csv    逐组合 obs / cf 配对值
#   selection_counterfactual_patients.csv  患者层面聚合（推断用）
#
# 注意：组合之间不独立（同一患者多个程序），组合层面的符号检验只作描述；
#       推断以患者层面的配对 Wilcoxon 为准。

suppressPackageStartupMessages(library(data.table))
tab <- "results/tables"

## ---------------- 1. 读三队列的形状表 ----------------
coh <- list(
  nc2025    = fread(file.path(tab, "nc2025_fig5_shape_calls.csv")),
  Bassez    = fread(file.path(tab, "bassez_v2_profile.csv")),
  NeoLetExe = fread(file.path(tab, "neoletexe_shape.csv"))
)
coh$Bassez[, `:=`(z50 = d50,
                  ex_trail = fifelse(d50 > 0, d25 - d50, (d75 - d50) * -1))]

stopifnot(all(vapply(coh, function(d) all(c("z50", "ex_trail") %in% names(d)), logical(1))))

## ---------------- 2. 阈值扫描 ----------------
scan <- rbindlist(lapply(c(0, 0.3, 0.5, 0.7, 0.9), function(th) {
  rbindlist(lapply(names(coh), function(nm) {
    u <- coh[[nm]][z50 > th]
    data.table(cohort = nm, threshold = th, n_up = nrow(u),
               trail_median = round(median(u$ex_trail), 3),
               trail_p90    = round(as.numeric(quantile(u$ex_trail, .9)), 3),
               trail_max    = round(max(u$ex_trail), 3),
               n_trail_gt_0.5 = sum(u$ex_trail > 0.5))
  }))
}))
fwrite(scan, file.path(tab, "trailing_excess_threshold_scan.csv"))

## ---------------- 3. 逐组合选择反事实 ----------------
# 反事实的 trailing 超额（上移组合）。z25 的**绝对**位移约 1.3 IQR，
# 但 z50 也一起上移，所以相对超额只有约 +0.24 IQR —— 这正是 Methods
# 旧版把"约 1 IQR"误安到 ex_trail 上的原因。
cf_one <- function(base, iqr) {
  sel <- base[base >= quantile(base, 0.80)]
  as.numeric((quantile(sel, .25) - quantile(base, .25)) / iqr -
             (quantile(sel, .50) - quantile(base, .50)) / iqr)
}

# 3a. nc2025：density 表 patient 是 01，形状表是 P01；phase = baseline / on-treatment
dn <- fread(file.path(tab, "nc2025_v3_density.csv"))
dn[, patient := sprintf("P%02d", as.integer(patient))][, score := as.numeric(score)]
r1 <- rbindlist(lapply(seq_len(nrow(coh$nc2025)), function(i) {
  d <- coh$nc2025[i]
  bb <- dn[patient == d$patient & mp == d$mp & phase == "baseline", score]
  if (length(bb) < 20) return(NULL)
  iqr <- as.numeric(quantile(bb, .75) - quantile(bb, .25)); if (iqr <= 0) return(NULL)
  data.table(cohort = "nc2025", unit = paste(d$patient, d$mp),
             z50 = d$z50, obs = d$ex_trail, cf = cf_one(bb, iqr))
}))

# 3b. Bassez：细胞级程序分数是宽表，Pre = 基线
bs <- fread(file.path(tab, "bassez_v2_cell_scores.csv"))
bp <- coh$Bassez
r2 <- rbindlist(lapply(seq_len(nrow(bp)), function(i) {
  prog <- bp$prog[i]; if (!prog %in% names(bs)) return(NULL)
  bb <- as.numeric(bs[patient == bp$patient[i] & timepoint == "Pre", get(prog)])
  bb <- bb[!is.na(bb)]; if (length(bb) < 20) return(NULL)
  iqr <- as.numeric(quantile(bb, .75) - quantile(bb, .25)); if (iqr <= 0) return(NULL)
  data.table(cohort = "Bassez", unit = paste(bp$patient[i], prog),
             z50 = bp$z50[i], obs = bp$ex_trail[i], cf = cf_one(bb, iqr))
}))

# 3c. NeoLetExe：tp = A/B/C，from = 该组合的基线端
ns <- fread(file.path(tab, "neoletexe_cell_scores.csv"))
nx <- coh$NeoLetExe
r3 <- rbindlist(lapply(seq_len(nrow(nx)), function(i) {
  prog <- nx$program[i]; if (!prog %in% names(ns)) return(NULL)
  bb <- as.numeric(ns[patient == nx$patient[i] & tp == nx$from[i], get(prog)])
  bb <- bb[!is.na(bb)]; if (length(bb) < 20) return(NULL)
  iqr <- as.numeric(quantile(bb, .75) - quantile(bb, .25)); if (iqr <= 0) return(NULL)
  data.table(cohort = "NeoLetExe",
             unit = paste(nx$patient[i], nx$pair[i], prog),
             z50 = nx$z50[i], obs = nx$ex_trail[i], cf = cf_one(bb, iqr))
}))

paired <- rbindlist(list(r1, r2, r3))
paired[, lower := obs < cf]
fwrite(paired, file.path(tab, "selection_counterfactual_paired.csv"))

## ---------------- 4. 患者层面聚合（推断用） ----------------
paired[, pat := sub(" .*", "", unit)]
pat <- paired[z50 > 0.5, .(n = .N, obs = median(obs), cf = median(cf)),
              by = .(cohort, pat)]
pat[, lower := obs < cf]
fwrite(pat, file.path(tab, "selection_counterfactual_patients.csv"))

## ---------------- 5. 汇总输出 ----------------
cat("\n=== 阈值扫描（上移侧后缘超额） ===\n"); print(scan)
cat("\n=== z50 > 0.5 的配对检验 ===\n")
for (nm in c("nc2025", "Bassez", "NeoLetExe")) {
  u <- paired[cohort == nm & z50 > 0.5]
  cat(sprintf("%-9s 组合 %2d  obs<cf %2d/%2d  Wilcoxon p=%.3g\n",
              nm, nrow(u), sum(u$lower), nrow(u),
              suppressWarnings(wilcox.test(u$cf, u$obs, paired = TRUE)$p.value)))
}
u <- paired[z50 > 0.5]
cat(sprintf("合并     组合 %2d  obs<cf %2d/%2d (%.1f%%)\n",
            nrow(u), sum(u$lower), nrow(u), 100 * mean(u$lower)))
cat(sprintf("合并     患者 %2d  obs<cf %2d/%2d (%.1f%%)  配对 Wilcoxon p=%.3g\n",
            nrow(pat), sum(pat$lower), nrow(pat), 100 * mean(pat$lower),
            wilcox.test(pat$cf, pat$obs, paired = TRUE)$p.value))
cat("\n-> 三张表已写入 results/tables/\n")
