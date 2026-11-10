# 38：置换检验——上端展宽是否超出"无位移"的噪声
#
# 目的：Fig 5 正文写了 "every decile moves ... the upper deciles move furthest"。
#       这句话是分位数之间的差，可能被分位数估计噪声污染。本脚本给零分布。
#
# 设计：
#   对每个 patient x programme 组合，把 baseline 与 on-treatment 细胞合并，
#   随机重新分成原来的两组大小（即打乱时点标签），完整重算 q10/q25/q50/q75/q90、
#   iqr（用新的 baseline 组）与 z。这样零分布包含了整个定量流程的噪声。
#
#   统计量 asym = z90 + z10 - 2*z50
#     刚性平移 -> 0；上端相对展宽 -> >0；上端被削 -> <0。
#     该量与 z50 的平移无关（z50 被减两次），是纯形状量。
#
# 同时输出：零假设下 |z50|>0.5 的比例。
#   若接近 0，则"按 |z50|>0.5 筛选后再算 spread 与 |z50| 的相关"在零假设下无法定义，
#   该检验设计不成立——这本身是必须报告的结论。

suppressPackageStartupMessages(library(data.table))

tab <- "results/tables"
dn <- fread(file.path(tab, "nc2025_v3_density.csv"))
dn[, patient := sprintf("P%02d", as.integer(patient))]
dn[, score := as.numeric(score)]

QS <- c(.10, .25, .50, .75, .90)
QN <- c("q10", "q25", "q50", "q75", "q90")

# ---------- 实测剖面（与 R/28 同口径：type = 7，以基线 IQR 为单位） ----------
qb <- dn[phase == "baseline",
         setNames(as.list(quantile(score, QS, type = 7)), QN),
         by = .(patient, mp)]
qa <- dn[phase == "on-treatment",
         setNames(as.list(quantile(score, QS, type = 7)), QN),
         by = .(patient, mp)]
setnames(qb, QN, paste0("b_", QN))
setnames(qa, QN, paste0("a_", QN))
P <- merge(qb, qa, by = c("patient", "mp"))
P[, iqr := b_q75 - b_q25]
for (k in QN) P[, (paste0("z_", k)) := (get(paste0("a_", k)) - get(paste0("b_", k))) / iqr]

P[, asym   := z_q90 + z_q10 - 2 * z_q50]
P[, spread := z_q90 - z_q10]

# ---------- 与既有形状表对账（确认同一口径） ----------
sh <- fread(file.path(tab, "nc2025_v3_shape.csv"))
sh[, patient := sprintf("P%02d", as.integer(patient))]
chk <- merge(P[, .(patient, mp, z50_here = z_q50)],
             sh[, .(patient, mp, z50_ref = as.numeric(z50))], by = c("patient", "mp"))
cat("== 对账：z50 here vs shape table ==\n")
cat("  n =", nrow(chk), " max abs diff =",
    sprintf("%.4f", max(abs(chk$z50_here - chk$z50_ref))), "\n\n")

# ---------- 置换 ----------
NPERM <- 1000
set.seed(20261109)

n_comb <- nrow(P)
res <- vector("list", n_comb)
t0 <- Sys.time()

for (i in seq_len(n_comb)) {
  pt <- P$patient[i]; pg <- P$mp[i]
  b <- dn[patient == pt & mp == pg & phase == "baseline", score]
  a <- dn[patient == pt & mp == pg & phase == "on-treatment", score]
  nb <- length(b); allv <- c(b, a); tot <- length(allv)

  zz <- matrix(NA_real_, NPERM, 3)      # 列：z10, z50, z90
  for (k in seq_len(NPERM)) {
    idx <- sample.int(tot, nb)
    qbb <- quantile(allv[idx],  QS, type = 7)
    qaa <- quantile(allv[-idx], QS, type = 7)
    z <- (qaa - qbb) / (qbb[4] - qbb[2])
    zz[k, ] <- c(z[1], z[3], z[5])
  }
  res[[i]] <- data.table(
    patient = pt, mp = pg,
    null_z10 = zz[, 1], null_z50 = zz[, 2], null_z90 = zz[, 3],
    null_asym = zz[, 3] + zz[, 1] - 2 * zz[, 2],
    null_spread = zz[, 3] - zz[, 1]
  )
  if (i %% 10 == 0) cat("  ", i, "/", n_comb, "\n")
}
null <- rbindlist(res)
cat("置换完成，用时", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "分钟\n\n")

# ---------- 零假设下能筛出多少个位移组合 ----------
cat("== 零假设下 |z50| > 0.5 的比例 ==\n")
cat("  逐组合中位：",
    sprintf("%.4f", median(null[, mean(abs(null_z50) > 0.5), by = .(patient, mp)]$V1)), "\n")
cat("  全体行占比：", sprintf("%.4f", mean(abs(null$null_z50) > 0.5)), "\n")
cat("  零分布 z50 的 SD（逐组合中位）：",
    sprintf("%.4f", median(null[, sd(null_z50), by = .(patient, mp)]$V1)), "\n\n")

# ---------- 逐组合 p 值 ----------
summ <- null[, .(null_asym_q50 = median(null_asym),
                 null_asym_q95 = quantile(null_asym, .95),
                 null_asym_q05 = quantile(null_asym, .05),
                 null_spread_q50 = median(null_spread),
                 null_z50_sd = sd(null_z50)),
             by = .(patient, mp)]
S <- merge(P[, .(patient, mp, obs_z50 = z_q50, obs_asym = asym, obs_spread = spread)],
           summ, by = c("patient", "mp"))

S[, p_asym := {
  nd <- null[patient == .BY$patient & mp == .BY$mp, null_asym]
  (1 + sum(nd >= obs_asym)) / (1 + length(nd))
}, by = .(patient, mp)]
S[, direction := fifelse(obs_z50 > 0.5, "up", fifelse(obs_z50 < -0.5, "down", "flat"))]

fwrite(S, file.path(tab, "nc2025_permutation_asym.csv"))

# ---------- 汇总：上移组合 ----------
up <- S[direction == "up"]
cat("== 上移组合（obs z50 > 0.5），n =", nrow(up), "==\n")
print(up[order(-obs_asym), .(patient, mp, obs_z50 = round(obs_z50, 2),
                             obs_asym = round(obs_asym, 2),
                             null_q95 = round(null_asym_q95, 2),
                             p = round(p_asym, 3))])

cat("\n  实测 asym > 零分布 95% 分位的组合数：",
    sum(up$obs_asym > up$null_asym_q95), "/", nrow(up), "\n")
cat("  实测 asym 中位：", sprintf("%.3f", median(up$obs_asym)), "\n")
cat("  零分布 asym 中位（各组合取中位后再中位）：",
    sprintf("%.3f", median(up$null_asym_q50)), "\n")

# Fisher 合并（注意：同患者多程序，p 非独立，仅作参考）
chi <- -2 * sum(log(up$p_asym))
cat("  Fisher 合并（伪独立，仅供参考）：chi2 =", sprintf("%.2f", chi),
    " df =", 2 * nrow(up), "\n")

cat("\n== 全部组合概览 ==\n")
print(S[order(-obs_asym), .(patient, mp, dir = direction,
                            z50 = round(obs_z50, 2), asym = round(obs_asym, 2),
                            p = round(p_asym, 3))])
cat("\n-> nc2025_permutation_asym.csv\n")
