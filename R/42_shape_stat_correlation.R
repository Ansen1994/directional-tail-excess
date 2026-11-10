#!/usr/bin/env Rscript
# ==============================================================================
# 42_shape_stat_correlation.R
#
# 目的（对应审稿意见 Moderate-1）：
#   正文把 z50、ex_trail、asym、A 四个形状量当四条并列证据使用（Results S4-S6），
#   但它们在数学上同源：
#     asym = z75 - z25
#     A    = z10 + z90 - 2*z50
#     ex_trail = (z25 - z50)*sign(z50)            (z50 > 0)
#              = (z75 - z50)*sign(z50)            (z50 < 0)
#   读者无法判断证据是否重叠。本脚本给出三队列各自的 4x4 Spearman 矩阵。
#
# 口径：
#   nc2025    分母 = nc2025_v3_shape.csv 的 iqr_base（原值表给定值）
#   Bassez    分母 = bassez_v2_profile.csv 的 iqr（原值表给定值；该列无法从仓库复现，
#             但四个统计量对分母是线性的，同一组合内同除一个常数，秩不变）
#   NeoLetExe 分母 = 基线 q75 - q25（从细胞全量算，已验证与形状表逐位一致）
#
# 【术语警告】本脚本里 `asym` 一律指 z75 - z25（与 Fig 5/6/7 一致），
#   不是 R/38 置换检验里的 z90 + z10 - 2*z50（那个量在本脚本叫 A）。
#
# 输出：
#   results/tables/shape_stat_correlation.csv   长表：cohort, stat1, stat2, n, rho, p
# ==============================================================================

suppressPackageStartupMessages(library(data.table))

tab <- "results/tables"
QS  <- c(.10, .25, .50, .75, .90)
STATS <- c("z50", "ex_trail", "asym", "A")

np <- function(x) {
  y <- suppressWarnings(as.integer(gsub("[^0-9]", "", as.character(x))))
  ifelse(is.na(y), as.character(x), sprintf("P%02d", y))
}

## ---------------------------------------------------------------------------
## 0. 四个统计量（输入：五个分位 z，均以基线 IQR 为单位）
## ---------------------------------------------------------------------------
four_stats <- function(z) {
  z10 <- z[1]; z25 <- z[2]; z50 <- z[3]; z75 <- z[4]; z90 <- z[5]
  list(z50      = z50,
       ex_trail = if (z50 > 0) z25 - z50 else (z75 - z50),
       asym     = z75 - z25,
       A        = z10 + z90 - 2 * z50)
}

## ---------------------------------------------------------------------------
## 1. 三队列逐组合重算
## ---------------------------------------------------------------------------
per <- list()

# --- nc2025：长表，baseline 恒 T1，post 因患者而异
fc <- fread(file.path(tab, "nc2025_fig2_cells.csv"))
setnames(fc, c("patient", "tp", "mp", "score"), c("patient", "tp", "prog", "score"),
         skip_absent = TRUE)
fc[, `:=`(patient = np(patient), prog = as.character(prog),
          tp = as.character(tp), score = as.numeric(score))]

sh_nc <- fread(file.path(tab, "nc2025_v3_shape.csv"))
sh_nc[, `:=`(patient = np(patient), prog = as.character(mp), to = as.character(tp))]
nc_meta <- unique(sh_nc[, .(patient, prog, to = to, iqr_base = iqr_base)])
nc_meta[, `:=`(from = "T1", pair = paste0("T1>", to))]

d <- merge(fc, nc_meta[, .(patient, prog, from, to, pair, iqr_base)],
           by = c("patient", "prog"), allow.cartesian = TRUE)
d <- d[!is.na(score) & (tp == from | tp == to)]
d[, phase := fifelse(tp == from, "baseline", "post")]

per$nc2025 <- d[, {
  b <- score[phase == "baseline"]; a <- score[phase == "post"]
  iq <- iqr_base[1]
  if (length(b) < 10 || length(a) < 10 || !is.finite(iq) || iq <= 0) NULL
  else {
    z <- (quantile(a, QS, type = 7, names = FALSE) -
          quantile(b, QS, type = 7, names = FALSE)) / iq
    as.data.table(four_stats(z))
  }
}, by = .(patient, prog, pair)]

# --- Bassez：分位数与分母都取自 profile（与原值一致）
pf <- fread(file.path(tab, "bassez_v2_profile.csv"))
per$Bassez <- pf[, .(patient = as.character(patient), prog = as.character(prog),
                     pair = "Pre>On",
                     z50 = d50,
                     ex_trail = fifelse(d50 > 0, d25 - d50, d75 - d50),
                     asym = d75 - d25,
                     A = d10 + d90 - 2 * d50)]

# --- NeoLetExe：全部取自形状表
ns <- fread(file.path(tab, "neoletexe_shape.csv"))
per$NeoLetExe <- ns[, .(patient = as.character(patient), prog = as.character(program),
                        pair = pair,
                        z50 = as.numeric(z50),
                        ex_trail = as.numeric(ex_trail),
                        asym = as.numeric(asym),
                        A = as.numeric(z10) + as.numeric(z90) - 2 * as.numeric(z50))]

cat("== 组合数 ==\n")
for (nm in names(per)) {
  x <- per[[nm]]
  cat(sprintf("  %-10s n=%3d  患者 %2d  ex_trail 缺失 %d\n",
              nm, nrow(x), uniqueN(x$patient), sum(!is.finite(x$ex_trail))))
}

## ---------------------------------------------------------------------------
## 2. 4x4 Spearman（队列内）
## ---------------------------------------------------------------------------
rho_mat <- function(x) {
  m <- as.matrix(x[, ..STATS])
  cor(m, method = "spearman", use = "complete.obs")
}

out <- list()
for (nm in names(per)) {
  x <- per[[nm]]
  r <- rho_mat(x)
  cat(sprintf("\n== %s（n=%d）==\n", nm, nrow(x)))
  print(round(r, 3))
  for (i in 1:(length(STATS) - 1)) for (j in (i + 1):length(STATS)) {
    s1 <- STATS[i]; s2 <- STATS[j]
    ct <- suppressWarnings(cor.test(x[[s1]], x[[s2]], method = "spearman",
                                    use = "complete.obs"))
    out[[length(out) + 1]] <- data.table(cohort = nm, stat1 = s1, stat2 = s2,
                                         n = sum(complete.cases(x[, ..STATS])),
                                         rho = unname(ct$estimate),
                                         p = ct$p.value)
  }
}

## ---------------------------------------------------------------------------
## 2b. 上移组合内（z50 > 0）：正文对这四条的并列引用主要发生在这一子集
##     ex_trail 在这里只有一支定义（z25 - z50），与 asym 的不对称来自分支切换
## ---------------------------------------------------------------------------
cat("\n== 上移组合内（z50 > 0）==\n")
for (nm in names(per)) {
  x <- per[[nm]][z50 > 0]
  if (nrow(x) < 8) { cat(sprintf("  %-10s n=%d，样本太小，跳过\n", nm, nrow(x))); next }
  cat(sprintf("  %-10s n=%d\n", nm, nrow(x)))
  print(round(rho_mat(x), 3))
  for (i in 1:(length(STATS) - 1)) for (j in (i + 1):length(STATS)) {
    s1 <- STATS[i]; s2 <- STATS[j]
    ct <- suppressWarnings(cor.test(x[[s1]], x[[s2]], method = "spearman",
                                    use = "complete.obs"))
    out[[length(out) + 1]] <- data.table(cohort = paste0(nm, "_up"), stat1 = s1, stat2 = s2,
                                         n = nrow(x), rho = unname(ct$estimate),
                                         p = ct$p.value)
  }
}

## ---------------------------------------------------------------------------
## 3. 队列内秩标准化后合并（跨队列不混幅度，只混秩）
## ---------------------------------------------------------------------------
pool <- rbindlist(lapply(names(per), function(nm) {
  x <- copy(per[[nm]])
  for (s in STATS) x[, (s) := as.numeric(rank(get(s), na.last = "keep")) / .N]
  x[, cohort := nm]
  x
}), use.names = TRUE)
rp <- rho_mat(pool)
cat(sprintf("\n== 合并（队列内秩标准化，n=%d）==\n", nrow(pool)))
print(round(rp, 3))
for (i in 1:(length(STATS) - 1)) for (j in (i + 1):length(STATS)) {
  s1 <- STATS[i]; s2 <- STATS[j]
  ct <- suppressWarnings(cor.test(pool[[s1]], pool[[s2]], method = "spearman"))
  out[[length(out) + 1]] <- data.table(cohort = "pooled_rank", stat1 = s1, stat2 = s2,
                                       n = nrow(pool), rho = unname(ct$estimate),
                                       p = ct$p.value)
}

res <- rbindlist(out)
res[, `:=`(rho = round(rho, 4), p = signif(p, 3))]
setorder(res, cohort, -rho)
fwrite(res, file.path(tab, "shape_stat_correlation.csv"))

cat("\n== 高相关对（|rho| >= 0.7，各队列并列）==\n")
hi <- res[abs(rho) >= 0.7]
for (i in seq_len(nrow(hi)))
  cat(sprintf("  %-11s %-9s ~ %-9s  rho=%+.3f  p=%.3g  n=%d\n",
              hi$cohort[i], hi$stat1[i], hi$stat2[i], hi$rho[i], hi$p[i], hi$n[i]))

cat("\n写出：shape_stat_correlation.csv\n")
