# 40：形状判据校准（Supplementary Figure 4）
#
# 目的：ex_trail 一直被当作"能区分选择与整体平移"的判据使用，但此前没有正对照。
#       本脚本用真实基线细胞构造两种已知真值：
#         (a) 纯选择  —— 只保留基线端 top p 的细胞，其余丢弃（真值：100% 选择）
#         (b) 整群平移 —— 每个细胞同幅右移 delta（真值：0% 选择）
#       再与 16 个实测上移组合对比，给出假阴性率、假阳性率与 AUC。
#
# 口径：与 R/28、R/38、R/39 完全一致
#       z = (q_post - q_base) / IQR_base,  IQR_base = q75_base - q25_base
#       ex_trail = (z25 - z50) * sign(z50)  当 z50 > 0
#                  (z75 - z50) * sign(z50)  当 z50 < 0
#       本脚本只造上移，故恒用第一式。
#
# 输出：results/figures/figS4_calibration.png / .pdf
#       results/tables/shape_calibration.csv

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
})

tab <- "results/tables"; fig <- "results/figures"

theme_set(theme_minimal(base_size = 10) +
            theme(plot.title = element_text(face = "bold", size = 10.2),
                  plot.subtitle = element_text(size = 7.3, colour = "grey30",
                                               lineheight = 1.06),
                  panel.grid.minor = element_blank(),
                  plot.title.position = "plot"))

COL_B <- "#7F8C8D"; COL_T <- "#C0392B"; COL_U <- "#2E5FA3"; COL_G <- "#7D8B99"

pat_key <- function(x) {
  y <- suppressWarnings(as.integer(x))
  ifelse(is.na(y), as.character(x), sprintf("P%02d", y))
}

dn <- fread(file.path(tab, "nc2025_v3_density.csv"))
sh <- fread(file.path(tab, "nc2025_v3_shape.csv"))
dn[, patient := pat_key(patient)]; sh[, patient := pat_key(patient)]
dn[, score := as.numeric(score)]
for (v in c("z25", "z50", "z75", "asym", "iqr_base")) sh[, (v) := as.numeric(get(v))]

QS <- c(.25, .50, .75)
etrail <- function(z25, z50, z75) ifelse(z50 > 0, z25 - z50, (z75 - z50) * -1)

# ---------- 实测：16 个上移组合 ----------
obs <- sh[z50 > 0.5, .(patient, mp, z50, et = etrail(z25, z50, z75))]
cat("实测上移组合 n =", nrow(obs), "  z50 范围",
    sprintf("%.2f-%.2f", min(obs$z50), max(obs$z50)),
    "  ex_trail 中位", sprintf("%.3f", median(obs$et)), "\n\n")

# ---------- 模拟 ----------
combos <- dn[phase == "baseline", .(n = .N), by = .(patient, mp)][n >= 100]
cat("参与模拟的组合 n =", nrow(combos), "（基线细胞 >= 100）\n")

set.seed(20261109)
B <- 200
t0 <- Sys.time()

sim <- rbindlist(lapply(seq_len(nrow(combos)), function(i) {
  pt <- combos$patient[i]; pg <- combos$mp[i]
  base <- dn[patient == pt & mp == pg & phase == "baseline", score]
  qb <- quantile(base, QS, type = 7); iqr_b <- qb[3] - qb[1]
  if (!is.finite(iqr_b) || iqr_b <= 0) return(NULL)

  sel_rep <- function(p) {
    rbindlist(lapply(seq_len(B), function(b) {
      s <- sample(base, length(base), replace = TRUE)     # 抽样噪声
      thr <- quantile(s, 1 - p, type = 7)
      post <- s[s >= thr]
      if (length(post) < 30) return(NULL)
      z <- (quantile(post, QS, type = 7) - qb) / iqr_b
      data.table(kind = "selection", param = p, z50 = z[2],
                 et = etrail(z[1], z[2], z[3]),
                 n_post = length(post))
    }), fill = TRUE)
  }
  sh_rep <- function(t) {
    delta <- t * iqr_b
    rbindlist(lapply(seq_len(B), function(b) {
      s <- sample(base, length(base), replace = TRUE)
      z <- (quantile(s + delta, QS, type = 7) - qb) / iqr_b
      data.table(kind = "shift", param = t, z50 = z[2],
                 et = etrail(z[1], z[2], z[3]), n_post = length(s))
    }), fill = TRUE)
  }

  r <- rbind(
    rbindlist(lapply(seq(0.10, 0.50, by = 0.05), sel_rep)),   # 逐档调用，p 必须是标量
    rbindlist(lapply(seq(0.2, 1.4, by = 0.1), sh_rep))
  )
  if (nrow(r) == 0) return(NULL)
  r[, `:=`(patient = pt, mp = pg)]
  r
}))

cat("模拟完成，用时", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2),
    "分钟，共", format(nrow(sim), big.mark = ","), "行\n\n")

# ---------- 校准指标 ----------
auc <- function(pos, neg) {
  if (!length(pos) || !length(neg)) return(NA_real_)
  r <- rank(c(pos, neg))            # Mann-Whitney U，避免 outer() 生成 O(n^2) 矩阵
  n1 <- as.numeric(length(pos)); n2 <- as.numeric(length(neg))
  (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * n2)
}

cal <- sim[, .(n = .N,
               et_median = median(et),
               et_q25 = quantile(et, .25),
               et_q75 = quantile(et, .75),
               frac_et_pos = mean(et > 0),
               frac_et_gt05 = mean(et > 0.5),
               z50_median = median(z50)),
           by = .(kind, param)]

# 同一 z50 区间内的公平比较
band <- sim[z50 >= 0.5 & z50 <= 1.0]
auc_band <- auc(band[kind == "selection", et], band[kind == "shift", et])
auc_all  <- auc(sim[kind == "selection", et], sim[kind == "shift", et])

cat("== 绝对分布（跨组合） ==\n")
cat("  纯选择下 ex_trail > 0 的比例：",
    sprintf("%.3f", sim[kind == "selection", mean(et > 0)]), "\n")
cat("  纯选择下 ex_trail 中位 / 上界：",
    sprintf("%.3f / %.3f", sim[kind == "selection", median(et)],
            sim[kind == "selection", max(et)]), "\n")
cat("  整群平移下 ex_trail 中位 / 95% 分位：",
    sprintf("%.3f / %.3f", sim[kind == "shift", median(et)],
            sim[kind == "shift", quantile(et, .95)]), "\n")
cat("  整群平移下 ex_trail > 0.5 的比例：",
    sprintf("%.4f", sim[kind == "shift", mean(et > 0.5)]), "（假阳性率）\n")
cat("  AUC（选择 vs 平移，全范围）：", sprintf("%.3f", auc_all), "\n")
cat("  AUC（z50 0.5-1.0 内）：", sprintf("%.3f", auc_band), "\n\n")

# ---------- 配对比较（与 R/39 同一口径：逐组合自身 top 20%，确定性） ----------
det <- rbindlist(lapply(seq_len(nrow(combos)), function(i) {
  pt <- combos$patient[i]; pg <- combos$mp[i]
  base <- dn[patient == pt & mp == pg & phase == "baseline", score]
  qb <- quantile(base, QS, type = 7); iqr_b <- qb[3] - qb[1]
  sel <- base[base >= quantile(base, 0.80, type = 7)]
  cf <- (quantile(sel, .25, type = 7) - qb[1]) / iqr_b -
        (quantile(sel, .50, type = 7) - qb[2]) / iqr_b
  data.table(patient = pt, mp = pg, cf_det = as.numeric(cf))
}))

sel20 <- sim[kind == "selection" & abs(param - 0.20) < 1e-9]
summ <- sel20[, .(sel_med = median(et),
                  sel_q05 = quantile(et, .05),
                  sel_q95 = quantile(et, .95)),
              by = .(patient, mp)]
M <- merge(merge(obs, summ, by = c("patient", "mp")), det, by = c("patient", "mp"))
M[, diff := cf_det - et]
M[, below := et < cf_det]
M[, below_q05 := et < sel_q05]

n_below <- sum(M$below)
cat("== 配对比较（16 个实测上移组合 vs 自身确定性 top-20% 反事实） ==\n")
cat("  实测低于自身反事实：", n_below, "/", nrow(M), "\n")
cat("  实测低于自身重抽样分布 5% 分位：", sum(M$below_q05), "/", nrow(M), "\n")
cat("  中位差（反事实 - 实测）：", sprintf("%.3f", median(M$diff)), "IQR\n")
cat("  中位差范围：", sprintf("%.3f ~ %.3f", min(M$diff), max(M$diff)), "\n\n")
print(M[order(-diff), .(patient, mp, z50 = round(z50, 2), obs = round(et, 3),
                        cf = round(cf_det, 3), diff = round(diff, 3), below)])
cat("\n")

fwrite(cal, file.path(tab, "shape_calibration.csv"))

# ---------- 图：配对斜率图 ----------
# 每个上移组合一条线：左端 = 若为纯选择时的 ex_trail，右端 = 实测。
# 15/16 向下倾斜 = 实测低于选择，与 R/39 的 62/66 逐组合计数同源。
Md <- M[order(-cf_det)]

p <- ggplot(Md) +
  geom_hline(yintercept = 0, linetype = "11", colour = "grey82", linewidth = .3) +
  geom_segment(aes(x = 1, xend = 2, y = cf_det, yend = et),
               colour = "grey68", linewidth = .42) +
  geom_point(aes(x = 1, y = cf_det), colour = COL_U, size = 1.85) +
  geom_point(aes(x = 2, y = et), colour = COL_T, size = 1.85) +
  scale_x_continuous(breaks = c(1, 2),
                     labels = c("if that pair\nwere selection", "observed"),
                     limits = c(0.94, 2.06)) +
  coord_cartesian(ylim = c(min(Md$et) - .3, 0.5)) +
  labs(title = "The observed excess sits below\nits own selection counterfactual",
       subtitle = sprintf("%d of %d upward pairs fall below; median gap %+.2f IQR",
                          n_below, nrow(M), median(M$diff)),
       x = NULL, y = "trailing excess (ex_trail, baseline IQR)") +
  theme(panel.grid.major.x = element_blank(),
        plot.title = element_text(lineheight = 1.05, margin = margin(b = 5)),
        axis.text.x = element_text(size = 8.2, lineheight = .95))

ggsave(file.path(fig, "figS4_calibration.png"), p, width = 4.3, height = 4.0, dpi = 340)
ggsave(file.path(fig, "figS4_calibration.pdf"), p, width = 4.3, height = 4.0)

cat("\n-> results/figures/figS4_calibration.png\n")
cat("-> results/tables/shape_calibration.csv\n")
