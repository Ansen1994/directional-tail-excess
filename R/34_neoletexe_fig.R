#!/usr/bin/env Rscript
# ==============================================================================
# 34_neoletexe_fig.R —— Fig 7：第三队列（内分泌治疗）的形态与克隆层检验
#
# 一句话：内分泌治疗下的程序位移同样是整条分布平移，且位移发生在克隆内部，
#         克隆组成变化解释不了它。
#
# A  分位位移散点（Q25 vs Q75），三点对分面；点落在对角线上 = 上下端一起动
# B  定向尾部超额地图：紫色虚线 = 29 个组合各自的选择反事实中位（+0.14 IQR），
#    上移组合无一达到；红色空心圈 = 4 个越过 0.5 IQR 描述性判据的下移 Proliferation
#    组合（"上端削平"、中位不动）。旧灰区及其 "selection would land here" 标注已删除
# C  每个程序在三个时点对上的位移：ER / 腔面分泌下降，IFN / EMT 上升，B 之后趋平
# D  克隆层：患者整体位移 vs 克隆内位移中位，几乎重合；组成反事实接近 0
#
# 排版：图内一律 ASCII；颜色沿用前六图（红=上移，蓝=下移，灰=不变）。
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
})

TAB <- "results/tables"; FIG <- "results/figures"
COL_UP <- "#C0392B"; COL_DN <- "#2E5FA3"; COL_NO <- "#B9C2CB"
COL_A <- "#7F8C8D"; COL_B <- "#E67E22"; COL_C <- "#16A085"
COL_CF <- "#6C3483"   # 选择反事实参照（与 Fig 5D/6A 同方法；避开红=上移、蓝=下移）

theme_set(theme_minimal(base_size = 10) +
            theme(plot.title = element_text(face = "bold", size = 10.2),
                  plot.subtitle = element_text(size = 7.3, colour = "grey30",
                                               lineheight = 1.06),
                  panel.grid.minor = element_blank(),
                  plot.title.position = "plot"))

sh <- fread(file.path(TAB, "neoletexe_shape.csv"))
lay <- fread(file.path(TAB, "neoletexe_clone_layer.csv"))
sh[, pair := factor(pair, levels = c("A>B", "B>C", "A>C"))]
sh[, dir := fifelse(z50 > 0, "up", "down")]
sh[, prog_lab := gsub("_", " ", program)]

CALL_COL <- c("whole shift up" = COL_UP, "whole shift down" = COL_DN,
              "dispersion up" = "#16A085", "dispersion down" = COL_A,
              "right-tail expansion" = "#E67E22", "left truncation" = "#8E44AD",
              "no change" = COL_NO)

## ---------------- A ----------------
## 用分位剖面而不是 Q25-Q75 散点：下移组合天然落在对角线下方，
## 散点读起来会像是"尾部扩张"，剖面则一眼看出整条分布是否同步平移。
prof <- melt(sh[abs(z50) > 0.5, .(patient, pair, program, dir, z10, z25, z50, z75, z90)],
             id.vars = c("patient", "pair", "program", "dir"),
             variable.name = "q", value.name = "d")
prof[, q := as.numeric(sub("z", "", q))]
med <- prof[, .(d = median(d)), by = .(q, dir)]
n_sub <- uniqueN(prof[, paste(patient, pair, program)])

sub_a <- paste0("Only the ", n_sub, " patient x programme x interval combinations ",
                "with |median shift| > 0.5 IQR are shown.\n",
                "Thin lines = one combination; thick lines = median profile. ",
                "Every quantile moves and none stays anchored; the trailing end\n",
                "lags rather than leads, which is the pattern a selection model ",
                "cannot produce. Y axis clipped.")

pA <- ggplot(prof, aes(q, d, group = interaction(patient, pair, program))) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey55") +
  geom_line(aes(colour = dir), linewidth = 0.25, alpha = 0.16) +
  geom_line(data = med, aes(q, d, colour = dir), linewidth = 1.15,
            inherit.aes = FALSE) +
  geom_point(data = med, aes(q, d, colour = dir), size = 1.1,
             inherit.aes = FALSE) +
  scale_x_continuous(breaks = c(10, 25, 50, 75, 90),
                     labels = c("Q10", "Q25", "Q50", "Q75", "Q90")) +
  scale_colour_manual(values = c(up = COL_UP, down = COL_DN), name = NULL,
                      labels = c(up = "shift up", down = "shift down")) +
  coord_cartesian(ylim = c(-3.2, 3.4)) +
  labs(x = "Quantile of the programme score",
       y = "Shift in quantile\n(baseline IQR)",
       title = "A  Every quantile moves, in the same direction",
       subtitle = sub_a)

## ---------------- B ----------------
# ---- 旧灰区（x>0.3）与 "selection signal would land here" 已删除 ----
# 两个问题：① 按 Q10 直觉画的固定界与 ex_trail 的 Q25 口径不符（NeoLetExe
# 反事实中位仅 +0.14）；② 灰区实际圈进 23 个下移点、9 个 Proliferation，
# 而副标题只说 4 个。现与 Fig 5D/6A 同口径：直读 R/39 落盘的 29 个上移组合
# 各自的 top-20% 反事实，画中位竖线；4 个越过 S7 正文 0.5 IQR 描述性判据的
# 下移 Proliferation 组合用红色空心圈标出。
# 该队列反事实的 10-90% 分位跨 -0.99~+0.36 IQR（基线分数离散，最小 -8.18），
# 画成带会盖住半张图且无信息量，故此处只画中位，分位范围写进图注。
cf7 <- fread(file.path(TAB, "selection_counterfactual_paired.csv"))
cf7 <- cf7[cohort == "NeoLetExe" & z50 > 0.5]
cf7_med  <- median(cf7$cf)
cf7_band <- as.numeric(quantile(cf7$cf, c(.10, .90)))
n_cf7    <- nrow(cf7)
sh[, trunc_ev := dir == "down" & ex_trail > 0.5]
n_trunc7 <- sum(sh$trunc_ev)
stopifnot(n_cf7 == nrow(sh[z50 > 0.5]),
          sum(sh[z50 > 0.5]$ex_trail < cf7_med) == n_cf7)
cat("\n[7B] NeoLetExe 反事实中位 ", sprintf("%.3f", cf7_med),
    " 10-90% ", sprintf("%.3f", cf7_band[1]), " ~ ", sprintf("%.3f", cf7_band[2]),
    "（带宽过大，图上只画中位）  上移低于中位 ", n_cf7, "/", n_cf7,
    "  下移越过 0.5 判据 ", n_trunc7, "\n")

pB <- ggplot(sh, aes(ex_trail, ex_lead, colour = dir)) +
  geom_vline(xintercept = cf7_med, colour = COL_CF,
             linetype = "22", linewidth = 0.55) +
  annotate("text", x = cf7_med + 0.05, y = 3.55,
           label = paste0("selection\ncounterfactual\n+",
                          sprintf("%.2f", cf7_med), " IQR"),
           size = 2.05, colour = COL_CF, hjust = 0, vjust = 0.5,
           lineheight = 0.95) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey70") +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey70") +
  geom_point(size = 1.1, alpha = 0.8) +
  geom_point(data = sh[trunc_ev == TRUE], aes(ex_trail, ex_lead),
             shape = 1, size = 3.0, colour = "#8B1A1A", stroke = 0.7,
             inherit.aes = FALSE) +
  scale_colour_manual(values = c(up = COL_UP, down = COL_DN), name = NULL,
                      labels = c(up = "shift up", down = "shift down")) +
  labs(x = "Excess at the trailing end (IQR)",
       y = "Excess at the leading end (IQR)",
       title = "B  Upward shifts move every quantile; downward ones are cut from above",
       subtitle = paste0(
         "Trailing excess = the quantile selection must inflate. Upward shifts:\n",
         "max ", sprintf("%.2f", max(sh[z50 > 0.5]$ex_trail)),
         " IQR, while leading excess reaches ", sprintf("%.2f", max(sh$ex_lead)),
         " IQR.\n",
         "Purple dashed = median +", sprintf("%.2f", cf7_med), " IQR of ",
         n_cf7, " pair-specific selection\n",
         "counterfactuals; no upward pair reaches it. Open red circles = ",
         n_trunc7, "\ndownward Proliferation shifts above the 0.5 IQR ",
         "descriptive cut, all with |z50| < 0.5.")) +
  coord_cartesian(xlim = c(-1.5, 1.5), ylim = c(-1.5, 4.8))

## ---------------- C ----------------
ord <- c("ER_signal", "Luminal_secretory", "Proliferation", "IFN_response",
         "EMT_matrix")
sh[, prog_lab := factor(gsub("_", " ", program), levels = gsub("_", " ", ord))]
pC <- ggplot(sh, aes(x = pair, y = z50)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey70") +
  geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.3,
               colour = "grey45", fill = "grey96") +
  geom_jitter(aes(colour = dir), width = 0.14, height = 0, size = 1.0, alpha = 0.75) +
  facet_wrap(~ prog_lab, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = c(up = COL_UP, down = COL_DN), guide = "none") +
  labs(x = "Interval", y = "Median shift (baseline IQR)",
       title = "C  Endocrine therapy lowers ER and luminal output, raises immune and matrix programmes",
       subtitle = paste0("A>B = first aromatase inhibitor (3 months), B>C = after the switch, ",
                         "A>C = whole course.\nBoxes = all patients. Cohort medians are ",
                         "largest in A>B and close to zero after the switch,\nthough ",
                         "individual patients still move in B>C."))

## ---------------- D ----------------
big <- lay[abs(patient_shift) > 0.5]
pD <- ggplot(big, aes(patient_shift, med_clone_shift)) +
  geom_abline(slope = 1, intercept = 0, linetype = "22", linewidth = 0.35,
              colour = "grey45") +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey70") +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey70") +
  geom_point(aes(colour = cf_frac), size = 2.0, alpha = 0.9) +
  scale_colour_gradient2(low = "#2E5FA3", mid = "grey85", high = COL_UP,
                         midpoint = 0.1, limits = c(-0.25, 0.25),
                         name = "share of shift\nexplained by\nclone composition") +
  labs(x = "Whole-patient shift (IQR)", y = "Median within-clone shift (IQR)",
       title = "D  The shift happens inside clones, not between them",
       subtitle = paste0("One point per patient x programme with |shift| > 0.5 ",
                         "(n = ", nrow(big), ", ", uniqueN(big$patient), " patients).\n",
                         sprintf("%.1f%% of these have every clone moving the same way; ",
                                 100 * mean(big$all_same)),
                         "clone composition\nexplains a median of ",
                         sprintf("%.2f%%", 100 * median(big$cf_frac, na.rm = TRUE)),
                         " of the shift.")) +
  coord_equal(xlim = c(-3.4, 3.4), ylim = c(-3.4, 3.4))

## ---------------- 组装 ----------------
## 版面顺序 = 阅读顺序 = 正文顺序：A 形状 → B 尾部 → C 方向 → D 克隆层。
## （旧版是 A / C / B|D，读者先看到方向再回看尾部，容易误以为 B 是 C 的补充。）
## B 与 D 都是散点、天然偏方，各占约 2/3 与 1/2 宽，避免整幅宽度把点压扁、
## 也避免纵轴标题被裁。
p <- pA /
  ((pB | plot_spacer()) + plot_layout(widths = c(2, 1))) /
  pC /
  (pD | plot_spacer()) +
  plot_layout(heights = c(1.00, 0.98, 0.86, 1.00))

ggsave(file.path(FIG, "fig32_fig7_neoletexe.png"), p, width = 13, height = 12.8,
       dpi = 300, bg = "white")
cat("-> fig32_fig7_neoletexe.png\n")
