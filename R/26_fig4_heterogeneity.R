# fig28：Fig 4 异质性——跨患者程序（A-D）+ 患者内轨迹（E-G）
#   全部使用统一口径：state_z / z50 = 位移 / 基线 IQR（与 fig21/fig27 一致）
#   数据源：nc2025_v3_shape.csv（末次时点 63 组合）、nc2025_v3_dominant.csv、
#           nc2025_v3_density.csv、nc2025_v2_l1_trajectory.csv（全时点 state_z）
#
# 一句话：位移最大的程序只在个别患者出现，出现频率最高的程序位移低于阈值
#         ——程序层面不存在可跨患者推广的共同靶点。
#
# 排版约定：图内一律 ASCII，避免 locale 把多字节符号渲染成乱码；
#           半宽 panel 的 subtitle 每行 <= 84 字符，否则横向溢入相邻 panel。

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

COL_B <- "#7F8C8D"; COL_T <- "#C0392B"; COL_C <- "#2C3E50"
COL_U <- "#2E5FA3"

MPNAME <- c(MP2 = "MP2 Prolif", MP8 = "MP8 Lum-sec", MP7 = "MP7 ER-I",
            MP3 = "MP3 Stromal", MP11 = "MP11 EMT-III", MP6 = "MP6 TGFb",
            MP9 = "MP9 Neural", MP4 = "MP4 IFN", MP5 = "MP5 IEG")
MPORD  <- c("MP2", "MP8", "MP7", "MP3", "MP11", "MP6", "MP9", "MP4", "MP5")
PATORD <- c("P01", "P09", "P12", "P13", "P16", "P18", "P23")   # 臂：chemo 先，ICI 后
YLAB   <- "Shift in median (baseline IQR)"
YSHORT <- "Shift (baseline IQR)"   # 竖排轴标题，半宽 panel 用短版

sh <- fread(file.path(tab, "nc2025_v3_shape.csv"))
dom <- fread(file.path(tab, "nc2025_v3_dominant.csv"))
dn  <- fread(file.path(tab, "nc2025_v3_density.csv"))
tr  <- fread(file.path(tab, "nc2025_v2_l1_trajectory.csv"))
for (d in list(sh, dom, dn, tr))
  d[, patient := sprintf("P%02d", as.integer(patient))]

sh[, z := as.numeric(z50)]
sh[, mp := factor(mp, levels = rev(MPORD))]
sh[, patient := factor(patient, levels = PATORD)]
TP <- unique(sh[, .(patient, tp, rcb, arm)])
TP[, lab := sprintf("%s\n%s %s", patient, tp, ifelse(arm == "Chemo", "chemo", "ICI"))]
TP[, lab := factor(lab, levels = TP[order(match(patient, PATORD))]$lab)]

# ================= A：9 x 7 位移热图 =================
shA <- merge(sh[, .(patient, mp, z)], TP[, .(patient, lab)], by = "patient")
shA[, is_top := FALSE]
for (p in PATORD) {
  idx <- which(shA$patient == p)
  if (!length(idx)) next
  shA$is_top[idx[which.max(abs(shA$z[idx]))]] <- TRUE
}

pA <- ggplot(shA, aes(x = lab, y = mp)) +
  geom_tile(aes(fill = pmax(pmin(z, 4), -4)), colour = "white", linewidth = 0.7) +
  geom_tile(data = shA[is_top == TRUE], colour = "grey10", linewidth = 0.85, fill = NA) +
  geom_text(aes(label = sprintf("%+.1f", z), colour = abs(z) > 2.2),
            size = 2.75, fontface = "bold") +
  scale_fill_gradient2(low = COL_U, mid = "white", high = COL_T, midpoint = 0,
                       limits = c(-4, 4), name = "Shift in median\n(baseline IQR)",
                       breaks = c(-4, -2, 0, 2, 4)) +
  scale_colour_manual(values = c(`FALSE` = "grey15", `TRUE` = "white"), guide = "none") +
  scale_y_discrete(labels = MPNAME[levels(shA$mp)]) +
  labs(title = "A  Every patient mobilises a different programme",
       subtitle = paste0("Shift in median programme score, standardised by that patient's baseline IQR. ",
                         "Black box = dominant programme of that patient.\n",
                         "Fill is censored at |z| > 4, printed values are exact. Last on-treatment ",
                         "timepoint per patient (P23: T2, four cases: T3, two: T4)."),
       x = NULL, y = NULL) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 7.6, lineheight = 0.92),
        axis.text.y = element_text(size = 7.8),
        legend.position = "right",
        legend.key.height = unit(0.95, "cm"),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6.6))

# ================= B：每例主导 / 次位程序 =================
dB <- copy(dom)
dB[, arm_lab := ifelse(arm == "Chemo", "chemo", "ICI")]
dB[, ylab := factor(sprintf("%s  (%s)", patient, arm_lab),
                    levels = rev(sprintf("%s  (%s)", PATORD,
                                         ifelse(TP[match(PATORD, patient)]$arm == "Chemo",
                                                "chemo", "ICI"))))]

pB <- ggplot(dB) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_vline(xintercept = c(-0.5, 0.5), colour = "grey87", linetype = "dashed",
             linewidth = 0.32) +
  geom_segment(aes(x = 0, xend = as.numeric(second_z50), y = ylab, yend = ylab),
               colour = "grey78", linewidth = 0.7) +
  geom_segment(aes(x = 0, xend = as.numeric(top_z50), y = ylab, yend = ylab),
               colour = COL_T, linewidth = 1.0) +
  geom_point(aes(x = as.numeric(top_z50), y = ylab), colour = COL_T, size = 2.6) +
  geom_point(aes(x = as.numeric(second_z50), y = ylab), colour = "grey55", size = 1.7) +
  geom_text(aes(x = as.numeric(top_z50), y = ylab, label = paste0("  ", top_mp)),
            hjust = 0, size = 2.2, colour = "#8B1A1A", vjust = -0.8) +
  scale_x_continuous(limits = c(-1.0, 9.6), breaks = c(0, 2, 4, 6, 8)) +
  labs(title = "B  Seven patients, five different programmes",
       subtitle = paste0("Red = dominant programme, grey = second. Five distinct programmes ",
                         "(MP5, MP7, MP4, MP3, MP9)\ntake the top slot across seven patients; none ",
                         "leads in more than two. Dashed: 0.5 threshold."),
       x = YLAB, y = NULL) +
  theme(panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 7.4))

# ================= C：一致性与幅度不可兼得 =================
conc <- sh[, .(n_pos = sum(z > 0), n_neg = sum(z < 0),
               med_abs = median(abs(z)), max_abs = max(abs(z)),
               med_signed = median(z)), by = mp]
conc[, n_cons := pmax(n_pos, n_neg)]
conc[, dir := ifelse(n_pos >= n_neg, "up", "down")]
setorder(conc, n_cons, med_abs)
MPO <- as.character(conc$mp)

shC <- merge(sh[, .(patient, mp, z)], conc[, .(mp, dir)], by = "mp")
shC[, mp := factor(as.character(mp), levels = MPO)]
shC[, agree := ifelse(dir == "up", z > 0, z < 0)]

# 一致性直接写进 y 轴标签，右侧不再留文字列，横轴可用满宽
YLAB_C <- setNames(sprintf("%s  (%d/7 %s)", MPNAME[MPO], conc$n_cons,
                           ifelse(conc$dir == "up", "up", "down")), MPO)

pC <- ggplot(shC, aes(x = z, y = mp)) +
  geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.35) +
  geom_vline(xintercept = c(-0.5, 0.5), colour = "grey88", linetype = "dashed",
             linewidth = 0.32) +
  geom_segment(data = conc, aes(x = med_signed, xend = med_signed,
                                y = as.integer(factor(as.character(mp), levels = MPO)) - 0.33,
                                yend = as.integer(factor(as.character(mp), levels = MPO)) + 0.33),
               inherit.aes = FALSE, colour = COL_C, linewidth = 0.65) +
  geom_point(aes(colour = agree, shape = agree), size = 1.95, alpha = 0.9) +
  scale_colour_manual(values = c(`TRUE` = COL_T, `FALSE` = COL_U), guide = "none") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 17), guide = "none") +
  scale_x_continuous(breaks = c(-1, 0, 1, 2, 4, 6, 8)) +
  scale_y_discrete(labels = YLAB_C) +
  labs(title = "C  Consistency and magnitude never coincide",
       subtitle = paste0("Each point = one patient; red circle = shift in the majority direction, ",
                         "blue triangle = against it.\nVertical bar = median. The three most ",
                         "consistent programmes (6 of 7) move by 0.29-0.43 IQR - below the\n0.5 ",
                         "threshold - while the three largest shifts (6.8-7.9 IQR) each occur in a ",
                         "single patient."),
       x = YLAB, y = NULL) +
  theme(panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 7.2))

# ================= D：最大几个位移的单细胞分布 =================
CASES <- data.table(patient = c("P12", "P16", "P23", "P01"),
                    mp      = c("MP4", "MP3", "MP7", "MP5"))
CASES[, zv := sh$z[match(paste(patient, mp), paste(sh$patient, sh$mp))]]
CASES[, tag := factor(sprintf("%s / %s  (%+.1f)", patient, mp, zv),
                      levels = sprintf("%s / %s  (%+.1f)", patient, mp, zv))]

dD <- rbindlist(lapply(seq_len(nrow(CASES)), function(i) {
  cs <- CASES[i]
  r <- dn[patient == cs$patient & mp == cs$mp]
  if (!nrow(r)) return(NULL)
  r <- r[score <= quantile(score, 0.995)]
  r[, .(score, phase, tag = cs$tag)]
}))
dD[, phase := factor(phase, levels = c("baseline", "on-treatment"),
                     labels = c("Baseline", "Last on-treatment"))]

pD <- ggplot(dD, aes(x = score, fill = phase, colour = phase)) +
  geom_density(alpha = 0.26, linewidth = 0.45, adjust = 1.15) +
  facet_wrap(~ tag, nrow = 2, scales = "free") +
  scale_fill_manual(values = c(Baseline = COL_B, `Last on-treatment` = COL_T), name = NULL) +
  scale_colour_manual(values = c(Baseline = "#5D6D7E", `Last on-treatment` = COL_T), name = NULL) +
  labs(title = "D  The largest shifts are whole-population moves",
       subtitle = paste0("Single-cell programme scores, four of the largest displacements, one per ",
                         "patient. The whole\ndistribution slides; no shoulder and no emerging ",
                         "subpopulation."),
       x = "Programme score per malignant cell", y = "Density") +
  theme(strip.text = element_text(size = 7),
        legend.position = "bottom", legend.key.size = unit(0.32, "cm"),
        legend.text = element_text(size = 6.4),
        panel.grid.major.y = element_blank(),
        axis.text = element_text(size = 6.8),
        axis.text.y = element_blank())

# ================= E/F：同一程序在不同患者体内方向相反 =================
tr[, z := as.numeric(state_z)]
tr[, tpn := as.integer(sub("T", "", tp))]
tr[, patient := factor(patient, levels = PATORD)]
ENDTP <- c("T1", "T2", "T3", "T4", "T5")

traj <- function(mps, keep, ttl, sub, ylim = NULL, nudge = NULL) {
  d <- tr[mp %in% mps & patient %in% keep][!is.na(z)]
  lastp <- d[, .SD[which.max(tpn)], by = patient]
  lastp[, yl := z + if (is.null(nudge)) 0 else
    ifelse(as.character(patient) %in% names(nudge), nudge[as.character(patient)], 0)]
  g <- ggplot(d, aes(x = tpn, y = z)) +
    geom_hline(yintercept = 0, colour = "grey45", linewidth = 0.35) +
    geom_hline(yintercept = c(-0.5, 0.5), colour = "grey88", linetype = "dashed",
               linewidth = 0.3) +
    geom_line(aes(group = patient, colour = patient), linewidth = 0.62) +
    geom_point(aes(colour = patient), size = 1.9) +
    geom_text(data = lastp,
              aes(y = yl, label = sprintf("%s %+.1f", tp, z), colour = patient),
              size = 2.0, hjust = 0, vjust = 0.5, show.legend = FALSE) +
    scale_x_continuous(breaks = 1:5, labels = ENDTP, limits = c(0.85, 5.65)) +
    scale_colour_manual(values = c(P01 = COL_U, P09 = "#16A085", P12 = "grey55",
                                   P13 = "#E67E22", P16 = COL_T, P18 = "#8E44AD",
                                   P23 = "#B7950B"), name = NULL) +
    labs(title = ttl, subtitle = sub, x = NULL, y = YSHORT) +
    theme(panel.grid.major.x = element_blank(),
          legend.position = "none",
          axis.text.x = element_text(size = 7.2),
          axis.text.y = element_text(size = 7.2))
  if (!is.null(ylim)) g <- g + coord_cartesian(ylim = ylim)
  g
}

pE <- traj("MP7", c("P01", "P09", "P16", "P23"),
           "E  One programme, four directions",
           paste0("MP7 (ER-I). P23 reaches +2.20 by T2, P09 climbs to +0.82, P01 falls to -0.31, ",
                  "and P16\nis -0.66 at T4 before rebounding. Lines end at each patient's last ",
                  "timepoint;\na group mean would cancel all four."),
           ylim = c(-1.35, 3.30),
           nudge = c(P01 = -0.42, P09 = 0.32, P16 = 0.30, P23 = 0.30))

pF <- traj("MP3", c("P01", "P09", "P13", "P16"),
           "F  Same direction, incomparable size",
           paste0("MP3 (stromal/collagen). Three of four rise, but by +0.42, +0.90 and +7.33 IQR ",
                  "- an 18-fold\nspread - while P09 falls to -0.51. Direction alone is not a ",
                  "finding; size is patient-specific."),
           ylim = c(-1.35, 8.75),
           nudge = c(P01 = 0.55, P09 = -0.50, P13 = 1.55)) +
  annotate("text", x = 3.92, y = 7.60, label = "+7.3", size = 2.1,
           colour = "#8B1A1A", hjust = 0.5, vjust = 0)

# ================= G：P16 的手术标本让轨迹断裂 =================
pg <- tr[patient == "P16" & !is.na(z)][, x := match(tp, c("T1", "T4", "T5"))]
pg[, hl := mp %in% c("MP11", "MP3", "MP7")]
pg[, labx := factor(MPNAME[mp])]
glv <- c("MP11 EMT-III" = COL_T, "MP3 Stromal" = "#B7950B", "MP7 ER-I" = COL_U)
lg <- pg[x == 3]

pG <- ggplot(pg, aes(x = x, y = z)) +
  geom_hline(yintercept = 0, colour = "grey45", linewidth = 0.35) +
  geom_line(data = pg[hl == FALSE], aes(group = mp), colour = "grey80", linewidth = 0.45) +
  geom_point(data = pg[hl == FALSE], colour = "grey72", size = 1.2) +
  geom_line(data = pg[hl == TRUE], aes(group = mp, colour = labx), linewidth = 0.85) +
  geom_point(data = pg[hl == TRUE], aes(colour = labx), size = 2.1) +
  annotate("text", x = 3.06, y = 2.95, label = "+2.8", size = 2.1,
           colour = "#B7950B", hjust = 0, vjust = 0) +
  annotate("text", x = 3.06, y = 0.55, label = "+0.4", size = 2.1,
           colour = COL_T, hjust = 0, vjust = -0.3) +
  annotate("text", x = 3.0, y = -1.15, label = "surgery", size = 2.2,
           colour = "grey35", hjust = 0.5, vjust = 1) +
  scale_x_continuous(breaks = 1:3, labels = c("T1\nbaseline", "T4\npre-surgery", "T5\nsurgery"),
                     limits = c(0.85, 3.6)) +
  scale_colour_manual(values = glv, name = NULL) +
  scale_y_continuous(limits = c(-1.5, 8.6)) +
  labs(title = "G  The surgical specimen breaks the trajectory",
       subtitle = paste0("P16 is the only patient with a surgical specimen, and the only one with ",
                         "five timepoints.\nFrom T4 to T5 five programmes fall back toward ",
                         "baseline (MP11 -6.4, MP3 -4.5 IQR) while four\nmove further out (MP7 +1.4, ",
                         "MP8 +1.4); two change sign. The direction is not uniform, so T5\nmeasures ",
                         "a different specimen rather than a later state, and the analysis stops at T4."),
       x = NULL, y = YSHORT) +
  theme(panel.grid.major.x = element_blank(),
        legend.position = "bottom", legend.key.size = unit(0.30, "cm"),
        legend.text = element_text(size = 6.4),
        legend.box.margin = margin(-6, 0, 0, 0),
        axis.text.x = element_text(size = 7.2),
        axis.text.y = element_text(size = 7.2),
        plot.margin = margin(2, 3, 10, 2)) +
  guides(colour = guide_legend(nrow = 1, title = NULL))

# ================= 拼图 =================
## 图级标题与图内解释段于 2026-11-10 移除（改由图注承担），与其后主图一致。
design <- "
AABB
AACC
DDEE
FFGG
"
f4 <- wrap_plots(A = pA, B = pB, C = pC, D = pD, E = pE, F = pF, G = pG,
                 design = design) +
  plot_layout(heights = c(1.10, 1.10, 1.00, 1.00))

ggsave(file.path(fig, "fig28_fig4_heterogeneity.png"), f4,
       width = 11.2, height = 10.4, dpi = 300, bg = "white")
cat("-> fig28_fig4_heterogeneity.png\n")

# ================= 导出 C 的核对表 =================
out <- conc[, .(mp = as.character(mp), n_up = n_pos, n_down = n_neg,
                n_consistent = n_cons, direction = dir,
                median_abs_z = round(med_abs, 3), max_abs_z = round(max_abs, 3),
                median_signed_z = round(med_signed, 3))]
fwrite(out, file.path(tab, "nc2025_fig4_concordance.csv"))
cat("\n-- C 核对 --\n"); print(out[order(-n_consistent, median_abs_z)])

cat("\n-- 关键数字 --\n")
cat(sprintf("主导程序：%s\n", paste(sprintf("%s=%s", dom$patient, dom$top_mp), collapse = " ")))
cat(sprintf("不同主导程序 %d 个，重复主导：%s\n", length(unique(dom$top_mp)),
            paste(names(table(dom$top_mp))[table(dom$top_mp) > 1], collapse = " ")))
big <- sh[abs(z) > 1.5][order(-abs(z))]
cat(sprintf("|z50|>1.5：%d 个组合，覆盖 %d 例（%s）\n", nrow(big),
            length(unique(big$patient)), paste(unique(big$patient), collapse = " ")))
print(big[, .(patient, mp, z = round(z, 2))])
cat(sprintf("\nP16 T4->T5：回落 %d/9，上升 %d/9，最大变化 %.2f\n",
            sum(abs(pg[x == 3]$z) < abs(pg[x == 2]$z)),
            sum(abs(pg[x == 3]$z) > abs(pg[x == 2]$z)),
            max(abs(pg[x == 3]$z - pg[x == 2]$z))))
