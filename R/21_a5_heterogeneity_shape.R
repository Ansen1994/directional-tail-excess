# fig21：程序响应的患者异质性（Fig 5）
#   A  9 个程序 x 7 例的标准化中位位移 z50 热图（每列标出主导程序）
#   B  每例的主导程序与次位程序（lollipop）
#
# fig22：分布形状判别（Fig 7）
#   A  z50 vs asym 散点：区分整体位移 / 尾部变化 / 展宽
#   B  代表性案例的分位数曲线（整体平移 vs 尾部扩张）
#   C  密度曲线：P16-MP11 与 P23-MP7

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(ggrepel); library(scales)
})

tab <- "results/tables"; fig <- "results/figures"

shape <- fread(file.path(tab, "nc2025_v3_shape.csv"))
dom   <- fread(file.path(tab, "nc2025_v3_dominant.csv"))
dens  <- fread(file.path(tab, "nc2025_v3_density.csv"))

# ---- 统一标签 ----
MPNAME <- c(MP2 = "MP2  Proliferation", MP8 = "MP8  Luminal-secretory",
            MP7 = "MP7  ER-I", MP3 = "MP3  Collagen/stroma",
            MP11 = "MP11  EMT-III", MP6 = "MP6  TGFb/EMT",
            MP9 = "MP9  Neural-like", MP4 = "MP4  IFN response",
            MP5 = "MP5  Immediate-early")
MPORD <- c("MP2", "MP8", "MP7", "MP3", "MP11", "MP6", "MP9", "MP4", "MP5")
PATORD <- c("01", "13", "16", "09", "12", "18", "23")   # 臂分组：Chemo | ICI

key <- function(x) sprintf("P%02d", as.integer(x))
for (d in list(shape, dom, dens)) d[, patient := key(patient)]

shape[, mp := factor(mp, levels = rev(MPORD))]
shape[, patient := factor(patient, levels = key(PATORD))]
dom[, patient := factor(patient, levels = key(PATORD))]

grp <- unique(dom[, .(patient, rcb, arm)])
grp[, arm_lab := ifelse(arm == "Chemo", "Chemo", "ICI")]
grp[, col_lab := paste0(patient, "\n", ifelse(rcb == "favorable", "Fav", "Unfav"),
                        " / ", arm_lab)]

# ================= Fig 5A 热图 =================
sh <- copy(shape)
sh[, z_plot := pmin(pmax(as.numeric(z50), -4), 4)]
sh[, lab := ifelse(abs(as.numeric(z50)) >= 4,
                   sprintf("%+.1f", as.numeric(z50)),
                   sprintf("%+.1f", as.numeric(z50)))]
sh[, is_top := FALSE]
for (p in levels(sh$patient)) {
  idx <- which(sh$patient == p)
  if (length(idx) == 0) next
  top <- idx[which.max(abs(as.numeric(sh$z50[idx])))]
  sh$is_top[top] <- TRUE
}
sh <- merge(sh, grp[, .(patient, col_lab)], by = "patient")
sh[, col_lab := factor(col_lab, levels = grp$col_lab[order(match(grp$patient, key(PATORD)))])]

pa <- ggplot(sh, aes(x = col_lab, y = mp, fill = z_plot)) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_tile(data = sh[is_top == TRUE], colour = "black", linewidth = 0.9, fill = NA) +
  geom_text(aes(label = lab,
                colour = abs(z_plot) > 2.2), size = 2.9, fontface = "bold") +
  scale_fill_gradient2(low = "#2E5FA3", mid = "white", high = "#C0392B",
                       midpoint = 0, limits = c(-4, 4),
                       name = "Shift in median\n(IQR units)",
                       breaks = c(-4, -2, 0, 2, 4),
                       labels = c("-4", "-2", "0", "+2", "+4")) +
  scale_colour_manual(values = c(`FALSE` = "grey15", `TRUE` = "white"), guide = "none") +
  scale_y_discrete(labels = MPNAME[levels(sh$mp)]) +
  labs(x = NULL, y = NULL,
       title = "A  Each patient mobilises a different programme",
       subtitle = "Shift in median score, standardised by baseline IQR. Black box = dominant programme of that patient. Colour is censored at |z| > 4; printed values are exact.") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 8.6, lineheight = 0.95),
        axis.text.y = element_text(size = 9),
        plot.title = element_text(face = "bold", size = 11.5),
        plot.subtitle = element_text(size = 8.2, colour = "grey30"),
        legend.position = "right",
        legend.key.height = unit(1.1, "cm"),
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8))
pa

# ================= Fig 5B 每例的位移谱 =================
sB <- copy(shape)[, .(patient, mp, z50n = as.numeric(z50))]
sB <- merge(sB, unique(dom[, .(patient, arm)]), by = "patient")
sB[, pat_lab := paste0(patient, "  (", ifelse(arm == "Chemo", "chemo", "ICI"), ")")]
sB[, pat_lab := factor(pat_lab,
                       levels = rev(unique(sB[order(match(patient, key(PATORD)))]$pat_lab)))]
sB[, is_top := FALSE]
for (p in unique(sB$patient)) {
  idx <- which(sB$patient == p)
  sB$is_top[idx[which.max(abs(sB$z50n[idx]))]] <- TRUE
}
toplab <- sB[is_top == TRUE]

pb <- ggplot(sB, aes(x = z50n, y = pat_lab)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_vline(xintercept = c(-0.5, 0.5), colour = "grey87", linetype = "dashed", linewidth = 0.35) +
  geom_point(data = sB[is_top == FALSE], colour = "grey62", size = 1.6, alpha = 0.8) +
  geom_point(data = toplab, colour = "#C0392B", size = 3.3) +
  geom_text_repel(data = toplab, aes(label = mp), size = 2.9, colour = "#8B1A1A",
                  nudge_y = 0.32, direction = "x", hjust = 0.5,
                  segment.size = 0.2, segment.colour = "grey60",
                  min.segment.length = 0, max.overlaps = 20, seed = 11) +
  scale_x_continuous(limits = c(-1.6, 9)) +
  labs(x = "Standardised shift in median score (IQR units)", y = NULL,
       title = "B  Full displacement profile of each patient",
       subtitle = "Grey = the nine programmes; red = the dominant one. Dashed lines: 0.5 effect threshold.") +
  theme_bw(base_size = 10.5) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8.2, colour = "grey30"),
        axis.text.y = element_text(size = 9))

fig21 <- pa / pb + plot_layout(heights = c(1.55, 1))
ggsave(file.path(fig, "fig21_a5_program_heterogeneity.png"), fig21,
       width = 11, height = 8.6, dpi = 300)
cat("fig21 saved\n")

# ================= Fig 7A 形状散点 =================
CALL_ORD <- c("whole shift up", "whole shift down", "right-tail expansion",
              "left truncation", "dispersion up", "dispersion down", "no change")
CALL_COL <- c("whole shift up" = "#C0392B", "whole shift down" = "#2E5FA3",
              "right-tail expansion" = "#27AE60", "left truncation" = "#16A085",
              "dispersion up" = "#E67E22", "dispersion down" = "#8E44AD",
              "no change" = "#BDC3C7")
shape[, call := factor(call, levels = CALL_ORD)]
shape[, z50n := as.numeric(z50)]
shape[, asymn := as.numeric(asym)]

lab_pts <- shape[call != "no change" & (abs(z50n) > 1.5 | abs(asymn) > 1.2)]
lab_pts[, lab := paste0(patient, " ", mp)]

qa <- ggplot(shape, aes(x = z50n, y = asymn)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_hline(yintercept = c(-0.5, 0.5), colour = "grey85", linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = c(-0.5, 0.5), colour = "grey85", linetype = "dashed", linewidth = 0.3) +
  geom_point(aes(colour = call), size = 2.4, alpha = 0.9) +
  geom_text_repel(data = lab_pts, aes(label = lab), size = 2.7, colour = "grey20",
                  min.segment.length = 0, box.padding = 0.3, max.overlaps = 30,
                  segment.colour = "grey70", segment.size = 0.22, seed = 7) +
  scale_colour_manual(values = CALL_COL, name = NULL) +
  annotate("text", x = -1.4, y = 7.5, label = "tail expands,\nmedian flat", size = 2.6,
           colour = "grey35", hjust = 0, vjust = 1) +
  annotate("text", x = -1.4, y = -3.2, label = "left end rises,\nmedian falls", size = 2.6,
           colour = "grey35", hjust = 0, vjust = 1) +
  labs(x = "Shift in median (IQR units)", y = "Shift in spread (z75 - z25)",
       title = "A  Shape of the change: where in the distribution it moves",
       subtitle = "Each point = one patient x one programme (n = 63). Dashed lines mark the 0.5 effect threshold.") +
  theme_bw(base_size = 10.5) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        legend.position = "bottom", legend.text = element_text(size = 7.6),
        legend.key.size = unit(0.35, "cm")) +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE))

# ================= Fig 7B 分位数曲线 =================
CASES <- data.table(
  patient = c("P16", "P12", "P23", "P01"),
  mp      = c("MP11", "MP4", "MP7", "MP5"),
  tag     = c("P16 / MP11 (EMT-III)", "P12 / MP4 (IFN)",
              "P23 / MP7 (ER-I)", "P01 / MP5 (imm.-early)"))
qs <- shape[mp %in% CASES$mp]
# shape 表只存了 Q25/Q50/Q75，故分位数曲线用这三点
pdat <- rbindlist(lapply(seq_len(nrow(CASES)), function(i) {
  cs <- CASES[i]
  r <- qs[patient == cs$patient & mp == cs$mp]
  if (nrow(r) == 0) return(NULL)
  data.table(tag = factor(cs$tag, levels = CASES$tag),
             q = rep(c(25, 50, 75), 2),
             val = c(as.numeric(r$q25_base), as.numeric(r$q50_base), as.numeric(r$q75_base),
                     as.numeric(r$q25_post), as.numeric(r$q50_post), as.numeric(r$q75_post)),
             phase = factor(rep(c("Baseline", "On-treatment"), each = 3),
                            levels = c("Baseline", "On-treatment")))
}))

qb <- ggplot(pdat, aes(x = q, y = val, colour = phase, group = phase)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.9) +
  facet_wrap(~ tag, nrow = 2, scales = "free_y") +
  scale_x_continuous(breaks = c(25, 50, 75), labels = c("Q25", "Q50", "Q75")) +
  scale_colour_manual(values = c(Baseline = "#7F8C8D", `On-treatment` = "#C0392B"), name = NULL) +
  labs(x = NULL, y = "Programme score",
       title = "B  Whole-distribution displacement in the largest responses",
       subtitle = "Parallel upward shift of Q25/Q50/Q75 = the whole population moves, not just a tail") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(size = 8.2),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        legend.position = "bottom", legend.key.size = unit(0.4, "cm"),
        axis.text = element_text(size = 8))

# ================= Fig 7C 密度曲线 =================
dc <- dens[(mp == "MP11" & patient %in% c("P16", "P18")) |
             (mp == "MP7" & patient == "P23") |
             (mp == "MP4" & patient == "P12")]
dc[, phase := factor(phase, levels = c("baseline", "on-treatment"),
                     labels = c("Baseline", "On-treatment"))]
dc[, tag := paste0(patient, " / ", mp)]
dc[, tag := factor(tag, levels = c("P16 / MP11", "P18 / MP11", "P23 / MP7", "P12 / MP4"))]
dc[, score := as.numeric(score)]
# 去掉极端右尾以便看清形状
dc <- dc[score <= quantile(score, 0.995)]

qc <- ggplot(dc, aes(x = score, fill = phase, colour = phase)) +
  geom_density(alpha = 0.28, linewidth = 0.5, adjust = 1.1) +
  facet_wrap(~ tag, nrow = 1, scales = "free") +
  scale_fill_manual(values = c(Baseline = "#7F8C8D", `On-treatment` = "#C0392B"), name = NULL) +
  scale_colour_manual(values = c(Baseline = "#5D6D7E", `On-treatment` = "#C0392B"), name = NULL) +
  labs(x = "Programme score per malignant cell", y = "Density",
       title = "C  Single-cell distributions behind those shifts",
       subtitle = "P16-MP11 and P23-MP7 move as a block; P18-MP11 stretches its upper tail") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(size = 8.6),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        legend.position = "bottom", legend.key.size = unit(0.4, "cm"),
        axis.text = element_text(size = 8))

fig22 <- qa / (qb | qc) + plot_layout(heights = c(1.05, 1))
ggsave(file.path(fig, "fig22_a5_shape_discriminate.png"), fig22,
       width = 13.2, height = 9.6, dpi = 300)
cat("fig22 saved\n")
