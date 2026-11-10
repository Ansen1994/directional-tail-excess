# fig26：Fig 1 队列与设计（6 图方案的第一张）
#   A  研究设计与时点：两臂、五个时点
#   B  patient x timepoint 肿瘤细胞数矩阵（配对 7 例与末次时点标注）
#   C  配对 vs 未配对 的 RCB 构成（软肋：配对组 1 好 / 6 差）
#   D  9 个元程序的定义（top 基因）
#   E  分析漏斗：249,379 -> 93,621 -> 65,206 -> 24,166（克隆层 6 例，从表里数）
#
# 口径与 py/13_clone_layering_all.py 完全一致：
#   末次时点 LAST = {01,09,12,13: T3; 16,18: T4; 23: T2}，不是统一 T4。

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
})

tab <- "results/tables"; fig <- "results/figures"; raw <- "data/raw/nc2025"

theme_set(theme_minimal(base_size = 10) +
            theme(plot.title = element_text(face = "bold", size = 11),
                  plot.subtitle = element_text(size = 8, colour = "grey30"),
                  panel.grid.minor = element_blank()))

COL_B <- "#7F8C8D"   # 基线 / 中性
COL_T <- "#C0392B"   # 随访 / 强调
COL_C <- "#2C3E50"   # 深色文字块

co <- fread(file.path(tab, "nc2025_fig1_cohort.csv"))
co[, tp := factor(tp, levels = c("T1", "T2", "T3", "T4", "T5"))]

# 行序：按臂分组，组内按患者号
ARMORD <- c("Chemo->Combo", "ICI->Combo", "not reported")
co[, arm := factor(arm, levels = ARMORD)]
PORD <- co[, .(arm = arm[1]), by = patient][order(arm, patient)]$patient
co[, patient := factor(patient, levels = rev(PORD))]

# ---------------- A：设计与时点 ----------------
tp_lab <- data.table(
  x = 1:5, tp = c("T1", "T2", "T3", "T4", "T5"),
  lab = c("Baseline", "W3D1\nafter mono", "W7D1\non combo",
          "Pre-surgery", "Surgery"))
lane <- data.table(
  arm = c("Chemo -> Combo", "ICI -> Combo"), y = c(2.35, 1.15))
seg <- lane[, .(arm, x = 1, xend = 5, y, yend = y)]

pA <- ggplot() +
  geom_segment(data = seg, aes(x = x, xend = xend, y = y, yend = yend),
               colour = "grey60", linewidth = 0.8,
               arrow = arrow(length = unit(0.16, "cm"), type = "closed",
                             ends = "last")) +
  geom_point(data = lane[, .(arm, y)], aes(x = 1, y = y),
             size = 3.6, colour = COL_B, shape = 21, fill = "white", stroke = 1.2) +
  geom_point(data = lane[, .(arm, y)], aes(x = 2, y = y),
             size = 3.0, colour = "grey45", shape = 21, fill = "grey90") +
  geom_point(data = lane[, .(arm, y)], aes(x = 3, y = y),
             size = 3.0, colour = "grey45", shape = 21, fill = "grey90") +
  geom_point(data = lane[, .(arm, y)], aes(x = 4, y = y),
             size = 3.0, colour = COL_T, shape = 21, fill = "white", stroke = 1.0) +
  geom_point(data = lane[, .(arm, y)], aes(x = 5, y = y),
             size = 3.0, colour = COL_T, shape = 23, fill = "white", stroke = 1.0) +
  geom_text(data = lane[, .(arm, y)], aes(x = 0.62, y = y, label = arm),
            hjust = 1, size = 2.9, colour = "grey20") +
  geom_text(data = tp_lab, aes(x = x, y = 0.52, label = lab),
            size = 2.5, colour = "grey30", lineheight = 0.92, vjust = 1) +
  geom_text(data = tp_lab, aes(x = x, y = 3.02, label = tp),
            size = 3.1, fontface = "bold", colour = "grey15") +
  annotate("rect", xmin = 1.0, xmax = 2.0, ymin = 0.82, ymax = 2.68,
           fill = "grey70", alpha = 0.10) +
  annotate("rect", xmin = 2.0, xmax = 3.0, ymin = 0.82, ymax = 2.68,
           fill = COL_T, alpha = 0.06) +
  coord_cartesian(xlim = c(-1.6, 5.35), ylim = c(0.30, 3.15), clip = "off") +
  labs(title = "A  Design and sampling timepoints") +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank(), plot.title = element_text(size = 10.5))

# ---------------- B：patient x timepoint 肿瘤细胞数 ----------------
hb <- co[, .(patient, tp, arm, n_tumor, n_total, is_paired, is_last)]
hb[, lab := ifelse(is_paired == "yes", paste0(patient, " *"), as.character(patient))]
hb[, lab := factor(lab, levels = rev(unique(
  hb[order(patient)][, .(lab = lab[1]), by = patient][order(
    match(patient, PORD))]$lab)))]
hb[, has := n_total > 0]

pB <- ggplot(hb, aes(x = tp, y = lab)) +
  geom_tile(aes(fill = ifelse(has, log10(n_tumor + 1), NA)),
            colour = "white", linewidth = 0.55) +
  geom_tile(data = hb[is_last == "yes"], fill = NA, colour = "grey10",
            linewidth = 1.05) +
  geom_text(aes(label = ifelse(has,
                               ifelse(n_tumor >= 1000,
                                      sprintf("%.1fk", n_tumor / 1000),
                                      as.character(n_tumor)),
                               "-")),
            size = 2.0,
            colour = ifelse(hb$has & hb$n_tumor >= 400, "white", "grey25")) +
  facet_grid(arm ~ ., scales = "free_y", space = "free_y", switch = "y",
             labeller = as_labeller(c(`Chemo->Combo` = "Chemo -> Combo",
                                      `ICI->Combo` = "ICI -> Combo",
                                      `not reported` = "not reported"))) +
  scale_fill_gradient(low = "#EAF2F8", high = COL_C, na.value = "grey94",
                      name = "cells", limits = c(0.9, 4.0),
                      breaks = c(1, 2, 3), labels = c("10", "100", "1000")) +
  labs(title = "B  Tumour cells per sample",
       subtitle = "* paired case. Box = last analysable biopsy.",
       x = NULL, y = NULL) +
  theme(panel.grid = element_blank(),
        strip.placement = "outside",
        strip.text.y.left = element_text(size = 7.0, angle = 90),
        axis.text.y = element_text(size = 7.0, family = "mono"),
        axis.text.x = element_text(size = 7.4),
        legend.key.height = unit(0.50, "cm"),
        legend.key.width = unit(0.26, "cm"),
        legend.title = element_text(size = 6.4),
        legend.text = element_text(size = 6.0),
        legend.margin = margin(0, 0, 0, 0),
        plot.title = element_text(size = 10.5))

# ---------------- C：配对 vs 未配对 的 RCB 构成 ----------------
cc <- co[!duplicated(patient)]
cc[, grp := ifelse(is_paired == "yes", "Paired (n = 7)", "Not paired (n = 13)")]
cc[, rcb3 := fifelse(rcb == "RCB 0-I", "RCB 0-I (favourable)",
             fifelse(rcb == "RCB II-III", "RCB II-III (unfavourable)",
                     "not reported"))]
cc[, rcb3 := factor(rcb3, levels = c("RCB 0-I (favourable)",
                                     "RCB II-III (unfavourable)",
                                     "not reported"))]
tab_c <- cc[, .N, by = .(grp, rcb3)]
tot_c <- tab_c[, .(N = sum(N)), by = grp]
tab_c <- merge(tab_c, tot_c, by = "grp", suffixes = c("", "_tot"))

pC <- ggplot(tab_c, aes(x = grp, y = N, fill = rcb3)) +
  geom_col(width = 0.62, colour = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(N > 0, N, "")),
            position = position_stack(vjust = 0.5), size = 2.9,
            colour = "white", fontface = "bold") +
  geom_text(data = unique(tab_c[, .(grp, N_tot)]), aes(x = grp, y = N_tot + 0.55,
                                                       label = paste0("n = ", N_tot)),
            inherit.aes = FALSE, size = 2.8, colour = "grey30") +
  scale_fill_manual(values = c("RCB 0-I (favourable)" = "#5DADE2",
                               "RCB II-III (unfavourable)" = COL_T,
                               "not reported" = "grey80"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(title = "C  Response composition",
       subtitle = "Only 1 of the 7 paired cases achieved RCB 0-I",
       x = NULL, y = "patients") +
  theme(legend.position = "right", legend.key.size = unit(0.32, "cm"),
        legend.text = element_text(size = 6.4),
        axis.text.x = element_text(size = 8),
        plot.title = element_text(size = 10.5))

# ---------------- D：9 个元程序的 top 基因 ----------------
gmt <- readLines(file.path(raw, "MP_Programs.gmt"))
mpg <- rbindlist(lapply(gmt, function(ln) {
  f <- strsplit(ln, "\t")[[1]]
  data.table(mp = f[1], gene = f[3:min(length(f), 9)])
}))
mpg[, rk := seq_len(.N), by = mp]

MPORD <- c("MP2", "MP8", "MP7", "MP3", "MP11", "MP6", "MP9", "MP4", "MP5")
MPNAME <- c(MP2 = "Prolif", MP8 = "Lum-sec", MP7 = "ER-I", MP3 = "Stromal-like",
            MP11 = "EMT-III", MP6 = "TGFb", MP9 = "Neural", MP4 = "IFN",
            MP5 = "IEG")
mpg[, mp := factor(mp, levels = MPORD)]
mpg[, rk := factor(rk, levels = rev(sort(unique(rk))))]

pD <- ggplot(mpg, aes(x = mp, y = rk)) +
  geom_tile(fill = "grey97", colour = "white", linewidth = 0.6) +
  geom_text(aes(label = gene), size = 2.15, colour = "grey20") +
  scale_x_discrete(labels = setNames(paste0(MPORD, "\n", MPNAME[MPORD]), MPORD)) +
  labs(title = "D  Nine meta-programmes (top 7 genes each)",
       subtitle = "Programme identity is author-defined; no single programme dominates the response",
       x = NULL, y = "rank") +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 6.8, lineheight = 0.85),
        axis.text.y = element_text(size = 6.4),
        axis.title.y = element_text(size = 7),
        plot.title = element_text(size = 10.5))

# ---------------- E：分析漏斗 ----------------
## 末级口径 = 真正进入克隆层分析的细胞：有亚克隆标签、且两时点可比（P12 基线 0 个亚克隆，
## 不在其中，见 nc2025_v4_clone_verdict.csv）。从表里数，不再硬编码。
cl <- fread(file.path(tab, "nc2025_v4_clone_cells.csv"))
N_CL <- nrow(cl); N_CLP <- uniqueN(cl$patient)
LAB_CL <- sprintf("Clone-resolved\n(%d patients)", N_CLP)
fun <- data.table(
  step = factor(c("All cells", "Epithelial\ncompartment", "Tumour cells", LAB_CL),
                levels = rev(c("All cells", "Epithelial\ncompartment",
                               "Tumour cells", LAB_CL))),
  n = c(249379, 93621, 65206, N_CL))

pE <- ggplot(fun, aes(x = n, y = step)) +
  geom_col(fill = c(COL_C, "#5D6D7E", "#85929E", COL_T), width = 0.60) +
  geom_text(aes(label = sprintf("%s  (%.1f%%)", comma(n), 100 * n / 249379)),
            hjust = -0.06, size = 2.6, colour = "grey15") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.16)),
                     labels = comma) +
  labs(title = "E  Analysis funnel", x = "cells", y = NULL) +
  theme(axis.text.y = element_text(size = 6.8, lineheight = 0.85),
        axis.text.x = element_text(size = 6.4),
        axis.title.x = element_text(size = 7),
        plot.title = element_text(size = 10.5))

# ---------------- 拼图 ----------------
design <- "
AAAB
CCCB
DDDB
EEEE
"

fig1 <- wrap_plots(A = pA, B = pB, C = pC, D = pD, E = pE, design = design) +
  plot_layout(heights = c(1.15, 1.05, 1.15, 0.95))

## Figure-level title and the in-figure explanatory block were removed on
## 2026-11-10: the legend carries the title and the panel-by-panel reading, and
## production sets figure titles in the legend, not inside the artwork. Panel
## tags and the one-line panel titles are kept.
ggsave(file.path(fig, "fig26_fig1_cohort.png"), fig1,
       width = 10.0, height = 8.5, dpi = 300, bg = "white")
cat("-> fig26_fig1_cohort.png\n")

# 数字核对
cat("\n-- 核对 --\n")
cat("配对例数:", uniqueN(co[is_paired == "yes"]$patient), "\n")
print(cc[, .N, by = .(grp, rcb3)][order(grp, rcb3)])
cat("末次时点:\n"); print(co[is_last == "yes", .(patient, tp, n_tumor)])
