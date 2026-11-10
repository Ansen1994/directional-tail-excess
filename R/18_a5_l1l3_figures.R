# A5 主臂：诱导 vs 选择 的第一轮结果图
# 输入 results/tables/nc2025_v2_*.csv
# 输出 results/figures/fig18_a5_mp7_trajectory.png, fig19_a5_clone_evidence.png

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

tab <- "results/tables"
fig <- "results/figures"
dir.create(fig, showWarnings = FALSE, recursive = TRUE)

# 患者编号补零为两位字符，按数值排序定因子水平（避免 "01" 与 "1" 不匹配）
pat_key <- function(x) sprintf("%02d", as.integer(x))
pat_levels <- function(x) {
  u <- unique(pat_key(x))
  u[order(as.integer(u))]
}

l1   <- fread(file.path(tab, "nc2025_v2_l1_trajectory.csv"))
l3   <- fread(file.path(tab, "nc2025_v2_l3_within.csv"))
lost <- fread(file.path(tab, "nc2025_v2_lostclone.csv"))
call <- fread(file.path(tab, "nc2025_v2_call.csv"))
setnames(call, "call", "call_lab")
call[, patient := pat_key(patient)]

l1[, tp := factor(tp, levels = c("T1","T2","T3","T4","T5"))]
l1[, patient := factor(pat_key(patient), levels = pat_levels(patient))]

dir_col <- c("01"="#C0392B","09"="#E67E22","12"="#7F8C8D","13"="#95A5A6",
             "16"="#8E44AD","18"="#2980B9","23"="#16A085")

traj <- function(mp_val, ycol, ylab, title) {
  d <- l1[mp == mp_val]
  ggplot(d, aes(x = tp, y = get(ycol), group = patient, colour = patient)) +
    geom_line(linewidth = 0.8) +
    geom_point(aes(size = n), alpha = 0.9) +
    scale_size_continuous(range = c(1.2, 4), guide = "none") +
    scale_colour_manual(values = dir_col, name = "Patient", drop = FALSE) +
    labs(x = "Timepoint (T1 = baseline; T5 = surgery)", y = ylab, title = title) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 11))
}

p1 <- traj("MP7",  "med",     "Median MP7 score (ER-I)", "A  MP7 / ER-I: median")
p2 <- traj("MP7",  "hi_frac", "Cells above baseline Q75 (%)", "B  MP7 / ER-I: high-cell fraction")
p3 <- traj("MP11", "med",     "Median MP11 score (EMT-III)", "C  MP11 / EMT-III: median")
p4 <- traj("MP11", "hi_frac", "Cells above baseline Q75 (%)", "D  MP11 / EMT-III: high-cell fraction")

fig18 <- (p1 | p2) / (p3 | p4) +
  plot_annotation(
    title = "Program dynamics in paired HR+ breast cancer (7 patients; Nat Commun 2025 public data)",
    subtitle = "Point size = malignant cells available at that timepoint",
  )
ggsave(file.path(fig, "fig18_a5_mp7_trajectory.png"), fig18,
       width = 11, height = 7.5, dpi = 300)

# ---------------- fig19 ----------------
# 只用"判定所用的最后一个时间点"，保证三个面板口径一致
wm <- l3[subclone == "ALL", .(patient = pat_key(patient), tp, mp, within_med)]
wm[, patient := as.character(patient)]
cc <- copy(call)
cc <- merge(cc, wm, by = c("patient", "tp", "mp"), all.x = TRUE)
cc[, patient := factor(patient, levels = pat_levels(patient))]
cc[, mp := factor(mp, levels = c("MP7", "MP11"),
                  labels = c("MP7 / ER-I", "MP11 / EMT-III"))]

# A：哑铃图 —— 总体 Δ中位 vs 克隆内 Δ中位；两点重合即为"诱导"
dd <- melt(cc[, .(patient, mp, `All malignant cells` = d_med, `Within clone` = within_med)],
           id.vars = c("patient", "mp"),
           variable.name = "scope", value.name = "delta")
dd <- dd[!is.na(delta)]

pa <- ggplot(dd, aes(y = patient, x = delta)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_line(aes(group = patient), colour = "grey70", linewidth = 0.6) +
  geom_point(aes(colour = scope, shape = scope), size = 2.8) +
  scale_colour_manual(values = c("All malignant cells" = "#7F8C8D",
                                 "Within clone" = "#C0392B"), name = NULL) +
  scale_shape_manual(values = c("All malignant cells" = 1, "Within clone" = 16), name = NULL) +
  facet_wrap(~ mp, nrow = 1) +
  labs(x = "Change in median program score (last timepoint vs baseline)",
       y = "Patient",
       title = "A  Overall shift vs shift inside surviving clones",
       subtitle = "Points coincide = the shift happens inside cells (induction), not by clone replacement") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        legend.key.size = unit(0.4, "cm"),
        plot.title = element_text(face = "bold", size = 11))

# B：高表达细胞占比变化（pp）
pb <- ggplot(cc, aes(x = patient, y = d_hi, fill = d_hi > 0)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%+.0f", d_hi),
                vjust = ifelse(d_hi > 0, -0.4, 1.2)), size = 2.8) +
  scale_fill_manual(values = c(`TRUE` = "#C0392B", `FALSE` = "#2980B9"), guide = "none") +
  facet_wrap(~ mp, nrow = 1) +
  labs(x = "Patient", y = "Change in high-cell fraction (pp)",
       title = "B  Expansion of the high-program tail") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 11))

# C：克隆命运 vs 基线 MP7
lost[, patient := factor(pat_key(patient), levels = pat_levels(patient))]
lost[, fate := factor(fate, levels = c("消失", "扩张", "维持"),
                      labels = c("Lost", "Expanded", "Stable"))]

pc <- ggplot(lost, aes(x = patient, y = med_MP7_base, colour = fate, size = n_base)) +
  geom_hline(yintercept = median(lost$med_MP7_base), linetype = "dotted", colour = "grey50") +
  geom_point(alpha = 0.85, position = position_jitter(width = 0.08, height = 0)) +
  scale_size_continuous(range = c(1.5, 6), name = "Baseline cells") +
  scale_colour_manual(values = c("Lost" = "#C0392B", "Expanded" = "#27AE60",
                                 "Stable" = "#95A5A6"), name = "Clone fate",
                      drop = FALSE) +
  labs(x = "Patient", y = "Baseline median MP7 of the clone",
       title = "C  Lost and retained clones do not differ in baseline ER program") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        legend.key.size = unit(0.4, "cm"),
        plot.title = element_text(face = "bold", size = 11))

fig19 <- (pa / pb / pc) +
  plot_annotation(
    title = "Clonal selection is real, but it does not explain the ER-program rise",
    subtitle = "7 paired patients. Panel C: dotted line = median of retained clones; patient 12 has no baseline clone label"
  )
ggsave(file.path(fig, "fig19_a5_clone_evidence.png"), fig19,
       width = 10.5, height = 10.5, dpi = 300)

cat("done\n")
cat("patients plotted:", nlevels(l1$patient), "\n")
print(table(l1$patient, l1$tp))
