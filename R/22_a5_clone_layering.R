# fig23：克隆分层检验（全 6 例）
#   A  每例主导程序的克隆内分布：基线 vs 随访，同一克隆并列箱线
#   B  哑铃图：各克隆中位数的位移，灰带=基线克隆间范围，红带=随访克隆间范围
#
# 判据说明：separation = 随访最低克隆中位数 - 基线最高克隆中位数。
# > 0 表示治疗后最差的克隆都已高于治疗前最好的克隆，位移越出基线克隆异质性。

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
})

tab <- "results/tables"; fig <- "results/figures"

# 全局主题：使 patchwork 的 plot_annotation 标题自动采用统一字号
theme_set(theme_minimal(base_size = 10) +
            theme(plot.title = element_text(face = "bold", size = 12),
                  plot.subtitle = element_text(size = 8.2, colour = "grey30")))

cells <- fread(file.path(tab, "nc2025_v4_clone_cells.csv"))
layer <- fread(file.path(tab, "nc2025_v4_clone_layer.csv"))
summ  <- fread(file.path(tab, "nc2025_v4_clone_summary.csv"))
ver   <- fread(file.path(tab, "nc2025_v4_clone_verdict.csv"))

MPNAME <- c(MP2 = "Proliferation", MP8 = "Luminal-secretory", MP7 = "ER-I",
            MP3 = "Collagen/stroma", MP11 = "EMT-III", MP6 = "TGFb/EMT",
            MP9 = "Neural-like", MP4 = "IFN response", MP5 = "Immediate-early")

PATORD <- c("01", "09", "13", "16", "18", "23")
key <- function(x) sprintf("P%02d", as.integer(x))
for (d in list(cells, layer, summ, ver)) d[, patient := key(patient)]

dom <- summ[is_dominant == "yes"][order(match(patient, key(PATORD)))]
SEP <- setNames(dom$sep_z, dom$patient)

cells[, clone_f := factor(clone, levels = sort(unique(clone)))]
cells[, phase := factor(phase, levels = c("baseline", "on-treatment"))]
layer[, `:=`(mb = as.numeric(med_base), mp_ = suppressWarnings(as.numeric(med_post)))]
layer[, clone_f := factor(clone, levels = sort(unique(clone)))]

COL_B <- "#7F8C8D"   # 基线
COL_T <- "#C0392B"   # 随访

# ================= fig23：并列箱线 =================
mk_box <- function(p) {
  d <- cells[patient == p]
  s <- dom[patient == p]
  lb <- layer[patient == p & mp == s$mp]
  bmax <- max(lb$mb, na.rm = TRUE)
  pmin <- min(lb[status == "retained"]$mp_, na.rm = TRUE)
  flag <- grepl("not_separated", ver[patient == p]$verdict)
  # 仅在极端离群明显压缩箱体时裁剪视图（P13 的 MP3 有 5 倍于 Q99 的离群）
  q99 <- as.numeric(quantile(d$score, 0.99))
  top <- max(d$score)
  ytop <- if (top > 2.5 * q99) as.numeric(quantile(d$score, 0.995)) else top

  ggplot(d, aes(x = clone_f, y = score, fill = phase)) +
    geom_boxplot(position = position_dodge(0.74), width = 0.6,
                 outlier.size = 0.2, outlier.alpha = 0.22,
                 linewidth = 0.3, colour = "grey25") +
    stat_summary(fun = median, geom = "point", shape = 18, size = 1.5,
                 colour = "grey10",
                 position = position_dodge(0.74), show.legend = FALSE) +
    geom_hline(yintercept = bmax, colour = COL_B, linetype = "22", linewidth = 0.45) +
    geom_hline(yintercept = pmin, colour = COL_T, linetype = "22", linewidth = 0.45) +
    coord_cartesian(ylim = c(min(d$score), ytop)) +
    scale_fill_manual(values = c(baseline = COL_B, `on-treatment` = COL_T)) +
    labs(x = "Subclone", y = "Programme score",
         title = paste0("P", substr(p, 2, 3), "  ", s$mp, " (", MPNAME[s$mp], ")"),
         subtitle = sprintf("median shift %+.2f IQR\nseparation %+.2f IQR%s",
                            as.numeric(s$z50), as.numeric(s$sep_z),
                            ifelse(flag, "  (not separated)", ""))) +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.border = element_rect(colour = ifelse(flag, "#E67E22", "grey80"),
                                      fill = NA, linewidth = ifelse(flag, 1.1, 0.5)),
          plot.title = element_text(face = "bold", size = 10.2),
          plot.subtitle = element_text(size = 7.8, lineheight = 1.05,
                                       colour = ifelse(flag, "#B9770E", "grey35")),
          legend.position = "none",
          axis.text = element_text(size = 8))
}

ps <- lapply(key(PATORD), mk_box)
p23 <- wrap_plots(ps, nrow = 2) +
  plot_annotation(
    title = "Clonal layering of the treatment-induced programme shift",
    subtitle = paste0("Each panel: the dominant programme of one patient. Grey = baseline, red = last on-treatment biopsy before surgery.\n",
                      "Dashed lines mark the highest baseline clone (grey) and the lowest surviving clone on treatment (red).\n",
                      "In 5 of 6 patients the red line sits above the grey line, i.e. no pre-existing clone reached the on-treatment level. P09 is the exception.\n",
                      "Diamonds mark clone medians; separation is computed between clone medians, while boxes show cell-level spread within each clone."))

ggsave(file.path(fig, "fig23_a5_clone_layering.png"), p23,
       width = 9.6, height = 6.4, dpi = 300, bg = "white")
cat("-> fig23_a5_clone_layering.png\n")

# ================= fig24：哑铃图 + 区间带 =================
dl <- merge(layer[status != "new"], dom[, .(patient, mp, sep_z)],
            by = c("patient", "mp"))
band <- dl[, .(bmin = min(mb, na.rm = TRUE), bmax = max(mb, na.rm = TRUE),
               tmin = min(mp_, na.rm = TRUE), tmax = max(mp_, na.rm = TRUE)),
           by = .(patient, mp, sep_z)]
band[, patient := factor(patient, levels = key(PATORD))]

dl[, patient := factor(patient, levels = key(PATORD))]
dl[, status_lab := ifelse(status == "retained", "Retained", "Lost")]
dl_ret <- dl[!is.na(mp_)]
dl_drop <- dl[is.na(mp_)]

p24 <- ggplot(dl) +
  geom_rect(data = band, aes(xmin = bmin, xmax = bmax, ymin = -Inf, ymax = Inf),
            fill = COL_B, alpha = 0.13) +
  geom_rect(data = band, aes(xmin = tmin, xmax = tmax, ymin = -Inf, ymax = Inf),
            fill = COL_T, alpha = 0.13) +
  geom_segment(data = dl_ret,
               aes(x = mb, xend = mp_, y = clone_f, yend = clone_f),
               arrow = arrow(length = unit(0.13, "cm"), type = "closed"),
               colour = "grey30", linewidth = 0.5) +
  geom_point(data = dl, aes(x = mb, y = clone_f), colour = COL_B, size = 2.1) +
  geom_point(data = dl_ret, aes(x = mp_, y = clone_f), colour = COL_T,
             size = 2.1, shape = 17) +
  geom_text(data = dl_drop, aes(x = mb, y = clone_f, label = "lost"),
            hjust = -0.3, size = 2.3, colour = "grey45") +
  facet_wrap(~ patient, scales = "free", nrow = 2) +
  labs(x = "Programme score", y = "Subclone",
       title = "All surviving subclones shift, and in 5 of 6 patients the shift clears the baseline range",
       subtitle = paste0("Circle = baseline clone median, triangle = on-treatment median.\n",
                         "Grey band = span of all baseline clone medians; red band = span of surviving clone medians.\n",
                         "Grey 'lost' labels mark clones absent at the on-treatment timepoint.\n",
                         "P09 is the exception: its red band still overlaps the grey band, so clonal selection cannot be excluded.")) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.text = element_text(face = "bold", size = 9.5),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 8.2, colour = "grey30"),
        legend.position = "top",
        axis.text = element_text(size = 8))

ggsave(file.path(fig, "fig24_a5_clone_dumbbell.png"), p24,
       width = 9.6, height = 5.8, dpi = 300, bg = "white")
cat("-> fig24_a5_clone_dumbbell.png\n")

# ---- 控制台核对 ----
cat("\n分离度（IQR 单位）：\n")
print(dom[, .(patient, mp, z50, separation, sep_z, icc_post)])