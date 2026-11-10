# fig27：Fig 2 测量效度——为什么必须用单细胞，以及批次威胁已排除
#   A  pseudo-bulk 之间高度相似（20 例两两余弦）
#   B  方差分解：大部分变异在患者内的细胞之间
#   C  批次与患者高度共线（6 批次 x 20 患者）
#   D  批次效应量级：batch SD 与 patient SD、真实位移的比较
#   E  批次中心化后 63 个位移的方向稳健性（仅 2 个翻转）
#
# 数据：py/16_fig2_bulk.py + py/17_batch_sensitivity.py

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
})

tab <- "results/tables"; fig <- "results/figures"

theme_set(theme_minimal(base_size = 10) +
            theme(plot.title = element_text(face = "bold", size = 10.5),
                  plot.subtitle = element_text(size = 8, colour = "grey30"),
                  panel.grid.minor = element_blank()))

COL_B <- "#7F8C8D"; COL_T <- "#C0392B"; COL_C <- "#2C3E50"
MPNAME <- c(MP2 = "Prolif", MP8 = "Lum-sec", MP7 = "ER-I", MP3 = "Stromal-like",
            MP11 = "EMT-III", MP6 = "TGFb", MP9 = "Neural", MP4 = "IFN",
            MP5 = "IEG")
PAIRED <- c("P01", "P09", "P12", "P13", "P16", "P18", "P23")

icc  <- fread(file.path(tab, "nc2025_fig2_icc.csv"))
sim  <- fread(file.path(tab, "nc2025_fig2_pseudobulk_sim.csv"))
ks   <- fread(file.path(tab, "nc2025_fig2_pair_ks.csv"))
bmap <- fread(file.path(tab, "nc2025_fig2_batch_map.csv"))
beff <- fread(file.path(tab, "nc2025_fig2_batch_effect.csv"))
bdel <- fread(file.path(tab, "nc2025_fig2_batch_delta.csv"))
sv   <- fread(file.path(tab, "nc2025_fig2_sim_vs_ks.csv"))

# ================= A：pseudo-bulk 相似性 vs 细胞级差异 =================
r_sv <- cor(sv$cosine, sv$ks_median)
hi   <- sv[cosine > 0.985]
n_hi <- nrow(hi)
rng_hi <- sprintf("%.2f-%.2f", min(hi$ks_median), max(hi$ks_median))

pA <- ggplot(sv, aes(x = cosine, y = ks_median)) +
  annotate("rect", xmin = 0.985, xmax = 1.0, ymin = -Inf, ymax = Inf,
           fill = COL_T, alpha = 0.07) +
  geom_point(colour = COL_B, size = 1.4, alpha = 0.70) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              colour = COL_C, linewidth = 0.5, linetype = "22") +
  ggrepel::geom_text_repel(
    data = hi[which.max(ks_median)],
    aes(label = sprintf("%s / %s\ncosine %.3f, KS %.2f",
                        p1, p2, cosine, ks_median)),
    size = 1.9, colour = COL_T, direction = "y", nudge_y = 0.26,
    lineheight = 0.9, min.segment.length = 0, segment.size = 0.2, seed = 2) +
  scale_y_continuous(limits = c(0, 1.0)) +
  labs(title = "A  Similar pseudo-bulk, different cell distributions",
       subtitle = sprintf(paste0(
         "Patient pairs, n = %d. Pearson r = %.2f: bulk carries some information.\n",
         "But at cosine > 0.985 (n = %d pairs) the median KS still spans %s."),
         nrow(sv), r_sv, n_hi, rng_hi),
       x = "pseudo-bulk cosine similarity",
       y = "median KS across 9 programmes") +
  theme(plot.title = element_text(size = 10.2))

# ================= B：方差分解 =================
icc[, icc_num := as.numeric(icc)]
MPORD_ICC <- icc[order(icc_num)]$mp
icc[, mp := factor(mp, levels = MPORD_ICC)]
il <- melt(icc[, .(mp, between = pct_between, within = pct_within)],
           id.vars = "mp", variable.name = "part", value.name = "pct")
il[, part := factor(part, levels = c("within", "between"),
                    labels = c("between cells of one patient",
                               "between patients"))]

pB <- ggplot(il, aes(x = mp, y = pct, fill = part)) +
  geom_col(width = 0.70, colour = "white", linewidth = 0.4) +
  geom_text(data = icc, aes(x = mp, y = pct_between, label = pct_between),
            inherit.aes = FALSE, hjust = 1.25, size = 2.3, colour = "white") +
  geom_hline(yintercept = 50, linetype = "12", colour = "grey40",
             linewidth = 0.35) +
  scale_fill_manual(values = c("between cells of one patient" = COL_B,
                               "between patients" = COL_C), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)),
                     breaks = c(0, 25, 50, 75, 100)) +
  coord_flip() +
  labs(title = "B  Variance is at least as large within as between patients",
       subtitle = "Variance of programme score decomposed by patient (%).",
       x = NULL, y = "% of total variance") +
  theme(legend.position = "bottom", legend.key.size = unit(0.30, "cm"),
        legend.text = element_text(size = 6.4),
        axis.text.y = element_text(size = 7.4),
        plot.title = element_text(size = 10.2))

# ================= C：配对病例的批次轨迹 =================
pl <- unique(rbind(
  bdel[, .(patient, stage = "Baseline", batch = batch_base)],
  bdel[, .(patient, stage = "Last on-treatment", batch = batch_last)]))
pl[, stage := factor(stage, levels = c("Baseline", "Last on-treatment"))]
pl[, batch := factor(batch, levels = rev(sort(unique(c(bdel$batch_base,
                                                       bdel$batch_last)))))]
pl[, moved := patient %in% bdel[batch_base != batch_last, patient]]
pl[, stage_x := as.integer(stage)]
# 抖动避免同一批次上多条线重叠
pl[, yj := as.integer(batch) + (as.integer(factor(patient)) %% 5 - 2) * 0.075]

pC <- ggplot(pl, aes(x = stage_x, y = yj)) +
  geom_line(aes(group = patient, colour = moved), linewidth = 0.55) +
  geom_point(aes(colour = moved), size = 1.9) +
  ggrepel::geom_text_repel(data = pl[stage_x == 2 & moved == TRUE],
                           aes(label = patient, colour = moved),
                           size = 2.3, nudge_x = 0.25, segment.size = 0.2,
                           min.segment.length = 0, seed = 3) +
  scale_x_continuous(breaks = 1:2,
                     labels = c("Baseline", "Last on-\ntreatment"),
                     limits = c(0.85, 2.55)) +
  scale_y_continuous(breaks = 1:6,
                     labels = rev(sort(unique(c(bdel$batch_base,
                                                bdel$batch_last))))) +
  scale_colour_manual(values = c(`FALSE` = COL_B, `TRUE` = COL_T),
                      labels = c("same batch", "batch changes"),
                      name = NULL) +
  labs(title = "C  Only one paired case crosses batches",
       subtitle = "Six of seven paired cases stay in one batch; other six: P01 P09 P12 P13 P18 P23.",
       x = NULL, y = "batch") +
  theme(legend.position = "bottom", legend.key.size = unit(0.30, "cm"),
        legend.text = element_text(size = 6.4),
        panel.grid.major.x = element_blank(),
        axis.text.y = element_text(size = 7.0),
        axis.text.x = element_text(size = 7.4),
        plot.title = element_text(size = 10.2))

# ================= D：批次效应量级 =================
bd <- beff[, .(mp, batch_sd, patient_sd,
               delta = as.numeric(delta_max_abs))]
bl <- melt(bd, id.vars = "mp",
           variable.name = "kind", value.name = "val")
bl[, kind := factor(kind, levels = c("batch_sd", "patient_sd", "delta"),
                    labels = c("batch effect (SD of 6 batch means)",
                               "between-patient (SD)",
                               "largest observed shift"))]
bl[, mp := factor(mp, levels = bd[order(patient_sd)]$mp)]

pD <- ggplot(bl, aes(x = val, y = mp, colour = kind, shape = kind)) +
  geom_line(aes(group = mp), colour = "grey80", linewidth = 0.4) +
  geom_point(size = 2.0) +
  scale_colour_manual(values = c("batch effect (SD of 6 batch means)" = COL_T,
                                 "between-patient (SD)" = COL_B,
                                 "largest observed shift" = COL_C), name = NULL) +
  scale_shape_manual(values = c(15, 16, 17), name = NULL) +
  scale_x_log10() +
  labs(title = "D  Batch effect is smaller than the effect being measured",
       subtitle = "Programme score units. Every programme shifts further than its batch spread.",
       x = "SD or absolute shift (log scale)", y = NULL) +
  theme(legend.position = "bottom", legend.key.size = unit(0.30, "cm"),
        legend.text = element_text(size = 5.9),
        legend.box.margin = margin(0, 0, 0, 0),
        axis.text.y = element_text(size = 7.4),
        plot.title = element_text(size = 10.2))

# ================= E：批次校正稳健性 =================
bdel[, raw := as.numeric(delta_raw)]
bdel[, cor := as.numeric(delta_batchcorr)]
bdel[, flip := raw * cor < 0]
n_flip <- sum(bdel$flip)
r_rho <- cor(bdel$raw, bdel$cor, method = "spearman")

pE <- ggplot(bdel, aes(x = raw, y = cor)) +
  geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "grey55",
              linewidth = 0.4) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_point(aes(colour = flip), size = 1.7, alpha = 0.85) +
  ggrepel::geom_text_repel(data = bdel[flip == TRUE],
                           aes(label = paste0(patient, " ", mp)),
                           size = 2.2, colour = COL_T, min.segment.length = 0,
                           segment.size = 0.25, seed = 1) +
  scale_colour_manual(values = c(`FALSE` = COL_B, `TRUE` = COL_T),
                      labels = c("direction kept", "direction flipped"),
                      name = NULL) +
  labs(title = "E  Shifts survive batch centring",
       subtitle = sprintf(paste("Each point is one programme in one patient (n = %d).",
                                "Spearman rho = %.3f; %d of %d flip sign."),
                          nrow(bdel), r_rho, n_flip, nrow(bdel)),
       x = "shift, raw (score units)", y = "shift, batch-centred") +
  theme(legend.position = "bottom", legend.key.size = unit(0.30, "cm"),
        legend.text = element_text(size = 6.4),
        plot.title = element_text(size = 10.2))

# ================= 拼图 =================
design <- "
AABB
CCDD
EEEE
"
fig2 <- wrap_plots(A = pA, B = pB, C = pC, D = pD, E = pE, design = design) +
  plot_layout(heights = c(1.15, 1.20, 0.90))

## Figure-level title and in-figure explanatory block removed 2026-11-10: the
## legend carries both. Panel tags and one-line panel titles are kept.
ggsave(file.path(fig, "fig27_fig2_validity.png"), fig2,
       width = 10.0, height = 8.6, dpi = 300, bg = "white")
cat("-> fig27_fig2_validity.png\n")

cat(sprintf("\n-- 核对 --\n余弦 vs KS: r = %.3f；余弦中位 %.3f\n",
            r_sv, median(sv$cosine)))
cat(sprintf("余弦>0.985 的 %d 对：KS 中位 %.3f - %.3f\n",
            n_hi, min(hi$ks_median), max(hi$ks_median)))
cat(sprintf("ICC 中位 %.3f（范围 %.3f-%.3f）\n",
            median(icc$icc), min(icc$icc), max(icc$icc)))
cat(sprintf("批次校正：%d/%d 翻转，rho=%.4f\n", n_flip, nrow(bdel), r_rho))
print(bdel[flip == TRUE, .(patient, mp, raw, cor)])
