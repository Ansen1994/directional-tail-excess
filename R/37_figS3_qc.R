#!/usr/bin/env Rscript
# ==============================================================================
# 37_figS3_qc.R
#
# 补充图 S3：质控与细胞数漏斗。
#
# 要交代三件事：
#   1. 472,737 个细胞经过 QC、去双细胞、上皮门控、CNV 判定后剩多少（漏斗）。
#   2. 有没有样本因细胞过少而使配对分析失真（7 个样本上皮 < 30）。
#   3. 三时点的 QC 指标是否可比（若手术样本 nFeature 系统性偏低，
#      位移可能部分来自技术差异）。
#
# 输出：
#   results/tables/neoletexe_qc_funnel.csv
#   results/figures/figS3_qc.png
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

TAB <- "results/tables"; FIG <- "results/figures"

qc <- fread(file.path(TAB, "neoletexe_qc_per_sample.csv"))
ct <- fread(file.path(TAB, "neoletexe_celltype.csv"))
cv <- fread(file.path(TAB, "neoletexe_cnv_cells.csv"))

## ---------------- 漏斗 ----------------
n_epi <- ct[celltype == "Epi", sum(N)]
n_all_ct <- ct[, sum(N)]
n_mal <- cv[malignant == TRUE, .N]

fn <- data.table(
  step = c("Sequenced cells", "Pass QC", "Singlets",
           "Epithelial gate", "Malignant (CNV)"),
  n = c(sum(qc$n_raw), sum(qc$n_qc), sum(qc$n_singlet), n_epi, n_mal))
fn[, pct := 100 * n / n[1]]
fwrite(fn, file.path(TAB, "neoletexe_qc_funnel.csv"))
cat("===== 漏斗 =====\n"); print(fn)

## ---------------- 样本级 ----------------
sm <- merge(qc, ct[celltype == "Epi", .(n_epi = sum(N)), by = "sample"],
            by = "sample", all.x = TRUE)
sm[is.na(n_epi), n_epi := 0]
sm <- merge(sm, cv[malignant == TRUE, .(n_mal = .N), by = "sample"],
            by = "sample", all.x = TRUE)
sm[is.na(n_mal), n_mal := 0L]
sm[, `:=`(dbl = 1 - n_singlet / n_qc,
          epi_frac = n_epi / n_singlet,
          mal_frac = n_mal / pmax(n_epi, 1))]
sm[, tp := factor(tp, levels = c("A","B","C"))]

cat("\n===== 每样本上皮细胞数 =====\n")
print(quantile(sm$n_epi, c(0,.1,.25,.5,.75,.9,1)))
cat("上皮 < 30 的样本:", sm[n_epi < 30, .N], "\n")
print(sm[n_epi < 30, .(sample, patient, tp, n_singlet, n_epi)][order(n_epi)])

cat("\n===== 三时点 QC（样本中位的中位）=====\n")
print(sm[, .(n = .N, feat = round(median(med_feature)), mt = round(median(med_mt), 2),
             dbl = round(100 * median(dbl), 1), epi = round(100 * median(epi_frac), 1),
             mal = round(100 * median(mal_frac), 1)), by = tp][order(tp)])

cat("\n===== 时点间配对检验（nFeature）=====\n")
w <- dcast(sm, patient ~ tp, value.var = "med_feature")
w <- w[complete.cases(w)]
for (pr in list(c("A","B"), c("B","C"), c("A","C"))) {
  if (nrow(w) < 5) next
  tt <- wilcox.test(w[[pr[1]]], w[[pr[2]]], paired = TRUE)
  cat(sprintf("  %s -> %s : 中位变化 %+.0f  P = %.3f\n",
              pr[1], pr[2], median(w[[pr[2]]] - w[[pr[1]]]), tt$p.value))
}

## ---------------- 出图 ----------------
BL <- "#4477AA"; RED <- "#B22222"
TPC <- c(A = "#4477AA", B = "#EE9944", C = "#66AA55")

fn[, step := factor(step, levels = rev(step))]
pA <- ggplot(fn, aes(step, n)) +
  geom_col(width = 0.62, fill = BL) +
  geom_text(aes(label = sprintf("%s (%.1f%%)", format(n, big.mark = ","), pct)),
            hjust = -0.05, size = 1.95) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 7.2e5),
                     labels = function(x) format(x, big.mark = ",", scientific = FALSE),
                     expand = expansion(mult = c(0, 0.01))) +
  labs(x = NULL, y = "Cells",
       title = "A  From raw matrix to malignant cells",
       subtitle = paste0("The deposited matrix is already filtered, so our cut-offs\n",
                         "remove nothing further; the epithelial gate drives the loss.")) +
  theme_classic(base_size = 7.5) +
  theme(plot.title = element_text(size = 8, face = "bold"),
        plot.subtitle = element_text(size = 6.4, colour = "grey35"),
        axis.text.y = element_text(size = 6.2))

pB <- ggplot(sm, aes(n_singlet, pmax(n_epi, 1), colour = tp)) +
  geom_abline(slope = 1, intercept = 0, linetype = "22",
              linewidth = 0.3, colour = "grey60") +
  geom_hline(yintercept = 30, linetype = "11", linewidth = 0.3, colour = RED) +
  geom_point(size = 1.15, alpha = 0.85) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = TPC, name = NULL) +
  annotate("text", x = 60, y = 22, label = "30 epithelial cells",
           size = 1.95, colour = RED, hjust = 0) +
  labs(x = "Singlets per biopsy (log)", y = "Epithelial cells (log)",
       title = "B  Seven biopsies fall below the pairing floor",
       subtitle = paste0("Dashed diagonal = all singlets are epithelial.\n",
                         "Those biopsies drop out of paired contrasts, not of QC.")) +
  theme_classic(base_size = 7.5) +
  theme(plot.title = element_text(size = 8, face = "bold"),
        plot.subtitle = element_text(size = 6.4, colour = "grey35"),
        legend.position = c(0.86, 0.16), legend.key.size = unit(0.28, "cm"),
        legend.text = element_text(size = 5.8), legend.background = element_blank())

mk <- melt(sm[, .(sample, tp, `Genes per cell` = med_feature, `Mitochondrial %` = med_mt)],
           id.vars = c("sample", "tp"))
pC1 <- ggplot(mk[variable == "Genes per cell"], aes(tp, value)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, fill = "grey93",
               colour = "grey40", linewidth = 0.3) +
  geom_jitter(aes(colour = tp), width = 0.13, size = 0.95, alpha = 0.9) +
  scale_colour_manual(values = TPC, guide = "none") +
  labs(x = NULL, y = "Median genes per cell",
       title = "C  Depth is comparable across biopsies") +
  theme_classic(base_size = 7.5) +
  theme(plot.title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 6))

pC2 <- ggplot(mk[variable == "Mitochondrial %"], aes(tp, value)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, fill = "grey93",
               colour = "grey40", linewidth = 0.3) +
  geom_jitter(aes(colour = tp), width = 0.13, size = 0.95, alpha = 0.9) +
  scale_colour_manual(values = TPC, guide = "none") +
  labs(x = NULL, y = "Median mitochondrial %") +
  theme_classic(base_size = 7.5) +
  theme(axis.text = element_text(size = 6))

pC <- (pC1 | pC2)

pD <- ggplot(sm, aes(tp, 100 * mal_frac)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, fill = "grey93",
               colour = "grey40", linewidth = 0.3) +
  geom_line(aes(group = patient), colour = "grey55", linewidth = 0.22,
            alpha = 0.55) +
  geom_point(aes(colour = tp), size = 1.15, alpha = 0.9) +
  scale_colour_manual(values = TPC, guide = "none") +
  labs(x = NULL, y = "Malignant fraction of epithelial cells (%)",
       title = "D  Tumour content declines",
       subtitle = paste0("Each line is one patient.\n",
                         "Most fall from A to C, as expected\n",
                         "under neoadjuvant endocrine therapy.")) +
  theme_classic(base_size = 7.5) +
  theme(plot.title = element_text(size = 8, face = "bold"),
        plot.subtitle = element_text(size = 6.4, colour = "grey35"))

pTop <- ggplot() + theme_void() +
  annotate("text", x = 0, y = 0.98, hjust = 0, vjust = 1, size = 3.3,
           fontface = "bold",
           label = "Figure S3  Quality control and the cell-count cascade") +
  annotate("text", x = 0, y = 0.55, hjust = 0, vjust = 1, size = 2.2,
           colour = "grey30", lineheight = 1.05,
           label = paste0("NeoLetExe, 73 biopsies from 25 patients (24 with three ",
                          "timepoints). Cells were filtered at 300 genes,\n",
                          "500 counts and 20% mitochondrial reads; doublets were ",
                          "removed per biopsy before the epithelial gate.")) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))

p <- pTop / (pA | pB) / (pC | pD) + plot_layout(heights = c(0.22, 1, 1))

ggsave(file.path(FIG, "figS3_qc.png"), p, width = 7.2, height = 6.2,
       dpi = 300, bg = "white")
cat("\n-> figS3_qc.png\n")
