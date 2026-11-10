## =============================================================================
## 08_nmf_diagnostics.R -- can the NMF programs carry a cohort-level claim?
##
## A program is only usable across patients if its share varies more between
## cells within a patient than between patients. If a program is basically a
## patient label, any "subtype X is enriched for program Y" statement is really
## "subtype X happens to contain the two patients that carry program Y".
##
## Three confounds are tested here:
##   1. sample identity      -> eta-squared of share explained by sample
##   2. dissociation stress  -> immediate-early / heat-shock genes induced while
##                              the tissue is being digested (van den Brink 2017)
##   3. technical covariates -> nFeature (library complexity), percent_mito
##
## Writes: results/tables/nmf_diagnostics.csv
##         results/figures/fig13_nmf_diagnostics.png
## =============================================================================

source("config/00_config.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

nm <- readRDS(file.path(PROC_DIR, "07_epi_nmf.rds"))
share <- nm$usage                       # cells x programs, rows sum to 1
so <- readRDS(file.path(PROC_DIR, "06_epi_clustered.rds"))
stopifnot(identical(rownames(share), colnames(so)))

progs <- colnames(share)

## ------------------------------------------------------------ 1. variance split
## eta^2 = between-group sum of squares / total sum of squares
eta2 <- function(y, g) {
  ss <- summary(aov(y ~ g))[[1]][["Sum Sq"]]
  ss[1] / sum(ss)
}

## ---------------------------------------------------------- 2. stress signature
STRESS <- c("FOS", "FOSB", "JUN", "JUNB", "JUND", "EGR1", "ATF3", "HSPA1A",
            "HSPA1B", "HSP90AA1", "DNAJB1", "ZFP36", "IER2", "NR4A1",
            "PPP1R15A", "DUSP1", "SOCS3", "KLF6", "BTG2")
gs <- intersect(STRESS, rownames(so))
m <- as.matrix(GetAssayData(so, layer = "data")[gs, , drop = FALSE])
z <- (m - rowMeans(m)) / sqrt(rowSums((m - rowMeans(m))^2) / (ncol(m) - 1))
stress <- colMeans(z)
rm(m, z); invisible(gc(FALSE))
msg("dissociation-stress score from %d/%d genes; range %.2f..%.2f",
    length(gs), length(STRESS), min(stress), max(stress))

## ------------------------------------------------------------- 3. diagnostics
## Subtype is nested in sample (each patient contributes one subtype), so
## eta2_sample is mechanically >= eta2_subtype and comparing the two proves
## nothing. The informative quantity is the sample effect measured *inside* each
## subtype: if patients of the same subtype still disagree strongly, the program
## is a patient label rather than a recurrent cell state.
SUBS <- c("ER+", "HER2+", "TNBC")
w_sub <- as.numeric(table(factor(so$ref_subtype, levels = SUBS)))

eta2_within <- function(y) {
  v <- sapply(SUBS, function(s) {
    i <- so$ref_subtype == s
    if (length(unique(so$sample[i])) < 2) return(NA_real_)
    eta2(y[i], so$sample[i])
  })
  stats::weighted.mean(v, w_sub, na.rm = TRUE)
}

diag_tab <- data.frame(
  program        = progs,
  eta2_sample    = round(sapply(progs, function(p) eta2(share[, p], so$sample)), 3),
  eta2_subtype   = round(sapply(progs, function(p) eta2(share[, p], so$ref_subtype)), 3),
  eta2_smp_in_sub = round(sapply(progs, function(p) eta2_within(share[, p])), 3),
  rho_stress     = round(sapply(progs, function(p)
    cor(share[, p], stress, method = "spearman")), 3),
  rho_nfeature   = round(sapply(progs, function(p)
    cor(share[, p], so$nFeature_RNA, method = "spearman")), 3),
  rho_mito       = round(sapply(progs, function(p)
    cor(share[, p], so$percent_mito, method = "spearman")), 3),
  mean_share     = round(colMeans(share), 4),
  stringsAsFactors = FALSE)
write_csv_safe(diag_tab, file.path(TAB_DIR, "nmf_diagnostics.csv"))
print(diag_tab)
msg("dissociation-stress score itself: eta2(sample)=%.3f eta2(subtype)=%.3f eta2(smp in subtype)=%.3f",
    eta2(stress, so$sample), eta2(stress, so$ref_subtype), eta2_within(stress))

## --------------------------------------------------------------- 4. per sample
smp <- sort(unique(so$sample))
smp_sub <- sapply(split(so$ref_subtype, so$sample), function(x) names(which.max(table(x))))
smp_mean <- sapply(smp, function(s) colMeans(share[so$sample == s, , drop = FALSE]))
smp_n <- as.integer(table(so$sample)[smp])

ord <- order(factor(smp_sub[smp], levels = c("ER+", "HER2+", "TNBC")), -smp_mean["MP1", ])
smp <- smp[ord]; smp_sub <- smp_sub[ord]; smp_mean <- smp_mean[, ord]; smp_n <- smp_n[ord]

## long format: as.vector() fills column-wise, i.e. one sample's programs in a
## block, so `sample` must repeat each= and `program` must repeat times=.
sdf <- data.frame(
  sample  = factor(rep(smp, each = length(progs)), levels = smp),
  subtype = factor(rep(smp_sub, each = length(progs)), levels = c("ER+", "HER2+", "TNBC")),
  program = factor(rep(progs, times = length(smp)), levels = progs),
  share   = as.vector(smp_mean),
  n       = rep(smp_n, each = length(progs)))
stopifnot(all(abs(tapply(sdf$share, sdf$sample, sum) - 1) < 1e-8))

p1 <- ggplot(sdf, aes(sample, share, fill = program)) +
  geom_col(width = .85, position = "fill") +
  facet_grid(~ subtype, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c(MP1 = "#4A6FA5", MP2 = "#C0504D",
                               MP3 = "#9BBB59", MP4 = "#8064A2")) +
  scale_y_continuous(labels = function(z) sprintf("%.0f%%", 100 * z), expand = c(0, 0)) +
  labs(x = NULL, y = "mean share per sample", fill = NULL,
       title = "Program composition tracks the sample, not the subtype") +
  theme_bw(base_size = 8) +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, size = 6),
        legend.position = "top", legend.key.size = unit(.35, "cm"),
        strip.background = element_rect(fill = "grey95", colour = "grey80"),
        strip.text = element_text(size = 7))

## eta^2 comparison: the sample effect measured inside subtypes is the one that
## decides whether a program travels across patients
e2 <- rbind(
  data.frame(program = diag_tab$program, group = "subtype",
             eta2 = diag_tab$eta2_subtype),
  data.frame(program = diag_tab$program, group = "sample (within subtype)",
             eta2 = diag_tab$eta2_smp_in_sub))
e2$group <- factor(e2$group, levels = c("subtype", "sample (within subtype)"))
p2 <- ggplot(e2, aes(program, eta2, fill = group)) +
  geom_col(position = position_dodge(width = .75), width = .7) +
  geom_hline(yintercept = .5, linetype = 2, colour = "grey40", linewidth = .3) +
  scale_fill_manual(values = c(subtype = "#4A6FA5", `sample (within subtype)` = "#C0504D")) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(x = NULL, y = expression(eta^2 ~ "of program share"), fill = NULL,
       title = "Variance explained") +
  theme_bw(base_size = 8) +
  theme(panel.grid.major.x = element_blank(), legend.position = "top",
        legend.key.size = unit(.35, "cm"))

## stress decile profile
q <- cut(stress, breaks = unique(quantile(stress, 0:10 / 10)), include.lowest = TRUE,
         labels = FALSE)
qm <- aggregate(share, by = list(decile = q), FUN = mean)
ql <- data.frame(decile = rep(qm$decile, times = length(progs)),
                 program = factor(rep(progs, each = nrow(qm)), levels = progs),
                 share = as.vector(as.matrix(qm[, progs])))
p3 <- ggplot(ql, aes(decile, share, colour = program)) +
  geom_line(linewidth = .5) + geom_point(size = 1.2) +
  scale_colour_manual(values = c(MP1 = "#4A6FA5", MP2 = "#C0504D",
                                 MP3 = "#9BBB59", MP4 = "#8064A2")) +
  scale_x_continuous(breaks = 1:10) +
  labs(x = "dissociation-stress score (decile)", y = "mean share", colour = NULL,
       title = "Program share vs stress score") +
  theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(), legend.position = "top",
        legend.key.size = unit(.35, "cm"))

## technical correlation heatmap
ct <- rbind(
  data.frame(program = diag_tab$program, cov = "stress",  rho = diag_tab$rho_stress),
  data.frame(program = diag_tab$program, cov = "nFeature", rho = diag_tab$rho_nfeature),
  data.frame(program = diag_tab$program, cov = "percent_mito", rho = diag_tab$rho_mito))
p4 <- ggplot(ct, aes(cov, program, fill = rho)) +
  geom_tile(colour = "white", linewidth = .4) +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 2.4, colour = "grey10") +
  scale_fill_gradient2(low = "#4A6FA5", mid = "#F2F5F7", high = "#C0504D",
                       midpoint = 0, limits = c(-1, 1)) +
  scale_y_discrete(limits = rev(progs)) +
  labs(x = NULL, y = NULL, fill = "Spearman\nrho",
       title = "Correlation with technical covariates") +
  theme_bw(base_size = 8) +
  theme(panel.grid = element_blank())

fig <- (p1 | p3) / (p2 | p4) + plot_layout(heights = c(1.15, 1))
ggsave(file.path(FIG_DIR, "fig13_nmf_diagnostics.png"), fig, width = 11, height = 7, dpi = 250)
ggsave(file.path(FIG_DIR, "fig13_nmf_diagnostics.pdf"), fig, width = 11, height = 7)

## ------------------------------------------- 5. does a program survive per sample?
## for each program, fraction of samples where it is the dominant program in a
## meaningful share of that sample's cells (>= 20%)
dom <- progs[max.col(share, ties.method = "first")]
dom_tab <- table(so$sample, dom)
dom_frac <- sweep(dom_tab, 1, rowSums(dom_tab), "/")
carry <- data.frame(
  program   = progs,
  n_samples_ge20pct = sapply(progs, function(p) sum(dom_frac[, p] >= .20)),
  max_sample_share  = round(sapply(progs, function(p) max(dom_frac[, p])), 3),
  top_sample        = sapply(progs, function(p) rownames(dom_frac)[which.max(dom_frac[, p])]))
write_csv_safe(carry, file.path(TAB_DIR, "nmf_program_carry.csv"))
print(carry)

writeLines(capture.output(sessionInfo()), file.path(LOG_DIR, "sessionInfo_08_diag.txt"))
msg("done")
