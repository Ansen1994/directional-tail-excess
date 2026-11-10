## ============================================================================
## 11_bulk_dualaxis_check.R
##
## 前置判断 3：双轴结构能否在 bulk 中复现。
##
## 单细胞里测到的是细胞层面（轴耦合 rho、象限组成）；bulk 只有患者层面一个数。
## 这里要回答的是：两臂在患者层面的关系，bulk 与单细胞是否一致，
## 以及两臂是否在各队列中独立预后 —— 若只在 METABRIC 成立，分解就没有外部基础。
##
## 依赖 ../rnaseq-analysis 的数据与结果（只读，不改动对方任何文件）。
## ============================================================================

source("config/00_config.R")
suppressPackageStartupMessages({
  library(survival); library(ggplot2); library(patchwork)
})

msg("== 11. bulk dual-axis reproducibility ==")

BLK <- normalizePath("../rnaseq-analysis", mustWork = FALSE)
if (!dir.exists(BLK)) stop("bulk project not found: ", BLK)

SIG <- read.csv(file.path(BLK, "results", "geo", "signature_genes.csv"),
                stringsAsFactors = FALSE)
UP <- SIG$gene[SIG$direction == "Up"]
DN <- SIG$gene[SIG$direction == "Down"]
msg("signature: %d up / %d down", length(UP), length(DN))

arm_from_matrix <- function(expr, up, dn, what) {
  ## samples x genes matrix; z-score each gene across samples, then arm mean
  u <- intersect(up, colnames(expr)); d <- intersect(dn, colnames(expr))
  if (length(u) < 5 || length(d) < 5) {
    msg("  %s: too few genes (up %d, down %d), skipped", what, length(u), length(d))
    return(NULL)
  }
  z <- function(g) {
    m <- expr[, g, drop = FALSE]
    zz <- scale(m)
    zz[!is.finite(zz)] <- NA
    rowMeans(zz, na.rm = TRUE)
  }
  data.frame(sample = rownames(expr), RS_up = z(u), RS_dn = z(d),
             n_up = length(u), n_dn = length(d), stringsAsFactors = FALSE)
}

cox_both <- function(df, time, event, label) {
  df <- df[is.finite(df[[time]]) & df[[time]] > 0 & !is.na(df[[event]]), ]
  df$t <- df[[time]]; df$e <- as.integer(df[[event]])
  if (nrow(df) < 30 || sum(df$e) < 10) {
    msg("  %s: too few events (%d), cox skipped", label, sum(df$e))
    return(NULL)
  }
  df$u <- as.numeric(scale(df$RS_up)); df$d <- as.numeric(scale(df$RS_dn))
  f1 <- coxph(Surv(t, e) ~ u, data = df)
  f2 <- coxph(Surv(t, e) ~ d, data = df)
  f3 <- coxph(Surv(t, e) ~ u + d, data = df)
  s1 <- summary(f1); s2 <- summary(f2); s3 <- summary(f3)
  data.frame(cohort = label, n = nrow(df), events = sum(df$e),
             rho_arms = round(cor(df$RS_up, df$RS_dn, method = "spearman"), 3),
             HR_up_uni = round(s1$conf.int[1, 1], 3), p_up_uni = signif(s1$coefficients[1, 5], 3),
             HR_dn_uni = round(s2$conf.int[1, 1], 3), p_dn_uni = signif(s2$coefficients[1, 5], 3),
             HR_up_adj = round(s3$conf.int[1, 1], 3), p_up_adj = signif(s3$coefficients[1, 5], 3),
             HR_dn_adj = round(s3$conf.int[2, 1], 3), p_dn_adj = signif(s3$coefficients[2, 5], 3),
             stringsAsFactors = FALSE)
}

out <- list()

## ------------------------------------------------------------- METABRIC
mb <- read.csv(file.path(BLK, "results", "metabric", "resistance_score_rfs.csv"),
               stringsAsFactors = FALSE)
mb <- mb[!is.na(mb$rfs_m) & mb$rfs_m > 0, ]
out$METABRIC <- cox_both(mb, "rfs_m", "rfs_ev", "METABRIC (RFS)")
msg("METABRIC: n=%d, rho(arms)=%.3f", nrow(mb),
    cor(mb$RS_up, mb$RS_dn, method = "spearman"))

## ------------------------------------------------------- GEO cohorts
load_geo <- function(tag, expr_file, clin_file, er_col, er_pos, time_col, event_col,
                     time_scale, label) {
  E <- readRDS(file.path(BLK, "data", "geo", "processed", expr_file))
  C <- read.csv(file.path(BLK, "data", "geo", "processed", clin_file),
                stringsAsFactors = FALSE, check.names = FALSE)
  samples <- intersect(colnames(E), C$sample)
  if (!length(samples)) { msg("  %s: no overlap", label); return(NULL) }
  M <- t(E[, samples, drop = FALSE])                 # samples x genes
  a <- arm_from_matrix(M, UP, DN, label)
  C <- C[match(samples, C$sample), ]
  df <- cbind(a, C[, c(er_col, time_col, event_col), drop = FALSE])
  names(df)[(ncol(df) - 2):ncol(df)] <- c("er", "tt", "ev")
  df$er <- as.character(df$er)
  keep <- if (is.null(er_pos)) rep(TRUE, nrow(df)) else df$er %in% er_pos
  msg("  %s: %d samples, ER+ %d", label, nrow(df), sum(keep, na.rm = TRUE))
  df$tt <- suppressWarnings(as.numeric(df$tt)) * time_scale
  df$ev <- suppressWarnings(as.integer(df$ev))
  cox_both(df[which(keep), ], "tt", "ev", label)
}

out$GSE17705 <- load_geo("gse17705", "gse17705_expr.rds", "gse17705_clin.csv",
                         "estrogen receptor (er) status", c("ER+", "ER +", "Positive", "pos"),
                         "event time (years)", "distant relapse (1=dr, 0 censored)",
                         12, "GSE17705 (ER+, DMFS)")
out$GSE6532 <- load_geo("gse6532", "gse6532_expr.rds", "gse6532_clin.csv",
                        "er", c("1"), "t.rfs", "e.rfs", 1 / 30.44, "GSE6532 (ER+, RFS)")
out$GSE9195 <- load_geo("gse9195", "gse9195_expr.rds", "gse9195_clin.csv",
                        "er", c("1"), "t.rfs", "e.rfs", 1 / 30.44, "GSE9195 (ER+, RFS)")

res <- do.call(rbind, out[!sapply(out, is.null)])
write_csv_safe(res, file.path(TAB_DIR, "fx_F1_bulk_dualaxis.csv"))
print(res, row.names = FALSE)

## single-cell patient-level rho, for direct comparison with the bulk column
fx <- readRDS(file.path(PROC_DIR, "10_feasibility.rds"))
pc <- aggregate(cbind(fx_up, fx_dn) ~ sample, data = data.frame(
  sample = fx$pat, fx_up = fx$u, fx_dn = fx$d), FUN = mean)
rho_sc_patient <- cor(pc$fx_up, pc$fx_dn, method = "spearman")
rho_sc_cell    <- cor(fx$u, fx$d, method = "spearman")
msg("single cell: rho(arms) at cell level %.3f, at patient level (n=%d) %.3f",
    rho_sc_cell, nrow(pc), rho_sc_patient)

cmp <- data.frame(
  level = c("ER+ scRNA, cell level", "ER+ scRNA, patient level", "random-arm null (cell)",
            res$cohort),
  n = c(length(fx$u), nrow(pc), length(fx$null_rho), res$n),
  rho = round(c(rho_sc_cell, rho_sc_patient, mean(fx$null_rho), res$rho_arms), 3),
  stringsAsFactors = FALSE)
write_csv_safe(cmp, file.path(TAB_DIR, "fx_F2_arm_correlation_compare.csv"))
print(cmp, row.names = FALSE)

## ------------------------------------------------------------- figure
pl <- rbind(
  data.frame(cohort = res$cohort, arm = "up arm (proliferation)",
             HR = res$HR_up_adj, lo = NA, hi = NA, p = res$p_up_adj),
  data.frame(cohort = res$cohort, arm = "down arm (luminal)",
             HR = res$HR_dn_adj, lo = NA, hi = NA, p = res$p_dn_adj))

p1 <- ggplot(cmp, aes(reorder(level, rho), rho)) +
  geom_col(fill = "#4A6FA5", width = .62) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = .3) +
  coord_flip() +
  labs(x = NULL, y = expression("Spearman " * rho * " between arms"),
       title = "Arm coupling: cell level vs patient level vs bulk") +
  theme_bw(base_size = 8) +
  theme(panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 6.5))

p2 <- ggplot(pl, aes(HR, reorder(cohort, HR), colour = arm, shape = arm)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey60", linewidth = .3) +
  geom_point(size = 2, position = position_dodge(width = .5)) +
  scale_colour_manual(values = c(`up arm (proliferation)` = "#C0504D",
                                 `down arm (luminal)` = "#4A6FA5"), name = NULL) +
  scale_shape_manual(values = c(3, 16), name = NULL) +
  scale_x_log10() +
  labs(x = "adjusted HR per SD (both arms in one model)", y = NULL,
       title = "Independent prognostic value of each arm") +
  theme_bw(base_size = 8) +
  theme(panel.grid.major.y = element_blank(), legend.position = "top",
        legend.key.size = unit(.3, "cm"), legend.text = element_text(size = 6.5))

fig <- p1 | p2
ggsave(file.path(FIG_DIR, "fig17_bulk_dualaxis.png"), fig, width = 9.2, height = 3.4, dpi = 250)
ggsave(file.path(FIG_DIR, "fig17_bulk_dualaxis.pdf"), fig, width = 9.2, height = 3.4)

msg("done")
