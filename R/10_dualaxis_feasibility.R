## ============================================================================
## 10_dualaxis_feasibility.R
##
## 课题 v2（评分分解）的三项前置判断，动笔前必须通过。
##
##   A0 臂随机化零假设：两臂的负相关是不是 LogNormalize 的成分约束造成的
##      （若不是，全部结论作废；这是最容易犯的错误）
##   A1 轴耦合在患者间是否一致（连续量，避免中位数二分的任意性）
##   B  四象限组成是否稳定（二项抽样噪声 vs 真实患者异质）
##   C  象限是否为跨患者共享的细胞状态（留一患者判别 + 随机基因集对照）
##   D  患者内 silhouette：象限结构是否超出患者身份
##   E  技术混淆（测序深度、线粒体、双细胞、解离应激、细胞周期）
##
## 输出：results/tables/fx_*.csv、results/figures/fig16_*.png
## ============================================================================

source("config/00_config.R")
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(presto)
  library(pROC); library(ggplot2); library(patchwork)
})
set.seed(WB_SEED)

msg("== 10. dual-axis feasibility ==")

SIG_CSV <- file.path("..", "rnaseq-analysis", "results", "geo", "signature_genes.csv")
if (!file.exists(SIG_CSV))
  stop("signature_genes.csv not found: ", SIG_CSV)

so_all <- readRDS(file.path(PROC_DIR, "06_epi_clustered.rds"))
so <- subset(so_all, cells = colnames(so_all)[so_all$ref_subtype == "ER+"])
rm(so_all); invisible(gc(FALSE))
msg("ER+ malignant epithelium: %d cells, %d patients", ncol(so), length(unique(so$sample)))

sig <- read.csv(SIG_CSV, stringsAsFactors = FALSE)
UP <- intersect(sig$gene[sig$direction == "Up"],   rownames(so))
DN <- intersect(sig$gene[sig$direction == "Down"], rownames(so))
msg("signature genes present: up %d/%d, down %d/%d",
    length(UP), sum(sig$direction == "Up"), length(DN), sum(sig$direction == "Down"))

expr <- GetAssayData(so, layer = "data")
cells <- colnames(so)
pat   <- as.character(so$sample)
pts   <- names(sort(table(pat), decreasing = TRUE))

## gene universe for marker selection / null: expressed in >=10% of ER+ cells
pct_in <- Matrix::rowMeans(expr > 0)
UNIV   <- rownames(expr)[pct_in >= 0.10]
msg("gene universe (pct >= 10%%): %d", length(UNIV))

arm_score <- function(g, cc = cells) {
  g <- intersect(g, rownames(expr))
  if (!length(g)) return(rep(NA_real_, length(cc)))
  as.numeric(Matrix::colMeans(expr[g, cc, drop = FALSE]))
}

u <- arm_score(UP); d <- arm_score(DN)
so$fx_up <- u; so$fx_dn <- d

## ---------------------------------------------------------------- A0 null
## LogNormalize 固定了每个细胞的总计数：一群基因占得多，另一群机械性地占得少，
## 因此任意两组表达量相近的基因都会呈现负相关。必须把观测相关放到这个零分布里看。
expr_mean <- Matrix::rowMeans(expr)
bins <- cut(rank(expr_mean[UNIV]), breaks = 10, labels = FALSE)
names(bins) <- UNIV

draw_arm <- function(n) {
  ## 按表达量十分位分层抽样，保证随机臂与真实臂的表达量分布匹配
  k <- table(bins[UP])                      # 用真实上调臂的分层结构做模板
  out <- unlist(lapply(names(k), function(b)
    sample(UNIV[bins == as.integer(b)], k[[b]], replace = FALSE)), use.names = FALSE)
  if (length(out) < n) out <- c(out, sample(setdiff(UNIV, out), n - length(out)))
  out[seq_len(n)]
}
dn_bins <- table(bins[DN])

draw_arm_dn <- function(n) {
  out <- unlist(lapply(names(dn_bins), function(b)
    sample(UNIV[bins == as.integer(b)], dn_bins[[b]], replace = FALSE)), use.names = FALSE)
  if (length(out) < n) out <- c(out, sample(setdiff(UNIV, out), n - length(out)))
  out[seq_len(n)]
}

N_NULL <- as.integer(Sys.getenv("WB_FX_NULL", unset = "200"))
set.seed(WB_SEED)
null_rho <- numeric(N_NULL)
for (i in seq_len(N_NULL)) {
  ur <- arm_score(draw_arm(length(UP)))
  dr <- arm_score(draw_arm_dn(length(DN)))
  null_rho[i] <- cor(ur, dr, method = "spearman")
}
obs_rho <- cor(u, d, method = "spearman")
msg("A0: observed rho(up,dn) = %.3f | random-arm null %.3f [%.3f, %.3f] | p = %.3f",
    obs_rho, mean(null_rho), quantile(null_rho, .025), quantile(null_rho, .975),
    (sum(null_rho <= obs_rho) + 1) / (N_NULL + 1))

write_csv_safe(data.frame(
  statistic = c("rho_observed", "null_mean", "null_q025", "null_q975", "p_empirical"),
  value = round(c(obs_rho, mean(null_rho), quantile(null_rho, .025),
                  quantile(null_rho, .975),
                  (sum(null_rho <= obs_rho) + 1) / (N_NULL + 1)), 4)),
  file.path(TAB_DIR, "fx_A0_arm_null.csv"))

## ------------------------------------------------------- A1 per-patient
boot_rho <- function(x, y, B = 500) {
  n <- length(x)
  v <- replicate(B, {
    i <- sample.int(n, n, replace = TRUE)
    if (sd(x[i]) == 0 || sd(y[i]) == 0) return(NA_real_)
    cor(x[i], y[i], method = "spearman")
  })
  quantile(v, c(.025, .975), na.rm = TRUE)
}

PROLIF <- c("MKI67","TOP2A","PCNA","CCNB1","CDK1","BUB1","AURKB","PLK1","CDC20",
            "NUSAP1","UBE2C","CCNA2","TYMS","RRM2","BIRC5","STMN1","TUBB","CENPF")
STRESS <- c("FOS","FOSB","JUN","JUNB","JUND","EGR1","ATF3","HSPA1A","HSPA1B",
            "HSP90AA1","DNAJB1","ZFP36","IER2","NR4A1","PPP1R15A","DUSP1","SOCS3",
            "KLF6","BTG2")
so$fx_prolif <- arm_score(PROLIF)
so$fx_stress <- arm_score(STRESS)

res_c <- lapply(pts, function(s) {
  i <- which(pat == s)
  ci <- boot_rho(u[i], d[i])
  ## partial rho given proliferation, to see whether the negative coupling is
  ## simply "proliferating cells cannot also be luminal"
  pr <- tryCatch({
    ru <- residuals(lm(rank(u[i]) ~ rank(so$fx_prolif[i])))
    rd <- residuals(lm(rank(d[i]) ~ rank(so$fx_prolif[i])))
    cor(ru, rd)
  }, error = function(e) NA_real_)
  data.frame(patient = s, n_cells = length(i),
             rho = round(cor(u[i], d[i], method = "spearman"), 3),
             lo = round(ci[1], 3), hi = round(ci[2], 3),
             rho_prolif_adj = round(pr, 3),
             mean_up = round(mean(u[i]), 4), mean_dn = round(mean(d[i]), 4),
             stringsAsFactors = FALSE)
})
res_c <- do.call(rbind, res_c)
res_c <- res_c[order(res_c$rho), ]
write_csv_safe(res_c, file.path(TAB_DIR, "fx_A1_patient_coupling.csv"))
print(res_c, row.names = FALSE)
msg("A1: negative coupling in %d/%d patients; pooled rho = %.3f",
    sum(res_c$rho < 0 & res_c$hi < 0), nrow(res_c),
    cor(u, d, method = "spearman"))

## ------------------------------------------------------------ B quadrants
## z 标准化后以 0 为界；同时给出患者内标准化（去掉患者均值差异）的版本，
## 因为若象限只是患者间差异的投影，患者内标准化会把它抹掉。
zu <- as.numeric(scale(u)); zd <- as.numeric(scale(d))
quad <- ifelse(zu > 0 & zd <= 0, "up-high/dn-low",
        ifelse(zu <= 0 & zd > 0, "up-low/dn-high",
        ifelse(zu > 0 & zd > 0, "both-high", "both-low")))
so$fx_quad <- quad

zu_w <- ave(u, pat, FUN = function(v) as.numeric(scale(v)))
zd_w <- ave(d, pat, FUN = function(v) as.numeric(scale(v)))
quad_w <- ifelse(zu_w > 0 & zd_w <= 0, "up-high/dn-low",
          ifelse(zu_w <= 0 & zd_w > 0, "up-low/dn-high",
          ifelse(zu_w > 0 & zd_w > 0, "both-high", "both-low")))
so$fx_quad_within <- quad_w

QL <- c("up-high/dn-low", "up-low/dn-high", "both-high", "both-low")
tab <- table(factor(pat, levels = pts), factor(quad, levels = QL))
frac <- prop.table(tab, 1)
n_pt <- rowSums(tab)

## binomial null for the spread of fractions across patients
p_bar <- colSums(tab) / sum(tab)
null_sd <- sapply(seq_along(QL), function(j) {
  v <- replicate(2000, {
    x <- rbinom(length(n_pt), n_pt, p_bar[j]) / n_pt
    sd(x)
  })
  mean(v)
})
obs_sd <- apply(frac, 2, sd)
disp <- data.frame(quadrant = QL,
                   overall_pct = round(100 * p_bar, 1),
                   obs_sd = round(obs_sd, 4),
                   binom_sd = round(null_sd, 4),
                   dispersion = round(obs_sd / null_sd, 2),
                   min_pct = round(100 * apply(frac, 2, min), 1),
                   max_pct = round(100 * apply(frac, 2, max), 1),
                   n_pat_ge5pct = apply(frac, 2, function(v) sum(v >= .05)),
                   stringsAsFactors = FALSE)
write_csv_safe(disp, file.path(TAB_DIR, "fx_B1_quadrant_dispersion.csv"))
print(disp, row.names = FALSE)

chisq_h <- suppressWarnings(chisq.test(tab))
msg("B: chi-square homogeneity X2 = %.1f, df = %d, phi = %.2f",
    chisq_h$statistic, chisq_h$parameter,
    as.numeric(chisq_h$statistic) / as.numeric(chisq_h$parameter))

## bootstrap reliability: how much of the between-patient spread is sampling noise
B_SUB <- as.integer(Sys.getenv("WB_FX_SUB", unset = "200"))
rel <- do.call(rbind, lapply(seq_along(pts), function(k) {
  idx <- which(pat == pts[k]); m <- min(300, length(idx))
  mm <- replicate(B_SUB, {
    ii <- sample(idx, m)
    as.numeric(prop.table(table(factor(quad[ii], levels = QL))))
  })
  data.frame(patient = pts[k], n_cells = length(idx),
             q = QL, frac = as.numeric(frac[k, ]),
             sub_sd = round(apply(mm, 1, sd), 4))
}))
write_csv_safe(rel, file.path(TAB_DIR, "fx_B2_quadrant_reliability.csv"))

## prop.table() keeps class "table": as.data.frame() on it returns a long
## Var1/Var2/Freq frame, so matrix dispatch is required here.
frac_df <- as.data.frame.matrix(frac, stringsAsFactors = FALSE)
names(frac_df) <- QL
frac_df$patient <- rownames(frac)
frac_df$n <- as.integer(n_pt)
write_csv_safe(frac_df, file.path(TAB_DIR, "fx_B0_quadrant_frac.csv"))

## patient-internal standardisation: does the structure survive removing the
## patient mean shift? If not, the quadrants are just patient identity in disguise.
tab_w <- table(factor(pat, levels = pts), factor(quad_w, levels = QL))
frac_w <- prop.table(tab_w, 1)
write_csv_safe(data.frame(patient = rownames(frac_w), round(as.data.frame.matrix(frac_w), 3),
                          check.names = FALSE),
               file.path(TAB_DIR, "fx_B3_quadrant_frac_within.csv"))
msg("B: within-patient split still gives all four quadrants in %d/%d patients",
    sum(apply(frac_w, 1, function(v) all(v > .02))), nrow(frac_w))

## ------------------------------------- C leave-one-patient-out recognition
## 训练：其余患者的象限 marker；测试：留出患者。随机基因集作对照。
MIN_TRAIN <- as.integer(Sys.getenv("WB_FX_MIN_CELL", unset = "60"))
TOP_N     <- as.integer(Sys.getenv("WB_FX_TOP", unset = "40"))

markers_in <- function(in_cells, out_cells, top = TOP_N) {
  if (length(in_cells) < MIN_TRAIN || length(out_cells) < MIN_TRAIN) return(character(0))
  X <- expr[, c(in_cells, out_cells), drop = FALSE]
  y <- c(rep("in", length(in_cells)), rep("out", length(out_cells)))
  a <- tryCatch(suppressWarnings(presto::wilcoxauc(X, y)), error = function(e) NULL)
  if (is.null(a)) return(character(0))
  a <- a[a$group == "in" & a$pct_in > 0.15 & a$auc > 0.55, ]
  if (!nrow(a)) return(character(0))
  a <- a[order(-a$auc), ]
  head(a$feature, top)
}

score_set <- function(cc, g) {
  g <- intersect(g, rownames(expr))
  if (!length(g)) return(rep(NA_real_, length(cc)))
  m <- as.matrix(expr[g, cc, drop = FALSE])
  s <- apply(m, 1, sd)
  z <- (m - rowMeans(m)) / pmax(s, 1e-6)
  as.numeric(colMeans(z))
}

auc_of <- function(lab, sc) {
  ok <- !is.na(sc)
  if (sum(ok) < 20 || length(unique(lab[ok])) < 2) return(NA_real_)
  r <- suppressMessages(pROC::roc(lab[ok], sc[ok], quiet = TRUE, direction = "<"))
  as.numeric(pROC::auc(r))
}

set.seed(WB_SEED)
lop <- list(); ctrl <- list()
for (k in seq_along(pts)) {
  test_s <- pts[k]
  test_i <- which(pat == test_s)
  train_i <- which(pat != test_s)
  if (length(test_i) < MIN_TRAIN) next

  mk <- lapply(QL, function(q)
    markers_in(train_i[quad[train_i] == q], train_i[quad[train_i] != q]))
  names(mk) <- QL

  for (q in QL) {
    lab <- quad[test_i] == q
    if (sum(lab) < 10 || sum(!lab) < 10) next
    a_own <- auc_of(lab, score_set(cells[test_i], mk[[q]]))

    ## wrong-quadrant control: markers of another quadrant applied here
    oth <- setdiff(QL, q)
    a_oth <- max(sapply(oth, function(q2)
      auc_of(lab, score_set(cells[test_i], mk[[q2]]))), na.rm = TRUE)

    ## random gene sets, matched for size
    rg <- replicate(50, auc_of(lab, score_set(cells[test_i],
             sample(UNIV, max(5, length(mk[[q]]))))))
    lop[[length(lop) + 1]] <- data.frame(
      patient = test_s, quadrant = q, n_in = sum(lab), n_out = sum(!lab),
      n_genes = length(mk[[q]]),
      auc_lopo = round(a_own, 3),
      auc_null_mean = round(mean(rg, na.rm = TRUE), 3),
      auc_null_q95 = round(quantile(rg, .95, na.rm = TRUE), 3),
      auc_lopo_minus_null = round(a_own - mean(rg, na.rm = TRUE), 3),
      auc_wrong_quadrant = round(a_oth, 3),
      stringsAsFactors = FALSE)
  }
  msg("  C: held out %s (%d cells)", test_s, length(test_i))
}
lop <- do.call(rbind, lop)
write_csv_safe(lop, file.path(TAB_DIR, "fx_C1_lopo_recognition.csv"))

lop_sum <- do.call(rbind, lapply(QL, function(q) {
  s <- lop[lop$quadrant == q, ]
  data.frame(quadrant = q, n_patients = nrow(s),
             auc_median = round(median(s$auc_lopo, na.rm = TRUE), 3),
             auc_min = round(min(s$auc_lopo, na.rm = TRUE), 3),
             auc_null_median = round(median(s$auc_null_mean, na.rm = TRUE), 3),
             delta_median = round(median(s$auc_lopo_minus_null, na.rm = TRUE), 3),
             n_pat_auc_ge65 = sum(s$auc_lopo >= .65, na.rm = TRUE),
             stringsAsFactors = FALSE)
}))
write_csv_safe(lop_sum, file.path(TAB_DIR, "fx_C2_lopo_summary.csv"))
print(lop_sum, row.names = FALSE)

## ------------------------------------------- D within-patient silhouette
## 若象限是真实细胞状态，同一患者内同象限细胞应比跨象限细胞更相似。
## 用 PCA 前 15 维的欧氏距离计算，并以打乱象限标签作零分布。
red <- if ("pca" %in% names(so@reductions)) "pca" else names(so@reductions)[1]
emb <- Embeddings(so, red)[, 1:min(15, ncol(Embeddings(so, red))), drop = FALSE]
msg("D: using %s (%d dims)", red, ncol(emb))

sil_calc <- function(D, lab) {
  n <- length(lab); s <- numeric(n)
  for (i in seq_len(n)) {
    same <- which(lab == lab[i]); same <- setdiff(same, i)
    oth <- setdiff(unique(lab), lab[i])
    if (!length(same) || !length(oth)) { s[i] <- NA; next }
    a <- mean(D[i, same])
    b <- min(sapply(oth, function(g) mean(D[i, lab == g])))
    s[i] <- (b - a) / max(a, b)
  }
  mean(s, na.rm = TRUE)
}

CAP <- as.integer(Sys.getenv("WB_FX_CAP", unset = "700"))
sil <- do.call(rbind, lapply(pts, function(s) {
  i <- which(pat == s)
  if (length(i) < 100) return(NULL)
  if (length(i) > CAP) i <- sample(i, CAP)
  D <- as.matrix(dist(emb[i, , drop = FALSE]))
  lab <- quad[i]
  q_all <- length(unique(lab)) == 4 && min(table(lab)) >= 10
  obs <- if (q_all) sil_calc(D, lab) else NA_real_
  ## sample-size matched null: same quadrant sizes, randomly assigned
  nullv <- if (q_all) replicate(100, sil_calc(D, sample(lab))) else NA_real_
  ## reference: patient-internal split, in case the pooled split is patient-driven
  labw <- quad_w[i]
  obsw <- if (length(unique(labw)) == 4 && min(table(labw)) >= 10) sil_calc(D, labw) else NA_real_
  data.frame(patient = s, n_used = length(i), silhouette = round(obs, 4),
             null_mean = round(mean(nullv, na.rm = TRUE), 4),
             null_q95 = round(quantile(nullv, .95, na.rm = TRUE), 4),
             silhouette_within = round(obsw, 4), stringsAsFactors = FALSE)
}))
write_csv_safe(sil, file.path(TAB_DIR, "fx_D1_silhouette.csv"))
print(sil, row.names = FALSE)
msg("D: median silhouette %.3f (null %.3f), %d/%d patients above null 95%%",
    median(sil$silhouette, na.rm = TRUE), median(sil$null_mean, na.rm = TRUE),
    sum(sil$silhouette > sil$null_q95, na.rm = TRUE), sum(!is.na(sil$silhouette)))

## ------------------------------------------------------- E confounders
CC <- intersect(cc.genes.updated.2019$s.genes, rownames(so))
GG <- intersect(cc.genes.updated.2019$g2m.genes, rownames(so))
so$fx_cc <- arm_score(c(CC, GG))

dbl <- if ("scDblFinder.score" %in% colnames(so@meta.data)) so$scDblFinder.score else NA_real_

e_tab <- do.call(rbind, lapply(seq_along(pts), function(k) {
  i <- which(pat == pts[k])
  data.frame(
    patient = pts[k], n_cells = length(i),
    med_nFeature = median(so$nFeature_RNA[i]),
    med_mito = round(median(so$percent_mito[i]), 2),
    med_doublet = if (all(is.na(dbl))) NA_real_ else round(median(dbl[i], na.rm = TRUE), 3),
    mean_stress = round(mean(so$fx_stress[i]), 4),
    mean_prolif = round(mean(so$fx_prolif[i]), 4),
    mean_cc = round(mean(so$fx_cc[i]), 4),
    pct_up_high = round(mean(zu[i] > 0), 3),
    pct_dn_high = round(mean(zd[i] > 0), 3),
    stringsAsFactors = FALSE)
}))
write_csv_safe(e_tab, file.path(TAB_DIR, "fx_E1_confounders.csv"))
print(e_tab, row.names = FALSE)

e_cor <- data.frame(
  contrast = c("pct_up_high", "pct_dn_high"),
  vs = rep(c("n_cells", "med_nFeature", "med_doublet", "mean_stress", "mean_prolif"), each = 2),
  rho = NA_real_, p = NA_real_, stringsAsFactors = FALSE)
for (r in seq_len(nrow(e_cor))) {
  x <- e_tab[[e_cor$contrast[r]]]; y <- e_tab[[e_cor$vs[r]]]
  if (all(is.na(y))) next
  ct <- suppressWarnings(cor.test(x, y, method = "spearman"))
  e_cor$rho[r] <- round(unname(ct$estimate), 3)
  e_cor$p[r]   <- signif(ct$p.value, 3)
}
write_csv_safe(e_cor, file.path(TAB_DIR, "fx_E2_confounder_cor.csv"))
print(e_cor, row.names = FALSE)

## within-patient technical dependence: if the up arm tracks library complexity
## inside a patient, the arm is partly a depth statistic, not a cell state.
e_within <- do.call(rbind, lapply(pts, function(s) {
  i <- which(pat == s)
  if (length(i) < 100) return(NULL)
  data.frame(patient = s, n_cells = length(i),
             rho_up_feature = round(cor(u[i], so$nFeature_RNA[i], method = "spearman"), 3),
             rho_up_count   = round(cor(u[i], so$nCount_RNA[i],   method = "spearman"), 3),
             rho_dn_feature = round(cor(d[i], so$nFeature_RNA[i], method = "spearman"), 3),
             rho_dn_count   = round(cor(d[i], so$nCount_RNA[i],   method = "spearman"), 3),
             stringsAsFactors = FALSE)
}))
write_csv_safe(e_within, file.path(TAB_DIR, "fx_E3_within_patient_technical.csv"))
print(e_within, row.names = FALSE)
msg("E: median within-patient rho(up, nFeature) = %.3f | rho(dn, nFeature) = %.3f",
    median(e_within$rho_up_feature), median(e_within$rho_dn_feature))

## ------------------------------------------------------------- figures
pal4 <- c("up-high/dn-low" = "#C0504D", "up-low/dn-high" = "#4A6FA5",
          "both-high" = "#9BBB59", "both-low" = "#8C8C8C")

pA <- ggplot(res_c, aes(reorder(patient, rho), rho)) +
  geom_hline(yintercept = mean(null_rho), linetype = 2, colour = "grey45", linewidth = .3) +
  geom_hline(yintercept = 0, colour = "grey75", linewidth = .25) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = .18, linewidth = .35, colour = "grey55") +
  geom_point(size = 1.9, colour = "#4A6FA5") +
  labs(x = NULL, y = expression("Spearman " * rho * " (up vs down arm)"),
       title = "A  Arm coupling per patient",
       subtitle = sprintf("pooled %.2f | random-arm null %.2f", obs_rho, mean(null_rho))) +
  theme_bw(base_size = 8) +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 6),
        plot.subtitle = element_text(size = 6.5, colour = "grey35"))

sdf <- frac_df
sdf$patient <- factor(sdf$patient, levels = pts)
sl <- reshape(sdf, varying = QL, v.names = "frac", timevar = "q",
              times = QL, direction = "long")
sl$q <- factor(sl$q, levels = QL)
sl$lab <- factor(paste0(sl$patient, " (", sl$n, ")"),
                 levels = paste0(pts, " (", n_pt, ")"))
pB <- ggplot(sl, aes(lab, frac, fill = q)) +
  geom_col(width = .82) +
  scale_fill_manual(values = pal4) +
  scale_y_continuous(expand = c(0, 0), labels = function(z) sprintf("%.0f%%", 100 * z)) +
  labs(x = NULL, y = "fraction of cells", fill = NULL,
       title = "B  Quadrant composition per patient") +
  theme_bw(base_size = 8) +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 6),
        legend.position = "top", legend.key.size = unit(.3, "cm"))

ld <- lop
ld$key <- paste(ld$patient, ld$quadrant)
pC <- ggplot(ld, aes(auc_lopo, reorder(key, auc_lopo), colour = quadrant)) +
  geom_segment(aes(x = auc_null_mean, xend = auc_lopo, yend = key),
               colour = "grey80", linewidth = .3) +
  geom_point(aes(x = auc_null_mean), shape = 124, size = 2.6, colour = "grey45") +
  geom_point(size = 1.7) +
  geom_vline(xintercept = .5, linetype = 2, colour = "grey60", linewidth = .25) +
  scale_colour_manual(values = pal4, guide = "none") +
  labs(x = "AUC in held-out patient (| = random-gene null)",
       y = NULL, title = "C  Leave-one-patient-out recognition") +
  theme_bw(base_size = 8) +
  theme(panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 5))

sd2 <- rbind(
  data.frame(patient = sil$patient, value = sil$silhouette, type = "observed"),
  data.frame(patient = sil$patient, value = sil$null_q95, type = "null 95%"),
  data.frame(patient = sil$patient, value = sil$silhouette_within, type = "within-patient split"))
pD <- ggplot(sd2[!is.na(sd2$value), ], aes(value, reorder(patient, value), colour = type, shape = type)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey60", linewidth = .25) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = c(observed = "#C0504D", `null 95%` = "grey55",
                                 `within-patient split` = "#4A6FA5"), name = NULL) +
  scale_shape_manual(values = c(observed = 16, `null 95%` = 4, `within-patient split` = 17), name = NULL) +
  labs(x = "mean silhouette (PCA space)", y = NULL, title = "D  Within-patient quadrant structure") +
  theme_bw(base_size = 8) +
  theme(panel.grid.major.y = element_blank(), legend.position = "top",
        legend.key.size = unit(.3, "cm"), legend.text = element_text(size = 6))

fig <- (pA | pB) / (pC | pD)
ggsave(file.path(FIG_DIR, "fig16_feasibility.png"), fig, width = 11, height = 7.6, dpi = 250)
ggsave(file.path(FIG_DIR, "fig16_feasibility.pdf"), fig, width = 11, height = 7.6)

saveRDS(list(u = u, d = d, quad = quad, quad_w = quad_w, pat = pat,
             up = UP, dn = DN, null_rho = null_rho,
             coupling = res_c, dispersion = disp, lopo = lop,
             lopo_sum = lop_sum, silhouette = sil),
        file.path(PROC_DIR, "10_feasibility.rds"))

msg("done: fig16 + %d tables", length(list.files(TAB_DIR, pattern = "^fx_")))
