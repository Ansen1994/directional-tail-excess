## =============================================================================
## 09_epi_nmf_per_sample.R -- programs that recur across patients
##
## Why this replaces 07: a single NMF fitted on the pooled malignant epithelium
## learns basis vectors that are dominated by between-patient differences
## (within-subtype eta^2 of program share was 0.50-0.89). A program that only
## exists in one patient cannot support a cohort-level claim, and pooling makes
## that invisible.
##
## The fix is the cNMF idea (Kotliar et al. 2019): fit each patient separately,
## align the resulting programs across patients, then keep only the ones that
## recur in most patients. Recurrence becomes a filter rather than an assumption.
##
## A cluster is a *clean* consensus program only if no patient contributes more
## than one member to it. A cluster that absorbs several programs from the same
## patient is a KMeans catch-all, not a matched program, and is reported as
## mixed rather than counted as recurrent.
##
## Scope: ER+ only. HER2+ has 3 usable patients and TNBC is left for a later
## run; a 5/7 threshold is meaningless with fewer than ~5 patients.
## =============================================================================

source("config/00_config.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(RcppML)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

set.seed(WB_SEED)
t0 <- Sys.time()

so <- readRDS(file.path(PROC_DIR, "06_epi_clustered.rds"))
nm <- readRDS(file.path(PROC_DIR, "07_epi_nmf.rds"))
feat <- nm$genes                       # the 788 genes of the global fit
rm(nm); invisible(gc(FALSE))

SUBTYPE <- "ER+"
n_by_smp <- table(so$sample[so$ref_subtype == SUBTYPE])
patients <- sort(names(n_by_smp)[n_by_smp >= NMF_PS_MIN_CELL])
msg("%s patients with >= %d malignant cells: %d", SUBTYPE, NMF_PS_MIN_CELL, length(patients))
msg("cells: %s", paste(sprintf("%s=%d", patients, n_by_smp[patients]), collapse = " "))

cells <- colnames(so)[so$ref_subtype == SUBTYPE & so$sample %in% patients]
V <- GetAssayData(so, layer = "data")[feat, cells, drop = FALSE]
V <- as(as(V, "CsparseMatrix"), "dgCMatrix")
vmax <- max(V)
V <- V / vmax                          # one scale for all patients, so the
invisible(gc(FALSE))                   # projection below is on the same scale
msg("matrix: %d genes x %d cells (max %.2f)", nrow(V), ncol(V), vmax)

## ------------------------------------------------------- 1. per-patient fits
.ari <- function(tb) {
  n <- sum(tb); a <- rowSums(tb); b <- colSums(tb)
  sij <- sum(choose(tb, 2))
  sa <- sum(choose(a, 2)); sb <- sum(choose(b, 2))
  exp <- sa * sb / choose(n, 2)
  mx <- (sa + sb) / 2
  if (mx == exp) return(1)
  (sij - exp) / (mx - exp)
}

RANKS <- NMF_RANKS
ps <- list()
for (s in patients) {
  idx <- which(so$sample[cells] == s)
  v <- V[, idx, drop = FALSE]
  ps[[s]] <- list()
  for (k in RANKS) {
    runs <- lapply(seq_len(NMF_PS_NRUN), function(i)
      nmf(v, k = k, seed = WB_SEED + i, tol = NMF_TOL, maxit = NMF_MAXIT, verbose = FALSE))
    labs <- lapply(runs, function(m) max.col(t(m$h), ties.method = "first"))
    aris <- numeric(0)
    for (a in seq_len(NMF_PS_NRUN - 1)) for (b in (a + 1):NMF_PS_NRUN)
      aris <- c(aris, .ari(table(labs[[a]], labs[[b]])))
    ari_mat <- matrix(NA_real_, NMF_PS_NRUN, NMF_PS_NRUN); j <- 1
    for (a in seq_len(NMF_PS_NRUN - 1)) for (b in (a + 1):NMF_PS_NRUN) {
      ari_mat[a, b] <- ari_mat[b, a] <- aris[j]; j <- j + 1
    }
    medoid <- which.max(rowMeans(ari_mat, na.rm = TRUE))
    W <- runs[[medoid]]$w
    rownames(W) <- feat
    colnames(W) <- paste0(s, "_P", seq_len(ncol(W)))
    W <- sweep(W, 2, pmax(colSums(W), 1e-12), "/")   # L1 so cosine is scale-free
    ps[[s]][[as.character(k)]] <- list(W = W, stab = mean(aris), medoid = medoid)
  }
  msg("  %s (%d cells): stability %s", s, length(idx),
      paste(sprintf("k%d=%.2f", RANKS,
                    sapply(RANKS, function(k) ps[[s]][[as.character(k)]]$stab)), collapse = " "))
}
saveRDS(ps, file.path(PROC_DIR, "09_ps_fits.rds"))

stab_tab <- expand.grid(patient = patients, k = RANKS, stringsAsFactors = FALSE)
stab_tab$stability <- mapply(function(p, k) ps[[p]][[as.character(k)]]$stab,
                             stab_tab$patient, stab_tab$k)
stab_tab$n_cells <- n_by_smp[stab_tab$patient]
write_csv_safe(stab_tab, file.path(TAB_DIR, "ps_nmf_stability.csv"))
print(stab_tab, row.names = FALSE)

## ------------------------------------------- 2. align programs within each k
cosine <- function(a, b) sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2)))

match_at_k <- function(k) {
  X <- do.call(cbind, lapply(patients, function(s) ps[[s]][[as.character(k)]]$W))
  Xn <- sweep(X, 2, sqrt(colSums(X^2)), "/")
  km <- kmeans(t(Xn), centers = k, nstart = 200, iter.max = 300)
  cl <- km$cluster

  ## members stays a list of index vectors; putting it into a data.frame via
  ## rbind would split it into V1..Vk columns of uneven length.
  members <- lapply(seq_len(k), function(i) which(cl == i))
  rec <- do.call(rbind, lapply(seq_len(k), function(i) {
    cols <- members[[i]]
    if (!length(cols)) return(NULL)
    smp <- sub("_P[0-9]+$", "", colnames(X)[cols])
    pw <- if (length(cols) > 1) {
      cmb <- combn(cols, 2)
      mean(apply(cmb, 2, function(pr) cosine(Xn[, pr[1]], Xn[, pr[2]])))
    } else NA_real_
    per_smp <- table(smp)
    ## centroid of the cluster, for gene identity
    ctr <- rowMeans(X[, cols, drop = FALSE])
    data.frame(cluster = i, n_patients = length(unique(smp)), n_members = length(cols),
               max_per_patient = as.integer(max(per_smp)),
               mean_cosine = pw, patients = paste(sort(unique(smp)), collapse = ","),
               k = k, stringsAsFactors = FALSE)
  }))
  list(rec = rec, X = X, members = members)
}

matched <- lapply(RANKS, match_at_k)
names(matched) <- as.character(RANKS)

match_tab <- do.call(rbind, lapply(matched, function(x) {
  d <- x$rec[, c("k", "cluster", "n_patients", "n_members", "max_per_patient",
                 "mean_cosine", "patients")]
  d$mean_cosine <- round(d$mean_cosine, 3)
  d$clean <- d$max_per_patient == 1
  d
}))
write_csv_safe(match_tab, file.path(TAB_DIR, "ps_nmf_matching.csv"))

## top genes per consensus cluster, so each cluster can be named and judged
cg <- do.call(rbind, lapply(names(matched), function(kk) {
  x <- matched[[kk]]
  do.call(rbind, lapply(x$rec$cluster, function(i) {
    cols <- x$members[[i]]
    ctr <- rowMeans(x$X[, cols, drop = FALSE])
    g <- names(sort(ctr, decreasing = TRUE))[1:15]
    data.frame(k = as.integer(kk), cluster = i,
               top_genes = paste(g, collapse = ","), stringsAsFactors = FALSE)
  }))
}))
write_csv_safe(cg, file.path(TAB_DIR, "ps_nmf_cluster_genes.csv"))

## -------------------------------------------- 3. recurrence vs cosine cutoff
## the cutoff is a judgement call, so report how the count moves with it
cut_grid <- c(0.50, 0.60, 0.65, 0.70, 0.75)
sel_long <- expand.grid(k = RANKS, cutoff = cut_grid, stringsAsFactors = FALSE)
sel_long$n <- mapply(function(kk, th) {
  r <- match_tab[match_tab$k == kk, ]
  sum(r$clean & r$n_patients >= NMF_PS_MIN_SMP &
        !is.na(r$mean_cosine) & r$mean_cosine >= th)
}, sel_long$k, sel_long$cutoff)
write_csv_safe(sel_long, file.path(TAB_DIR, "ps_nmf_cutoff_sensitivity.csv"))

sel <- data.frame(
  k = RANKS,
  n_clean_5of7 = mapply(function(kk) {
    r <- match_tab[match_tab$k == kk, ]
    sum(r$clean & r$n_patients >= NMF_PS_MIN_SMP &
          !is.na(r$mean_cosine) & r$mean_cosine >= NMF_PS_MIN_COS)
  }, RANKS),
  mean_cos_all = round(sapply(RANKS, function(kk)
    mean(match_tab$mean_cosine[match_tab$k == kk], na.rm = TRUE)), 3))
write_csv_safe(sel, file.path(TAB_DIR, "ps_nmf_k_selection.csv"))
print(sel, row.names = FALSE)

best_k <- sel$k[order(-sel$n_clean_5of7, sel$k)][1]
msg("k=%d has %d clean recurrent programs; all k are reported below rather than",
    best_k, sel$n_clean_5of7[sel$k == best_k])
msg("picking one, because different k keep different programs (stress at k=4, luminal at k=6)")

## --------------------------------- 4. evaluate every candidate, not one k
## Choosing a single k silently decides which biology survives: at k=4 the only
## clean recurrent program is dissociation stress, at k=6 it is the luminal
## secretory program. Both are reported with their credibility tags.
cand <- match_tab[match_tab$clean & match_tab$n_patients >= NMF_PS_MIN_SMP, ]
msg("%d clean candidates covering >= %d/%d patients", nrow(cand),
    NMF_PS_MIN_SMP, length(patients))

cand_ctr <- lapply(seq_len(nrow(cand)), function(i) {
  x <- matched[[as.character(cand$k[i])]]
  cols <- x$members[[cand$cluster[i]]]
  ctr <- rowMeans(x$X[, cols, drop = FALSE])
  ctr / sum(ctr)
})
C <- do.call(cbind, cand_ctr)
rownames(C) <- feat
colnames(C) <- sprintf("k%dc%d", cand$k, cand$cluster)

## collapse candidates that are the same program reached at two different k
if (ncol(C) > 1) {
  Cn <- sweep(C, 2, sqrt(colSums(C^2)), "/")
  cs <- crossprod(Cn)
  keep <- rep(TRUE, ncol(C))
  for (i in seq_len(ncol(C))) {
    if (!keep[i]) next
    for (j in seq_len(ncol(C))) {
      if (j <= i || !keep[j]) next
      if (cs[i, j] > 0.85) {
        ## keep whichever covers more patients, then the higher cosine
        better <- if (cand$n_patients[j] != cand$n_patients[i])
          cand$n_patients[j] > cand$n_patients[i]
        else cand$mean_cosine[j] > cand$mean_cosine[i]
        if (better) keep[i] <- FALSE else keep[j] <- FALSE
      }
    }
  }
  msg("deduplicated %d candidates to %d distinct programs (cosine > 0.85)", ncol(C), sum(keep))
  C <- C[, keep, drop = FALSE]; cand <- cand[keep, , drop = FALSE]
}

## project every ER+ cell onto each candidate program separately
STRESS <- c("FOS", "FOSB", "JUN", "JUNB", "JUND", "EGR1", "ATF3", "HSPA1A",
            "HSPA1B", "HSP90AA1", "DNAJB1", "ZFP36", "IER2", "NR4A1",
            "PPP1R15A", "DUSP1", "SOCS3", "KLF6", "BTG2")
gs <- intersect(STRESS, rownames(so))
ms <- as.matrix(GetAssayData(so, layer = "data")[gs, cells, drop = FALSE])
zs <- (ms - rowMeans(ms)) / sqrt(rowSums((ms - rowMeans(ms))^2) / (ncol(ms) - 1))
stress <- colMeans(zs)
eta2 <- function(y, g) { ss <- summary(aov(y ~ g))[[1]][["Sum Sq"]]; ss[1] / sum(ss) }

## The data layer is already log-normalised, so the projection weight IS the
## activity; dividing by nCount again would manufacture a negative correlation
## with sequencing depth and mislabel every program as technical.
abs_w <- sapply(seq_len(ncol(C)), function(i) {
  as.numeric(project(V, w = C[, i, drop = FALSE]))
})
dimnames(abs_w) <- list(cells, colnames(C))

tag <- data.frame(
  program = colnames(C),
  k = cand$k, cluster = cand$cluster,
  n_patients = cand$n_patients, mean_cosine = round(cand$mean_cosine, 3),
  mean_weight = round(colMeans(abs_w), 3),
  rho_stress = round(sapply(seq_len(ncol(abs_w)), function(i)
    cor(abs_w[, i], stress, method = "spearman")), 3),
  rho_nfeature = round(sapply(seq_len(ncol(abs_w)), function(i)
    cor(abs_w[, i], so$nFeature_RNA[cells], method = "spearman")), 3),
  eta2_patient = round(sapply(seq_len(ncol(abs_w)), function(i)
    eta2(abs_w[, i], so$sample[cells])), 3),
  stringsAsFactors = FALSE)
tag$identity <- ifelse(tag$rho_nfeature < -0.5, "technical (complexity-driven)",
                ifelse(tag$rho_stress > 0.6, "dissociation stress", "biological"))
write_csv_safe(tag, file.path(TAB_DIR, "ps_nmf_candidates.csv"))
print(tag, row.names = FALSE)

## per-patient mean: a recurrent program must be present in each patient it was
## called recurrent for, otherwise the alignment matched noise
act_smp <- sapply(patients, function(s)
  colMeans(abs_w[so$sample[cells] == s, , drop = FALSE]))
act_smp <- matrix(act_smp, nrow = ncol(abs_w),
                  dimnames = list(colnames(abs_w), patients))
smp_out <- data.frame(program = rownames(act_smp), round(act_smp, 3), check.names = FALSE)
write_csv_safe(smp_out, file.path(TAB_DIR, "ps_nmf_activity_by_patient.csv"))
print(smp_out, row.names = FALSE)

## top genes per candidate
pg <- lapply(colnames(C), function(p) {
  g <- names(sort(C[, p], decreasing = TRUE))[1:NMF_TOP_GENES]
  data.frame(program = p, rank = seq_along(g), gene = g,
             weight = round(C[g, p], 5), stringsAsFactors = FALSE)
})
pg <- do.call(rbind, pg)
write_csv_safe(pg, file.path(TAB_DIR, "ps_nmf_genes.csv"))
for (p in colnames(C))
  msg("  %-8s [%s] %s", p, tag$identity[tag$program == p],
      paste(pg$gene[pg$program == p][1:15], collapse = ", "))

## ------------------------------------------------------------------ figures
## every cluster from every k, so the reader sees the whole field rather than
## only the survivors; the cutoff only moves the count between 0 and 1 anyway
mt <- match_tab
mt$clean_lab <- ifelse(mt$clean, "clean", "mixed (same patient twice)")
## no labels here on purpose: k6c1 and k8c5 sit on top of each other at 5/0.74,
## and the right-hand panel already names them
p1 <- ggplot(mt, aes(mean_cosine, n_patients, colour = clean_lab)) +
  geom_point(size = 1.7, alpha = .85) +
  geom_hline(yintercept = NMF_PS_MIN_SMP, linetype = 2, colour = "grey50", linewidth = .3) +
  geom_vline(xintercept = NMF_PS_MIN_COS, linetype = 2, colour = "grey50", linewidth = .3) +
  scale_colour_manual(values = c(clean = "#4A6FA5", `mixed (same patient twice)` = "#B9C4CE")) +
  scale_y_continuous(limits = c(0, length(patients) + 1.2)) +
  labs(x = "mean cosine within cluster", y = "patients covered", colour = NULL,
       title = "Every consensus cluster, all k") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        legend.key.size = unit(.3, "cm"))

## which k keeps which program: the point of the figure
kp <- data.frame(k = cand$k, n_patients = cand$n_patients,
                 label = colnames(C), identity = tag$identity)
p2 <- ggplot(kp, aes(factor(k), n_patients, colour = identity)) +
  geom_hline(yintercept = NMF_PS_MIN_SMP, linetype = 2, colour = "grey50", linewidth = .3) +
  geom_point(size = 2.4) +
  geom_text(aes(label = label), size = 2.1, vjust = -1, show.legend = FALSE) +
  scale_y_continuous(limits = c(0, length(patients) + 1)) +
  scale_colour_manual(values = c(`technical (complexity-driven)` = "#C0504D",
                                 `dissociation stress` = "#8064A2",
                                 biological = "#4A6FA5")) +
  labs(x = "k", y = "patients covered", colour = NULL,
       title = "What survives at each k") +
  theme_bw(base_size = 9) + theme(panel.grid.major.x = element_blank())

rt <- tag[, c("program", "rho_stress", "rho_nfeature", "eta2_patient")]
rl <- rbind(
  data.frame(program = rt$program, metric = "rho nFeature", value = rt$rho_nfeature),
  data.frame(program = rt$program, metric = "rho stress", value = rt$rho_stress),
  data.frame(program = rt$program, metric = "eta2 patient", value = rt$eta2_patient))
p3 <- ggplot(rl, aes(program, value, fill = metric)) +
  geom_col(position = position_dodge(width = .75), width = .7) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = .3) +
  scale_fill_manual(values = c(`rho nFeature` = "#C0504D", `rho stress` = "#8064A2",
                               `eta2 patient` = "#4A6FA5")) +
  labs(x = NULL, y = "value", fill = NULL, title = "Credibility tags") +
  theme_bw(base_size = 9) +
  theme(panel.grid.major.x = element_blank(), legend.position = "top",
        legend.key.size = unit(.32, "cm"))

fig <- p1 | (p2 / p3)
ggsave(file.path(FIG_DIR, "fig14_ps_nmf_recurrence.png"), fig, width = 10.5, height = 5, dpi = 250)
ggsave(file.path(FIG_DIR, "fig14_ps_nmf_recurrence.pdf"), fig, width = 10.5, height = 5)

emb <- Embeddings(so, "umap")[cells, , drop = FALSE]
ul <- data.frame(
  umap_1   = rep(emb[, 1], times = ncol(abs_w)),
  umap_2   = rep(emb[, 2], times = ncol(abs_w)),
  program  = factor(rep(colnames(C), each = nrow(abs_w)),
                    levels = colnames(C)),
  weighted = as.vector(as.matrix(abs_w)))
## each program is rescaled to its own 95th percentile: raw weights differ by
## ~2x between programs, which would make one panel look uniformly pale
ul$rel <- ave(ul$weighted, ul$program, FUN = function(x) x / quantile(x, .95))
p4 <- ggplot(ul, aes(umap_1, umap_2, colour = rel)) +
  geom_point(size = .12) +
  facet_wrap(~ program, nrow = 1) +
  scale_colour_gradient(low = "#F2F5F7", high = "#4A6FA5", limits = c(0, 1.1),
                        oob = scales::squish) +
  labs(x = "UMAP 1", y = "UMAP 2", colour = "weight\n(rel. p95)",
       title = sprintf("Programs recurring in >=5/7 %s patients (n=%d cells)",
                       SUBTYPE, length(cells))) +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank(), legend.position = "right")
ggsave(file.path(FIG_DIR, "fig15_ps_nmf_umap.png"), p4,
       width = 2.5 * ncol(abs_w) + 1.1, height = 3.1, dpi = 250)
ggsave(file.path(FIG_DIR, "fig15_ps_nmf_umap.pdf"), p4,
       width = 2.5 * ncol(abs_w) + 1.1, height = 3.1)

saveRDS(list(W = C, activity = abs_w, tags = tag, matching = match_tab,
             genes = feat, patients = patients),
        file.path(PROC_DIR, "09_ps_nmf.rds"))
writeLines(capture.output(sessionInfo()), file.path(LOG_DIR, "sessionInfo_09_ps_nmf.txt"))
msg("done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))
