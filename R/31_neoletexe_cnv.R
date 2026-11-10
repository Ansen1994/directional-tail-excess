#!/usr/bin/env Rscript
# ==============================================================================
# 31_neoletexe_cnv.R
#
# NeoLetExe 上皮细胞 -> 分箱 CNV 剖面 -> 恶性判定
#
# 为什么自己实现而不用现成包：
#   SCEVAN 的 GitHub 仓库已不可访问；inferCNV 未安装且依赖较重。
#   本项目在 nc2025 上已验证过「分箱 + 留一样本基线 + 平滑 + 剖面离散度」
#   这套口径（R/05_cnv.R），这里沿用并细化到分箱级，因为下游还要用它做亚克隆。
#
# 与 R/05 的两点差别：
#   1. 05 是臂级（整条染色体臂一个值），只能判恶/非恶；
#      这里用 100 基因一箱，箱内取细胞均值，剖面上有几十个点，可供聚类。
#   2. 基线用**留一样本**的上皮均值（LOSO）：某一样本的基线由其余 72 个样本
#      的上皮细胞给出。不用同一个样本自己的免疫细胞作基线——05 里实测过，
#      免疫基线带着谱系差异，AUC 只有 0.34。
#
# 输出：
#   data/proc/neoletexe/epi_combined.rds     全部上皮细胞 counts + meta（复用）
#   data/proc/neoletexe/cnv_profile.rds      分箱剖面（箱 x 细胞，已减基线、平滑）
#   results/tables/neoletexe_cnv_cells.csv   每细胞 CNV 分数与恶性判定
#   results/tables/neoletexe_cnv_bins.csv    箱坐标
# ==============================================================================

suppressPackageStartupMessages({
  library(Matrix); library(data.table)
})

RAW  <- "data/raw/neoletexe/mat"
PROC <- "data/proc/neoletexe"
SAMP <- file.path(PROC, "samples")
TAB  <- "results/tables"; FIG <- "results/figures"

NBIN_TARGET <- 100L     # 每箱基因数
SMOOTH_W    <- 11L      # 平滑窗口（箱）
CLIP        <- 0.35     # 残差截断，与 R/05 一致
MIN_CELLS   <- 30L      # 每样本进入 CNV 的最少上皮细胞

dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

## ---------------- 1. 汇总上皮细胞 ----------------
CB <- file.path(PROC, "epi_combined.rds")
if (file.exists(CB)) {
  ep <- readRDS(CB)
  cat(sprintf("复用已汇总上皮细胞：%d 个\n", ncol(ep$counts)))
} else {
  sb <- fread(file.path(RAW, "sample_blocks.csv"))[order(sample)]
  mats <- list(); metas <- list()
  for (s in sb$sample) {
    f <- file.path(SAMP, sprintf("%s.rds", s))
    if (!file.exists(f)) { cat("缺样本 ", s, "\n"); next }
    o <- readRDS(f)
    sel <- o$meta$celltype == "Epi"
    if (!any(sel)) next
    mats[[s]]  <- o$counts[, sel, drop = FALSE]
    metas[[s]] <- o$meta[sel]
    rm(o); invisible(gc(FALSE))
  }
  ## 统一基因顺序后按列拼接
  g <- rownames(mats[[1]])
  mats <- lapply(mats, function(m) m[g, , drop = FALSE])
  counts <- do.call(cbind, mats)
  meta <- rbindlist(metas)
  stopifnot(identical(colnames(counts), meta$cell))
  ep <- list(counts = counts, meta = meta)
  saveRDS(ep, CB, compress = FALSE)
  cat(sprintf("上皮细胞汇总：%d 个 / %d 样本\n", ncol(counts), uniqueN(meta$sample)))
}
cnt <- ep$counts; meta <- copy(ep$meta)

## ---------------- 2. 分箱 ----------------
gmap <- fread(file.path(PROC, "gene_position.csv"))
gmap <- gmap[gene %in% rownames(cnt)]
## 没有坐标的基因（CHRLOC 缺失）不能定位，直接剔除，否则分组里全是 NA
gmap <- gmap[!is.na(start) & is.finite(start)]
gmap <- gmap[!duplicated(gene)]
setorder(gmap, chr_num, start)
## 染色体内按坐标排序；缺失坐标的基因丢弃（无法定位就不能算 CNV）
keep_g <- gmap$gene
cnt <- cnt[keep_g, , drop = FALSE]

gmap[, tmp := ceiling(seq_len(.N) / NBIN_TARGET), by = chr_num]
gmap[, bin := .GRP, by = .(chr_num, tmp)]   # 箱不跨染色体边界
gmap[, tmp := NULL]
bin_info <- gmap[, .(n_gene = .N, chr = paste(unique(chr), collapse = "/"),
                     chr_num = as.integer(chr_num[1]),
                     start = as.numeric(min(start)),
                     end = as.numeric(max(start)),
                     genes = paste(head(gene, 3), collapse = ",")),
                 by = .(bin)]
setorder(bin_info, chr_num, start)
cat(sprintf("分箱：%d 个基因 -> %d 箱（每箱约 %d 基因）\n",
            nrow(gmap), nrow(bin_info), NBIN_TARGET))
fwrite(bin_info, file.path(TAB, "neoletexe_cnv_bins.csv"))

## 箱指示矩阵（基因 x 箱），用于快速求箱均值
bm <- sparseMatrix(i = seq_len(nrow(gmap)), j = gmap$bin, x = 1,
                   dims = c(nrow(gmap), nrow(bin_info)))
rownames(bm) <- gmap$gene

## ---------------- 3. 归一化 + 箱均值 ----------------
lib <- Matrix::colSums(cnt)
ln <- cnt
ln@x <- log1p(ln@x / rep.int(lib, diff(ln@p)) * 1e4)
binsum <- as.matrix(Matrix::t(bm) %*% ln)         # 箱 x 细胞
bincnt <- as.numeric(Matrix::colSums(bm))
B <- binsum / bincnt                               # 箱内基因均值
rm(binsum, ln, cnt); invisible(gc(FALSE))
cat(sprintf("箱均值矩阵：%d 箱 x %d 细胞\n", nrow(B), ncol(B)))

## ---------------- 4. 留一样本基线 ----------------
sidx  <- match(meta$sample, sort(unique(meta$sample)))
n_s   <- tabulate(sidx, nbins = max(sidx))
S     <- sapply(seq_len(max(sidx)), function(k) {
  if (n_s[k] == 0) return(rep(0, nrow(B)))
  rowSums(B[, sidx == k, drop = FALSE])
})
tot   <- rowSums(S); n_tot <- sum(n_s)
base  <- (tot - S[, sidx, drop = FALSE]) / (n_tot - n_s[sidx])
R     <- B - base
rm(B, S, base); invisible(gc(FALSE))
R[R >  CLIP] <-  CLIP
R[R < -CLIP] <- -CLIP

## 平滑：沿箱方向滑动平均（边界窗口自动收窄）
if (SMOOTH_W > 1L) {
  h <- SMOOTH_W %/% 2L
  nb <- nrow(R)
  ii <- integer(0); jj <- integer(0); xx <- numeric(0)
  for (d in -h:h) {
    from <- max(1, 1 - d); to <- min(nb, nb - d)
    if (from > to) next
    ii <- c(ii, from:to); jj <- c(jj, (from + d):(to + d)); xx <- c(xx, rep(1, to - from + 1))
  }
  Sm <- sparseMatrix(i = ii, j = jj, x = xx, dims = c(nb, nb))
  Sm <- Sm / Matrix::rowSums(Sm)
  R  <- as.matrix(Sm %*% R)
  rm(Sm, ii, jj, xx); invisible(gc(FALSE))
}

## ---------------- 5. CNV 分数与恶性判定 ----------------
## 分数 = 细胞剖面去掉自身中位后的标准差（与 R/05 同定义）：
## 整条染色体增减不改变剖面形状，多点上下起伏才拉高标准差。
med <- matrixStats::colMedians(R)
Rc  <- R - med          # 按列回收，等价于每列减去自身中位
score <- sqrt(colMeans(Rc^2))
meta[, cnv_score := score]

## 阈值：低分群体（底部 25%）视为二倍体核心，取其稳健离散度
core <- score <= quantile(score, 0.25, type = 7)
m0 <- median(score[core]); s0 <- mad(score[core])
THR <- max(m0 + 5 * s0, m0 + 0.02)
meta[, malignant := score > THR]
cat(sprintf("\nCNV 分数：二倍体核心 m0=%.4f  s0=%.4f  阈值=%.4f\n", m0, s0, THR))
cat(sprintf("判定恶性 %d / %d 上皮细胞 (%.1f%%)\n",
            sum(meta$malignant), nrow(meta), 100 * mean(meta$malignant)))

fwrite(meta[, .(cell, sample, patient, tp, nFeature, nCount, pct_mt,
                scEpi, cnv_score = round(cnv_score, 5), malignant)],
       file.path(TAB, "neoletexe_cnv_cells.csv"))

cat("\n=== 每样本恶性比例 ===\n")
ps <- meta[, .(n_epi = .N, n_mal = sum(malignant),
               frac = round(100 * mean(malignant), 1),
               med_score = round(median(cnv_score), 4)), by = .(sample, patient, tp)]
print(ps[order(patient, tp)])
cat("\n恶性比例分位：\n"); print(quantile(ps$frac, c(0, .1, .25, .5, .75, .9, 1)))

## ---------------- 6. 落盘剖面（供亚克隆用） ----------------
mal <- meta$malignant
saveRDS(list(R = R[, mal, drop = FALSE], cells = meta$cell[mal],
             bins = bin_info),
        file.path(PROC, "cnv_profile.rds"), compress = FALSE)
cat(sprintf("剖面已存：%d 箱 x %d 恶性细胞\n", nrow(R), sum(mal)))
