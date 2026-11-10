#!/usr/bin/env Rscript
# ==============================================================================
# 30_neoletexe_qc_gate.R
#
# NeoLetExe（Nat Genet 2026，25 例 ER+/HER2-，芳香化酶抑制剂新辅助，
# 三时点 A=基线 / B=换药 / C=术前）——QC、去双细胞、上皮门控。
#
# 为什么按样本逐个处理：23 GB 矩阵已拆成 73 个 per-sample RDS（R/23），
# 全程不需要把 47 万细胞放进内存，单个样本最大 1 万细胞。
#
# 处理链（每个样本独立）：
#   1. QC 指标：nFeature / nCount / 线粒体比例 / 核糖体比例
#   2. 过滤：nFeature >= 300、nCount >= 500、线粒体 <= 20%
#   3. scDblFinder 去双细胞（双细胞里的肿瘤+免疫融合会被 CNV 误判为恶性）
#   4. LogNormalize，按四套标记算谱系分数，取最高者作为细胞类型
#   5. 落盘：过滤后的 counts（稀疏）+ 元数据，供 R/31 做 CNV
#
# 【数值口径核查 2026-11-09】
#   文件名为 raw_counts，但实测**不是整数计数**：首个细胞 3438 个非零值里
#   只有 1.2% 是整数，最小正值 0.725，细胞求和 10,777（另一些细胞 17,449），
#   全表仅 1.6% 为整数；也不是对数尺度（最大值达 11,447，指数回代会溢出）。
#   最合理的解释是做过环境 RNA 校正 / 深度均衡的计数（会留下分数值）。
#   因此按"计数类"数据对待：LogNormalize 与 scDblFinder 均可用，
#   但**不能**把 nCount 当作测序深度绝对值来解读。这一条已写入分析记录。
#
# 输出：
#   data/proc/neoletexe/samples/<sample>.rds   list(counts, meta)
#   data/proc/neoletexe/gene_position.csv      基因坐标（CNV 分箱用）
#   results/tables/neoletexe_qc_per_sample.csv
#   results/tables/neoletexe_celltype.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(Matrix); library(data.table)
})

RAW   <- "data/raw/neoletexe/mat"
PROC  <- "data/proc/neoletexe"
SAMP  <- file.path(PROC, "samples")
TAB   <- "results/tables"
NCORE <- 4L

dir.create(SAMP, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB,  recursive = TRUE, showWarnings = FALSE)

MIN_FEATURE <- 300L
MIN_COUNT   <- 500L
MAX_MT      <- 20

## ---------------- 样本表 ----------------
sb <- fread(file.path(RAW, "sample_blocks.csv"))
sb[, patient := substr(sample, 1, 3)]
sb[, tp := substr(sample, 4, 4)]
sb <- sb[order(sample)]
## 便于先用少量样本试跑：NEOLET_ONLY="P24A,P25B"
only <- Sys.getenv("NEOLET_ONLY")
if (nzchar(only)) sb <- sb[sample %in% strsplit(only, ",")[[1]]]
cat(sprintf("样本 %d 个，患者 %d 例\n", nrow(sb), uniqueN(sb$patient)))
print(sb[, .N, by = .(patient)][, .N, by = N])

## ---------------- 基因坐标（CNV 分箱用） ----------------
MAP <- file.path(PROC, "gene_position.csv")
if (file.exists(MAP)) {
  gmap <- fread(MAP)
  cat("基因坐标缓存：", nrow(gmap), " 个\n")
} else {
  suppressPackageStartupMessages({
    library(AnnotationDbi); library(org.Hs.eg.db)
  })
  genes <- gsub('"', "", fread("data/raw/neoletexe/genes.csv", header = TRUE)[[1]])
  cat("构建基因坐标，共 ", length(genes), " 个 ...\n")
  ch <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = genes,
          keytype = "SYMBOL", columns = c("ENTREZID", "CHR")))
  ch <- as.data.table(ch)[!is.na(CHR)]
  ch <- ch[!duplicated(SYMBOL)]
  loc <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db,
           keys = unique(ch$ENTREZID), keytype = "ENTREZID", columns = "CHRLOC"))
  loc <- as.data.table(loc)[!is.na(CHRLOC)]
  loc <- loc[, .(ENTREZID, start = abs(CHRLOC))][!duplicated(ENTREZID)]
  gmap <- merge(ch, loc, by = "ENTREZID", all.x = TRUE)
  gmap <- gmap[CHR %in% c(as.character(1:22), "X", "Y")]
  gmap[, chr_num := fifelse(CHR == "X", 23L, fifelse(CHR == "Y", 24L, as.integer(CHR)))]
  setorder(gmap, chr_num, start)
  gmap <- gmap[, .(gene = SYMBOL, chr = CHR, chr_num, start)]
  fwrite(gmap, MAP)
  cat("基因坐标：", nrow(gmap), " 个可用\n")
}

## ---------------- 标记集 ----------------
# 注意：parLapply 的 worker 里函数的外围环境是 worker 的 globalenv，
# 不会带上主进程的 MK / RAW 等变量，故这些常量一律放进函数体内或显式传参。
MK <- list(
  Epi = c("EPCAM","KRT8","KRT18","KRT19","KRT7","CDH1","ELF3","FOXA1",
          "KRT15","KRT5","KRT14","AGR2","ESR1"),
  Imm = c("PTPRC","CD3D","CD3E","CD2","MS4A1","CD79A","MZB1","IGKC",
          "LYZ","CD68","CSF1R","NKG7","IL7R"),
  Str = c("COL1A1","COL1A2","COL3A1","DCN","LUM","PDGFRB","ACTA2","TAGLN","FBLN1"),
  End = c("PECAM1","VWF","CDH5","PLVAP")
)

## ---------------- 单样本处理 ----------------
process_one <- function(sn, RAW, SAMP, MIN_FEATURE, MIN_COUNT, MAX_MT, MK) {
  suppressPackageStartupMessages({
    library(Matrix); library(data.table); library(scDblFinder)
    library(SingleCellExperiment)
  })
  m <- readRDS(file.path(RAW, sprintf("%s.rds", sn)))
  g <- rownames(m)
  mt   <- grepl("^MT-", g)
  ribo <- grepl("^RP[SL]", g)

  nn <- Matrix::colSums(m > 0)
  nc <- Matrix::colSums(m)
  pm <- if (any(mt)) Matrix::colSums(m[mt, , drop = FALSE]) / pmax(nc, 1) * 100 else 0
  pr <- if (any(ribo)) Matrix::colSums(m[ribo, , drop = FALSE]) / pmax(nc, 1) * 100 else 0

  keep <- nn >= MIN_FEATURE & nc >= MIN_COUNT & pm <= MAX_MT
  n0 <- ncol(m)
  m <- m[, keep, drop = FALSE]
  nn <- nn[keep]; nc <- nc[keep]; pm <- pm[keep]; pr <- pr[keep]

  ## 去双细胞（细胞数过少则跳过）
  dbl <- rep("singlet", ncol(m))
  if (ncol(m) >= 400) {
    sce <- SingleCellExperiment(list(counts = m))
    sce <- scDblFinder(sce, verbose = FALSE)
    dbl <- as.character(sce$scDblFinder.class)
    rm(sce)
  }
  sing <- dbl == "singlet"
  m <- m[, sing, drop = FALSE]
  nn <- nn[sing]; nc <- nc[sing]; pm <- pm[sing]; pr <- pr[sing]; dbl <- dbl[sing]

  ## 归一化 + 谱系分数
  ls <- Matrix::colSums(m)
  ln <- m
  ln@x <- ln@x / rep.int(ls, diff(ln@p))
  ln@x <- log1p(ln@x * 1e4)

  sc <- sapply(MK, function(gs) {
    gs <- intersect(gs, g)
    if (!length(gs)) return(rep(NA_real_, ncol(ln)))
    as.numeric(Matrix::colMeans(ln[gs, , drop = FALSE]))
  })
  colnames(sc) <- names(MK)
  best <- max.col(sc, ties.method = "first")
  mx   <- sc[cbind(seq_len(nrow(sc)), best)]
  lab  <- names(MK)[best]
  lab[!is.finite(mx) | mx < 0.10] <- "Other"

  meta <- data.table(
    cell = colnames(m), sample = sn,
    patient = substr(sn, 1, 3), tp = substr(sn, 4, 4),
    nFeature = as.integer(nn), nCount = as.numeric(nc),
    pct_mt = round(pm, 2), pct_ribo = round(pr, 2),
    scEpi = sc[, "Epi"], scImm = sc[, "Imm"],
    scStr = sc[, "Str"], scEnd = sc[, "End"], celltype = lab)

  saveRDS(list(counts = m, meta = meta), file.path(SAMP, sprintf("%s.rds", sn)),
          compress = FALSE)
  rm(m, ln); invisible(gc(FALSE))

  list(qc = data.table(sample = sn, patient = substr(sn, 1, 3),
                       tp = substr(sn, 4, 4), n_raw = n0, n_qc = sum(keep),
                       n_singlet = nrow(meta),
                       med_feature = as.integer(median(meta$nFeature)),
                       med_count = as.integer(median(meta$nCount)),
                       med_mt = median(meta$pct_mt)),
       ct = meta[, .N, by = .(sample, patient, tp, celltype)])
}

## ---------------- 并行执行 ----------------
cat(sprintf("\n开始处理 %d 个样本，%d 进程 ...\n", nrow(sb), NCORE))
t0 <- Sys.time()
run_one <- function(s) process_one(s, RAW, SAMP, MIN_FEATURE, MIN_COUNT, MAX_MT, MK)
if (NCORE > 1L) {
  cl <- parallel::makeCluster(NCORE)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  ## worker 的 globalenv 是空的，函数体里引用的名字必须显式送过去
  parallel::clusterExport(cl, c("process_one", "MK", "RAW", "SAMP",
                                "MIN_FEATURE", "MIN_COUNT", "MAX_MT"),
                          envir = globalenv())
  res <- parallel::parLapplyLB(cl, sb$sample, function(s) {
    tryCatch(process_one(s, RAW, SAMP, MIN_FEATURE, MIN_COUNT, MAX_MT, MK),
             error = function(e) {
      list(qc = data.table(sample = s, n_raw = NA, n_qc = NA, n_singlet = NA,
                           med_feature = NA, med_count = NA, med_mt = NA,
                           patient = NA, tp = NA, err = conditionMessage(e)),
           ct = NULL)
    })
  })
} else {
  res <- lapply(sb$sample, function(s) tryCatch(run_one(s), error = function(e) {
    list(qc = data.table(sample = s, n_raw = NA, n_qc = NA, n_singlet = NA,
                         med_feature = NA, med_count = NA, med_mt = NA,
                         patient = NA, tp = NA, err = conditionMessage(e)),
         ct = NULL)
  }))
}

qc <- rbindlist(lapply(res, `[[`, "qc"), fill = TRUE)
ct <- rbindlist(lapply(res, `[[`, "ct"), fill = TRUE)
fwrite(qc, file.path(TAB, "neoletexe_qc_per_sample.csv"))
fwrite(ct, file.path(TAB, "neoletexe_celltype.csv"))

cat(sprintf("\n完成，耗时 %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat("\n=== QC 汇总 ===\n")
cat(sprintf("原始细胞 %d -> 过滤后 %d (%.1f%%)  双细胞去除平均 %.1f%%\n",
            sum(qc$n_raw), sum(qc$n_qc), 100 * sum(qc$n_qc) / sum(qc$n_raw),
            100 * mean(1 - qc$n_singlet / qc$n_qc)))
print(qc[, .(n_sample = .N, cells = sum(n_singlet),
             med_feature = as.integer(median(med_feature)),
             med_mt = round(median(med_mt), 2)), by = .(patient)][order(patient)])

cat("\n=== 细胞类型构成 ===\n")
tot <- ct[, .(N = sum(N)), by = celltype][order(-N)]
tot[, pct := round(100 * N / sum(N), 1)]
print(tot)
cat("\n=== 每个样本的上皮细胞数（分位数）===\n")
epi <- ct[celltype == "Epi", .(n_epi = sum(N)), by = .(sample, patient, tp)]
print(quantile(epi$n_epi, c(0, .1, .25, .5, .75, .9, 1)))
