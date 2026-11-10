#!/usr/bin/env Rscript
# ==============================================================================
# NeoLetExe raw_counts.mtx -> 按样本拆分的 dgCMatrix + 细胞级 QC 表
#
# 背景：raw_counts.mtx 23.25 GB 纯文本（MatrixMarket coordinate，
#       35461 基因 x 472737 细胞，768,497,319 非零），无法整块读入。
#       单遍流式扫描：一边解析、一边按 barcode 前缀分流到 73 个样本缓冲，
#       样本切换时落盘。
#
# 用法：
#   Rscript R/23_neoletexe_load.R [mtx] [输出目录] [每块行数]
# 默认：
#   mtx  = data/raw/neoletexe/raw_counts.mtx
#   out  = data/raw/neoletexe/mat
#
# 修订（2026-11-09）：
#   1. 原 `v[jj] <- v[jj] + 1L` 在 jj 重复时只有最后一条生效（每块每列只计 1 次），
#      QC 严重失真；而 tabulate(weights=) 在 R 4.4 中并不存在，直接报错。
#      改为**循环里不累加**：每个样本落盘后由稀疏矩阵自身算 QC——
#      nFeature = diff(m@p)（列指针差，精确的每列非零数），
#      nCount = colSums(m)，基因总量 = rowSums(m)。
#   2. 增加逐块进度与耗时估计。
#   3. 记录样本边界表，供后续按样本并行读取。
# ==============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
MTX   <- if (length(args) >= 1) args[1] else "data/raw/neoletexe/raw_counts.mtx"
OUTD  <- if (length(args) >= 2) args[2] else "data/raw/neoletexe/mat"
CHUNK <- if (length(args) >= 3) as.integer(args[3]) else 4000000L
MIN_FEATURE <- 200L

dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)
QCDIR <- file.path(OUTD, "qc")
dir.create(QCDIR, recursive = TRUE, showWarnings = FALSE)

## ---------------- 1. barcode -> 样本归属 ----------------
bc <- fread("data/raw/neoletexe/barcodes.csv", header = TRUE)
setnames(bc, 1L, "barcode")
bc[, barcode := trimws(gsub('"', "", barcode))]
bc[, sample := sub("_.*$", "", barcode)]
bc[, cell_idx := .I]
cat(sprintf("barcodes: %d 个细胞 / %d 个样本\n", nrow(bc), uniqueN(bc$sample)))

## ---------------- 2. 文件头 ----------------
con <- file(MTX, "r")
hdr <- readLines(con, n = 6L, warn = FALSE)
cat("--- 文件头 ---\n"); cat(paste(hdr, collapse = "\n"), "\n")
dimline <- hdr[!grepl("^%", hdr)][1]
dims    <- as.numeric(strsplit(trimws(dimline), "[[:space:]]+")[[1]])
NG <- dims[1]; NC <- dims[2]; NNZ <- dims[3]
cat(sprintf("--- 维度: %d 基因 x %d 细胞, 非零 %d ---\n", NG, NC, NNZ))
if (NC != nrow(bc)) stop(sprintf("列数 %d 与 barcode 数 %d 不一致", NC, nrow(bc)))

genes <- fread("data/raw/neoletexe/genes.csv", header = TRUE)[[1]]
genes <- gsub('"', "", genes)
if (length(genes) != NG) stop(sprintf("基因数 %d 与矩阵行数 %d 不一致", length(genes), NG))

## 样本边界（要求 barcode 按样本分块排列）
sb <- bc[, .(start = min(cell_idx), end = max(cell_idx), n = .N), by = sample][order(start)]
gap <- sum(sb$end - sb$start + 1) - sum(sb$n)
cat(sprintf("样本块数 %d；区间总长 %d vs 实际 %d（差 %d）\n",
            nrow(sb), sum(sb$end - sb$start + 1), sum(sb$n), gap))
if (gap != 0) stop("barcode 未按样本连续排列，边界假设失效")
fwrite(sb, file.path(OUTD, "sample_blocks.csv"))
brk <- sb$end

## ---------------- 3. 单遍流式扫描 ----------------
gene_tot  <- numeric(NG)
cell_nnz  <- numeric(NC)
cell_sum  <- numeric(NC)

cur_s <- NA_integer_
bi <- integer(0); bj <- integer(0); bx <- numeric(0)
n_written <- 0L

flush_sample <- function(sid) {
  if (is.na(sid) || length(bi) == 0L) return(invisible(NULL))
  sn  <- sb$sample[sid]
  j0  <- sb$start[sid] - 1L
  ncl <- sb$end[sid] - j0
  cols <- (j0 + 1L):(j0 + ncl)
  jrel <- bj - j0
  if (min(jrel) < 1L || max(jrel) > ncl)
    stop(sprintf("样本 %s 的列号越界: [%d, %d]，应在 [1, %d]",
                 sn, min(jrel), max(jrel), ncl))
  if (min(bi) < 1L || max(bi) > NG) stop(sprintf("样本 %s 的行号越界", sn))
  m <- sparseMatrix(i = bi, j = jrel, x = bx, dims = c(NG, ncl),
                    dimnames = list(genes, bc$barcode[cols]),
                    giveCsparse = TRUE)
  ## QC 直接取稀疏结构，避免循环内累加的重复索引问题
  cell_nnz[cols] <<- diff(m@p)
  cell_sum[cols] <<- Matrix::colSums(m)
  gene_tot       <<- gene_tot + Matrix::rowSums(m)

  saveRDS(m, file.path(OUTD, sprintf("%s.rds", sn)), compress = FALSE)
  n_written <<- n_written + 1L
  cat(sprintf("  [落盘 %2d/%d] %-6s %6d 细胞  nnz=%d\n",
              n_written, nrow(sb), sn, ncl, length(bi)))
  invisible(NULL)
}

t0 <- Sys.time()
n_read <- 0L
repeat {
  lines <- readLines(con, n = CHUNK, warn = FALSE)
  if (length(lines) == 0L) break
  n_read <- n_read + length(lines)

  sp <- scan(text = lines, what = list(integer(), integer(), numeric()),
             sep = " ", quiet = TRUE)
  ii <- sp[[1]]; jj <- sp[[2]]; xx <- sp[[3]]
  rm(lines, sp)

  ## 区间是 (brk[k-1], brk[k]]，故先减 1 再查；否则每个样本的最后一个细胞
  ## （j 恰好等于 brk[k]）会被判给下一个样本，导致 j - j0 = 0 报错
  sid <- findInterval(jj - 1L, brk) + 1L
  sid[sid > nrow(sb)] <- nrow(sb)

  r <- rle(sid)
  pos <- cumsum(c(1L, r$lengths))
  for (k in seq_along(r$values)) {
    a <- pos[k]; b <- pos[k + 1L] - 1L
    s <- r$values[k]
    if (!is.na(cur_s) && s != cur_s) {
      flush_sample(cur_s); bi <- integer(0); bj <- integer(0); bx <- numeric(0)
    }
    cur_s <- s
    bi <- c(bi, ii[a:b]); bj <- c(bj, jj[a:b]); bx <- c(bx, xx[a:b])
  }
  rm(ii, jj, xx, sid, r, pos); invisible(gc(FALSE))

  el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf("  %5.1f%%  %s / %s 行  已用 %.1f min  预计剩余 %.1f min\n",
              100 * n_read / NNZ, format(n_read, big.mark = ","),
              format(NNZ, big.mark = ","), el,
              el * (NNZ - n_read) / max(n_read, 1)))
}
close(con)
flush_sample(cur_s)

## ---------------- 4. QC 表 ----------------
qc <- data.table(barcode = bc$barcode, sample = bc$sample,
                 nCount = as.numeric(cell_sum), nFeature = as.integer(cell_nnz))
qc[, keep := nFeature >= MIN_FEATURE]
fwrite(qc, file.path(QCDIR, "cell_qc.csv"))
fwrite(data.table(gene = genes, total = gene_tot),
       file.path(QCDIR, "gene_totals.csv"))

cat(sprintf("\n完成。细胞 %d，其中 nFeature>=%d 的 %d 个 (%.1f%%)\n",
            nrow(qc), MIN_FEATURE, sum(qc$keep), 100 * mean(qc$keep)))
cat(sprintf("中位 nFeature = %d, 中位 nCount = %d\n",
            as.integer(median(qc$nFeature)), as.integer(median(qc$nCount))))
cat(sprintf("输出目录: %s\n", normalizePath(OUTD)))
