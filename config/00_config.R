## =============================================================================
## 00_config.R — 全项目唯一参数中心
##
## 用法：所有脚本假定工作目录为项目根目录，开头 source("config/00_config.R")。
## 任何环境变量均可覆盖同名参数（见下），便于抽样试跑而不改代码。
##
## 重要环境限制：本机 R 启动时 LC_CTYPE=C，无法打开含中文的路径。
## 因此数据通过 ASCII 目录联接 data/raw_sc 访问（指向原始目录），
## 所有输出也一律落在纯 ASCII 路径下。勿把 RAW_DIR 改回中文绝对路径。
## =============================================================================

WB_ROOT   <- normalizePath(".")
RAW_DIR   <- Sys.getenv("WB_SC_RAW_DIR",  unset = "data/raw_sc")
PROC_DIR  <- Sys.getenv("WB_SC_PROC_DIR", unset = "data/processed")
RES_DIR   <- Sys.getenv("WB_SC_RES_DIR",  unset = "results")
FIG_DIR   <- file.path(RES_DIR, "figures")
TAB_DIR   <- file.path(RES_DIR, "tables")
LOG_DIR   <- "logs"

for (d in c(PROC_DIR, RES_DIR, FIG_DIR, TAB_DIR, LOG_DIR)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

## ------------------------------------------------------------------ 随机种子
WB_SEED <- as.integer(Sys.getenv("WB_SEED", unset = "20261108"))

## ------------------------------------------------------------ 抽样试跑开关
## SC_SAMPLE_LIMIT    >0 时只取前 N 个样本
## SC_CELLS_PER_SAMPLE >0 时每个样本随机抽 N 个细胞（试跑用，正式分析必须为 0）
WB_SC_SAMPLE_LIMIT   <- as.integer(Sys.getenv("WB_SC_SAMPLE_LIMIT",   unset = "0"))
WB_SC_CELLS_PER_SAMPLE <- as.integer(Sys.getenv("WB_SC_CELLS_PER_SAMPLE", unset = "0"))

## ------------------------------------------------------------------ QC 阈值
## 硬阈值兜底 + 每样本自适应 MAD 收紧，取两者中更严格的一侧（见 01_qc.R）
QC_MIN_FEATURES   <- as.integer(Sys.getenv("WB_QC_MIN_FEATURES", unset = "200"))
QC_MIN_COUNTS     <- as.integer(Sys.getenv("WB_QC_MIN_COUNTS",   unset = "500"))
QC_MAX_MITO       <- as.numeric(Sys.getenv("WB_QC_MAX_MITO",     unset = "20"))
QC_MAX_FEATURES   <- as.integer(Sys.getenv("WB_QC_MAX_FEATURES", unset = "6000"))
QC_MAD_N          <- as.numeric(Sys.getenv("WB_QC_MAD_N",        unset = "3"))
QC_USE_MAD        <- Sys.getenv("WB_QC_USE_MAD", unset = "1") == "1"
QC_MIN_CELLS_GENE <- as.integer(Sys.getenv("WB_QC_MIN_CELLS_GENE", unset = "3"))

## ------------------------------------------------------------- 双细胞检测
## 10x Chromium 的经验装载量-双细胞率关系：约每 1000 个回收细胞 0.8%。
## 逐样本按各自细胞数推算先验，避免合并后先验被大样本主导。
DBL_RATE_PER_1K <- as.numeric(Sys.getenv("WB_DBL_RATE_PER_1K", unset = "0.008"))
DBL_MAX_RATE    <- as.numeric(Sys.getenv("WB_DBL_MAX_RATE",    unset = "0.10"))

## ------------------------------------------------------------ CNV 恶性判定
## copykat 参数：ngene.chr 每条染色体最少基因数；win.size 滑窗基因数；
## KS.cut 分段阈值（0.05 更敏感但碎片化，0.1 为作者推荐默认）。
## CNV_MAX_REF 限制参考细胞数，避免样本间参考量差异过大拖慢运行。
##
## win.size 取 101 而非默认 25：copykat 的耗时与窗口数成正比（逐细胞 dlm
## 平滑），实测 25 基因窗下单例 4,352 细胞需约 50 分钟，全量要 2 小时以上。
## 101 基因窗给出的是染色体臂级分辨率，分辨"恶性 vs 正常上皮"足够，但对
## 局灶扩增不敏感。本项目只用它做恶性判定，不做克隆级 CNV 结构分析，
## 因此接受这一取舍。若将来要做亚克隆，须改回 25 重跑。
CNV_NGENE_CHR <- as.integer(Sys.getenv("WB_CNV_NGENE_CHR", unset = "5"))
CNV_WIN_SIZE  <- as.integer(Sys.getenv("WB_CNV_WIN_SIZE",  unset = "101"))
CNV_KS_CUT    <- as.numeric(Sys.getenv("WB_CNV_KS_CUT",    unset = "0.1"))
CNV_MAX_REF   <- as.integer(Sys.getenv("WB_CNV_MAX_REF",   unset = "600"))
## 轻量 CNV 方案（05_cnv_light.R）参数
## CNV_WIN       基因组上平滑窗口（基因数），101 ≈ 染色体臂级
## CNV_CLIP      残差截断，抑制单个高表达基因的杠杆效应
## CNV_SD_MULT   阈值 = 参考免疫细胞分数均值 + N × 标准差
CNV_WIN        <- as.integer(Sys.getenv("WB_CNV_WIN",      unset = "101"))
CNV_CLIP       <- as.numeric(Sys.getenv("WB_CNV_CLIP",     unset = "0.35"))
CNV_SD_MULT    <- as.numeric(Sys.getenv("WB_CNV_SD_MULT",  unset = "3"))
CNV_MIN_EPI    <- as.integer(Sys.getenv("WB_CNV_MIN_EPI",  unset = "50"))
CNV_MIN_REF    <- as.integer(Sys.getenv("WB_CNV_MIN_REF",  unset = "100"))
## 阈值取参考（其他样本的正常上皮）分数的 95 分位
CNV_REF_Q      <- as.numeric(Sys.getenv("WB_CNV_REF_Q",    unset = "0.95"))
## 分片格式 "k/n"：只跑第 k 片（从 1 计），用于多进程并行；空字符串为不分片
WB_CNV_SHARD  <- Sys.getenv("WB_CNV_SHARD", unset = "")

## ------------------------------------------------------------ NMF 程序分解
## 学习集每患者抽 NMF_CELLS_PER_SAMPLE 个细胞，避免大样本主导基向量；
## 秩的候选集与重复次数（重复次数决定 cophenetic 与稳定性估计的可靠度）。
NMF_CELLS_PER_SAMPLE <- as.integer(Sys.getenv("WB_NMF_CELLS_PER_SAMPLE", unset = "150"))
NMF_MIN_PCT     <- as.numeric(Sys.getenv("WB_NMF_MIN_PCT", unset = "0.10"))
NMF_TOP_GENES   <- as.integer(Sys.getenv("WB_NMF_TOP_GENES", unset = "30"))
## 重复次数直接决定 cophenetic / 稳定性的估计误差，RcppML 单次拟合在
## 788 x 2647 上不足 1 秒，所以给到 20 次而不是 8 次。
NMF_NRUN        <- as.integer(Sys.getenv("WB_NMF_NRUN", unset = "20"))
NMF_RANKS       <- as.integer(strsplit(
  Sys.getenv("WB_NMF_RANKS", unset = "4,5,6,7,8"), ",")[[1]])
NMF_TOL         <- as.numeric(Sys.getenv("WB_NMF_TOL",   unset = "1e-4"))
NMF_MAXIT       <- as.integer(Sys.getenv("WB_NMF_MAXIT", unset = "300"))

## ------------------------------------------- 逐样本 NMF 与跨样本复现筛选
## 全局 NMF 学到的基向量优先捕获样本间差异（见 docs/NMF程序分解与诊断_v1.md），
## 所以改成每个样本单独拟合、再把各样本的程序对齐，只保留在多数样本中复现的。
NMF_PS_NRUN     <- as.integer(Sys.getenv("WB_NMF_PS_NRUN", unset = "10"))
NMF_PS_MIN_CELL <- as.integer(Sys.getenv("WB_NMF_PS_MIN_CELL", unset = "100"))
## 一个共识程序要在多少个样本里出现才算复现（分母是可匹配的样本数）
NMF_PS_MIN_SMP  <- as.integer(Sys.getenv("WB_NMF_PS_MIN_SMP", unset = "5"))
NMF_PS_MIN_COS  <- as.numeric(Sys.getenv("WB_NMF_PS_MIN_COS", unset = "0.70"))

## ------------------------------------------------------------- 整合与聚类
CL_N_HVG        <- as.integer(Sys.getenv("WB_CL_N_HVG",    unset = "2000"))
CL_N_PCS        <- as.integer(Sys.getenv("WB_CL_N_PCS",    unset = "30"))
CL_N_NEIGHBORS  <- as.integer(Sys.getenv("WB_CL_N_NEIGHBORS", unset = "30"))
CL_MIN_DIST     <- as.numeric(Sys.getenv("WB_CL_MIN_DIST", unset = "0.3"))
CL_BATCH_VAR    <- Sys.getenv("WB_CL_BATCH_VAR", unset = "sample")
CL_RESOLUTIONS  <- as.numeric(strsplit(
  Sys.getenv("WB_CL_RESOLUTIONS", unset = "0.1,0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.2"), ",")[[1]])
CL_MAIN_RES     <- as.numeric(Sys.getenv("WB_CL_MAIN_RES", unset = "0.5"))
CL_SKIP_HARMONY <- Sys.getenv("WB_CL_SKIP_HARMONY", unset = "0") == "1"

## ------------------------------------------------------------------ 差异/注释
DE_MIN_PCT  <- as.numeric(Sys.getenv("WB_DE_MIN_PCT",  unset = "0.25"))
DE_LOGFC    <- as.numeric(Sys.getenv("WB_DE_LOGFC",    unset = "0.25"))
DE_ONLY_POS <- Sys.getenv("WB_DE_ONLY_POS", unset = "1") == "1"

## ------------------------------------------------------------ 基因集与配色
## QC 指标用的基因前缀（按 symbol 匹配，大小写不敏感）
GENE_PAT_MITO <- "^MT-"
GENE_PAT_RIBO <- "^RP[SL][0-9]"
GENE_PAT_HB   <- "^HB[AB][0-9]?$"

## 经典谱系标记（用于人工核对独立聚类结果，不作为注释依据）
LINEAGE_MARKERS <- list(
  Epithelial    = c("EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"),
  Luminal       = c("ESR1", "PGR", "KRT7", "GATA3"),
  Basal         = c("KRT5", "KRT14", "KRT17", "TP63"),
  Proliferation = c("MKI67", "TOP2A", "PCNA"),
  T_cell        = c("CD3D", "CD3E", "CD2", "IL7R"),
  B_cell        = c("CD79A", "CD79B", "MS4A1", "IGHM"),
  Plasma        = c("MZB1", "JCHAIN", "DERL3"),
  Myeloid       = c("LYZ", "CD68", "CD14", "AIF1", "ITGAM"),
  Mast          = c("TPSAB1", "TPSB2", "CPA3", "KIT"),
  Endothelial   = c("PECAM1", "VWF", "EMCN", "ACKR1"),
  Fibroblast    = c("COL1A1", "COL1A2", "DCN", "LUM", "PDGFRB"),
  Perivascular  = c("ACTA2", "MYH11", "RGS5", "TAGLN"),
  Adipocyte     = c("ADIPOQ", "PLIN1", "LEP", "FABP4")
)

## 统一配色（与 rnaseq-analysis 一致的莫兰迪基调，避免两项目图件风格割裂）
PAL_MAJOR <- c(
  Epithelial    = "#7BA7BC", Luminal = "#A8C6D8", Basal = "#5B8CA8",
  T_cell        = "#E8A87C", B_cell = "#D98E73", Plasma = "#C97B63",
  Myeloid       = "#9CAF88", Mast = "#C4B7A6",
  Endothelial   = "#B8A9C9", Fibroblast = "#D4C4A8",
  Perivascular  = "#A99B8D", Adipocyte = "#E0D3C3", Other = "#BFBFBF"
)

## ------------------------------------------------------------------ 工具函数
## flush 必要：输出重定向到文件时 stdout 是块缓冲，不 flush 会导致日志长时间空白，
## 无法判断长任务是在运行还是卡死。
msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...), "\n", sep = "")
  flush(stdout())
}
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a

## 写 CSV 时统一用相对路径 + UTF-8，规避 locale 问题
write_csv_safe <- function(x, path, row.names = FALSE) {
  utils::write.csv(x, path, row.names = row.names, fileEncoding = "")
  invisible(path)
}
