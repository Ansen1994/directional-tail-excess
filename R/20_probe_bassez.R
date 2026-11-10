# 探查 Bassez 2021 处理后数据的结构
# 只读取对象并打印元数据摘要，不做分析

suppressPackageStartupMessages({
  library(Seurat); library(data.table)
})

f <- "data/raw/bassez/bassez_downsampled.rds"
cat("读取 ", f, " ...\n")
obj <- readRDS(f)

cat("\n=== 对象类别 ===\n"); print(class(obj))
if (is.list(obj) && !inherits(obj, "Seurat")) {
  cat("是 list，元素：\n"); print(names(obj))
  obj <- obj[[1]]
}

cat("\n=== 维度 ===\n")
print(dim(obj))
cat("\n=== assays ===\n"); print(Assays(obj))
cat("\n=== 默认 assay / slot ===\n")
print(DefaultAssay(obj))

md <- obj@meta.data
cat("\n=== meta.data 列 ===\n")
print(colnames(md))
cat("\n行数：", nrow(md), "\n")

cat("\n=== 各列唯一值数（<=30 的全部列出）===\n")
for (cn in colnames(md)) {
  v <- md[[cn]]
  u <- unique(v)
  if (length(u) <= 30) {
    cat(sprintf("\n[%s]  (%d 类)\n", cn, length(u)))
    print(sort(table(v, useNA = "ifany"), decreasing = TRUE))
  } else {
    cat(sprintf("\n[%s]  唯一值 %d 个，示例：%s\n", cn, length(u),
                paste(head(as.character(u), 5), collapse = ", ")))
  }
}

cat("\n=== 每个样本的细胞数（前 50）===\n")
cand <- intersect(c("sample", "Sample", "orig.ident", "sample_id", "Library",
                    "patient", "Patient", "patient_id"), colnames(md))
cat("候选样本列：", paste(cand, collapse = ", "), "\n")
if (length(cand)) {
  for (cn in cand) {
    cat("\n--", cn, "--\n")
    tb <- sort(table(md[[cn]]), decreasing = TRUE)
    print(head(tb, 50))
  }
}

cat("\n=== 细胞类型相关列 ===\n")
ct <- grep("cell|type|annot|cluster|ident|major", colnames(md),
           ignore.case = TRUE, value = TRUE)
print(ct)

cat("\n=== 内存 ===\n"); print(gc())
