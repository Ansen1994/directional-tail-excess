# 00_run_all.R —— 按依赖顺序跑本稿主线（§A），每个脚本在独立 R 进程里执行
#
# 用法（必须在项目根目录）：
#   Rscript R/00_run_all.R              # 全跑
#   Rscript R/00_run_all.R --dry-run    # 只打印计划，不执行
#   Rscript R/00_run_all.R --from=28    # 从 tag 为 28 的那一步开始
#
# 每一步跑完，自动把该脚本所在进程的 sessionInfo() 写到
#   logs/sessionInfo_<tag>.txt
# 标准输出/错误写到
#   logs/run_<tag>.log
#
# 为什么 R/28 出现两次：
#   R/28 画 Fig 5D 时要读 R/39 落盘的 selection_counterfactual_paired.csv，
#   而 R/39 要读 R/28 落盘的 nc2025_fig5_shape_calls.csv —— 两者互为输入。
#   仓库自带这两个表，所以常规情况一次即可；若完全从零重建（tables 为空），
#   第一遍 R/28 会在读反事实表处报错，属预期，跑完 R/39 再跑一遍即可。

args <- commandArgs(trailingOnly = TRUE)
dry  <- "--dry-run" %in% args
from_tag <- sub("^--from=", "", grep("^--from=", args, value = TRUE))
if (length(from_tag) == 0) from_tag <- NA_character_

if (!dir.exists("results/tables") || !dir.exists("R")) {
  stop("请在项目根目录运行：cd projects/scrna-analysis && Rscript R/00_run_all.R",
       call. = FALSE)
}
dir.create("logs", showWarnings = FALSE)

# ---- 执行计划 -------------------------------------------------------------
# phase A：产出/刷新细胞级与组合级表
# phase B：对照与检验
# phase C：全部图（消费 A、B 的产物）
plan <- data.frame(
  tag    = c("21", "29", "23", "30", "31", "32", "33",
             "38", "39", "40", "41", "42",
             "24", "25", "27", "26", "28", "29", "34", "35", "36", "37", "28"),
  script = c("R/21_bassez_validation.R",
             "R/29_fig6_external.R",
             "R/23_neoletexe_load.R",
             "R/30_neoletexe_qc_gate.R",
             "R/31_neoletexe_cnv.R",
             "R/32_neoletexe_shape.R",
             "R/33_neoletexe_clones.R",
             "R/38_permutation_shape.R",
             "R/39_selection_counterfactual.R",
             "R/40_shape_calibration.R",
             "R/41_depth_matched_downsample.R",
             "R/42_shape_stat_correlation.R",
             "R/24_fig1_cohort.R",
             "R/25_fig2_validity.R",
             "R/27_fig3_clones.R",
             "R/26_fig4_heterogeneity.R",
             "R/28_fig5_shape.R",
             "R/29_fig6_external.R",
             "R/34_neoletexe_fig.R",
             "R/35_figS1_stress.R",
             "R/36_figS2_sensitivity.R",
             "R/37_figS3_qc.R",
             "R/28_fig5_shape.R"),
  phase  = c(rep("A", 7), rep("B", 5), rep("C", 11)),
  note   = c("Bassez 形状表", "Bassez 细胞级分数（第一遍）",
             "载入 NeoLetExe", "QC 闸门", "CNV 分箱", "形状表 + 细胞级分数", "克隆层",
             "标签置换（Fig 5E 零分布）", "逐组合选择反事实（Fig 5D/6A/7B 参照线）",
             "形状判据校准（Fig S4）", "深度匹配下采样（Table S2）",
             "四个形状量相关矩阵（Table S3）",
             "Fig 1", "Fig 2", "Fig 3", "Fig 4", "Fig 5（第一遍）",
             "Fig 6（第二遍，用最终反事实）", "Fig 7", "Fig S1", "Fig S2", "Fig S3",
             "Fig 5（第二遍，用最终反事实）"),
  stringsAsFactors = FALSE
)
plan$order <- seq_len(nrow(plan))

if (!is.na(from_tag)) {
  hit <- which(plan$tag == from_tag)
  if (length(hit) == 0) stop("找不到 tag：", from_tag, call. = FALSE)
  plan <- plan[plan$order >= min(hit), ]
}

miss <- plan$script[!file.exists(plan$script)]
if (length(miss)) stop("缺少脚本：", paste(miss, collapse = ", "), call. = FALSE)

rscript <- file.path(R.home("bin"), "Rscript")

cat("执行计划：", nrow(plan), "步\n")
for (i in seq_len(nrow(plan))) {
  cat(sprintf("  %2d  [%s] %-32s %s\n", i, plan$phase[i],
              plan$script[i], plan$note[i]))
}
if (dry) {
  cat("\n--dry-run：未执行任何脚本。\n")
  quit(save = "no", status = 0)
}

# ---- 逐步执行 -------------------------------------------------------------
res <- data.frame(step = integer(0), tag = character(0), script = character(0),
                  status = character(0), secs = numeric(0))
for (i in seq_len(nrow(plan))) {
  sc <- plan$script[i]; tg <- plan$tag[i]
  info <- file.path("logs", sprintf("sessionInfo_%s.txt", tg))
  logf <- file.path("logs", sprintf("run_%s.log", tg))
  wrap <- file.path("logs", sprintf("_wrap_%s.R", tg))

  writeLines(c(
    sprintf('source("%s", echo = FALSE)', sc),
    sprintf('writeLines(capture.output(sessionInfo()), "%s")', info)
  ), wrap)

  cat(sprintf("\n[%d/%d] %s —— %s\n", i, nrow(plan), sc, plan$note[i]))
  t0 <- Sys.time()
  st <- try(system2(rscript, c("--vanilla", wrap),
                    stdout = logf, stderr = logf), silent = TRUE)
  secs <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  unlink(wrap)

  ok <- !inherits(st, "try-error") && identical(as.integer(st), 0L)
  cat(sprintf("      %s  (%.1f s, 日志 logs/run_%s.log)\n",
              if (ok) "OK" else "FAILED", secs, tg))

  res <- rbind(res, data.frame(step = i, tag = tg, script = sc,
                               status = if (ok) "ok" else "failed", secs = secs))
  if (!ok) {
    cat("\n已中止。修好后再用 --from=", tg, " 从这里继续。\n", sep = "")
    print(res, row.names = FALSE)
    quit(save = "no", status = 1)
  }
}

cat("\n全部完成。\n")
print(res, row.names = FALSE)
cat("\n每个脚本的 sessionInfo：logs/sessionInfo_<tag>.txt\n")
