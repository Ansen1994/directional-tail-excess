#!/usr/bin/env Rscript
# ==============================================================================
# 33_neoletexe_clones.R
#
# NeoLetExe：核型式亚克隆 -> 克隆内位移是否同向 / 组成能否解释位移
#
# 这是本队列相对原论文的增量：NeoLetExe 原文比较的是**克隆间**表达差异，
# 没有做同一克隆治疗前后的自身对比，也没做组成反事实。
#
# 为什么用"核型式"而不是层次聚类：
#   第一版用 227 箱剖面做 ward.D2 聚类，覆盖率只有 2-8%——100 基因一箱的
#   噪声盖过了 CNV 信号（箱内 SD 0.051，而染色体级差异只有 ±0.1）。
#   改为：每细胞先减去自身中位（去掉全局偏移），再取染色体级均值（24 维），
#   按 ±0.045 二值化为 增/减/不变，模式相同即同一克隆。
#   阈值 0.045 由覆盖度扫描定：0.025 覆盖 14-23%，0.035 覆盖 56-62%，
#   0.045 覆盖 82-86%，再往上克隆数骤减。等价于 SCEVAN 输出核型的粒度。
#
# 局限（写进结果）：染色体级分辨率，识别不了臂内或 focal 事件；
#   克隆边界依赖阈值，故同时报告 0.035 / 0.045 / 0.055 的敏感性。
#
# 判据（与 py/13_clone_layering_all.py 同口径）：
#   dir_agree   克隆内位移方向与患者整体位移方向一致的比例
#   sep_z       随访"最低"克隆中位 - 基线"最高"克隆中位（基线 IQR 为单位）
#               > 0 表示位移完全越出基线克隆异质性，预存克隆解释不了
#   cf_frac     组成反事实能解释的位移比例
#                 cf = Σ_c w_to,c·m_from,c - Σ_c w_from,c·m_from,c
#   icc_post    随访时点克隆间方差 / 总方差
# ==============================================================================

suppressPackageStartupMessages({
  library(Matrix); library(data.table)
})

PROC <- "data/proc/neoletexe"; TAB <- "results/tables"
THR_MAIN <- 0.045
THR_SENS <- c(0.035, 0.045, 0.055)
MIN_EPI   <- 50
MIN_CLONE <- 20
MIN_PAIR  <- 10

prof <- readRDS(file.path(PROC, "cnv_profile.rds"))
R0 <- prof$R; bins <- prof$bins; cells <- prof$cells
sc <- fread(file.path(TAB, "neoletexe_cell_scores.csv"))
sc <- sc[match(cells, cell)]
stopifnot(!anyNA(sc$cell))

## 每细胞减去自身中位后取染色体均值
Rc <- R0 - matrixStats::colMedians(R0)[rep(seq_len(ncol(R0)), each = nrow(R0))]
ch <- bins$chr_num; chs <- sort(unique(ch))
M <- sapply(chs, function(k) colMeans(Rc[ch == k, , drop = FALSE]))
colnames(M) <- paste0("chr", chs)
cat(sprintf("核型矩阵：%d 细胞 x %d 条染色体\n", nrow(M), ncol(M)))

make_clone <- function(thr) {
  B <- (M > thr) * 1L - (M < -thr) * 1L
  apply(B, 1, paste0, collapse = "")
}

PROGS <- c("ER_signal", "Proliferation", "IFN_response", "EMT_matrix",
           "Luminal_secretory")
PAIRS <- c("A>C", "A>B", "B>C")

run_layer <- function(pat, thr) {
  out <- list(); lay <- list()
  for (pt in sort(unique(sc$patient))) {
    i <- which(sc$patient == pt)
    tp_tab <- table(sc$tp[i])
    if (sum(tp_tab[c("A","B","C")] >= MIN_EPI, na.rm = TRUE) < 2) next
    pat_pt <- pat[i]
    tb <- table(pat_pt)
    keep <- names(tb)[tb >= MIN_CLONE]
    if (length(keep) < 2) next
    d <- data.table(cell = cells[i], patient = pt, tp = sc$tp[i],
                    clone = ifelse(pat_pt %in% keep, pat_pt, NA_character_))
    for (p in PROGS) d[[p]] <- sc[[p]][i]
    d <- d[!is.na(clone)]
    d[, clone := factor(clone)]

    for (pp in PAIRS) {
      fr <- strsplit(pp, ">")[[1]][1]; to <- strsplit(pp, ">")[[1]][2]
      a <- d[tp == fr]; b <- d[tp == to]
      if (nrow(a) < MIN_EPI || nrow(b) < MIN_EPI) next
      for (p in PROGS) {
        iqr <- IQR(a[[p]])
        if (!is.finite(iqr) || iqr <= 0) next
        pt_shift <- (median(b[[p]]) - median(a[[p]])) / iqr
        mb <- a[, .(mb = median(get(p)), nb = .N), by = clone]
        mo <- b[, .(mo = median(get(p)), no = .N), by = clone]
        cm <- merge(mb, mo, by = "clone")
        cm <- cm[nb >= MIN_PAIR & no >= MIN_PAIR]
        if (nrow(cm) < 2) next
        cm[, d_z := (mo - mb) / iqr]

        ## 组成反事实：克隆内分布不动，只让占比变
        wa <- a[, .N, by = clone]; wa[, w := N / sum(N)]
        wb <- b[, .N, by = clone]; wb[, w := N / sum(N)]
        mm <- merge(merge(cm[, .(clone, mb)], wa[, .(clone, w_from = w)], by = "clone"),
                    wb[, .(clone, w_to = w)], by = "clone")
        mm[is.na(w_to), w_to := 0]
        cf <- (sum(mm$w_to * mm$mb) - sum(mm$w_from * mm$mb)) / iqr

        b2 <- b[, .(v = get(p), clone)]
        gm <- mean(b2$v)
        mv <- b2[, .(m = mean(v), n = .N), by = clone]
        between <- sum(mv$n * (mv$m - gm)^2)
        total <- sum((b2$v - gm)^2)

        lay[[length(lay) + 1]] <- data.table(
          patient = pt, pair = pp, program = p, thr = thr,
          n_clone = nrow(cm), n_paired = nrow(cm),
          patient_shift = pt_shift, med_clone_shift = median(cm$d_z),
          dir_agree = mean(sign(cm$d_z) == sign(pt_shift)),
          all_same = all(sign(cm$d_z) == sign(pt_shift)),
          sep_z = (min(cm$mo) - max(cm$mb)) / iqr,
          base_spread = (max(cm$mb) - min(cm$mb)) / iqr,
          cf_shift = cf,
          cf_frac = if (abs(pt_shift) > 1e-6) cf / pt_shift else NA_real_,
          icc_post = if (total > 0) between / total else NA_real_)
        out[[length(out) + 1]] <- cm[, .(patient = pt, pair = pp, program = p,
                                         thr = thr, clone,
                                         n_from = nb, n_to = no,
                                         med_from = mb, med_to = mo, d_z = d_z,
                                         same_dir = sign(d_z) == sign(pt_shift))]
      }
    }
  }
  list(layer = rbindlist(lay), clone = rbindlist(out))
}

## ---------------- 主阈值 ----------------
pat_main <- make_clone(THR_MAIN)
res <- run_layer(pat_main, THR_MAIN)
fwrite(res$layer, file.path(TAB, "neoletexe_clone_layer.csv"))
fwrite(res$clone, file.path(TAB, "neoletexe_clone_cells.csv"))
lay <- res$layer

cat(sprintf("\n克隆 x 程序 x 时点对：%d 行；克隆层判定 %d 行；患者 %d 例\n",
            nrow(res$clone), nrow(lay), uniqueN(lay$patient)))

cat("\n=== 方向一致性（全部）===\n")
print(lay[, .(n = .N, med_agree = round(median(dir_agree), 3),
              frac_all_same = round(mean(all_same), 3),
              frac_ge_08 = round(mean(dir_agree >= 0.8), 3))])

cat("\n=== 只看患者层面实质位移（|shift|>0.5）===\n")
big <- lay[abs(patient_shift) > 0.5]
print(big[, .(n = .N, patients = uniqueN(patient),
              med_agree = round(median(dir_agree), 3),
              frac_all_same = round(mean(all_same), 3),
              sep_pos = sum(sep_z > 0),
              med_sep = round(median(sep_z), 2),
              med_cf_frac = round(median(cf_frac, na.rm = TRUE), 3),
              med_icc = round(median(icc_post, na.rm = TRUE), 2))])

cat("\n=== 组成反事实能解释多少位移 ===\n")
print(big[, .(program, patient, pair, shift = round(patient_shift, 2),
              cf = round(cf_shift, 2), cf_frac = round(cf_frac, 2),
              dir = round(dir_agree, 2), sep = round(sep_z, 2))][order(-abs(shift))][1:min(15, .N)])

## ---------------- 阈值敏感性 ----------------
sens <- rbindlist(lapply(THR_SENS, function(th) {
  r <- run_layer(make_clone(th), th)$layer
  r[abs(patient_shift) > 0.5,
    .(thr = th, n = .N, med_agree = round(median(dir_agree), 3),
      frac_all_same = round(mean(all_same), 3), sep_pos = sum(sep_z > 0),
      med_cf_frac = round(median(cf_frac, na.rm = TRUE), 3))]
}))
cat("\n=== 阈值敏感性 ===\n"); print(sens)
fwrite(sens, file.path(TAB, "neoletexe_clone_sensitivity.csv"))
cat("\n-> neoletexe_clone_layer.csv / neoletexe_clone_cells.csv / neoletexe_clone_sensitivity.csv\n")
