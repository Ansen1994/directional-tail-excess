## =============================================================================
## 01b_attach_labels.R - patch: attach the original curated labels to the
##                       already-written QC artefacts
##
## Why this exists: the first version of 01_qc.R did not carry the curated
## labels from metadata.csv into the meta table. This script reads only the 26
## metadata.csv files (milliseconds), never the 2.6 GB count matrices, so it can
## be re-run at any time after 01_qc.R without repeating the heavy step.
##
## Output: results/tables/labels_reference.csv   curated labels for every cell
##         data/processed/01_counts_qc.rds       columns added in place (backed up)
##
## NOTE: keep this file pure ASCII (see 02_integrate.R header).
## =============================================================================

source("config/00_config.R")
t0 <- Sys.time()

dirs <- sort(list.dirs(RAW_DIR, recursive = FALSE, full.names = TRUE))
dirs <- dirs[file.exists(file.path(dirs, "metadata.csv"))]

rows <- lapply(dirs, function(d) {
  md <- read.csv(file.path(d, "metadata.csv"), row.names = 1,
                 stringsAsFactors = FALSE, check.names = FALSE)
  get_col <- function(cl)
    if (cl %in% colnames(md)) as.character(md[[cl]]) else NA_character_
  data.frame(cell_id = rownames(md), sample = basename(d),
             ref_subtype = get_col("subtype"),
             ref_major   = get_col("celltype_major"),
             ref_minor   = get_col("celltype_minor"),
             ref_subset  = get_col("celltype_subset"),
             stringsAsFactors = FALSE)
})
lab <- do.call(rbind, rows); rownames(lab) <- NULL
msg("read %d samples, %d cells with curated labels", length(dirs), nrow(lab))
msg("missing labels: major %d / minor %d / subtype %d",
    sum(is.na(lab$ref_major)), sum(is.na(lab$ref_minor)), sum(is.na(lab$ref_subtype)))
write_csv_safe(lab, file.path(TAB_DIR, "labels_reference.csv"))

## ------------------------------------------- sync into 01_counts_qc.rds
rds <- file.path(PROC_DIR, "01_counts_qc.rds")
if (file.exists(rds)) {
  inp <- readRDS(rds)
  need <- c("ref_subtype", "ref_major", "ref_minor", "ref_subset")
  if (!all(need %in% colnames(inp$meta))) {
    i <- match(rownames(inp$meta), lab$cell_id)
    if (anyNA(i)) warning(sum(is.na(i)), " cells without a matched label")
    for (cl in need) inp$meta[[cl]] <- lab[[cl]][i]
    bak <- file.path(PROC_DIR, "01_counts_qc.bak_before_labels.rds")
    if (!file.exists(bak)) file.copy(rds, bak)
    saveRDS(inp, rds, compress = "gzip")
    msg("columns added, %s overwritten (backup: %s)", basename(rds), basename(bak))
  } else {
    msg("01_counts_qc.rds already has label columns, skipped")
  }
  tb <- table(inp$meta$ref_major)
  msg("label distribution among QC-passed cells: %s",
      paste(sprintf("%s=%d", names(tb), as.integer(tb)), collapse = "; "))
}

writeLines(capture.output(sessionInfo()), file.path(LOG_DIR, "sessionInfo_01b_labels.txt"))
msg("done in %.2f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))
