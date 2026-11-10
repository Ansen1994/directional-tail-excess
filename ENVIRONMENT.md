# Environment

Recorded from the `sessionInfo()` output of the runs that produced the manuscript's figures.

## R

```
R version 4.4.2 (2024-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64
```

## Packages

| Package | Version | Role |
|---|---|---|
| Seurat | 5.3.0 | single-cell objects; `AddModuleScore` for the five literature gene sets |
| SeuratObject | 5.1.0 | object class |
| scDblFinder | 1.20.2 | doublet detection |
| RcppML | 0.3.7.1 | per-patient NMF (cohort 1 metaprogram reproduction) |
| Matrix | 1.7-1 | sparse matrices |
| data.table | 1.16.4 | tabulation |
| matrixStats | 1.5.0 | quantile summaries |
| ggplot2 | 4.0.3 | figures |
| patchwork | 1.3.0 | figure assembly |
| SingleR | 2.8.0 | cell-type labels (cohort 3) |
| scran | 1.34.0 | normalisation helpers |
| scater | 1.34.1 | QC metrics |
| BiocParallel | 1.40.0 | parallel CNV binning |
| harmony | 1.2.3 | batch integration (cohort 3) |
| dplyr | 1.2.1 | data manipulation |
| RColorBrewer | 1.1-3 | figure palettes |

## Python

Python 3. The scripts import only the standard library, apart from `pandas` and `numpy`
in `py/16_fig2_bulk.py`.

## Dependency record

No `renv.lock` is shipped. The pipeline was run as a sequence of independent R processes
rather than against a single locked library, so a lockfile generated after the fact from the
current global library would not be an honest record of the versions that produced the
figures; we do not fabricate one. The versions tabulated above were transcribed from the
`sessionInfo()` outputs of those runs and match the versions stated in the manuscript.

**Shipped records.** `logs/sessionInfo_<tag>.txt` holds the nine `sessionInfo()` files that
were actually written during the early preprocessing stages (tags `01_qc`, `01b_labels`,
`03_annotate`, `04_doublets`, `05_cnv`, `06_epi`, `07_nmf`, `08_diag`, `09_ps_nmf`).
Run logs (`*.log`) are not shipped.

**Not shipped.** Steps `21`–`42` (shape, counterfactual, clone-layer and figure analyses)
were run as separate R processes at different times, and their per-step `sessionInfo()`
files were not archived. A full reproduction generates them: `R/00_run_all.R` writes
`logs/sessionInfo_<tag>.txt` and `logs/run_<tag>.log` after each step. The cell-level input
tables needed to start that run are likewise not shipped (see [README](README.md#data)).

## Conventions that must not drift

* Every quantile is R's `type = 7` (linear interpolation), set explicitly in each script.
* Seeds are fixed (`20261108`, `20261109`, `20261110`) at every resampling call site.
* Permutation p-values use `(1 + exceedances) / (1 + B)` with `B = 1000`; the smallest
  attainable value is `1/1001`.
* Scripts resolve paths relative to the repository root and must be run from there.
