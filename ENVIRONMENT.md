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

No `renv.lock` is shipped, because the pipeline was run as a sequence of independent R
processes rather than against a single locked library. The equivalent record is produced by
the runner: `R/00_run_all.R` writes **`logs/sessionInfo_<tag>.txt`** after every step, so a
full run generates one `sessionInfo()` per analysis stage. Those files are the authoritative
record of what each figure was produced with.

## Conventions that must not drift

* Every quantile is R's `type = 7` (linear interpolation), set explicitly in each script.
* Seeds are fixed (`20261108`, `20261109`, `20261110`) at every resampling call site.
* Permutation p-values use `(1 + exceedances) / (1 + B)` with `B = 1000`; the smallest
  attainable value is `1/1001`.
* Scripts resolve paths relative to the repository root and must be run from there.
