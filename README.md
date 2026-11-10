# directional-tail-excess

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24106083.svg)](https://doi.org/10.5281/zenodo.24106083)

Analysis code for the manuscript:

> **A distribution-level test separates subpopulation selection from a population-wide
> state transition** — Xu A, Chen C, Yin G, Luo J.

The repository contains the pipeline that produces every figure (Fig. 1–7, Fig. S1–S4) and
every number reported in the paper, from the processed single-cell objects of three
independent neoadjuvant cohorts.

Raw and processed sequencing data are **not** redistributed here; see [Data](#data).

---

## What the code computes

For every **patient × programme × interval**, cells are scored, quantiles of the score are
taken before and after treatment, and the post-treatment quantiles are divided by the
**baseline IQR of that same patient and programme**, so that all quantities are in units of
one baseline IQR and are comparable across patients and cohorts.

| Quantity | Definition | Used for |
|---|---|---|
| `z50` | (median_post − median_base) / IQR_base | whether the displacement is substantive; threshold **0.5** |
| `z25`, `z75` | same normalisation for the quartiles | shape |
| `asym` | `z75 − z25` | shape classification |
| `ex_lead` | (leading quartile − median) × sign(`z50`) | whether the leading end runs ahead |
| `ex_trail` | (trailing quartile − median) × sign(`z50`) | **the test statistic**: subpopulation selection predicts `ex_trail > 0` |
| `disp` | IQR_post / IQR_base | widening or narrowing; equals `1 + asym` only when the denominator is `q75 − q25` |

Both excesses are signed in the direction of the shift, so positive always means *that end
moved further than the median in the direction the distribution moved*. A **whole shift**
moves every quantile equally and gives `ex_lead = ex_trail = 0` by construction; a selection
event moves one end only and gives a positive excess at that end. Two thresholds are used
and are **not** interchangeable: **0.3** marks any whole shift, **0.5** marks a substantive
one.

Two further pieces sit on top of this:

* **selection counterfactual** (`R/39`) — resamples each combination's baseline cells under a
  pure-selection model and recomputes the trailing excess it would produce, giving a
  per-combination null that does not assume a functional form.
* **label permutation** (`R/38`) — shuffles the treatment label within patient to obtain the
  null distribution of `asym`, used for the `P ≤ 0.01` calls in Fig. 5.

Deciles are used only for the separate index `A = z10 + z90 − 2·z50`, which appears in the
permutation test and in Fig. 6C. `A` and `asym` are **not** interchangeable.

## Data

| | Cohort 1 (discovery) | Cohort 2 | Cohort 3 |
|---|---|---|---|
| Study | Fu et al., *Nat Commun* 2025 | Bassez et al., *Nat Med* 2021 | Altinönder et al., *Nat Genet* 2026 |
| Regimen | nab-paclitaxel + pembrolizumab, randomised 2-week single-agent lead-in | single dose anti-PD1 | aromatase inhibition |
| Access | dbGaP `phs002419.v2` (controlled); processed matrices at Zenodo `10.5281/zenodo.15799565` | EGA `EGAS00001004809` / `EGAD00001006608` (controlled) | EGA `EGAS50000001362` (controlled) |
| Cell-level tables expected by this code | `results/tables/nc2025_*.csv` | `results/tables/bassez_*.csv` | `results/tables/neoletexe_*.csv` |

None of the three cohorts can be redistributed. Access is by application to the original
data access committees; the code here starts from the processed objects that those
applications provide. Paths are relative (`results/tables`, `results/figures`), so the
pipeline runs from the repository root without editing any absolute path. Most parameters
can be overridden by environment variable (e.g. `WB_SEED`, `WB_SC_RES_DIR`).

## Requirements

* **R 4.4.2** — Seurat 5.3.0, SeuratObject 5.1.0, scDblFinder 1.20.2, RcppML 0.3.7.1,
  Matrix 1.7-1, data.table 1.16.4, matrixStats 1.5.0, ggplot2 4.0.3, patchwork 1.3.0.
  Full list in [`ENVIRONMENT.md`](ENVIRONMENT.md).
* **Python 3** — standard library only; `pandas` and `numpy` are used by one script.

Every quantile in this project is R's `type = 7` (linear interpolation), set explicitly in
each script so the convention cannot drift between cohorts.

## Running

```bash
Rscript R/00_run_all.R --dry-run    # print the 23-step plan, execute nothing
Rscript R/00_run_all.R              # run the plan in dependency order
Rscript R/00_run_all.R --from=38    # resume from a given step
```

The runner executes each script in its own R process, writes `logs/run_<tag>.log` and
`logs/sessionInfo_<tag>.txt`, and stops at the first failure with the `--from=` resume hint.

Two things are worth knowing before the first run:

1. **The Python layer runs first.** `R/00_run_all.R` starts at `R/21` and assumes the
   cell-level tables listed above already exist. Those are produced by `py/`; see
   [`py/README.md`](py/README.md) for the order.
2. **`R/28` and `R/39` are mutually dependent.** `R/28` draws Fig. 5D from the per-combination
   counterfactual that `R/39` writes, and `R/39` reads the shape calls that `R/28` writes.
   With the tables present, one pass is enough; from an empty `results/tables`, the first
   `R/28` run stops at the counterfactual read — this is expected, and the runner places a
   second `R/28` at the end of the plan so that `--from=38` completes it.

## Code map

| Script | Produces |
|---|---|
| `py/01_fetch_nc2025.py` | download and unpack cohort 1 |
| `py/03_l1_l3_discriminate.py` | lineage-marker separation tables |
| `py/11_shape_discriminate.py` | `nc2025_v3_shape.csv`, `_density.csv`, `_dominant.csv` |
| `py/13_clone_layering_all.py` | `nc2025_v4_clone_{cells,layer,summary,verdict}.csv` |
| `py/14_counterfactual_composition.py` | composition counterfactual |
| `py/15_fig1_cohort.py` | `nc2025_fig1_cohort.csv`, RCB labels |
| `py/16_fig2_bulk.py`, `py/16b_sim_vs_ks.py`, `py/17_batch_sensitivity.py` | Fig. 2 inputs |
| `R/21_bassez_validation.R` | cohort 2 shape table |
| `R/23`, `R/30`–`R/33` | cohort 3 loading, QC gate, CNV binning, shape, clones |
| `R/24` | Fig. 1 |
| `R/25` | Fig. 2 |
| `R/27` | Fig. 3 |
| `R/26` | Fig. 4 |
| `R/28` | Fig. 5 |
| `R/29` | Fig. 6 |
| `R/34` | Fig. 7 |
| `R/35`–`R/37` | Fig. S1–S3 |
| `R/38_permutation_shape.R` | label permutation null |
| `R/39_selection_counterfactual.R` | per-combination selection counterfactual |
| `R/40_shape_calibration.R` | Fig. S4 calibration |
| `R/41_depth_matched_downsample.R` | depth-matched downsampling (Supplementary Table 2) |
| `R/42_shape_stat_correlation.R` | correlation among the four shape quantities (Supplementary Table 3) |

`R/01`–`R/11` and `R/18`–`R/22` belong to an earlier pilot reanalysis of GSE176078 that is
**superseded** and produces none of the figures in this paper. `config/00_config.R` is used
only by that layer. They are kept for lineage; nothing in the manuscript depends on them.

## Reproducibility notes

* Random seeds are fixed (`20261108`, `20261109`, `20261110`) and every resampling call site
  sets one; the runner records a `sessionInfo()` per step.
* Permutation p-values use the conservative `(1 + exceedances) / (1 + B)` form with
  `B = 1000`, so the smallest attainable value is `1/1001` and no p-value is ever reported
  as zero.
* **Cohort 2 denominator.** `bassez_v2_profile.csv` normalises every `z` value by the column
  `iqr`, which is *not* equal to the baseline `q75 − q25` (ratio 0.48–3.28, median 1.22).
  Signs, directions and patient-level agreement are unaffected by a positive rescaling, and
  the column is retained so that the deposited tables match the published figures. The
  consequences are measurable and stated in the paper: `disp = 1 + asym` does not hold row-wise
  for cohort 2, and recomputing with the standard IQR would move the `|z50| > 0.5` count from
  38/145 to 29/145. Changing it requires editing three scripts in step and redrawing Fig. 6–7.
* Cohort 3's CNV caller is written in base R (`R/31_neoletexe_cnv.R`) from the arm-level
  procedure rather than taken from a package, so its behaviour is documented in Methods
  rather than by a version number.

## Citation

Citation metadata is in [`CITATION.cff`](CITATION.cff). If you use this code, please cite the
manuscript and the archived release.

**This version (v1.0.0)** — use this when you need to refer to exactly the code behind the
manuscript:

> Xu A, Chen C, Yin G, Luo J. *directional-tail-excess: a distribution-level test separating
> subpopulation selection from a population-wide state transition* (v1.0.0). Zenodo. 2026.
> doi:10.5281/zenodo.24106083

**Concept DOI** — `10.5281/zenodo.24106082`. This always resolves to the latest release, so use
it in a "code availability" statement that should stay valid as the repository is updated.

Archived snapshot: <https://zenodo.org/records/24106083>

## License

MIT — see [`LICENSE`](LICENSE).
