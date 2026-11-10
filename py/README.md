# Python layer

This layer produces the cohort-1 cell-level and combination-level tables that the R scripts
in `../R/` consume. It uses only the standard library apart from `pandas` and `numpy` in one
script, and it expects to be run from the **repository root** (paths are relative to
`results/tables`).

Raw cohort-1 data are controlled-access (dbGaP `phs002419.v2`). `01_fetch_nc2025.py`
assumes the processed objects released with the original publication are already present;
nothing here downloads or redistributes them.

## Order

```
python py/01_fetch_nc2025.py                 # unpack the processed objects
python py/03_l1_l3_discriminate.py           # nc2025_v2_*.csv (lineage / clone layers)
python py/11_shape_discriminate.py           # nc2025_v3_shape.csv, _density.csv, _dominant.csv
python py/13_clone_layering_all.py           # nc2025_v4_clone_{cells,layer,summary,verdict}.csv
python py/14_counterfactual_composition.py   # nc2025_v4_composition_counterfactual.csv
python py/15_fig1_cohort.py                  # nc2025_fig1_cohort.csv, nc2025_rcb_labels.csv
python py/16_fig2_bulk.py                    # nc2025_fig2_cells.csv, _icc.csv, _pair_ks.csv, ...
python py/16b_sim_vs_ks.py                   # nc2025_fig2_pseudobulk_sim.csv, _sim_vs_ks.csv
python py/17_batch_sensitivity.py            # nc2025_fig2_batch_{delta,effect,map}.csv
```

`02`, `04`–`10` and `12` are exploratory probes (metaprogram identity, overlap, confounding
checks) that produce no table used by the manuscript and are not required for the figures.
`07` recomputes the lineage sensitivity of `03` and is only needed if the trajectory
thresholds are changed.

## Notes

* `11_shape_discriminate.py` is where the shape quantities are first computed; the R scripts
  later re-derive them from cell-level scores. Both use R's `type = 7` quantile convention.
* `13_clone_layering_all.py` writes the six-patient clone layer (`n = 6`, 24,166 cells) that
  Fig. 1E and the clone funnel refer to. Patient P12 has no subclone at baseline and is
  excluded at this step, not later.
