# Non-Linear Trajectories of HIV Incidence in Africa and South-East Asia, 1990-2024

GAM analysis of WHO regional HIV incidence estimates. Full write-up: see `REPORT.md`.

**Headline:** female incidence peaked in Africa in 1993-94 (about 4.7 per 1,000) and fell about 86% by 2024;
in South-East Asia it peaked in 1996 and fell about 77% to the reporting floor by 2008.
Four regions sit at or near a 0.10 reporting floor and were not modelled.

## How to run
Run from the project root, in order: `R/01_clean.R` to `R/05_results.R`.
Requires R with `mgcv` and `ggplot2` (plus their dependencies). Raw data in `data/raw/` is never modified.
Each script appends one line to this log when it runs, so delete duplicate lines after re-running.

## Final figures
`outputs/figures/results_curves_*.png` (from `R/05_results.R`). The `gam_*` figures come from the
earlier all-years model and are kept only to document the reporting-floor problem.

## Reporting floor
Many series are published as exactly 0.10 per 1,000 (apparently a reporting floor). Fitting a smooth
through them produces artefact curves (see Western Pacific in `outputs/figures/western_pacific_k_sensitivity.png`).
`R/05_results.R` therefore fits only years above the floor. Scripts `03` and `04` are kept as the record of how
the problem was found.

## Decision log

- **2026-10-15 - `R/01_clean.R`:** Split `data/raw/RELAY_WHS.csv` on `DIM_GEO_CODE_TYPE` into country (`COUNTRY`), WHO region (`WHOREGION`), and global (`GLOBAL`) tables. Kept year (`DIM_TIME`), region name (`GEO_NAME_SHORT`), sex, rate, and lower/upper bounds (`RATE_PER_1000_N/NL/NU`). Added `rate_na` rather than dropping missing rates so later models can subset explicitly. Wrote CSVs to `data/clean/`.
- **2026-10-15 - `R/02_eda.R`:** WHO-region EDA uses the observed year range 1990-2024 as the complete grid (6 regions x 3 sexes x 35 years = 630 cells). Missing region-sex-year cells are left as NA so ribbons and lines gap instead of interpolating. Peak year is the earliest year at the maximum rate; `n_years_at_peak` records plateaus. Linear and log10 figures saved at 300 dpi; peak table in `outputs/tables/who_region_peaks.csv`.
- **2026-10-16 - `R/03_models.R`:** Gamma(log) GAMs, mgcv REML. Primary analysis: FEMALE series (complete for all 6 regions). Compared 3 models by AIC; best = m3_by_region_k20 (`models/gam_female.rds`). Secondary MALE/TOTAL fitted only for series with >= 15 obs. Excluded: Eastern Mediterranean (TOTAL, n = 3); Western Pacific (TOTAL, n = 4). Superseded for headline results by `R/05_results.R` (reporting floor, see above).
- **2026-10-16 - `R/04_diagnostics.R`:** Fitted-curve plots, Western Pacific k-sensitivity (k = 10, 15, 20), residual plots, ACF/Ljung-Box, k.check, concurvity. Showed Western Pacific (female) is a flat series at the 0.10 floor that the model was over-fitting. Outputs in `outputs/figures/` and `outputs/tables/`.
- **2026-10-16 - `R/05_results.R`:** Reporting floor of 0.10 per 1,000 found (many series are published as exactly 0.10). Models refit on above-floor years only, one Gamma(log) GAM per eligible region x sex (>= 15 years above floor and max/min >= 2). Peak year, % change and falling periods from 2000 posterior simulations. Earlier by-region models (`R/03`) are superseded for headline results.
