# Non-Linear Trajectories of HIV Incidence Across WHO Regions

## Decision log

- **2026-10-15 — `R/01_clean.R`:** Split `data/raw/RELAY_WHS.csv` on `DIM_GEO_CODE_TYPE` into country (`COUNTRY`), WHO region (`WHOREGION`), and global (`GLOBAL`) tables. Kept year (`DIM_TIME`), region name (`GEO_NAME_SHORT`), sex, rate, and lower/upper bounds (`RATE_PER_1000_N/NL/NU`). Added `rate_na` rather than dropping missing rates so later models can subset explicitly. Wrote CSVs to `data/clean/`.
- **2026-10-15 — `R/02_eda.R`:** WHO-region EDA uses the observed year range 1990–2024 as the complete grid (6 regions × 3 sexes × 35 years = 630 cells). Missing region–sex–year cells are left as NA so ribbons and lines gap instead of interpolating. Peak year is the earliest year at the maximum rate; `n_years_at_peak` records plateaus. Linear and log10 figures saved at 300 dpi; peak table in `outputs/tables/who_region_peaks.csv`.

- **2026-10-16 - `R/03_models.R`:** Gamma(log) GAMs, mgcv REML. Primary analysis: FEMALE series (complete for all 6 regions). Compared 3 models by AIC; best = m3_by_region_k20 (models/gam_female.rds). Secondary MALE/TOTAL fitted only for series with >= 15 obs. Excluded: Eastern Mediterranean (TOTAL, n = 3); Western Pacific (TOTAL, n = 4).

- **2026-10-16 - `R/03_models.R`:** Gamma(log) GAMs, mgcv REML. Primary analysis: FEMALE series (complete for all 6 regions). Compared 3 models by AIC; best = m3_by_region_k20 (models/gam_female.rds). Secondary MALE/TOTAL fitted only for series with >= 15 obs. Excluded: Eastern Mediterranean (TOTAL, n = 3); Western Pacific (TOTAL, n = 4).

- **2026-10-16 - `R/03_models.R`:** Gamma(log) GAMs, mgcv REML. Primary analysis: FEMALE series (complete for all 6 regions). Compared 3 models by AIC; best = m3_by_region_k20 (models/gam_female.rds). Secondary MALE/TOTAL fitted only for series with >= 15 obs. Excluded: Eastern Mediterranean (TOTAL, n = 3); Western Pacific (TOTAL, n = 4).
