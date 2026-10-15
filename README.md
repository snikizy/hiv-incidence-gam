# Non-Linear Trajectories of HIV Incidence Across WHO Regions

## Decision log

- **2026-10-15 — `R/01_clean.R`:** Split `data/raw/RELAY_WHS.csv` on `DIM_GEO_CODE_TYPE` into country (`COUNTRY`), WHO region (`WHOREGION`), and global (`GLOBAL`) tables. Kept year (`DIM_TIME`), region name (`GEO_NAME_SHORT`), sex, rate, and lower/upper bounds (`RATE_PER_1000_N/NL/NU`). Added `rate_na` rather than dropping missing rates so later models can subset explicitly. Wrote CSVs to `data/clean/`.
