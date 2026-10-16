# Non-Linear Trajectories of HIV Incidence in Africa and South-East Asia, 1990-2024
*A Generalized Additive Model analysis of WHO regional estimates*

## Question
How has HIV incidence changed over time across WHO regions, and is the change linear or does it follow a non-linear trajectory (rise, peak, decline)?

## Data
- WHO Global Health Observatory indicator `SDGHIV` (HIV infections per 1,000 population), annual, 1990-2024, by WHO region and sex (female, male, total).
- Each value comes with published lower and upper bounds. These are modelled estimates, not raw surveillance counts.

## A data problem found before modelling
Many series are published as exactly **0.10 per 1,000**, which appears to be a reporting floor (values below 0.1 shown as 0.10; confirm against the indicator metadata). A first model fitted to all values produced artefact curves: Western Pacific (female) is 0.01 for 1990-92 and exactly 0.10 for 1993-2024, yet the model used 18.4 effective degrees of freedom to chase it.

Years at the floor, female series (35 years each): Eastern Mediterranean 35, Western Pacific 32 (plus 3 at 0.01), Europe 30, Americas 20, South-East Asia 16, Africa 0.

## Method
- Models were fitted **only to years above the floor**, one Gamma-family GAM with log link per region and sex (`mgcv`, penalised smooth `s(year)`, REML).
- A series was modelled only if it had at least 15 years above the floor and a max/min ratio of at least 2 (rates are rounded to 2 decimals, so a narrow range just above 0.10 is mostly rounding).
- Peak year, percent change from peak, and the years in which the curve is falling were taken from 2,000 simulations of the fitted smooth.

## Results (female, primary analysis)
| Region | Peak year | Peak rate | Last year above floor | Rate then | Change from peak | Avg change per year | Falling since |
|---|---|---|---|---|---|---|---|
| Africa | 1994 | 4.71 | 2024 | 0.67 | -85.7% | -6.3% | 1994 (every year) |
| South-East Asia | 1996 | 0.48 | 2008 | 0.11 | -77.3% | -11.6% | 1996 (every year) |

**Like-for-like comparison, 1996-2008 (published values):** Africa fell from 4.50 to 2.30 (-49%, about -5.4% a year). South-East Asia fell from 0.50 to 0.11 (-78%, about -11.9% a year). South-East Asia declined roughly twice as fast as Africa over the same window, from a far lower level.

**Shape:** both regions show a non-linear rise, peak, and decline. Africa's decline is long and steady, and on the log scale appears to steepen after about 2015. South-East Asia's is a sharp hump that reaches the reporting floor by 2008, so nothing can be said about what happens after.

## Secondary results (male and total)
- Africa: male peak 1994, -89.1% to 2024; total peak 1994, -87.2%.
- South-East Asia: total peak 1996, -80.2% to 2010. Male series was already falling at the start of observation (peak at or before 1990), -83% to 2012.
- Europe (male): a different pattern, peaking late (2013, 0.26) and falling only slightly (-7.8% to 2024).

## Not modelled
Americas, Eastern Mediterranean and Western Pacific (all sexes, except where noted in `outputs/tables/excluded_series_floor.csv`), Europe (female and total): at or near the reporting floor for most or all years, or too narrow a range above it. These are reported descriptively only.

## Limitations
1. The data are modelled estimates, already smooth, so the GAM reproduces them almost exactly (deviance explained about 1.0). The intervals describe the smooth's fit, not real-world sampling uncertainty, and are much narrower than the true uncertainty.
2. The 0.10 floor hides the later trajectory of every region that reaches it, including South-East Asia after 2008.
3. Percent changes are measured over different windows for different regions and are not directly comparable except in the common-window comparison above.
4. Regional aggregates can hide large differences between countries.
5. No covariates were included, so these results describe trends and do not explain them (no causal claims about treatment, prevention, or other drivers).
6. The eligibility thresholds (15 years, 2-fold range) are judgement calls, documented here and adjustable in `R/05_results.R`.

## Reproducibility
`R/01_clean.R` to `R/05_results.R`, run in order from the project root. Raw data in `data/raw/` is never modified. Scripts `03` and `04` document how the floor problem was found; `05` produces the headline results.
