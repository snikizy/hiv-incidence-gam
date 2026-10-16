# 05_results.R
# Block 5: headline results, fitted only to observations ABOVE the reporting floor.
# Why: many series are published as exactly 0.10 per 1,000 (a floor), which is not a
# real measurement. Fitting a smooth through it produces artefact wiggles (see
# Western Pacific in 04_diagnostics.R). Here each eligible region x sex series gets its
# own Gamma(log) GAM on the above-floor years only.
#
# Eligibility: >= 15 years above the floor AND max/min >= 2 among those years
# (rates are rounded to 2 decimals, so a narrow range just above 0.10 is mostly rounding).
#
# Outputs: results_summary.csv, results_derivatives.csv, excluded_series_floor.csv,
#          results_curves_<sex>.png
# Run from the project root (RStudio: setwd(...) in the Console, then source("R/05_results.R")).
# Reads data/clean only. Never touches data/raw/.

library(mgcv)
set.seed(2026)

floor_val <- 0.10
min_obs   <- 15
min_ratio <- 2
n_sim     <- 2000

fig_dir <- file.path("outputs", "figures")
tab_dir <- file.path("outputs", "tables")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

d <- read.csv(file.path("data", "clean", "who_region.csv"), stringsAsFactors = FALSE)
d <- d[!is.na(d$rate), ]
combos <- unique(d[, c("sex", "region_name")])
combos <- combos[order(combos$sex, combos$region_name), ]

summ <- list(); deriv <- list(); excl <- list(); res <- list()

for (i in seq_len(nrow(combos))) {
  sx <- combos$sex[i]; rg <- combos$region_name[i]
  x  <- d[d$sex == sx & d$region_name == rg, ]
  x  <- x[order(x$year), ]
  xa <- x[x$rate > floor_val + 1e-9, ]
  n_above <- nrow(xa)
  ratio <- if (n_above > 0) max(xa$rate) / min(xa$rate) else NA

  if (n_above < min_obs || is.na(ratio) || ratio < min_ratio) {
    why <- if (n_above < min_obs) paste0("fewer than ", min_obs, " years above floor")
           else "range above floor under 2-fold (mostly rounding)"
    excl[[length(excl) + 1]] <- data.frame(sex = sx, region = rg, n_total = nrow(x),
                                           n_above_floor = n_above, range_ratio = round(ratio, 2),
                                           reason = why)
    next
  }

  # ---- Fit -----------------------------------------------------------------
  k <- min(10, floor(n_above / 2))
  m <- gam(rate ~ s(year, k = k), family = Gamma(link = "log"), data = xa, method = "REML")
  grid <- data.frame(year = seq(min(xa$year), max(xa$year)))
  Xp   <- predict(m, newdata = grid, type = "lpmatrix")
  beta <- coef(m); Vb <- vcov(m)
  fit_link <- as.numeric(Xp %*% beta)
  S <- Xp %*% t(rmvn(n_sim, beta, Vb))           # grid x n_sim, link scale
  last <- nrow(S)

  # ---- Peak year (with simulation interval) ---------------------------------
  pk_pt   <- grid$year[which.max(fit_link)]
  pk_draw <- grid$year[apply(S, 2, which.max)]
  pk_q    <- quantile(pk_draw, c(0.5, 0.025, 0.975), type = 1)

  # ---- Percent change, peak -> last year above floor ------------------------
  pct   <- 100 * (exp(S[last, ] - apply(S, 2, max)) - 1)
  pct_q <- quantile(pct, c(0.5, 0.025, 0.975))

  # ---- Year-to-year change in log rate (derivative) -------------------------
  D     <- S[-1, , drop = FALSE] - S[-last, , drop = FALSE]
  d_med <- apply(D, 1, median)
  d_lo  <- apply(D, 1, quantile, 0.025)
  d_hi  <- apply(D, 1, quantile, 0.975)
  cls   <- ifelse(d_hi < 0, "falling", ifelse(d_lo > 0, "rising", "not distinguishable"))
  deriv[[length(deriv) + 1]] <- data.frame(
    sex = sx, region = rg,
    from_year = grid$year[-last], to_year = grid$year[-1],
    pct_change_per_year = round(100 * (exp(d_med) - 1), 2),
    pct_lo = round(100 * (exp(d_lo) - 1), 2),
    pct_hi = round(100 * (exp(d_hi) - 1), 2),
    status = cls)

  yrs_after <- grid$year[last] - pk_pt
  avg_ann <- if (yrs_after > 0) 100 * (exp((fit_link[last] - max(fit_link)) / yrs_after) - 1) else NA

  summ[[length(summ) + 1]] <- data.frame(
    sex = sx, region = rg, n_years_above_floor = n_above,
    first_year = min(xa$year), last_year_above_floor = max(xa$year), k = k,
    edf = round(summary(m)$s.table[1, "edf"], 2),
    dev_explained = round(summary(m)$dev.expl, 3),
    peak_year_fit = pk_pt, peak_year_median = pk_q[1], peak_year_lo = pk_q[2], peak_year_hi = pk_q[3],
    peak_at_boundary = pk_pt == min(grid$year),
    peak_rate_fit = round(exp(max(fit_link)), 3),
    last_rate_fit = round(exp(fit_link[last]), 3),
    pct_change_peak_to_last = round(pct_q[1], 1),
    pct_lo = round(pct_q[2], 1), pct_hi = round(pct_q[3], 1),
    avg_annual_pct_after_peak = round(avg_ann, 1),
    first_year_falling = if (any(cls == "falling")) min(grid$year[-last][cls == "falling"]) else NA,
    n_intervals_falling = sum(cls == "falling"))

  res[[length(res) + 1]] <- list(sex = sx, region = rg, x = x, xa = xa, grid = grid,
                                 fit = exp(fit_link),
                                 lo = exp(apply(S, 1, quantile, 0.025)),
                                 hi = exp(apply(S, 1, quantile, 0.975)),
                                 pk = pk_q[1])
}

summ  <- do.call(rbind, summ)
deriv <- do.call(rbind, deriv)
excl  <- if (length(excl)) do.call(rbind, excl) else data.frame()

write.csv(summ,  file.path(tab_dir, "results_summary.csv"),        row.names = FALSE)
write.csv(deriv, file.path(tab_dir, "results_derivatives.csv"),    row.names = FALSE)
write.csv(excl,  file.path(tab_dir, "excluded_series_floor.csv"),  row.names = FALSE)

cat("\n== Modelled series (above-floor years only) ==\n"); print(summ, row.names = FALSE)
cat("\n== Not modelled (floor-dominated or too narrow) ==\n"); print(excl, row.names = FALSE)

# ---- Figures: one per sex, one panel per modelled region --------------------
for (sx in unique(summ$sex)) {
  items <- Filter(function(z) z$sex == sx, res)
  png(file.path(fig_dir, paste0("results_curves_", tolower(sx), ".png")),
      width = 1700 * length(items), height = 1600, res = 300)
  par(mfrow = c(1, length(items)), mar = c(4, 4, 3, 1))
  for (it in items) {
    plot(it$x$year, it$x$rate, type = "n", log = "y",
         ylim = range(c(it$lo, it$hi, it$x$rate)),
         xlab = "Year", ylab = "Rate per 1,000 (log scale)",
         main = paste0(it$region, " - ", tolower(sx)))
    polygon(c(it$grid$year, rev(it$grid$year)), c(it$lo, rev(it$hi)),
            col = adjustcolor("steelblue", 0.25), border = NA)
    lines(it$grid$year, it$fit, col = "steelblue", lwd = 2)
    points(it$xa$year, it$xa$rate, pch = 16, cex = 0.7, col = "grey30")
    fl <- it$x[it$x$rate <= floor_val + 1e-9, ]
    if (nrow(fl) > 0) points(fl$year, fl$rate, pch = 1, cex = 0.7, col = "grey60")
    abline(h = floor_val, lty = 3, col = "grey50")
    abline(v = it$pk, lty = 2, col = "firebrick")
    legend("topright", bty = "n", cex = 0.7,
           legend = c("fitted curve + 95% band", "observed (above floor)",
                      "at reporting floor (not fitted)", "estimated peak year"),
           col = c("steelblue", "grey30", "grey60", "firebrick"),
           lty = c(1, NA, NA, 2), pch = c(NA, 16, 1, NA), lwd = c(2, NA, NA, 1))
  }
  dev.off()
}

# ---- Log -------------------------------------------------------------------
cat(paste0("\n- **", Sys.Date(), " - `R/05_results.R`:** Reporting floor of 0.10 per 1,000 found ",
           "(many series are published as exactly 0.10). Models refit on above-floor years only, ",
           "one Gamma(log) GAM per eligible region x sex (>= ", min_obs, " years above floor and ",
           "max/min >= ", min_ratio, "). Peak year, % change and falling periods from ", n_sim,
           " posterior simulations. Earlier by-region models (R/03) are superseded for headline results.\n"),
    file = "README.md", append = TRUE)

cat("\nDone.\n")
