# 04_diagnostics.R
# Block 4: diagnostics for the FEMALE by-region GAM (models/gam_female.rds)
#  1. Fitted curves + 95% bands over observed data, every region
#  2. Western Pacific sensitivity: does the curve change with basis size k?
#  3. Deviance residuals by year, per region
#  4. Residual autocorrelation (ACF, lag-1, Ljung-Box), per region
#  5. k.check table (is k big enough?) and concurvity
# Run from the project root (RStudio: setwd(...) in the Console, then source("R/04_diagnostics.R")).
# Reads data/clean and models/ only. Never touches data/raw/.

library(mgcv)

fig_dir <- file.path("outputs", "figures")
tab_dir <- file.path("outputs", "tables")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load ------------------------------------------------------------------
best <- readRDS(file.path("models", "gam_female.rds"))
regions <- levels(best$model$region)

d <- read.csv(file.path("data", "clean", "who_region.csv"), stringsAsFactors = FALSE)
d <- d[!is.na(d$rate) & d$sex == "FEMALE", ]
d$region <- factor(d$region_name, levels = regions)
years <- seq(min(d$year), max(d$year))

# Predictions on a full year grid (link scale -> back-transform)
newd <- expand.grid(region = factor(regions, levels = regions), year = years)
pr <- predict(best, newdata = newd, type = "link", se.fit = TRUE)
newd$fit <- exp(pr$fit)
newd$lo  <- exp(pr$fit - 1.96 * pr$se.fit)
newd$hi  <- exp(pr$fit + 1.96 * pr$se.fit)

draw_panel <- function(reg, nd, ob, main = reg, ...) {
  plot(ob$year, ob$rate, pch = 16, cex = 0.6, col = "grey35",
       ylim = range(c(nd$lo, nd$hi, ob$rate)),
       xlab = "Year", ylab = "Rate per 1,000", main = main, ...)
  polygon(c(nd$year, rev(nd$year)), c(nd$lo, rev(nd$hi)),
          col = adjustcolor("steelblue", 0.25), border = NA)
  lines(nd$year, nd$fit, col = "steelblue", lwd = 2)
  points(ob$year, ob$rate, pch = 16, cex = 0.6, col = "grey35")
}

# ---- 1. Fitted curves, all regions ----------------------------------------
png(file.path(fig_dir, "gam_fitted_by_region.png"), width = 3300, height = 2000, res = 300)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
for (r in regions) {
  draw_panel(r, newd[newd$region == r, ], d[d$region == r, ])
}
dev.off()

# ---- 2. Western Pacific sensitivity to k -----------------------------------
wp <- "Western Pacific"
ks <- c(10, 15, 20)
cols <- c("#1b9e77", "#d95f02", "#7570b3")
nd_wp <- newd[newd$region == wp, ]
ob_wp <- d[d$region == wp, ]

sens <- data.frame()
curves <- list()
for (i in seq_along(ks)) {
  m <- gam(rate ~ region + s(year, by = region, k = ks[i]),
           family = Gamma(link = "log"), data = d, method = "REML")
  st <- summary(m)$s.table
  edf_wp <- st[grep(wp, rownames(st), fixed = TRUE), "edf"]
  curves[[i]] <- exp(predict(m, newdata = nd_wp, type = "link"))
  sens <- rbind(sens, data.frame(k = ks[i], edf_western_pacific = unname(edf_wp), AIC = AIC(m)))
}
cat("\n== Western Pacific: sensitivity to basis size k ==\n"); print(sens)
write.csv(sens, file.path(tab_dir, "wp_k_sensitivity.csv"), row.names = FALSE)

png(file.path(fig_dir, "western_pacific_k_sensitivity.png"), width = 3300, height = 1500, res = 300)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
for (lg in c("", "y")) {
  plot(ob_wp$year, ob_wp$rate, pch = 16, cex = 0.7, col = "grey35", log = lg,
       xlab = "Year", ylab = "Rate per 1,000",
       main = paste("Western Pacific", if (lg == "y") "(log scale)" else "(raw scale)"))
  for (i in seq_along(ks)) lines(nd_wp$year, curves[[i]], col = cols[i], lwd = 2)
  legend("topright", bty = "n", lwd = 2, col = cols, cex = 0.8,
         legend = paste0("k = ", ks, " (edf ", round(sens$edf_western_pacific, 1), ")"))
}
dev.off()

# ---- 3. Residuals by year --------------------------------------------------
rd <- data.frame(region = best$model$region, year = best$model$year,
                 res = residuals(best, type = "deviance"))

png(file.path(fig_dir, "gam_residuals_by_year.png"), width = 3300, height = 2000, res = 300)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
for (r in regions) {
  x <- rd[rd$region == r, ]; x <- x[order(x$year), ]
  plot(x$year, x$res, type = "b", pch = 16, cex = 0.6, xlab = "Year",
       ylab = "Deviance residual", main = r)
  abline(h = 0, lty = 2, col = "grey50")
}
dev.off()

# ---- 4. Residual autocorrelation ------------------------------------------
png(file.path(fig_dir, "gam_residual_acf.png"), width = 3300, height = 2000, res = 300)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
ac <- data.frame()
for (r in regions) {
  x <- rd[rd$region == r, ]; x <- x[order(x$year), ]
  a <- acf(x$res, lag.max = 15, main = r)
  lb <- Box.test(x$res, lag = 5, type = "Ljung-Box")
  ac <- rbind(ac, data.frame(region = r, n = nrow(x),
                             lag1_autocorr = a$acf[2],
                             ljung_box_p_lag5 = lb$p.value))
}
dev.off()
cat("\n== Residual autocorrelation per region ==\n"); print(ac)
write.csv(ac, file.path(tab_dir, "residual_autocorrelation.csv"), row.names = FALSE)

# ---- 5. k.check and concurvity ---------------------------------------------
set.seed(1)
kc <- as.data.frame(k.check(best))
kc$smooth <- rownames(kc)
cat("\n== k.check (edf close to k', or low k-index with small p, suggests k too small) ==\n")
print(kc)
write.csv(kc, file.path(tab_dir, "k_check_female.csv"), row.names = FALSE)

cat("\n== Concurvity (worst case) ==\n")
cc <- tryCatch(concurvity(best, full = TRUE), error = function(e) NULL)
if (!is.null(cc)) print(round(cc["worst", ], 2))

# ---- Log -------------------------------------------------------------------
cat(paste0("\n- **", Sys.Date(), " - `R/04_diagnostics.R`:** Fitted-curve plots, Western Pacific ",
           "k-sensitivity (k = 10, 15, 20), residual plots, ACF/Ljung-Box, k.check, concurvity. ",
           "Outputs in outputs/figures/ and outputs/tables/.\n"),
    file = "README.md", append = TRUE)

cat("\nDone.\n")
