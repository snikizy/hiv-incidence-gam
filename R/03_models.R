# 03_models.R
# Block 3: GAMs of HIV incidence by WHO region, 1990-2024
# Primary analysis: FEMALE series. Secondary: MALE and TOTAL (series with >= 15 obs).
# Run from the project root:  Rscript R/03_models.R
# Reads data/clean/who_region.csv only. Never touches data/raw/.

library(mgcv)

in_file   <- file.path("data", "clean", "who_region.csv")
model_dir <- "models"
tab_dir   <- file.path("outputs", "tables")
dir.create(model_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir,   showWarnings = FALSE, recursive = TRUE)

# ---- Load ------------------------------------------------------------------
d <- read.csv(in_file, stringsAsFactors = FALSE)
stopifnot(all(c("region_name", "sex", "year", "rate") %in% names(d)))
d$region <- factor(d$region_name)
d <- d[!is.na(d$rate), ]
stopifnot(all(d$rate > 0))          # Gamma family needs strictly positive values
fam <- Gamma(link = "log")

# ---- Primary analysis: FEMALE ----------------------------------------------
f <- d[d$sex == "FEMALE", ]
cat("FEMALE rows:", nrow(f), "| regions:", nlevels(droplevels(f$region)), "\n")

m1 <- gam(rate ~ region + s(year),
          family = fam, data = f, method = "REML")
m2 <- gam(rate ~ region + s(year, by = region),
          family = fam, data = f, method = "REML")
m3 <- gam(rate ~ region + s(year, by = region, k = 20),
          family = fam, data = f, method = "REML")

models <- list(m1_shared = m1, m2_by_region = m2, m3_by_region_k20 = m3)
aic_tab <- data.frame(
  model = names(models),
  edf   = sapply(models, function(m) sum(m$edf)),
  AIC   = sapply(models, AIC),
  dev_explained = sapply(models, function(m) summary(m)$dev.expl),
  row.names = NULL
)
aic_tab$delta_AIC <- aic_tab$AIC - min(aic_tab$AIC)
cat("\n== AIC comparison (FEMALE) ==\n"); print(aic_tab)

best_name <- aic_tab$model[which.min(aic_tab$AIC)]
best <- models[[best_name]]
cat("\nBest by AIC:", best_name, "\n")

edf_best <- as.data.frame(summary(best)$s.table)
edf_best$smooth <- rownames(edf_best)
cat("\n== Effective df of each smooth (best model) ==\n"); print(edf_best)

saveRDS(best, file.path(model_dir, "gam_female.rds"))
write.csv(aic_tab,  file.path(tab_dir, "aic_female.csv"), row.names = FALSE)
write.csv(edf_best, file.path(tab_dir, "edf_female.csv"), row.names = FALSE)

cat("\n== gam.check (best FEMALE model) ==\n")
gam.check(best)

# ---- Secondary analysis: MALE and TOTAL ------------------------------------
min_obs  <- 15
included <- data.frame(); excluded <- data.frame(); edf_sec <- data.frame()

for (sx in c("MALE", "TOTAL")) {
  s <- d[d$sex == sx, ]
  n_by <- table(droplevels(s$region))
  keep <- names(n_by)[n_by >= min_obs]
  drop <- names(n_by)[n_by <  min_obs]
  included <- rbind(included, data.frame(sex = sx, region = keep,
                                         n_obs = as.integer(n_by[keep])))
  if (length(drop) > 0)
    excluded <- rbind(excluded, data.frame(sex = sx, region = drop,
                                           n_obs = as.integer(n_by[drop])))
  if (length(keep) < 2) { cat("Skipping", sx, "\n"); next }

  s <- s[s$region %in% keep, ]
  s$region <- droplevels(s$region)
  k_use <- min(10, min(n_by[keep]) - 1)      # k must stay below the shortest series
  m <- gam(rate ~ region + s(year, by = region, k = k_use),
           family = fam, data = s, method = "REML")
  saveRDS(m, file.path(model_dir, paste0("gam_", tolower(sx), ".rds")))

  e <- as.data.frame(summary(m)$s.table)
  e$smooth <- rownames(e); e$sex <- sx; e$k_used <- k_use
  edf_sec <- rbind(edf_sec, e)
  cat("\n== Effective df,", sx, "(k =", k_use, ") ==\n"); print(e)
}

cat("\n== Series excluded (fewer than", min_obs, "observations) ==\n")
print(excluded)

write.csv(included, file.path(tab_dir, "included_series.csv"), row.names = FALSE)
write.csv(excluded, file.path(tab_dir, "excluded_series.csv"), row.names = FALSE)
write.csv(edf_sec,  file.path(tab_dir, "edf_secondary.csv"),   row.names = FALSE)

# ---- Log to README ----------------------------------------------------------
excl_txt <- if (nrow(excluded) == 0) "none" else
  paste0(excluded$region, " (", excluded$sex, ", n = ", excluded$n_obs, ")", collapse = "; ")
cat(paste0(
  "\n- **", Sys.Date(), " - `R/03_models.R`:** Gamma(log) GAMs, mgcv REML. ",
  "Primary analysis: FEMALE series (complete for all 6 regions). ",
  "Compared 3 models by AIC; best = ", best_name, " (models/gam_female.rds). ",
  "Secondary MALE/TOTAL fitted only for series with >= ", min_obs,
  " obs. Excluded: ", excl_txt, ".\n"),
  file = "README.md", append = TRUE)

cat("\nDone.\n")
