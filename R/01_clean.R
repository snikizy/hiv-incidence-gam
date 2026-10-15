# 01_clean.R
# Split WHO HIV incidence rows into country, WHO-region, and global tables.

raw_path <- file.path("data", "raw", "RELAY_WHS.csv")
out_dir <- file.path("data", "clean")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

raw <- read.csv(raw_path, stringsAsFactors = FALSE, na.strings = c("", "NA"))

keep_and_flag <- function(df) {
  out <- data.frame(
    year = as.integer(df$DIM_TIME),
    region_name = df$GEO_NAME_SHORT,
    sex = df$DIM_SEX,
    rate = as.numeric(df$RATE_PER_1000_N),
    rate_lower = as.numeric(df$RATE_PER_1000_NL),
    rate_upper = as.numeric(df$RATE_PER_1000_NU),
    stringsAsFactors = FALSE
  )
  out$rate_na <- is.na(out$rate)
  out
}

country <- keep_and_flag(raw[raw$DIM_GEO_CODE_TYPE == "COUNTRY", , drop = FALSE])
who_region <- keep_and_flag(raw[raw$DIM_GEO_CODE_TYPE == "WHOREGION", , drop = FALSE])
global <- keep_and_flag(raw[raw$DIM_GEO_CODE_TYPE == "GLOBAL", , drop = FALSE])

write.csv(country, file.path(out_dir, "country.csv"), row.names = FALSE)
write.csv(who_region, file.path(out_dir, "who_region.csv"), row.names = FALSE)
write.csv(global, file.path(out_dir, "global.csv"), row.names = FALSE)

tables <- list(country = country, who_region = who_region, global = global)
for (nm in names(tables)) {
  df <- tables[[nm]]
  cat(
    sprintf(
      "%s: dim = %d x %d; NA rates = %d\n",
      nm, nrow(df), ncol(df), sum(df$rate_na)
    )
  )
}
