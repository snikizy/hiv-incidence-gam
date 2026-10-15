# 02_eda.R
# WHO-region HIV incidence: missing grid, trajectories with uncertainty, peak years.

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", repos = "https://cloud.r-project.org")
}
library(ggplot2)

df <- read.csv(file.path("data", "clean", "who_region.csv"), stringsAsFactors = FALSE)

region_levels <- c(
  "Africa",
  "Americas",
  "South-East Asia",
  "Europe",
  "Eastern Mediterranean",
  "Western Pacific"
)
sex_levels <- c("FEMALE", "MALE", "TOTAL")
sex_labels <- c(FEMALE = "Female", MALE = "Male", TOTAL = "Total")

df$region_name <- factor(df$region_name, levels = region_levels)
df$sex <- factor(df$sex, levels = sex_levels)

years <- seq(min(df$year), max(df$year))
grid <- expand.grid(
  region_name = factor(region_levels, levels = region_levels),
  sex = factor(sex_levels, levels = sex_levels),
  year = years,
  stringsAsFactors = FALSE
)
plot_df <- merge(grid, df, by = c("region_name", "sex", "year"), all.x = TRUE)
plot_df <- plot_df[order(plot_df$region_name, plot_df$sex, plot_df$year), ]

cat("Regions:\n")
print(region_levels)
cat("\nSex:", paste(sex_levels, collapse = ", "), "\n")
cat("Years:", min(years), "-", max(years), "(", length(years), "years )\n", sep = " ")

missing <- plot_df[is.na(plot_df$rate), c("region_name", "sex", "year")]
cat("\nMissing region x sex x year combinations:", nrow(missing), "of", nrow(grid), "\n")
print(missing, row.names = FALSE)

fig_dir <- file.path("outputs", "figures")
tab_dir <- file.path("outputs", "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

region_cols <- c(
  "Africa" = "#D55E00",
  "Americas" = "#0072B2",
  "South-East Asia" = "#009E73",
  "Europe" = "#CC79A7",
  "Eastern Mediterranean" = "#E69F00",
  "Western Pacific" = "#56B4E9"
)

base_plot <- ggplot(
  plot_df,
  aes(x = year, y = rate, colour = region_name, fill = region_name)
) +
  geom_ribbon(aes(ymin = rate_lower, ymax = rate_upper), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.7, na.rm = FALSE) +
  facet_wrap(~ sex, labeller = labeller(sex = sex_labels)) +
  scale_colour_manual(values = region_cols, name = "WHO region") +
  scale_fill_manual(values = region_cols, name = "WHO region") +
  scale_x_continuous(breaks = seq(1990, 2024, by = 5)) +
  labs(
    x = "Year",
    y = "HIV infections per 1,000 population",
    title = "HIV incidence rate by WHO region",
    subtitle = "Shaded band is the published lower–upper uncertainty interval; gaps are missing years"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey95", colour = "grey80"),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(fig_dir, "who_region_rate_by_year.png"),
  plot = base_plot,
  width = 11,
  height = 6.5,
  dpi = 300
)

log_plot <- base_plot +
  scale_y_log10() +
  labs(
    y = "HIV infections per 1,000 population (log10)",
    title = "HIV incidence rate by WHO region (log scale)"
  )

ggsave(
  file.path(fig_dir, "who_region_rate_by_year_log.png"),
  plot = log_plot,
  width = 11,
  height = 6.5,
  dpi = 300
)

obs <- df[!is.na(df$rate), ]
split_keys <- interaction(obs$region_name, obs$sex, drop = TRUE, lex.order = TRUE)
peak_list <- lapply(split(obs, split_keys), function(x) {
  max_rate <- max(x$rate)
  peak_years <- sort(x$year[x$rate == max_rate])
  data.frame(
    region_name = as.character(x$region_name[1]),
    sex = as.character(x$sex[1]),
    peak_year = peak_years[1],
    peak_rate = max_rate,
    n_years_at_peak = length(peak_years),
    stringsAsFactors = FALSE
  )
})
peaks <- do.call(rbind, peak_list)
rownames(peaks) <- NULL
peaks$region_name <- factor(peaks$region_name, levels = region_levels)
peaks$sex <- factor(peaks$sex, levels = sex_levels)
peaks <- peaks[order(peaks$region_name, peaks$sex), ]

write.csv(peaks, file.path(tab_dir, "who_region_peaks.csv"), row.names = FALSE)
cat("\nPeak year and peak rate per region and sex:\n")
print(peaks, row.names = FALSE)
cat("\nWrote figures to", fig_dir, "and peak table to", tab_dir, "\n")
