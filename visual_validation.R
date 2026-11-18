library(tidyverse)

WARMING_LEVELS_FILE <- "warming_levels.csv"
INTERMEDIATE_DIR <- "intermediate_tables"
OUTPUT_DIR <- "visual_validation"

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

plot_dataset_comparison <- function(data) {
  data |>
    select(Year, DS_Global, GCM_Global, DS_Land_ext, GCM_Land) |>
    pivot_longer(-Year, names_to = "dataset") |>
    ggplot(aes(x = Year, y = value, color = dataset, group = dataset)) +
    geom_line() +
    geom_point() +
    scale_colour_manual(
      name = NULL,
      values = c(
        "GCM_Global" = "dodgerblue3",
        "GCM_Land" = "tomato",
        "DS_Land_ext" = "orange",
        "DS_Global" = "deepskyblue1"
      )
    ) +
    labs(
      title = "DS_Global = DS_Land_ext with exponential correction"
    )
}

plot_breaching_years <- function(data, wl_model, baseline) {
  breaches <-
    wl_model |>
    filter(!is.na(year)) |>
    rename(Year = year) |>
    left_join(data |> select(Year, DS_Global_anom_w), by = "Year") |>
    mutate(label = sprintf("WL %.1f: %d", wl, Year))

  data |>
    select(Year, DS_Global_anom_w) |>
    ggplot(aes(x = Year, y = DS_Global_anom_w)) +
    geom_line() +
    geom_point(
      data = breaches,
      size = 3,
      color = "red"
    ) +
    geom_text(
      data = breaches,
      aes(label = label),
      hjust = 1.1,
      vjust = 0.5,
      size = 4,
      color = "red"
    ) +
    labs(y = sprintf("anomaly from baseline (%.2f K)", baseline))
}

wl_table <- read_csv(WARMING_LEVELS_FILE, show_col_types = FALSE)

for (model in unique(wl_table$gcm)) {
  png_paths <- file.path(
    OUTPUT_DIR,
    c(
      paste0(model, "_dataset_comparison.png"),
      paste0(model, "_breaching_years.png")
    )
  )

  if (all(file.exists(png_paths))) {
    message("Figures for ", model, " already exist. Skipping.")
    next
  }

  smoothed_csv <- file.path(
    INTERMEDIATE_DIR,
    paste0("annual_temperatures_smoothed_", model, ".csv")
  )
  if (!file.exists(smoothed_csv)) {
    warning("Missing ", smoothed_csv, ". Skipping ", model, ".")
    next
  }

  data <- read_csv(smoothed_csv, show_col_types = FALSE)

  baseline <-
    data |>
    filter(Year >= 1850, Year <= 1900) |>
    pull(DS_Global) |>
    mean(na.rm = TRUE)

  wl_model <-
    wl_table |>
    filter(gcm == model) |>
    select(wl, year)

  message("Saving figures for ", model, "...")

  ggsave(
    png_paths[1],
    plot = plot_dataset_comparison(data),
    width = 8,
    height = 5
  )

  ggsave(
    png_paths[2],
    plot = plot_breaching_years(data, wl_model, baseline),
    width = 7,
    height = 5
  )
}
