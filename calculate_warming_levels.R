#!/usr/bin/env Rscript
#
# Orchestrates the warming-level workflow and writes the final results table.
#
# For each model in the model registry (fetched from GitHub), the script:
#   1. Generates the annual land/global temperature table with
#      process_gcm_temps.py, unless it already exists in intermediate_tables/.
#   2. Extends DS_Land before 1971 with GCM_Land, offset by the mean
#      DS_Land-minus-GCM_Land difference over the overlap period (1971-2100),
#      giving DS_Land_ext.
#   3. Builds DS_Global (filling the ocean gap) by adding a smooth exponential
#      fit of the GCM_Global - GCM_Land difference to DS_Land_ext.
#   4. Computes DS_Global_anom, the anomaly from the 1850-1900 mean of DS_Global,
#      and DS_Global_anom_w, its 21-year centred running mean (+/-10 years,
#      complete windows only).
#   5. Records the first year in which DS_Global_anom_w reaches each warming
#      level from 0.5 to 3.0 C in 0.5 C steps.
#
# Models already present in warming_levels.csv are skipped, and existing rows
# for a model are replaced when it is reprocessed.
#
# Outputs:
#   warming_levels.csv
#     One row per model and warming level, with columns:
#       gcm   CMIP6 model name (character)
#       wl    warming level above pre-industrial, in C (0.5, 1.0, ..., 3.0)
#       year  first year DS_Global_anom_w >= wl; NA if never reached
#   intermediate_tables/annual_temperatures_smoothed_<model>.csv
#     Per-model annual series: the inputs DS_Land, GCM_Global, and GCM_Land,
#     plus DS_Land_ext, DS_Global, DS_Global_anom, and DS_Global_anom_w.
#
# Requires ./process_gcm_temps.py to be executable and on the working directory,
# plus network access to the model registry and Google Cloud Storage.

library(tidyverse)
library(yaml)
library(slider)

REGISTRY_URL <- "https://raw.githubusercontent.com/carlosdobler/cmip6-downscale-0p1/refs/heads/main/model_registry.yaml"
OUTPUT_TABLE <- "warming_levels.csv"
INTERMEDIATE_DIR <- "intermediate_tables"
WARMING_LEVELS <- seq(0.5, 3.0, by = 0.5)

dir.create(INTERMEDIATE_DIR, showWarnings = FALSE, recursive = TRUE)

first_year_exceeding <- function(df, threshold) {
  res <- df |>
    filter(DS_Global_anom_w >= threshold) |>
    slice_min(Year, n = 1) |>
    pull(Year)
  if (length(res) == 0) NA_real_ else as.numeric(res[1])
}

message("Fetching model registry from GitHub...")
registry <- yaml::read_yaml(REGISTRY_URL)

# Read or initialize master warming levels table
if (file.exists(OUTPUT_TABLE)) {
  wl_master <- read_csv(OUTPUT_TABLE, show_col_types = FALSE)
  # Validate columns
  if (!all(c("gcm", "wl", "year") %in% names(wl_master))) {
    warning(
      "Existing warming_levels.csv does not have expected columns. Initializing fresh table."
    )
    wl_master <- tibble(gcm = character(), wl = double(), year = double())
  }
} else {
  wl_master <- tibble(gcm = character(), wl = double(), year = double())
}

for (model_key in names(registry)) {
  model_entry <- registry[[model_key]]
  model_name <- model_entry$model_name
  inst_hist <- model_entry$institution_hist
  inst_ssp <- model_entry$institution_ssp
  ens <- model_entry$ensemble_member
  grid <- model_entry$grid_label

  message("\n--------------------------------------------------")
  message(str_glue("Checking model: {model_name} ({model_key})"))

  # Tier 1 check: Is the model already processed in the final table?
  already_in_master <- model_name %in% wl_master$gcm
  if (already_in_master) {
    message(str_glue(
      "Model {model_name} already present in {OUTPUT_TABLE}. Skipping."
    ))
    next
  }

  raw_intermediate_csv <- file.path(
    INTERMEDIATE_DIR,
    str_glue("global_mean_temperatures_{model_name}.csv")
  )

  # Tier 2 check: Does the intermediate table already exist?
  if (!file.exists(raw_intermediate_csv)) {
    message(str_glue(
      "Intermediate raw table not found. Calling process_gcm_temps.py for {model_name}..."
    ))
    cmd_args <- c(
      "--model-name",
      model_name,
      "--institution-hist",
      inst_hist,
      "--institution-ssp",
      inst_ssp,
      "--ensemble-member",
      ens,
      "--grid-label",
      grid,
      "--output",
      raw_intermediate_csv
    )
    status <- system2("./process_gcm_temps.py", args = cmd_args)
    if (status != 0) {
      stop(str_glue(
        "Error running process_gcm_temps.py for model {model_name} (exit status: {status})"
      ))
    }
  } else {
    message(str_glue(
      "Found existing intermediate table: {raw_intermediate_csv}. Skipping GCS processing."
    ))
  }

  # Load raw annual means
  message(str_glue("Processing temperatures for {model_name}..."))
  data <- read_csv(raw_intermediate_csv, show_col_types = FALSE) |>
    arrange(Year)

  # Calculate bias between DS and GCM land-only
  diff_ds_gcm <- mean(data$GCM_Land - data$DS_Land, na.rm = TRUE) * -1

  data <- data |>
    mutate(
      DS_Land_ext = if_else(is.na(DS_Land), GCM_Land + diff_ds_gcm, DS_Land)
    )

  # Fit exponential model for GCM Global - Land difference
  year_center <- mean(data$Year, na.rm = TRUE)
  years_c <- data$Year - year_center
  diff_land_global <- data$GCM_Global - data$GCM_Land

  fit_exp <- tryCatch(
    {
      nls(
        delta ~ a * exp(b * yr) + c,
        data = data.frame(delta = diff_land_global, yr = years_c),
        start = list(a = -0.5, b = 0.01, c = -1.5),
        control = nls.control(maxiter = 500, warnOnly = TRUE)
      )
    },
    error = function(e) {
      nls(
        delta ~ a * exp(b * yr) + c,
        data = data.frame(delta = diff_land_global, yr = years_c),
        start = list(a = -1.0, b = 0.005, c = -1.0),
        control = nls.control(maxiter = 500, warnOnly = TRUE)
      )
    }
  )

  delta_smooth_exp <- predict(fit_exp)

  # Reconstruct DS Global, calculate baseline and smoothed anomaly
  data <- data |>
    mutate(
      DS_Global = DS_Land_ext + delta_smooth_exp
    )

  pre_ind <- data |>
    filter(Year >= 1850, Year <= 1900) |>
    pull(DS_Global) |>
    mean(na.rm = TRUE)

  data <- data |>
    mutate(
      DS_Global_anom = DS_Global - pre_ind,
      DS_Global_anom_w = slider::slide_dbl(
        DS_Global_anom,
        .f = mean,
        .before = 10,
        .after = 10,
        .complete = TRUE
      )
    )

  # Save per-model smoothed table
  smoothed_intermediate_csv <- file.path(
    INTERMEDIATE_DIR,
    str_glue("annual_temperatures_smoothed_{model_name}.csv")
  )
  write_csv(data, smoothed_intermediate_csv)
  message(str_glue("Saved smoothed series to {smoothed_intermediate_csv}."))

  # Compute breaching years
  model_wl <- tibble(
    gcm = model_name,
    wl = WARMING_LEVELS
  ) |>
    mutate(year = map_dbl(wl, ~ first_year_exceeding(data, .x)))

  # Update master table: replace any existing records for this model, keep all others
  wl_master <- bind_rows(
    wl_master |> filter(gcm != model_name),
    model_wl
  ) |>
    arrange(gcm, wl)

  write_csv(wl_master, OUTPUT_TABLE)
  message(str_glue("Updated {OUTPUT_TABLE} with results for {model_name}:"))
  print(model_wl)
}

message("\n==================================================")
message("Pipeline run complete.")
