# CMIP6 downscaled warming levels

This repo estimates the years in which global mean temperatures from downscaled climate models cross +0.5 to +3.0 °C above pre-industrial levels. The downscaled data comes from a 0.1° resolution product that covers land only and begins in 1971, so it cannot provide a pre-industrial baseline on its own.

The workflow fills the DS land series before 1971 from the driving GCM, reconstructs a DS global series by adding the GCM's land-to-global difference, and reads breaching years from the 21-year smoothed anomaly.

## Pipeline

1. `process_era5.py` (Python, xarray/dask) computes an ERA5-Land annual mean for checking the DS land series against observations.
2. `process_gcm_temps.py` (Python, xarray/dask): writes annual land-only DS, global GCM, and land-only GCM means to `intermediate_tables/global_mean_temperatures_<model>.csv`.
3. `calculate_warming_levels.R` (R): **main script**; it reads the model registry, runs step 2 for any model without an intermediate table, and writes `warming_levels.csv`. It also saves per-model smoothed series to `intermediate_tables/`.
4. `visual_validation.R` (R): saves comparison and breaching-year figures to `visual_validation/`.


## Outputs

- `warming_levels.csv`: columns `gcm`, `wl` (°C above pre-industrial), and `year` (first year the smoothed anomaly reaches `wl`; `NA` if never reached).
- `intermediate_tables/`: per-model annual series.
- `visual_validation/`: PNG figures.

## Requirements

- `uv` for the Python scripts (dependencies are declared inline).
- R with `tidyverse`, `yaml`, and `slider`, plus `ggplot2` for the figures.
- Network access to Google Cloud Storage and GitHub.

## Usage

```sh
Rscript calculate_warming_levels.R
Rscript visual_validation.R
```
