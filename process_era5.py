#!/usr/bin/env -S uv run --script
#
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "xarray",
#     "numpy<2",
#     "pandas",
#     "dask[complete]",
#     "zarr",
#     "fsspec",
#     "gcsfs",
#     "scipy",
# ]
# ///

import xarray as xr
import numpy as np
import pandas as pd
from dask.distributed import Client
import time
import warnings

# Suppress warnings
warnings.filterwarnings("ignore")


def compute_area_weighted_mean(da, lat_name="y", lon_name="x"):
    """Computes area-weighted spatial mean, ignoring 0 values (ocean)."""
    # Create mask for non-zero values
    mask = da != 0
    da_masked = da.where(mask)

    weights = np.cos(np.deg2rad(da[lat_name]))
    weights.name = "weights"

    # Apply weights and compute mean only over masked points
    return da_masked.weighted(weights).mean(dim=[lat_name, lon_name])


def main():
    print("Starting Dask client...")
    client = Client()

    t_start = time.time()
    era5_path = (
        "gs://clim_data_reg_useast1/era5_land/daily_aggregates/temperature_2m.zarr"
    )

    print(f"\nOpening ERA5 dataset: {era5_path}")
    ds = xr.open_dataset(era5_path, engine="zarr", chunks={})
    da = ds["temperature_2m"]

    print("Computing area-weighted annual mean (masking 0 values)...")
    # Spatial mean (daily)
    ts_daily = compute_area_weighted_mean(da, "y", "x")

    # Annual mean
    ts_annual = ts_daily.groupby("time.year").mean().compute()

    print("\nSaving results...")
    s_era5 = pd.Series(ts_annual.values, index=ts_annual.year.values, name="ERA5_Land")
    s_era5.index.name = "Year"

    output_file = "intermediate_tables/era5_land_annual_mean.csv"
    s_era5.to_csv(output_file)

    print(f"\nDone! Results saved to {output_file}")
    print(f"Total time: {time.time() - t_start:.1f} seconds")

    client.close()


if __name__ == "__main__":
    main()
