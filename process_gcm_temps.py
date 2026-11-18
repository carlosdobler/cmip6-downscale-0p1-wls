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
#     "cftime",
# ]
# ///

"""Compute annual land-only and global mean near-surface air temperature series.

For one CMIP6 model, this script produces three annual mean series of daily
`tas` (near-surface air temperature):

1. DS_Land: area-weighted land-only mean of the downscaled (DS) product,
   1971 onward. The DS product has no ocean cells, so its valid-data mask
   is used as the land mask.
2. GCM_Global: area-weighted global mean of the driving GCM, historical
   (from 1850) followed by ssp585 through 2100.
3. GCM_Land: the same GCM series restricted to the land mask. The DS mask
   is taken from the first time step, regridded to the GCM grid with
   nearest-neighbour interpolation, and thresholded at 0.5.

Method notes:
- Spatial means weight each grid cell by cos(latitude).
- DS longitudes are shifted from [-180, 180] to [0, 360] to match the GCM,
  after dropping the duplicate -180 meridian.
- Data are read lazily from Google Cloud Storage with xarray and dask on a
  local Dask cluster. The DS series is computed in batches of 10 years to 
  avoid out-of-memory issues.
- The latest date-stamped version (vYYYYMMDD) of each GCM store is selected
  automatically.

Usage:
    uv run process_gcm_temps.py --model-name MPI-ESM1-2-HR \
        --institution-hist MPI-M --institution-ssp DKRZ \
        --ensemble-member r1i1p1f1 --grid-label gn [--output PATH]

Output:
    A CSV with one row per year and an index column named "Year". Columns:

    | Column     | Years     | Content                                      |
    |------------|-----------|----------------------------------------------|
    | DS_Land    | 1971-end  | Land-only mean of the downscaled product     |
    | GCM_Global | 1850-2100 | Global mean of the driving GCM               |
    | GCM_Land   | 1850-2100 | Land-only mean of the GCM, using the DS mask |

    Values are annual means of daily `tas` in the source units (Kelvin for
    CMIP6 `tas`); no unit conversion is applied. The default output path is
    intermediate_tables/global_mean_temperatures_<model-name>.csv.
"""

import argparse
import os
import re
import time
import warnings
import dask
from dask.distributed import Client, LocalCluster
import fsspec
import numpy as np
import pandas as pd
import xarray as xr

# Suppress xarray/dask warnings for clean output
warnings.filterwarnings("ignore")

VERSION_PATTERN = re.compile(r"^v(\d{8})$")
_LATEST_VERSION_CACHE = {}


def get_latest_version(
    activity: str,
    institution: str,
    model_name: str,
    experiment: str,
    ensemble: str,
    mip_table: str,
    grid: str,
    var_name: str,
) -> str:
    """Scans the CMIP6 bucket for the available data versions of a given model
    combination and returns the latest one."""
    cache_key = (
        activity,
        institution,
        model_name,
        experiment,
        ensemble,
        mip_table,
        grid,
        var_name,
    )
    if cache_key in _LATEST_VERSION_CACHE:
        return _LATEST_VERSION_CACHE[cache_key]

    path = (
        f"gs://cmip6/CMIP6/{activity}/{institution}/{model_name}/{experiment}/"
        f"{ensemble}/{mip_table}/{var_name}/{grid}"
    )
    fs = fsspec.filesystem("gs")
    versions = [
        os.path.basename(entry.rstrip("/")) for entry in fs.ls(path, detail=False)
    ]
    versions = [v for v in versions if VERSION_PATTERN.match(v)]

    if not versions:
        raise RuntimeError(f"No date-stamped version directories found under '{path}'.")

    latest = sorted(versions)[-1]
    _LATEST_VERSION_CACHE[cache_key] = latest
    return latest


def compute_area_weighted_mean(da, lat_name="lat", lon_name="lon"):
    """Computes area-weighted spatial mean."""
    weights = np.cos(np.deg2rad(da[lat_name]))
    weights.name = "weights"
    return da.weighted(weights).mean(dim=[lat_name, lon_name])


def main():
    parser = argparse.ArgumentParser(
        description="Compute annual mean temperatures from GCM and Downscaled datasets."
    )
    parser.add_argument(
        "--model-name", required=True, help="Full CMIP6 model name (e.g. MPI-ESM1-2-HR)"
    )
    parser.add_argument(
        "--institution-hist",
        required=True,
        help="Institution for historical experiment",
    )
    parser.add_argument(
        "--institution-ssp", required=True, help="Institution for ssp585 experiment"
    )
    parser.add_argument(
        "--ensemble-member", required=True, help="Ensemble member (e.g. r1i1p1f1)"
    )
    parser.add_argument(
        "--grid-label", required=True, help="Grid label (e.g. gn, gr, gr1)"
    )
    parser.add_argument("--output", default=None, help="Output CSV path")

    args = parser.parse_args()

    model_name = args.model_name
    output_file = (
        args.output
        if args.output
        else f"intermediate_tables/global_mean_temperatures_{model_name}.csv"
    )

    output_dir = os.path.dirname(output_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    print(f"Resolving versions for {model_name}...")
    hist_version = get_latest_version(
        activity="CMIP",
        institution=args.institution_hist,
        model_name=model_name,
        experiment="historical",
        ensemble=args.ensemble_member,
        mip_table="day",
        grid=args.grid_label,
        var_name="tas",
    )
    ssp_version = get_latest_version(
        activity="ScenarioMIP",
        institution=args.institution_ssp,
        model_name=model_name,
        experiment="ssp585",
        ensemble=args.ensemble_member,
        mip_table="day",
        grid=args.grid_label,
        var_name="tas",
    )

    gcm_hist_path = (
        f"gs://cmip6/CMIP6/CMIP/{args.institution_hist}/{model_name}/historical/"
        f"{args.ensemble_member}/day/tas/{args.grid_label}/{hist_version}/"
    )
    gcm_ssp_path = (
        f"gs://cmip6/CMIP6/ScenarioMIP/{args.institution_ssp}/{model_name}/ssp585/"
        f"{args.ensemble_member}/day/tas/{args.grid_label}/{ssp_version}/"
    )
    ds_path = (
        f"gs://clim_data_reg_useast1/cmip6_downscaled_woodwell/daily/tas/"
        f"tas_{model_name}_ww-isimip_ssp585_day.zarr"
    )

    print(f"Historical path: {gcm_hist_path}")
    print(f"SSP585 path:     {gcm_ssp_path}")
    print(f"Downscaled path: {ds_path}")

    print("\nStarting Dask cluster for parallel computation...")
    dask.config.set(
        {
            "distributed.comm.timeouts.connect": "30s",
            "distributed.comm.retry.count": 3,
        }
    )
    cluster = LocalCluster(n_workers=4, threads_per_worker=2, memory_limit="8GB")
    client = Client(cluster)

    try:
        t_start = time.time()

        print("\nOpening datasets lazily...")
        open_kwargs = {"engine": "zarr", "chunks": {}}
        if hasattr(xr, "coders") and hasattr(xr.coders, "CFDatetimeCoder"):
            open_kwargs["decode_times"] = xr.coders.CFDatetimeCoder(use_cftime=True)
        else:
            open_kwargs["use_cftime"] = True

        ds = xr.open_dataset(ds_path, **open_kwargs)
        ds_data = ds["tas"]

        gcm_hist = xr.open_dataset(gcm_hist_path, **open_kwargs)["tas"]
        gcm_ssp = xr.open_dataset(gcm_ssp_path, **open_kwargs)["tas"]
        # Limit SSP to year 2100 to ensure consistency with other CMIP6 models and DS series
        gcm_ssp = gcm_ssp.sel(time=gcm_ssp["time.year"] <= 2100)
        gcm_data = xr.concat([gcm_hist, gcm_ssp], dim="time")

        print("\nGenerating and regridding land mask...")
        ds_mask = ds["tas"].isel(time=0).compute().notnull()

        ds_lat = "latitude" if "latitude" in ds_mask.coords else "lat"
        ds_lon = "longitude" if "longitude" in ds_mask.coords else "lon"

        # Resolve longitude overlap (-180 and 180 are both present)
        if -180.0 in ds_mask[ds_lon].values:
            ds_mask = ds_mask.drop_sel({ds_lon: -180.0})

        # Shift longitudes from [-180, 180] to [0, 360] to match GCM
        lon = ds_mask[ds_lon].values
        lon = np.where(lon < 0, lon + 360, lon)
        ds_mask = ds_mask.assign_coords({ds_lon: lon}).sortby(ds_lon)

        # Rename coordinates to match GCM before interpolation
        ds_mask = ds_mask.rename({ds_lat: "lat", ds_lon: "lon"})

        gcm_mask = (
            ds_mask.astype(float).interp(
                lat=gcm_data.lat, lon=gcm_data.lon, method="nearest"
            )
            > 0.5
        )

        print("\nSetting up computational graph for GCM...")
        # 1. GCM Global
        gcm_global_ts = compute_area_weighted_mean(gcm_data, "lat", "lon")
        gcm_global_ann = gcm_global_ts.groupby("time.year").mean()

        # 2. GCM Land-only
        gcm_land_data = gcm_data.where(gcm_mask)
        gcm_land_ts = compute_area_weighted_mean(gcm_land_data, "lat", "lon")
        gcm_land_ann = gcm_land_ts.groupby("time.year").mean()

        print("\nComputing GCM Time Series...")
        res_gcm_global, res_gcm_land = dask.compute(gcm_global_ann, gcm_land_ann)

        s_gcm_global = pd.Series(
            res_gcm_global.values, index=res_gcm_global.year.values, name="GCM_Global"
        )
        s_gcm_land = pd.Series(
            res_gcm_land.values, index=res_gcm_land.year.values, name="GCM_Land"
        )

        print("\nComputing DS Land Time Series incrementally (by decade)...")
        ds_years = np.unique(ds_data["time.year"].values)
        ds_results = []
        batch_size = 10
        for i in range(0, len(ds_years), batch_size):
            batch_years = ds_years[i : i + batch_size]
            print(f"Processing DS years {batch_years[0]} to {batch_years[-1]}...")
            ds_batch = ds_data.sel(time=ds_data["time.year"].isin(batch_years))
            ds_ts = compute_area_weighted_mean(ds_batch, ds_lat, ds_lon)
            ds_ann = ds_ts.groupby("time.year").mean().compute()
            ds_results.append(
                pd.Series(ds_ann.values, index=ds_ann.year.values, name="DS_Land")
            )

        print("\nAssembling final CSV...")
        s_ds_land = pd.concat(ds_results)

        # Combine into a single DataFrame
        df_final = pd.concat([s_ds_land, s_gcm_global, s_gcm_land], axis=1)
        df_final.index.name = "Year"

        df_final.to_csv(output_file)
        print(f"\nSuccess! Data exported to {output_file}")
        print(f"Total processing time: {time.time() - t_start:.1f} seconds")

    finally:
        print("Shutting down Dask cluster...")
        client.shutdown()
        cluster.close()


if __name__ == "__main__":
    main()
