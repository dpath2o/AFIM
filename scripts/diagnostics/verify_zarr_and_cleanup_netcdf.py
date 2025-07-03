#!/usr/bin/env python3
"""
verify_zarr_and_cleanup_netcdf.py

Verifies that each monthly Zarr directory (iceh_YYYY-MM.zarr) under a given
AFIM simulation archive faithfully represents all the expected daily NetCDF
files (iceh.YYYY-MM-DD.nc) in terms of:
  - Date coverage
  - Variable presence

If verified, and --delete is used, it optionally deletes the matching NetCDF files
(after user confirmation).

This script runs in parallel across months. Recommended usage:
  --> Run this script inside an interactive PBS job to utilize multiple CPUs.

Produces a persistent cleanup.log file under the simulation directory,
with a "last updated" timestamp and one entry per verified or skipped month.
"""

import xarray as xr
import argparse
import concurrent.futures
from pathlib import Path
from datetime import datetime, timedelta
import pandas as pd
import sys, os
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

# ------------- Argument Parsing -------------
def parse_args():
    parser = argparse.ArgumentParser(description="Verify and optionally clean up AFIM daily NetCDF files after Zarr conversion.")
    parser.add_argument("sim_name", type=str, help="Simulation name (e.g. FI-heavy)")
    parser.add_argument("--base", type=str, default=str(Path.home() / "AFIM_archive"), help="Base archive path [default: ~/AFIM_archive]")
    parser.add_argument("--dry-run", action="store_true", help="Perform checks only, do not delete NetCDF files")
    parser.add_argument("--delete", action="store_true", help="After verification, prompt to delete matching NetCDF files")
    parser.add_argument("--max-workers", type=int, default=4, help="Number of processes for parallel execution")
    return parser.parse_args()

# ------------- Helpers -------------
def get_month_range(year_month):
    dt0 = datetime.strptime(year_month + "-01", "%Y-%m-%d")
    dtN = (dt0.replace(day=28) + timedelta(days=4)).replace(day=1) - timedelta(days=1)
    return [dt0 + timedelta(days=i) for i in range((dtN - dt0).days + 1)]

def verify_month(args_tuple):
    zarr_path, nc_dir, done_marker, dry_run = args_tuple
    year_month = zarr_path.stem.split("_")[1]
    if done_marker.exists():
        return f"[SKIP] {year_month}: already verified (.done exists)"
    dt_list = get_month_range(year_month)
    nc_files = [nc_dir / f"iceh.{dt.strftime('%Y-%m-%d')}.nc" for dt in dt_list]
    existing_nc_files = [f for f in nc_files if f.exists()]
    if not existing_nc_files:
        return f"[SKIP] {year_month}: no NetCDFs remain"
    try:
        zarr = xr.open_zarr(zarr_path)
    except Exception as e:
        return f"[FAIL] {year_month}: cannot open Zarr: {e}"
    zarr_dates = pd.to_datetime(zarr["time"].values).normalize()
    expected_dates = pd.to_datetime([dt.date() for dt in dt_list])
    all_dates_present = set(expected_dates).issubset(set(zarr_dates))
    variables_expected = set()
    for f in existing_nc_files:
        try:
            ds = xr.open_dataset(f)
            variables_expected.update(ds.data_vars)
        except Exception as e:
            return f"[FAIL] {year_month}: cannot open NetCDF {f.name}: {e}"
    variables_missing = [v for v in variables_expected if v not in zarr.data_vars]
    if not all_dates_present:
        return f"[FAIL] {year_month}: missing time steps in Zarr"
    if variables_missing:
        return f"[FAIL] {year_month}: Zarr missing variables: {variables_missing}"
    if not dry_run:
        done_marker.touch()
    return f"[OK]   {year_month}: verified ({len(existing_nc_files)} files)"

# ------------- Main Execution -------------
def main():
    args = parse_args()
    sim_path = Path(args.base) / args.sim_name
    zarr_dir = sim_path / "zarr"
    daily_dir = sim_path / "history" / "daily"
    log_path = sim_path / "cleanup.log"
    zarr_months = sorted([p for p in zarr_dir.glob("iceh_????-??.zarr") if p.is_dir()])
    tasks = []
    for zarr_path in zarr_months:
        ym = zarr_path.stem.split("_")[1]
        done_marker = zarr_dir / f".done_{ym}"
        tasks.append((zarr_path, daily_dir, done_marker, args.dry_run))
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.max_workers) as executor:
        results = list(executor.map(verify_month, tasks))
    deletions = []
    for res in results:
        print(res)
    with open(log_path, "a") as logf:
        logf.write("\n# Last updated: " + datetime.now().strftime("%Y-%m-%d %H:%M:%S") + "\n")
        for res in results:
            logf.write(res + "\n")
    if args.delete:
        print("\nDeletion mode active.")
        verified_months = []
        total_files = []
        for zarr_path in zarr_months:
            ym = zarr_path.stem.split("_")[1]
            done_marker = zarr_dir / f".done_{ym}"
            if not done_marker.exists():
                continue  # skip unverified
            nc_files = list(daily_dir.glob(f"iceh.{ym}-??.nc"))
            if nc_files:
                verified_months.append((ym, nc_files))
                total_files.extend(nc_files)
        if not total_files:
            print("No deletable NetCDF files found.")
        else:
            print(f"\n🔍 {len(total_files)} NetCDF files across {len(verified_months)} verified months are eligible for deletion.")
            confirm = input("Confirm delete all these files? [y/N] ").strip().lower()
            if confirm == "y":
                for ym, files in verified_months:
                    for f in files:
                        try:
                            f.unlink()
                            print(f"[DELETED] {f.name}")
                            log_entries.append(f"[DELETED] {f}")
                        except Exception as e:
                            print(f"[ERROR] Could not delete {f.name}: {e}")
                            log_entries.append(f"[ERROR] Failed to delete {f}: {e}")
                log_entries.append(f"# Deletion complete: {len(total_files)} files removed")
            else:
                print("Deletion cancelled.")
                log_entries.append("# Deletion prompt declined — no files deleted")

if __name__ == "__main__":
    main()

