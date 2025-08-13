#!/usr/bin/env python3

import xarray as xr
from pathlib import Path
import sys

# Get sim_name from command line
sim_name = sys.argv[1]
base_path = Path.home() / "AFIM_archive" / sim_name / "zarr"
monthly_dir = base_path
daily_store = base_path / "iceh_daily.zarr"

daily_store.mkdir(parents=True, exist_ok=True)

monthly_zarrs = sorted(monthly_dir.glob("iceh_????-??.zarr"))
print(f"Found {len(monthly_zarrs)} monthly Zarrs for {sim_name}...")

for src in monthly_zarrs:
    group_name = src.stem.replace("iceh_", "")
    target_group_path = daily_store / group_name
    print(f"\n→ Processing {src.name} → writing group {group_name}")
    ds = xr.open_zarr(src, consolidated=False)
    ds.to_zarr(store=daily_store, group=group_name, consolidated=False, mode="w")
    print(f"Wrote group {group_name} to {daily_store}")

print(f"\nDone with {sim_name}.")
