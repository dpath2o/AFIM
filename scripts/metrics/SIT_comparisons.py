#!/usr/bin/env python3
import sys, os
from pathlib import Path
import xarray as xr
import numpy as np
import pandas as pd
from scipy.stats import binned_statistic_2d
from tqdm import tqdm

def compute_sea_ice_thickness(hi, aice, tarea):
    mask = aice > 0.15
    hi = hi.where(mask)
    aice = aice.where(mask)
    tarea = tarea.where(mask)
    sia = (aice * tarea).sum(dim=("nj", "ni"))
    siv = (hi * tarea).sum(dim=("nj", "ni"))
    return siv / sia

def main(year=None, sim_name=None):
    mod_path = '/home/581/da1339/AFIM/src/AFIM/src'
    sys.path.insert(0, mod_path)
    from sea_ice_toolbox import SeaIceToolbox

    # === Setup ===
    dt0_str = f"{year}-01-01"
    dtN_str = f"{year}-12-31"
    save_dir = Path.home() / "AFIM_archive" / sim_name
    save_dir.mkdir(parents=True, exist_ok=True)
    save_file = save_dir / f"SIT_{sim_name}_{year}.nc"

    # === AOM2-ERA5 ===
    SI_tools = SeaIceToolbox(sim_name="AOM2-ERA5", dt0_str=dt0_str, dtN_str=dtN_str)
    _, AOM2 = SI_tools.load_processed_cice(zarr_CICE=True)
    AOM2_SO = AOM2.isel(nj=SI_tools.hemisphere_dict['nj_slice'])
    sit_list_model = []
    for i in tqdm(range(len(AOM2_SO.time)), desc=f"Computing SIT for AOM2-ERA5 ({year})"):
        time_val = pd.to_datetime(AOM2_SO.time.values[i]).floor('D')
        hi_t = AOM2_SO['hi'].isel(time=i).load()
        aice_t = AOM2_SO['aice'].isel(time=i).load()
        tarea_t = AOM2_SO['tarea'].isel(time=i).load()
        sit = compute_sea_ice_thickness(hi_t, aice_t, tarea_t)
        sit_da = xr.DataArray(sit.expand_dims("time"), coords={"time": [time_val]}, name="SIT_AOM2")
        sit_list_model.append(sit_da)
    SIT_AOM2 = xr.concat(sit_list_model, dim="time")

    # === AFIM simulation ===
    SI_tools = SeaIceToolbox(sim_name=sim_name, dt0_str=dt0_str, dtN_str=dtN_str)
    _, CICE = SI_tools.load_processed_cice(zarr_CICE=True)
    CICE_SO = CICE.isel(nj=SI_tools.hemisphere_dict['nj_slice'])
    sit_list_model = []
    for i in tqdm(range(len(CICE_SO.time)), desc=f"Computing SIT for AOM2-ERA5 ({year})"):
        time_val = pd.to_datetime(CICE_SO.time.values[i]).floor('D')
        hi_t = CICE_SO['hi'].isel(time=i).load()
        aice_t = CICE_SO['aice'].isel(time=i).load()
        tarea_t = CICE_SO['tarea'].isel(time=i).load()
        sit = compute_sea_ice_thickness(hi_t, aice_t, tarea_t)
        sit_da = xr.DataArray(sit.expand_dims("time"), coords={"time": [time_val]}, name="SIT_CICE")
        sit_list_model.append(sit_da)
    SIT_CICE = xr.concat(sit_list_model, dim="time")   

    # === ESA-CCI ===
    lon_bins = np.arange(-180, 181, 0.25)
    lat_bins = np.arange(-90, -49.75, 0.25)
    lon_centers = 0.5 * (lon_bins[:-1] + lon_bins[1:])
    lat_centers = 0.5 * (lat_bins[:-1] + lat_bins[1:])
    D_search = Path.home() / "seaice" / "ESA_CCI" / "L2P" / "envisat" / "sh"
    F_search = f"ESACCI-SEAICE-L2P-SITHICK-RA2_ENVISAT-SH-{year}*.nc"
    P_found = sorted(D_search.rglob(F_search))
    sit_list_obs = []
    for f in tqdm(P_found, desc=f"Gridding ESA-CCI SIT ({year})"):
        ds = xr.open_dataset(f)
        lon_raw = ds["lon"].values
        lat_raw = ds["lat"].values
        sit_raw = ds["sea_ice_thickness"].values
        valid = np.isfinite(lon_raw) & np.isfinite(lat_raw) & np.isfinite(sit_raw)
        lon = lon_raw[valid]
        lat = lat_raw[valid]
        sit = sit_raw[valid]
        sit_binned, _, _, _ = binned_statistic_2d(lon, lat, sit, statistic="mean", bins=[lon_bins, lat_bins])
        sit_grid = np.transpose(sit_binned)
        try:
            time_val = pd.to_datetime(ds.time.values[0]).floor("D")
        except Exception:
            time_val = pd.to_datetime(f.name.split("-")[-2], format="%Y%m%d")
        sit_da = xr.DataArray(
            sit_grid[np.newaxis, :, :],
            dims=("time", "y", "x"),
            coords={"time": [time_val], "y": lat_centers, "x": lon_centers},
            name="SIT_obs_grid"
        )
        sit_list_obs.append(sit_da)
    ESA_SIT_reG = xr.concat(sit_list_obs, dim="time")
    ESA_SIT_reG = ESA_SIT_reG.sel(time=~ESA_SIT_reG.get_index("time").duplicated())
    SIT_obs = ESA_SIT_reG.mean(dim=("x", "y"), skipna=True).rename("SIT_obs")

    # === Save ===
    xr.Dataset({"SIT_obs": SIT_obs, "SIT_AOM2": SIT_AOM2, "SIT_CICE": SIT_CICE}).to_netcdf(save_file)
    print(f"Saved daily SIT data to {save_file}")

if __name__ == "__main__":
    import argparse
    from multiprocessing import freeze_support
    freeze_support()
    parser = argparse.ArgumentParser()
    parser.add_argument("--year", type=int, required=True, help="Year to process (e.g. 2002)")
    parser.add_argument("--sim_name", type=str, required=True, help="Simulation name (e.g. elps-min)")
    args = parser.parse_args()
    main(year=args.year, sim_name=args.sim_name)
