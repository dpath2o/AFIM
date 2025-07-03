#!/usr/bin/env python3

def bin_esa_thickness_to_cice_grid(ESA_CCI, CICE_SO):
    import numpy as np
    import xarray as xr
    from scipy.stats import binned_statistic_2d
    from tqdm import tqdm
    x_edges = np.linspace(CICE_SO['TLON'].min().item(), CICE_SO['TLON'].max().item(), CICE_SO['TLON'].shape[1] + 1)
    y_edges = np.linspace(CICE_SO['TLAT'].min().item(), CICE_SO['TLAT'].max().item(), CICE_SO['TLAT'].shape[0] + 1)
    unique_days = np.unique(ESA_CCI['time'].dt.floor('D').values)
    sit_list = []
    sigma_list = []
    time_list = []
    for day in tqdm(unique_days, desc="Binning ESA CCI by day"):
        ESA_day = ESA_CCI.where(ESA_CCI['time'].dt.floor('D') == day, drop=True)
        lon = ESA_day.lon.values
        lat = ESA_day.lat.values
        sit = ESA_day.sea_ice_thickness.values
        sit_sigma = ESA_day.sea_ice_thickness_uncertainty.values
        mask = np.isfinite(lon) & np.isfinite(lat) & np.isfinite(sit) & np.isfinite(sit_sigma)
        lon = lon[mask]
        lat = lat[mask]
        sit = sit[mask]
        sit_sigma = sit_sigma[mask]
        binned_sit, _, _, _ = binned_statistic_2d(lon, lat, sit, statistic="mean", bins=[x_edges, y_edges])
        binned_sigma, _, _, _ = binned_statistic_2d(lon, lat, sit_sigma, statistic="mean", bins=[x_edges, y_edges])
        sit_list.append(binned_sit.T[None, ...])
        sigma_list.append(binned_sigma.T[None, ...])
        time_list.append(np.datetime64(day))
    ESA_sit_reG = xr.DataArray(
        data=np.concatenate(sit_list, axis=0),
        dims=("time", "nj", "ni"),
        coords={"time": time_list, "TLON": (("nj", "ni"), CICE_SO.TLON), "TLAT": (("nj", "ni"), CICE_SO.TLAT)},
        name="ESA_sit"
    )
    ESA_sit_sigma_reG = xr.DataArray(
        data=np.concatenate(sigma_list, axis=0),
        dims=("time", "nj", "ni"),
        coords={"time": time_list, "TLON": (("nj", "ni"), CICE_SO.TLON), "TLAT": (("nj", "ni"), CICE_SO.TLAT)},
        name="ESA_sit_sigma"
    )
    return xr.Dataset({"ESA_sit": ESA_sit_reG, "ESA_sit_sigma": ESA_sit_sigma_reG})

def sort_time(ds):
    return ds.sortby("time")

def main():
    import sys, glob
    import xarray as xr
    from pathlib import Path
    # === Setup AFIM environment ===
    mod_path = '/home/581/da1339/AFIM/src/AFIM/src'
    sys.path.insert(0, mod_path)
    from sea_ice_toolbox import SeaIceToolbox
    # === Load model grid from SeaIceToolbox ===
    D_ESACCI = "/g/data/gv90/da1339/SeaIce/ESA_CCI/L2P/envisat/sh"
    SI_tools = SeaIceToolbox(sim_name='elps-min', dt0_str="2002-01-01", dtN_str="2012-12-31", P_log="/g/data/gv90/da1339/logs/ESACCI_SIT_reG.log")
    CICE_all = SI_tools.load_iceh_zarr()
    CICE_SO = CICE_all.isel(nj=SI_tools.hemisphere_dict['nj_slice'])
    # === Load ESA CCI L2P ===
    files = sorted(glob.glob(f"{D_ESACCI}/ESACCI-SEAICE-L2P-SITHICK-RA2_ENVISAT-SH-*.nc"))
    ESA_CCI = xr.open_mfdataset(files, combine="nested", concat_dim="time", preprocess=sort_time, parallel=False)
    # === Regrid and Save ===
    ds_out = bin_esa_thickness_to_cice_grid(ESA_CCI, CICE_SO)
    output_path = Path(f"{D_ESACCI}/reG/ESA_CCI_SIT_regridded.zarr")
    ds_out.to_zarr(output_path, mode="w")
    print(f"✅ Regridded ESA CCI saved to: {output_path}")

if __name__ == '__main__':
    from multiprocessing import freeze_support
    freeze_support()
    main()
