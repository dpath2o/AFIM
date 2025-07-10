from scipy.stats import binned_statistic_2d
from pathlib import Path
import numpy as np
import xarray as xr
import pandas as pd
import os

def compute_gridcell_area(lat_bins, lon_bins, radius=6371000):
    dlat = np.diff(np.radians(lat_bins))
    dlon = np.diff(np.radians(lon_bins))
    lat_rad = np.radians(0.5 * (lat_bins[:-1] + lat_bins[1:]))
    area = (radius ** 2) * dlon[np.newaxis, :] * dlat[:, np.newaxis] * np.cos(lat_rad[:, np.newaxis])
    return area  # shape: (lat, lon)

def main():
    lon_bins    = np.arange(-180, 181, 0.25)
    lat_bins    = np.arange(-90, -49.75, 0.25)
    lon_centers = 0.5 * (lon_bins[:-1] + lon_bins[1:])
    lat_centers = 0.5 * (lat_bins[:-1] + lat_bins[1:])
    area        = compute_gridcell_area(lat_bins, lon_bins)
    area_da     = xr.DataArray(area, coords={"y": lat_centers, "x": lon_centers}, dims=("y", "x"))
    D_ESA_CCI   = Path.home() / "seaice" / "ESA_CCI" / "L2P" / "envisat" / "sh"
    file_list   = sorted(D_ESA_CCI.glob("ESACCI-SEAICE-L2P-SITHICK-RA2_ENVISAT-SH-*.nc"))
    sit_list    = []
    siv_list    = []
    for f in file_list:
        try:
            print(f"working on {f}")
            ds = xr.open_dataset(f)
            lon_raw = ds['lon'].values
            lat_raw = ds['lat'].values
            sit_raw = ds['sea_ice_thickness'].values
            valid = np.isfinite(lon_raw) & np.isfinite(lat_raw) & np.isfinite(sit_raw)
            lon = lon_raw[valid]
            lat = lat_raw[valid]
            sit = sit_raw[valid]
            sit_binned, _, _, _ = binned_statistic_2d(lon, lat, sit, statistic='mean', bins=[lon_bins, lat_bins])
            sit_grid = np.transpose(sit_binned)  # shape: (lat, lon)
            time_val = pd.to_datetime(ds['time'].dt.floor('D').values[0])
            sit_da   = xr.DataArray(sit_grid[np.newaxis, :, :],
                                    dims=("time", "y", "x"),
                                    coords={"time": [time_val], "y": lat_centers, "x": lon_centers},
                                    name="SIT")
            sit_valid = sit_da.where(sit_da.notnull(), 0.0)
            siv       = (sit_valid * area_da).sum(dim=("y", "x"))
            sit_list.append(sit_da)
            siv_list.append(xr.DataArray(np.atleast_1d(siv_km3.values),  # ensures 1D array
                                        coords={"time": [time_val]},
                                        dims=["time"],
                                        name="SIV"))
        except Exception as e:
            print(f"Failed on {f.name}: {e}")
    sit_all = xr.concat(sit_list, dim="time")
    siv_all = xr.concat(siv_list, dim="time").rename("SIV")
    ESA_CCI_out = xr.merge([sit_all, siv_all])
    ESA_CCI_out['doy'] = ESA_CCI_out['time'].dt.dayofyear
    siv_clim = ESA_CCI_out['SIV'].groupby('doy').mean(dim='time').rename("SIV_clim")
    sit_clim = ESA_CCI_out['SIT'].groupby('doy').mean(dim='time').rename("SIT_clim")
    ESA_CCI_out_all = xr.Dataset({"SIT": ESA_CCI_out['SIT'],
                                "SIV": ESA_CCI_out['SIV'],
                                "SIT_clim": sit_clim,
                                "SIV_clim": siv_clim})
    ESA_CCI_out_all.to_netcdf( Path(Path.home(),"seaice","ESA_CCI","ESA_CCI_L2P_envisat_SH_SIT_SIV_daily.nc") )

if __name__ == '__main__':
    from multiprocessing import freeze_support
    freeze_support()
    main()
