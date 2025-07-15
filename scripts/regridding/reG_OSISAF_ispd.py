import sys, os, numpy as np, pandas as pd, xarray as xr
from datetime import datetime
from pathlib import Path
from tqdm import tqdm
import xesmf as xe

def main(year=None):
    mod_path = '/home/581/da1339/AFIM/src/AFIM/src'
    sys.path.insert(0, mod_path)
    from sea_ice_toolbox import SeaIceToolbox
    # --- Setup model and grid info ---
    afim_name        = "elps-min"
    dt0_str, dtN_str = f"{year}-01-01", f"{year}-12-31"
    SI_tools         = SeaIceToolbox(sim_name=afim_name, dt0_str=dt0_str, dtN_str=dtN_str)
    _, CICE          = SI_tools.load_processed_cice(zarr_CICE=True)
    CICE_ispd        = SI_tools.compute_ice_magnitude_from_ice_components_on_Bgrid(CICE, ivec_type="BT")
    CICE_SO          = CICE_ispd.isel(nj=SI_tools.hemisphere_dict['nj_slice'])
    G_CICE_SO        = {"lon": CICE_SO["TLON"].values, "lat": CICE_SO["TLAT"].values}
    aom2_name        = "AOM2-ERA5"
    SI_tools         = SeaIceToolbox(sim_name=aom2_name, dt0_str=dt0_str, dtN_str=dtN_str)
    _, AOM2          = SI_tools.load_processed_cice(zarr_CICE=True)
    AOM2_ispd        = SI_tools.compute_ice_magnitude_from_ice_components_on_Bgrid(AOM2, ivec_type="BT")
    AOM2_SO          = AOM2_ispd.isel(nj=SI_tools.hemisphere_dict['nj_slice'])
    # --- Locate relevant OSI-SAF files ---
    D_search = Path(Path.home(), "seaice", "OSI_SAF", "ice_drift_455m", "sh")
    F_search = "ice_drift_sh_ease2-750_cdr-v1p0_24h-*.nc"
    P_found  = sorted(D_search.rglob(F_search))
    P_year   = []
    for f in P_found:
        with xr.open_dataset(f) as ds:
            t = pd.to_datetime(ds["time"].values[0])
            if t.year == year:
                P_year.append(f)
    if not P_year:
        print(f"No OSI-SAF files found for year {year}.")
        return
    print(f"Found {len(P_year)} files for year {year}. Beginning regridding...")
    # --- Regridding Loop ---
    DS = []
    reG = None
    for i, F_ in tqdm(enumerate(P_year), total=len(P_year), desc=f"Regridding OSI-SAF for {year}"):
        osisaf = xr.open_dataset(F_)
        dt_obs = pd.to_datetime(osisaf["time"].values[0])
        last_cice_time = pd.to_datetime(CICE_SO["time"].values[-1])
        if dt_obs > last_cice_time:
            print(f"Stopping: reached {dt_obs}, beyond model end {last_cice_time}")
            break
        dt            = (osisaf["t1"] - osisaf["t0"]).astype("timedelta64[s]").astype(float)
        dt            = xr.where(dt == 0, np.nan, dt)
        disp_m        = np.sqrt(osisaf["dX"]**2 + osisaf["dY"]**2) * 1000
        ispd_obs      = (disp_m / dt).squeeze()
        ispd_obs.name = "ispd"
        ispd_obs.attrs["units"] = "m s-1"
        dt_str = dt_obs.strftime("%Y-%m-%d")
        if reG is None:
            reG = xe.Regridder(ispd_obs, G_CICE_SO, method="bilinear", periodic=True, reuse_weights=False)
        ispd_obs_reG = reG(ispd_obs).compute().rename({"y": "nj", "x": "ni"})
        ispd_cice = CICE_SO["ispd_BT"].sel(time=dt_str).compute()
        ispd_aom2 = AOM2_SO["ispd_BT"].sel(time=dt_str).compute()
        d_cice = ispd_obs_reG - ispd_cice
        d_aom2 = ispd_obs_reG - ispd_aom2
        d_sims = ispd_aom2 - ispd_cice
        rmse_cice = float(np.sqrt((d_cice**2).mean()))
        rmse_aom2 = float(np.sqrt((d_aom2**2).mean()))
        rmse_sims = float(np.sqrt((d_sims**2).mean()))
        ds_day = xr.Dataset(data_vars = dict(ispd_obs    = (("time", "nj", "ni"), ispd_obs_reG.expand_dims("time").data, {"units": "m/s", "long_name": "OSI-SAF regridded"}),
                                             ispd_CICE   = (("time", "nj", "ni"), ispd_cice.expand_dims("time").data, {"units": "m/s", "long_name": f"CICE6-SA {afim_name}"}),
                                             ispd_AOM2   = (("time", "nj", "ni"), ispd_aom2.expand_dims("time").data, {"units": "m/s", "long_name": "ACCESS-OM2-025 ERA5 forced"}),
                                             d_ispd_CICE = (("time", "nj", "ni"), d_cice.expand_dims("time").data, {"units": "m/s", "long_name": "obs - CICE"}),
                                             d_ispd_AOM2 = (("time", "nj", "ni"), d_aom2.expand_dims("time").data, {"units": "m/s", "long_name": "obs - AOM2"}),
                                             d_ispd_sims = (("time", "nj", "ni"), d_sims.expand_dims("time").data, {"units": "m/s", "long_name": "AOM2 - CICE"}),
                                             RMSE_CICE   = (("time",), [rmse_cice], {"units": "m/s", "long_name": "RMSE obs - CICE"}),
                                             RMSE_AOM2   = (("time",), [rmse_aom2], {"units": "m/s", "long_name": "RMSE obs - AOM2"}),
                                             RMSE_sims   = (("time",), [rmse_sims], {"units": "m/s", "long_name": "RMSE AOM2 - CICE"})),
                            coords = dict(time = ("time", [ispd_obs.time.values]),
                                          TLON = (("nj", "ni"), CICE_SO["TLON"].values),
                                          TLAT = (("nj", "ni"), CICE_SO["TLAT"].values)))
        DS.append(ds_day)
    # --- Save netCDF ---
    DS_ispds = xr.concat(DS, dim="time")
    out_path = Path.home() / "seaice" / "OSI_SAF" / "ice_drift_455m" / f"ice_speeds_and_differences_{afim_name}_{aom2_name}_{year}.nc"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    DS_ispds.to_netcdf(path=out_path, mode="w", format="NETCDF4", encoding={v: {"zlib": True, "complevel": 4} for v in DS_ispds.data_vars})
    print(f"Finished saving: {out_path}")

if __name__ == "__main__":
    import argparse
    from multiprocessing import freeze_support
    freeze_support()
    parser = argparse.ArgumentParser()
    parser.add_argument("--year", type=int, required=True, help="Year to process (e.g. 1994)")
    args = parser.parse_args()
    main(year=args.year)

