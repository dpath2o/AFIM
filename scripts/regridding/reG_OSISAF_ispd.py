import sys, os, argparse, pygmt
import numpy       as np
import pandas      as pd
import xarray      as xr
from datetime      import datetime
from pathlib       import Path
from glob          import glob
import warnings
warnings.filterwarnings("ignore",
                        message  = "Sending large graph of size",
                        category = UserWarning,
                        module   = "distributed.client")

def pygmt_regrid(da, lon, lat, grid_res=None, region=None, search_radius="200k"):
    mask = np.isfinite(da)
    df   = pd.DataFrame({"longitude": lon.values[mask].ravel(), 
                         "latitude" : lat.values[mask].ravel(),
                         "z"        : da.values[mask].ravel()})
    return pygmt.nearneighbor(data        = df,
                            spacing       = grid_res,
                            region        = region,
                            search_radius = search_radius)

def main(sim_name=None, year=None):
    # --- Setup ---
    mod_path = '/home/581/da1339/AFIM/src/AFIM/src'
    sys.path.insert(0, mod_path)
    from sea_ice_toolbox import SeaIceToolbox
    G_res      = "0.25/0.25"      #quarter degree (no tripole) only good for southern hemisphere
    pygmt_reg  = [-180,180,-90,0]
    afim_name  = sim_name
    dt0_str    = f"{year}-01-01"
    dtN_str    = f"{year}-12-31"
    # Load sim data
    SI_tools   = SeaIceToolbox(sim_name = afim_name,
                               dt0_str  = dt0_str,
                               dtN_str  = dtN_str,
                               P_log    = Path(Path.home(),"logs",f"reG_OSISAF_{afim_name}_{year}.log"))
    _, CICE    = SI_tools.load_processed_cice(zarr_CICE=True)
    CICE_ispd  = SI_tools.compute_ice_magnitude_from_ice_components_on_Bgrid(CICE, ivec_type="BT")
    CICE_SO    = CICE_ispd.isel(nj=SI_tools.hemisphere_dict['nj_slice'])
    TLAT       = CICE_SO["TLAT"]
    TLON       = CICE_SO["TLON"]
    aom2_name  = "AOM2-ERA5"
    SI_tools   = SeaIceToolbox(sim_name = aom2_name,
                               dt0_str  = dt0_str,
                               dtN_str  = dtN_str,
                               P_log    = Path(Path.home(),"logs",f"reG_OSISAF_{aom2_name}_{year}.log"))
    _, AOM2    = SI_tools.load_processed_cice(zarr_CICE=True)
    AOM2_ispd  = SI_tools.compute_ice_magnitude_from_ice_components_on_Bgrid(AOM2, ivec_type="BT")
    AOM2_SO    = AOM2_ispd.isel(nj=SI_tools.hemisphere_dict['nj_slice'])
    # Search OSI-SAF files
    D_search  = Path(Path.home(),"seaice","OSI_SAF","ice_drift_455m","sh")
    F_search  = f"ice_drift_sh_ease2-750_cdr-v1p0_24h-{year}*.nc"
    P_found   = sorted(D_search.rglob(F_search))
    DS        = []
    for F_ in P_found:
        SI_tools.logger.info(f"working on file {F_}")
        osisaf = xr.open_dataset(F_)
        dt_obs = pd.to_datetime(osisaf["time"].values[0])
        dt_str = dt_obs.strftime("%Y-%m-%d")
        # Compute observed speed
        SI_tools.logger.info(f"compute observed speeds in m/s")
        dt_sec   = (osisaf["t1"] - osisaf["t0"]).astype("timedelta64[s]").astype(float)
        dt_sec   = xr.where(dt_sec == 0, np.nan, dt_sec)
        disp_m   = np.sqrt(osisaf["dX"]**2 + osisaf["dY"]**2) * 1000  # km to m
        ispd_obs = (disp_m / dt_sec).squeeze()
        # Nearest-neighbour match
        SI_tools.logger.info(f"compute nearest neighbour")
        ispd_obs_py =  pygmt_regrid(ispd_obs, osisaf["lon"], osisaf["lat"], grid_res=G_res, region=pygmt_reg, search_radius="200k")
        # sim speeds for same day
        SI_tools.logger.info(f"re-grid sims")
        ispd_cice_py = pygmt_regrid(CICE_SO["ispd_BT"].sel(time=dt_str).load(), TLON, TLAT, grid_res=G_res, region=pygmt_reg, search_radius="100k")
        ispd_aom2_py = pygmt_regrid(AOM2_SO["ispd_BT"].sel(time=dt_str).load(), TLON, TLAT, grid_res=G_res, region=pygmt_reg, search_radius="100k")
        # Differences and RMSE
        SI_tools.logger.info(f"compute differences")
        d_cice = ispd_obs_py  - ispd_cice_py
        d_aom2 = ispd_obs_py  - ispd_aom2_py
        d_sims = ispd_aom2_py - ispd_cice_py
        SI_tools.logger.info(f"compute RMSE")
        rmse_cice = float(np.sqrt((d_cice ** 2).mean(skipna=True).compute()))
        rmse_aom2 = float(np.sqrt((d_aom2 ** 2).mean(skipna=True).compute()))
        rmse_sims = float(np.sqrt((d_sims ** 2).mean(skipna=True).compute()))
        # Package into daily dataset
        daily_ds = xr.Dataset(data_vars=dict(ispd_obs    = (("time", "ny", "nx"), ispd_obs_py.expand_dims("time").data),
                                             ispd_CICE   = (("time", "ny", "nx"), ispd_cice_py.expand_dims("time").data),
                                             ispd_AOM2   = (("time", "ny", "nx"), ispd_aom2_py.expand_dims("time").data),
                                             d_ispd_CICE = (("time", "ny", "nx"), d_cice.expand_dims("time").data),
                                             d_ispd_AOM2 = (("time", "ny", "nx"), d_aom2.expand_dims("time").data),
                                             d_ispd_sims = (("time", "ny", "nx"), d_sims.expand_dims("time").data),
                                             RMSE_CICE   = ("time", [rmse_cice]),
                                             RMSE_AOM2   = ("time", [rmse_aom2]),
                                             RMSE_sims   = ("time", [rmse_sims])),
                            coords = dict(time = ("time", [dt_obs]),
                                          lon  = (("ny", "nx"), ispd_cice_py.values),
                                          lat  = (("ny", "nx"), ispd_cice_py.values),))
        DS.append(daily_ds)
    DS_all = xr.concat(DS, dim="time")
    P_nc = Path(Path.home(),"seaice","OSI_SAF","ice_drift_455m",f"ispd_diffs_pygmt_nn_{afim_name}_{year}.nc")
    P_nc.parent.mkdir(parents=True, exist_ok=True)
    DS_all.to_netcdf(P_nc, mode="w", format="NETCDF4", encoding={var: {"zlib": True, "complevel": 4} for var in DS_all.data_vars})

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--sim" , type=str, required=False, help="Name of AFIM simulation", default="elps-min")
    parser.add_argument("--year", type=int, required=False, help="Year to process (e.g. 1994)", default=1994)
    args = parser.parse_args()
    main(sim_name=args.sim, year=args.year)