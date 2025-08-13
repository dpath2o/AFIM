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

def main(sim_name=None, year=None):
    # --- Setup ---
    mod_path = '/home/581/da1339/AFIM/src/AFIM/src'
    sys.path.insert(0, mod_path)
    from sea_ice_toolbox import SeaIceToolbox, SeaIceToolboxManager    
    G_res      = "0.25/0.25"      #quarter degree (no tripole) only good for southern hemisphere
    pygmt_reg  = [-180,180,-90,0]
    afim_name  = sim_name
    dt0_str    = f"{year}-01-01"
    dtN_str    = f"{year}-12-31"
    # manager
    P_log       = Path(Path.home(), "logs", "OSI-SAF_reG.log")
    SI_tool_mgr = SeaIceToolboxManager(P_log=P_log)
    # Load sim data
    sim_tools   = SI_tool_mgr.get_toolbox(sim_name = afim_name, dt0_str = dt0_str, dtN_str = dtN_str)
    sim_tools.load_bgrid(slice_hem=True)
    TLAT       = sim_tools.G_t['lat']
    TLON       = sim_tools.G_t['lon']
    CICE       = sim_tools.load_cice_zarr()
    CICE_ispd  = sim_tools.compute_ice_magnitude_from_ice_components_on_Bgrid(CICE, ivec_type="BT")
    CICE_SO    = sim_tools.slice_hemisphere(CICE_ispd)
    sim_tools  = SI_tool_mgr.get_toolbox(sim_name = "AOM2-ERA5", dt0_str = dt0_str, dtN_str = dtN_str)
    AOM2       = sim_tools.load_cice_zarr()
    AOM2_ispd  = sim_tools.compute_ice_magnitude_from_ice_components_on_Bgrid(AOM2, ivec_type="BT")
    AOM2_SO    = sim_tools.slice_hemisphere(AOM2_ispd)
    ORAS         = xr.open_zarr("/g/data/gv90/da1339/SeaIce/CMEMS/0p083/daily/reG/CMEMS-ORAS-SI.zarr")
    ORAS["ispd"] = xr.apply_ufunc(np.hypot, ORAS.usi, ORAS.vsi, dask='allowed')
    ORAS_SO      = sim_tools.slice_hemisphere(ORAS)
    # Search OSI-SAF files
    D_search  = Path(Path.home(),"seaice","OSI_SAF","ice_drift_455m","sh")
    F_search  = f"ice_drift_sh_ease2-750_cdr-v1p0_24h-{year}*.nc"
    P_found   = sorted(D_search.rglob(F_search))
    DS        = []
    for F_ in P_found:
        sim_tools.logger.info(f"working on file {F_}")
        osisaf = xr.open_dataset(F_)
        dt_obs = pd.to_datetime(osisaf["time"].values[0])
        dt_str = dt_obs.strftime("%Y-%m-%d")
        # Compute observed speed
        sim_tools.logger.info(f"compute observed speeds in m/s")
        dt_sec   = (osisaf["t1"] - osisaf["t0"]).astype("timedelta64[s]").astype(float)
        dt_sec   = xr.where(dt_sec == 0, np.nan, dt_sec)
        disp_m   = np.sqrt(osisaf["dX"]**2 + osisaf["dY"]**2) * 1000  # km to m
        ispd_obs = (disp_m / dt_sec).squeeze()
        # Observed vector components (m/s)
        u_obs = (osisaf["dX"] * 1000.0 / dt_sec).squeeze()  # dX, dY are in km on the OSI-SAF grid
        v_obs = (osisaf["dY"] * 1000.0 / dt_sec).squeeze()
        # Nearest-neighbour match
        sim_tools.logger.info(f"compute nearest neighbour for OSI-SAF")
        ispd_obs_py = sim_tools.pygmt_regrid(ispd_obs, osisaf["lon"], osisaf["lat"], grid_res=G_res, region=pygmt_reg, search_radius="200k")
        u_obs_py    = sim_tools.pygmt_regrid(u_obs, osisaf["lon"], osisaf["lat"], grid_res=G_res, region=pygmt_reg, search_radius="200k")
        v_obs_py    = sim_tools.pygmt_regrid(v_obs, osisaf["lon"], osisaf["lat"], grid_res=G_res, region=pygmt_reg, search_radius="200k")
        # sim speeds for same day
        sim_tools.logger.info(f"re-grid CICE")
        ispd_cice_py = sim_tools.pygmt_regrid(CICE_SO["ispd_BT"].sel(time=dt_str).load(), TLON, TLAT, grid_res=G_res, region=pygmt_reg, search_radius="100k")
        sim_tools.logger.info(f"re-grid AOM2")
        ispd_aom2_py = sim_tools.pygmt_regrid(AOM2_SO["ispd_BT"].sel(time=dt_str).load(), TLON, TLAT, grid_res=G_res, region=pygmt_reg, search_radius="100k")
        sim_tools.logger.info(f"re-grid ORAS")
        ispd_oras_py = sim_tools.pygmt_regrid(ORAS_SO["ispd"].sel(time=dt_str).load(), TLON, TLAT, grid_res=G_res, region=pygmt_reg, search_radius="100k")
        # sim components
        sim_tools.logger.info(f"load model velocity components")
        CICE_u = CICE_SO["uvel"].sel(time=dt_str).load()
        CICE_v = CICE_SO["vvel"].sel(time=dt_str).load()
        AOM2_u = AOM2_SO["uvel"].sel(time=dt_str).load()
        AOM2_v = AOM2_SO["vvel"].sel(time=dt_str).load()
        ORAS_u = ORAS_SO["usi"].sel(time=dt_str).load()
        ORAS_v = ORAS_SO["vsi"].sel(time=dt_str).load() 
        u_cice_py = sim_tools.pygmt_regrid(CICE_u, TLON, TLAT, grid_res=G_res, region=pygmt_reg, search_radius="100k")
        v_cice_py = sim_tools.pygmt_regrid(CICE_v, TLON, TLAT, grid_res=G_res, region=pygmt_reg, search_radius="100k")
        u_aom2_py = sim_tools.pygmt_regrid(AOM2_u, TLON, TLAT, grid_res=G_res, region=pygmt_reg, search_radius="100k")
        v_aom2_py = sim_tools.pygmt_regrid(AOM2_v, TLON, TLAT, grid_res=G_res, region=pygmt_reg, search_radius="100k")
        u_oras_py = sim_tools.pygmt_regrid(ORAS_u, ORAS_SO["lon"], ORAS_SO["lat"], grid_res=G_res, region=pygmt_reg, search_radius="100k")
        v_oras_py = sim_tools.pygmt_regrid(ORAS_v, ORAS_SO["lon"], ORAS_SO["lat"], grid_res=G_res, region=pygmt_reg, search_radius="100k")
        # compute cosine similarity
        cos_cice = sim_tools.cosine_vector_similarity(u_obs_py, v_obs_py, u_cice_py, v_cice_py)
        cos_aom2 = sim_tools.cosine_vector_similarity(u_obs_py, v_obs_py, u_aom2_py, v_aom2_py)
        cos_oras = sim_tools.cosine_vector_similarity(u_obs_py, v_obs_py, u_oras_py, v_oras_py)
        # Optional: angular difference in degrees
        ang_cice = np.degrees(sim_tools.vector_angle_diff(u_obs_py, v_obs_py, u_cice_py, v_cice_py))
        ang_aom2 = np.degrees(sim_tools.vector_angle_diff(u_obs_py, v_obs_py, u_aom2_py, v_aom2_py))
        ang_oras = np.degrees(sim_tools.vector_angle_diff(u_obs_py, v_obs_py, u_oras_py, v_oras_py))
        # Spatial means for quick daily scalar metrics
        cos_cice_mean = float(cos_cice.mean(skipna=True).compute())
        cos_aom2_mean = float(cos_aom2.mean(skipna=True).compute())
        cos_oras_mean = float(cos_oras.mean(skipna=True).compute())
        ang_cice_mean = float(np.nanmean(ang_cice))  # careful: not periodic-mean safe, but ok as a rough diagnostic
        ang_aom2_mean = float(np.nanmean(ang_aom2))
        ang_oras_mean = float(np.nanmean(ang_oras))
        # Differences and RMSE
        sim_tools.logger.info(f"compute differences")
        d_cice = ispd_cice_py - ispd_obs_py
        d_aom2 = ispd_aom2_py - ispd_obs_py 
        d_oras = ispd_oras_py - ispd_obs_py
        sim_tools.logger.info(f"compute RMSE")
        rmse_cice = float(np.sqrt((d_cice ** 2).mean(skipna=True).compute()))
        rmse_aom2 = float(np.sqrt((d_aom2 ** 2).mean(skipna=True).compute()))
        rmse_oras = float(np.sqrt((d_oras ** 2).mean(skipna=True).compute()))
        # Package into daily dataset
        daily_ds = xr.Dataset(data_vars=dict(ispd_obs       = (("time", "ny", "nx"), ispd_obs_py.expand_dims("time").data),
                                             ispd_CICE      = (("time", "ny", "nx"), ispd_cice_py.expand_dims("time").data),
                                             ispd_AOM2      = (("time", "ny", "nx"), ispd_aom2_py.expand_dims("time").data),
                                             ispd_ORAS      = (("time", "ny", "nx"), ispd_oras_py.expand_dims("time").data),
                                             d_ispd_CICE    = (("time", "ny", "nx"), d_cice.expand_dims("time").data),
                                             d_ispd_AOM2    = (("time", "ny", "nx"), d_aom2.expand_dims("time").data),
                                             d_ispd_ORAS    = (("time", "ny", "nx"), d_oras.expand_dims("time").data),
                                             cos_CICE       = (("time", "ny", "nx"), cos_cice.expand_dims("time").data),
                                             cos_AOM2       = (("time", "ny", "nx"), cos_aom2.expand_dims("time").data),
                                             cos_ORAS       = (("time", "ny", "nx"), cos_oras.expand_dims("time").data),
                                             ang_CICE_deg   = (("time", "ny", "nx"), ang_cice.expand_dims("time").data),
                                             ang_AOM2_deg   = (("time", "ny", "nx"), ang_aom2.expand_dims("time").data),
                                             ang_ORAS_deg   = (("time", "ny", "nx"), ang_oras.expand_dims("time").data),
                                             COS_CICE_mean  = ("time", [cos_cice_mean]),
                                             COS_AOM2_mean  = ("time", [cos_aom2_mean]),
                                             COS_ORAS_mean  = ("time", [cos_oras_mean]),
                                             ANG_CICE_mean  = ("time", [ang_cice_mean]),
                                             ANG_AOM2_mean  = ("time", [ang_aom2_mean]),
                                             ANG_ORAS_mean  = ("time", [ang_oras_mean]),
                                             RMSE_CICE      = ("time", [rmse_cice]),
                                             RMSE_AOM2      = ("time", [rmse_aom2]),
                                             RMSE_ORAS      = ("time", [rmse_oras])),
                            coords = dict(time = ("time", [dt_obs]),
                                          lon  = (("ny", "nx"), ispd_cice_py.values),
                                          lat  = (("ny", "nx"), ispd_cice_py.values),))
        DS.append(daily_ds)
    DS_all = xr.concat(DS, dim="time")
    P_nc   = Path(Path.home(),"seaice","OSI_SAF","ice_drift_455m",f"ispd_diffs_pygmt_nn_{afim_name}_{year}.nc")
    P_nc.parent.mkdir(parents=True, exist_ok=True)
    DS_all.to_netcdf(P_nc, mode="w", format="NETCDF4", encoding={var: {"zlib": True, "complevel": 4} for var in DS_all.data_vars})

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--sim" , type=str, required=False, help="Name of AFIM simulation", default="elps-min")
    parser.add_argument("--year", type=int, required=False, help="Year to process (e.g. 1994)", default=1994)
    args = parser.parse_args()
    main(sim_name=args.sim, year=args.year)