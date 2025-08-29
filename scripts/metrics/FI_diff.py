import os, sys, argparse, xarray as xr, pandas as pd, numpy as np
from pathlib import Path
from datetime import datetime
import xesmf as xe

def define_regular_grid(grid_res, region=[-180,180,-90,0]):
    lon_min, lon_max, lat_min, lat_max = region
    lon_regular = np.arange(lon_min, lon_max + grid_res, grid_res)
    lat_regular = np.arange(lat_min, lat_max + grid_res, grid_res)
    return xr.Dataset({"lon": (["lon"], lon_regular),
                       "lat": (["lat"], lat_regular)})


def main():
    mod_path = '/home/581/da1339/AFIM/src/AFIM/src'
    sys.path.insert(0, mod_path)
    from sea_ice_toolbox import SeaIceToolbox, SeaIceToolboxManager
    parser = argparse.ArgumentParser()
    parser.add_argument("year", type=int)
    parser.add_argument("--sim_name", type=str, default="elps-min")
    parser.add_argument("--grid_res", type=float, default=0.15)
    parser.add_argument("--search_radius", type=float, default=12)
    parser.add_argument("--xesmf", action="store_true", help="Use xESMF for regridding (default: False)")
    reG_method  = "bilinear"
    args        = parser.parse_args()
    year        = args.year
    srch_rds    = args.search_radius
    use_xesmf   = args.xesmf
    G_res       = f"{args.grid_res}/{args.grid_res}"
    G_region    = [-180,180,-90,0]
    G_regular   = define_regular_grid(args.grid_res, G_region)
    year_group  = f"{year}"
    sim_name    = args.sim_name
    dt0_str     = f"{year}-01-01"
    dtN_str     = f"{year}-12-31"
    P_log       = Path(Path.home(), "logs", f"FI_diff_{sim_name}_{year}.log")
    SI_tool_mgr = SeaIceToolboxManager(P_log=P_log)
    SI_tools    = SI_tool_mgr.get_toolbox(sim_name=sim_name, dt0_str=dt0_str, dtN_str=dtN_str)
    G_t         = SI_tools.define_cice_grid(grid_type='t', mask=True, build_grid_corners=True, slice_hem=True)
    P_zarr      = Path(SI_tools.D_zarr, "FI_diff.zarr")
    FI_bin      = SI_tools.load_classified_ice(bin_days=True)['FI_mask']
    CICE_SO     = SI_tools.load_cice_zarr(slice_hem=True, variables=['aice','tarea'])
    A_SO        = CICE_SO['tarea'].isel(time=0)
    FI_binly    = CICE_SO.where(FI_bin)
    FI_aice     = FI_binly['aice'].coarsen(time=15,boundary='trim').mean()
    FI_mask_mod = xr.where(FI_aice>0,1,0)
    P_AF2020    = Path(SI_tools.AF_FI_dict['D_AF2020_db_org'],f"FastIce_70_{year}.nc")
    SI_tools.logger.info(f"opening AF2020 FI dataset: {P_AF2020}")
    FI_obs      = xr.open_mfdataset(P_AF2020, engine='netcdf4', combine='by_coords')
    FI_mask     = xr.where(FI_obs['Fast_Ice_Time_series'] >= 4, 1.0, 0.0)
    lat_obs     = FI_obs.latitude
    lon_obs     = FI_obs.longitude
    if use_xesmf:
        #lat_b,lon_b = SI_tools.build_AF2020_grid_corners(lat_c.values, lon_c.values)
        G_obs       = xr.Dataset({"lat"   : (("Y", "X"), lat_obs.values),
                                "lon"   : (("Y", "X"), lon_obs.values),})
                                # "lat_b" : (("Y_b", "X_b"), lat_b),
                                # "lon_b" : (("Y_b", "X_b"), lon_b)})
        
        P_wghts_mod = Path(f"/g/data/gv90/da1339/grids/weights/CICE_to_regular-grid-{args.grid_res}_{reG_method}_extrap-ns2d.nc")
        if P_wghts_mod.exists(): 
            SI_tools.logger.info(f"re-using regridding weights for CICE: {P_wghts_mod}")
            reuse_weights = True
        else:
            SI_tools.logger.info(f"creating regridding weights for CICE: {P_wghts_mod}")
            reuse_weights = False
        SI_tools.logger.info("creating re-gridding object for FI-CICE")
        reG_mod = xe.Regridder(G_t, G_regular,
                               method            = reG_method,
                               periodic          = True,
                               ignore_degenerate = True,
                               extrap_method     = "nearest_s2d",
                               reuse_weights     = reuse_weights,
                               filename          = P_wghts_mod)
        P_wghts_obs = Path(f"/g/data/gv90/da1339/grids/weights/AF2020db_to_regular-grid-{args.grid_res}_{reG_method}.nc")
        if P_wghts_obs.exists(): 
            SI_tools.logger.info(f"re-using regridding weights for AF2020: {P_wghts_obs}")
            reuse_weights = True
        else:
            SI_tools.logger.info(f"creating regridding weights for AF2020: {P_wghts_obs}")
            reuse_weights = False
        SI_tools.logger.info("creating re-gridding object for FI-obs")
        reG_obs = xe.Regridder(G_obs, G_regular,
                               method            = reG_method,
                               periodic          = False,
                               ignore_degenerate = True,
                               reuse_weights     = reuse_weights,
                               filename          = P_wghts_obs)   
    FI_mod_reT  = FI_mask_mod.reindex(time=FI_obs.time, method="nearest")
    yr_slice    = []
    for i in range(FI_obs.dims['time']):
        ti = pd.to_datetime(FI_obs.isel(time=i).time.values)
        SI_tools.logger.info(f"working on date stamp: {ti}")
        fi_obs_slice = FI_mask.isel(time=i)
        fi_mod_slice = FI_mod_reT.isel(time=i)
        # --- Regrid observed and model ---
        if use_xesmf:
            FI_obs_reG = reG_obs(fi_obs_slice)
            FI_mod_reG = reG_mod(fi_mod_slice)
        else:
            FI_obs_reG = SI_tools.pygmt_regrid(fi_obs_slice, lon_obs, lat_obs, grid_res=G_res, region=G_region, search_radius=f"3k")
            FI_mod_reG = SI_tools.pygmt_regrid(fi_mod_slice, G_t['lon'], G_t['lat'], grid_res=G_res, region=G_region, search_radius=f"{srch_rds}k")
        # --- Mask invalid cells ---
        both_valid = xr.where(xr.ufuncs.isfinite(FI_obs_reG) & xr.ufuncs.isfinite(FI_mod_reG), 1, np.nan)
        FI_diff_i = (FI_mod_reG - FI_obs_reG) * both_valid
        FI_diff_i = FI_diff_i.assign_coords(time=ti)
        # Compute area and FIA
        dummy_da = xr.ones_like(FI_mod_reG)
        A_SO     = SI_tools.compute_regular_grid_area(dummy_da)
        FIA_obs_8, FIA_obs_tot_8 = SI_tools.compute_sector_FIA(FI_obs_reG, A_SO, SI_tools.Ant_8sectors)
        FIA_mod_8, FIA_mod_tot_8 = SI_tools.compute_sector_FIA(FI_mod_reG, A_SO, SI_tools.Ant_8sectors, GI_area=SI_tools.use_gi)
        SI_tools.logger.info(f"FIA obs total 8: {FIA_obs_tot_8}")
        SI_tools.logger.info(f"FIA mod total 8: {FIA_mod_tot_8}")
        ds_i = xr.Dataset(
            {
                "FI_diff": FI_diff_i,
                "FI_obs_reG": FI_obs_reG,
                "FI_mod_reG": FI_mod_reG,
                "FIA_obs_8sec": (["sector"], FIA_obs_8.values, {"sector": FIA_obs_8.sector.values.tolist()}),
                "FIA_mod_8sec": (["sector"], FIA_mod_8.values, {"sector": FIA_mod_8.sector.values.tolist()}),
                "FIA_obs_tot_8": ((), FIA_obs_tot_8),
                "FIA_mod_tot_8": ((), FIA_mod_tot_8),
            },
            coords={"time": [ti]},
        )
        yr_slice.append(ds_i)
    # Concatenate all time steps for the year
    ds_year = xr.concat(yr_slice, dim="time")
    P_zarr.parent.mkdir(parents=True, exist_ok=True)
    for k, v in list(ds_year.attrs.items()):
        if isinstance(v, (xr.DataArray, np.ndarray)):
            ds_year.attrs[k] = v.tolist()
        elif not isinstance(v, (str, int, float, list, dict, bool, type(None))):
            ds_year.attrs.pop(k)
    ds_year.to_zarr(P_zarr, group=year_group, mode="w")

if __name__ == '__main__':
    from multiprocessing import freeze_support
    freeze_support()  # safe on all platforms, even if not frozen
    main()  

