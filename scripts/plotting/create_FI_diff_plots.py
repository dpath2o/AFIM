import sys, os, glob, pygmt
import numpy as np
import pandas as pd
import xarray as xr
from pathlib import Path

def main():
    mod_path = '/home/581/da1339/AFIM/src/AFIM/src'
    sys.path.insert(0, mod_path)
    from sea_ice_toolbox import SeaIceToolbox
    sim_name = "elps-min"
    dt0_str  = "2000-01-01"
    dtN_str  = "2018-12-31"
    years    = range(2000, 2019)
    SI_tools = SeaIceToolbox(sim_name            = sim_name,
                            client               = None,
                            dt0_str              = dt0_str,
                            dtN_str              = dtN_str,
                            ice_speed_threshold  = 5e-4,
                            ice_speed_type       = "ispd_BT",
                            ice_type             = "FI_BT",
                            overwrite_zarr       = False,
                            save_new_figs        = True,
                            show_figs            = False,
                            overwrite_saved_figs = True)
    P_FI_diff    = Path(SI_tools.D_sim,f"FI-diff_obs-elps-min_{dt0_str[:3]}-{dtN_str[:3]}.nc")
    FI_raw, CICE = SI_tools.load_processed_cice( zarr_CICE = True )
    FI_bool      = SI_tools.boolean_fast_ice( FI_raw['FI_mask'] , window=11, min_count=9).compute()
    D_obs        = Path(SI_tools.AF_FI_dict['D_AF2020_db_org'])
    P_orgs       = sorted(D_obs.glob("FastIce_70_*.nc"))
    FI_obs       = xr.open_mfdataset(P_orgs, engine='netcdf4', combine='by_coords')
    FI_obs       = FI_obs.chunk({'time':1})
    def threshold_mask(block):
        return xr.where(block >= 4, 1, 0)
    FI_obs_mask              = FI_obs['Fast_Ice_Time_series'].map_blocks(threshold_mask)
    FI_obs_ts                = FI_obs_mask.astype(int)
    CICE_coar                = SI_tools.coarsen_and_align_simulated_FI_to_observed_FI(FI_bool, FI_obs_ts)
    FI_obs_algn, FI_sim_algn = xr.align(FI_obs, CICE_coar, join="inner")
    FI_diff_list = []
    time_list    = []
    for i in range(len(FI_obs_algn['time'].values)):
        dt     = pd.Timestamp(FI_obs_algn.isel(time=i).time.values)
        dt_str = f"{dt.year}-{dt.month:02d}-{dt.day:02d}"
        print(f"\nprocessing {dt_str}")
        print(f"slicing and re-gridding fast ice from observations")
        FI_obs_slc  = FI_obs_algn.isel(time=i)
        FI_obs_mask = xr.where(FI_obs_slc['Fast_Ice_Time_series'] >= 4, 1.0, 0.0).values.flatten()
        FI_obs_lat  = FI_obs_slc.latitude.values.flatten()
        FI_obs_lon  = FI_obs_slc.longitude.values.flatten()
        df_obs      = pd.DataFrame({"longitude": FI_obs_lon,
                                    "latitude" : FI_obs_lat,
                                    "z"        : FI_obs_mask})
        da_obs_reG  = pygmt.nearneighbor(data          = df_obs,
                                        region        = [0, 360, -90, -50],
                                        spacing       = "0.1/0.1",
                                        search_radius = "5k")
        LON,LAT     = np.meshgrid(da_obs_reG.lon, da_obs_reG.lat)
        LON         = (('lat','lon'), LON)
        LAT         = (('lat','lon'), LAT)
        da_obs_reG  = da_obs_reG.assign_coords({'LON':LON,'LAT':LAT})
        print(f"slicing and re-gridding fast ice from {sim_name}")
        FI_mask     = FI_sim_algn.chunk({'time': 1, 'nj': 540, 'ni': 1440})  # or smaller spatial chunks too
        FI_bool_slc = FI_mask.isel(time=i).compute()
        FI_bool_dat = FI_bool_slc.values.astype(int).flatten()
        FI_bool_lat = FI_bool.TLAT.values.flatten()
        FI_bool_lon = FI_bool.TLON.values.flatten()
        df_bool     = pd.DataFrame({"longitude": FI_bool_lon,
                                    "latitude" : FI_bool_lat,
                                    "z"        : FI_bool_dat})
        da_bool_reG = pygmt.nearneighbor(data          = df_bool,
                                        region        = [0, 360, -90, -50],
                                        spacing       = "0.1/0.1",
                                        search_radius = "50k")
        LON,LAT     = np.meshgrid(da_bool_reG.lon, da_bool_reG.lat)
        LON         = (('lat','lon'), LON)
        LAT         = (('lat','lon'), LAT)
        da_bool_reG = da_bool_reG.assign_coords({'LON':LON,'LAT':LAT})
        print(f"creating categorised 'difference' array")
        valid_mask  = (da_obs_reG  > 0) | (da_bool_reG >  0)
        only_obs    = (da_obs_reG  > 0) & (da_bool_reG <= 0)
        only_sim    = (da_bool_reG > 0) & (da_obs_reG  <= 0)
        agree       = (da_obs_reG  > 0) & (da_bool_reG >  0)
        FI_diff     = xr.full_like(da_obs_reG, np.nan, dtype=float)
        FI_diff     = xr.where(only_obs,  2, FI_diff)  # red
        FI_diff     = xr.where(only_sim,  1, FI_diff)  # blue
        FI_diff     = xr.where(agree,     0, FI_diff)  # green
        LON,LAT     = np.meshgrid(FI_diff.lon,FI_diff.lat)
        LON         = (('lat','lon'), LON)
        LAT         = (('lat','lon'), LAT)
        FI_diff     = FI_diff.assign_coords({'LON':LON,'LAT':LAT})
        print(f"plotting")
        SI_tools.pygmt_map_plot_one_var(FI_diff, 'FI_diff',
                                        plot_regions   = 8,
                                        time_stamp     = dt_str,
                                        tit_str        = dt_str,
                                        plot_GI        = False,
                                        cbar_label     = "fast ice difference (obs-sim)",
                                        cbar_units     = "",
                                        extend_cbar    = False,
                                        lon_coord_name = "LON",
                                        lat_coord_name = "LAT",
                                        var_sq_size    = 0.125,
                                        GI_sq_size     = 0.075,
                                        GI_fill_color  = "yellow",
                                        water_color='black',
                                        plot_bathymetry=False,
                                        show_fig=False)
        print(f"appending to dictionary")
        curr_time   = FI_obs_algn['time'].isel(time=i).values
        time_list.append(curr_time)
        FI_diff_expanded = FI_diff.expand_dims(time=[curr_time])
        FI_diff_list.append(FI_diff_expanded)
    print(f"concatenating")
    FI_diff_3D = xr.concat(FI_diff_list, dim="time")
    FI_diff_3D.name = "FI_diff"
    print(f"writing to disk")
    FI_diff_3D.to_netcdf(Path(SI_tools.D_sim,"FI_diff_time_series.nc"))

if __name__ == '__main__':
    from multiprocessing import freeze_support
    freeze_support()  # safe on all platforms, even if not frozen
    main()    