import sys, os, glob
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
    SI_tools = SeaIceToolbox(sim_name             = sim_name,
                            client               = None,
                            dt0_str              = dt0_str,
                            dtN_str              = dtN_str,
                            ice_speed_threshold  = 5e-4,
                            ice_speed_type       = "ispd_BT",
                            ice_type             = "FI_BT",
                            overwrite_zarr       = False,
                            save_new_figs        = True,
                            show_figs            = True,
                            overwrite_saved_figs = False)
    
    pattern = str(Path(SI_tools.D_sim, f"FI-sim-TS_{sim_name}_*.nc"))
    files = sorted(glob.glob(pattern))
    FI_sim = xr.open_mfdataset(files, engine="netcdf4")
    FI_obs_all = SI_tools.load_AF2020db().isel(nj=SI_tools.hemisphere_dict['nj_slice'])
    FI_diff_list = []
    for i in range(len(FI_obs_all.t_FI_obs)):
        t_idx_sim   = i+4   # FI_sim starts on "2000-01-01"
        t_idx_obs   = i     # FI_obs_all starts on "2000-03-01"
        dt_str      = pd.Timestamp(FI_sim.isel(t_FI_obs=t_idx_sim).t_FI_obs.values).strftime("%Y-%m-%d")
        FI_sim_aice = FI_sim.isel(t_FI_obs=t_idx_sim)['aice']
        FI_sim_mask = (FI_sim_aice > 0).astype(int)
        FI_obs_FI   = FI_obs_all.isel(t_FI_obs=t_idx_obs)['FI']
        FI_obs_mask = (FI_obs_FI > 0).astype(int)
        FI_diff     = FI_obs_mask - FI_sim_mask
        FI_diff = FI_diff.assign_coords({"lon"      : FI_sim["lon"],
                                        "lat"      : FI_sim["lat"],
                                        "t_FI_obs" : FI_sim.t_FI_obs.isel(t_FI_obs=t_idx_sim)})
        FI_diff_list.append(FI_diff)
        SI_tools.pygmt_map_plot_one_var(FI_diff, "FI_diff",
                                    plot_regions   = 2,
                                    time_stamp     = dt_str,
                                    tit_str        = dt_str,
                                    plot_GI        = True,
                                    cmap           = "cmocean/balance",
                                    series         = [-1,1],
                                    reverse        = False,
                                    cbar_label     = "binary fast ice grid cell difference (obs-sim)",
                                    cbar_units     = "unitless",
                                    lon_coord_name = 'lon',
                                    lat_coord_name = 'lat',
                                    var_sq_size    = 0.05,
                                    GI_sq_size     = 0.02,
                                    GI_fill_color  = "green",
                                    overwrite_fig  = False,
                                    show_fig       = False)
    FI_diff_all      = xr.concat(FI_diff_list, dim="t_FI_obs")
    FI_diff_all.name = "FI_diff"
    FI_diff_all.to_netcdf(Path(SI_tools.D_sim,f"FI_diff_with_obs_2000-2018.nc"))
    SI_tools.client.close()

if __name__ == '__main__':
    from multiprocessing import freeze_support
    freeze_support()  # safe on all platforms, even if not frozen
    main()    