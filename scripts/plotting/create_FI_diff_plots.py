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
    P_FI_diff = Path("/g/data/gv90/da1339/afim_output/elps-min/FI-diff_obs-elps-min_2000-2018.nc")
    if P_FI_diff.exists():
        FI_diff_computed = xr.open_dataset(P_FI_diff)
    else:    
        CICE_coar  = xr.open_mfdataset(f"/g/data/gv90/da1339/afim_output/elps-min/FI-sim-TS_elps-min_*.nc")
        FI_obs     = SI_tools.load_AF2020db()
        FI_obs_SO  = FI_obs.isel(nj=SI_tools.hemisphere_dict['nj_slice'])
        FI_obs_bin = FI_obs_SO['FI'].astype(bool).where(FI_obs_SO['FI'].notnull())
        FI_sim_bin = CICE_coar['FI_mask'].astype(bool).where(CICE_coar['FI_mask']==1)#FI_obs_SO['FI'].notnull())
        FI_obs_bin, FI_sim_bin = xr.align(FI_obs_bin, FI_sim_bin, join="inner")
        FI_diff  = xr.full_like(FI_obs_bin, np.nan, dtype=float)
        only_obs = (FI_obs_bin == True) & (FI_sim_bin != True)
        only_sim = (FI_sim_bin == True) & (FI_obs_bin != True)
        agree    = (FI_obs_bin == FI_sim_bin)
        FI_diff  = xr.where(only_obs,     2, FI_diff)  # observation-only
        FI_diff  = xr.where(only_sim,     1, FI_diff)  # simulation-only
        FI_diff  = xr.where(agree,        0, FI_diff)  # agreement
        FI_diff.name = "FI_diff"
        FI_diff_computed = FI_diff.chunk({"t_FI_obs": 13, "nj": 540, "ni": 1440}).compute()
        FI_diff_computed.to_netcdf(f"/g/data/gv90/da1339/afim_output/elps-min/FI-diff_obs-elps-min_2000-2018.nc")
    for i in range(len(FI_diff_computed.t_FI_obs)):
        dt  = pd.Timestamp(FI_diff_computed.isel(t_FI_obs=i)['t_FI_obs'].values)
        dt_str = f"{dt:%Y-%m-%d}"
        SI_tools.pygmt_map_plot_one_var(FI_diff_computed['FI_diff'].isel(t_FI_obs=i), 'FI_diff',
                                        plot_regions   = 8,
                                        time_stamp     = dt_str,
                                        tit_str        = dt_str,
                                        plot_GI        = True,
                                        cbar_label     = "fast ice mask difference (obs-sim)",
                                        cbar_units     = "",
                                        extend_cbar    = False,
                                        lon_coord_name = "TLON",
                                        lat_coord_name = "TLAT",
                                        var_sq_size    = 0.125,
                                        GI_sq_size     = 0.075,
                                        GI_fill_color  = "yellow")

if __name__ == '__main__':
    from multiprocessing import freeze_support
    freeze_support()  # safe on all platforms, even if not frozen
    main()    