import sys, os
from datetime import datetime
import pandas as pd
import numpy as np

def main():
    mod_path = '/home/581/da1339/AFIM/src/AFIM/src'
    sys.path.insert(0, mod_path)
    from sea_ice_toolbox import SeaIceToolbox# Extend path to AFIM modules
    sim_name = "elps-min"
    dt0_str  = "2000-01-01"
    dtN_str  = "2018-12-31"
    SI_tools = SeaIceToolbox(sim_name             = sim_name,
                             dt0_str              = dt0_str,
                             dtN_str              = dtN_str,
                             save_new_figs        = True,
                             show_figs            = False,
                             overwrite_saved_figs = False)
    DS, CICE = SI_tools.load_processed_cice(zarr_CICE=True)
    FI_obs   = SI_tools.load_AF2020db()
    DS_bool = SI_tools.boolean_fast_ice(DS['FI_mask'])
    for i in range(len(DS_bool.time)):
        dt     = pd.Timestamp(DS_bool.time[i].values)
        dt_str = dt.strftime('%Y-%m-%d')
        print(f"Plotting FI_mask for {dt_str} ...")
        SI_tools.pygmt_map_plot_one_var(DS_bool.isel(time=i), 'FI_mask',
                                        plot_regions   = 8,
                                        time_stamp     = dt_str,
                                        tit_str        = dt_str,
                                        plot_GI        = True,
                                        cmap           = "cmocean/amp",
                                        series         = [0, 1],
                                        reverse        = False,
                                        cbar_label     = "fast ice mask",
                                        cbar_units     = "",
                                        extend_cbar    = False,
                                        lon_coord_name = "TLON",
                                        lat_coord_name = "TLAT",
                                        var_sq_size    = 0.125,
                                        GI_sq_size     = 0.075,
                                        GI_fill_color  = "green",)

if __name__ == '__main__':
    from multiprocessing import freeze_support
    freeze_support()  # safe on all platforms, even if not frozen
    main()  
