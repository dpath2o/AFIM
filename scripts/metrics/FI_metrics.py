import sys, os, argparse, pygmt
sys.path.insert(0, '/home/581/da1339/AFIM/src/AFIM/src')
from sea_ice_toolbox import SeaIceToolbox, SeaIceToolboxManager
from pathlib import Path
import xarray as xr
import numpy as np
import pandas as pd

def main(sim_name, ispd_thresh, ice_type, BorC2T_type, dt0_str, dtN_str, overwrite_zarr, overwrite_png):
    print(f"ice_type passed is: {ice_type}")
    load_vars = ['aice','tarea','hi','dvidtt','daidtt','dvidtd','daidtd'] # 'strength'
    P_log     = Path(Path.home(), "logs", f"metrics_{sim_name}_ispd_thresh{ispd_thresh}.log")
    SI_mgr    = SeaIceToolboxManager(P_log=P_log)
    SI_tools  = SI_mgr.get_toolbox(sim_name             = sim_name,
                                   dt0_str              = dt0_str,
                                   dtN_str              = dtN_str,
                                   ice_type             = ice_type,
                                   list_of_BorC2T       = BorC2T_type,
                                   ice_speed_threshold  = ispd_thresh,
                                   overwrite_zarr       = overwrite_zarr,
                                   overwrite_saved_figs = overwrite_png)
    # Load classified ice masks and data
    FI_day  = SI_tools.load_classified_ice(bin_days=False)['FI_mask']
    FI_bin  = SI_tools.load_classified_ice(bin_days=True)['FI_mask']
    FI_rol  = SI_tools.load_classified_ice(bin_days=False, roll_mean=True)['FI_mask']
    CICE_SO = SI_tools.load_cice_zarr(slice_hem=True, variables=load_vars)
    A       = CICE_SO['tarea'].isel(time=0)
    # Apply the mask to the data
    FI_daily = CICE_SO.where(FI_day)
    FI_rolly = CICE_SO.where(FI_rol)
    FI_binly = CICE_SO.where(FI_bin)
    # Create FI dictionaries for dy, rl, and bn
    FI_dy = SI_tools.fast_ice_metrics_data_dict(FI_day, FI_daily, A)
    FI_rl = SI_tools.fast_ice_metrics_data_dict(FI_rol, FI_rolly, A)
    FI_bn = SI_tools.fast_ice_metrics_data_dict(FI_bin, FI_binly, A)
    # Compute metrics for each FI type
    SI_tools.define_fast_ice_class_name(BorC2T_type = BorC2T_type , fast_ice_class_method = 'raw')
    FI_dy_name = SI_tools.FI_class
    SI_tools.define_fast_ice_class_name(BorC2T_type = BorC2T_type , fast_ice_class_method = 'rolling-mean')
    FI_roll_name = SI_tools.FI_class
    SI_tools.define_fast_ice_class_name(BorC2T_type = BorC2T_type , fast_ice_class_method = 'binary-days')
    FI_bin_name = SI_tools.FI_class
    FI_types = [(FI_dy_name, FI_dy), (FI_roll_name, FI_rl), (FI_bin_name, FI_bn)]
    for FI_name, FI_data in FI_types:
        P_mets_zarr = Path(SI_tools.D_ispd_thresh, f"{FI_name}_{SI_tools.metrics_name}.zarr")
        SI_tools.compute_sea_ice_metrics(FI_data, 
                                         ice_type       = "FI",
                                         P_mets_zarr    = P_mets_zarr,
                                         ice_area_scale = SI_tools.FIC_scale)
    SI_tools.client.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute FIA and FIP metrics, apply boolean mask, and plot spatial + temporal outputs.")
    parser.add_argument("--sim_name", type=str, required=True)
    parser.add_argument("--ispd_thresh", type=float, required=True)
    parser.add_argument("--ice_type", default="FI", help="either FI, PI, SI or MI")
    parser.add_argument("--BorC2T_type", default="Tc", help="must be Tc, Ta, Tb, Tx, B or BT")
    parser.add_argument("--start_date", default="1994-01-01", help="Start date (YYYY-MM-DD), which is then added to FI_days as the first center-date")
    parser.add_argument("--end_date", default="1999-12-31", help="End date (YYYY-MM-DD), will stop processing when this end_date-FI_days")
    parser.add_argument("--overwrite_zarr", action="store_true")
    parser.add_argument("--overwrite_png", action="store_true")
    parser.set_defaults(compute_boolean=False, overwrite_zarr=False, overwrite_png=False)
    args = parser.parse_args()

    main(args.sim_name,
         args.ispd_thresh,
         args.ice_type,
         args.BorC2T_type,
         args.start_date,
         args.end_date,
         args.overwrite_zarr,
         args.overwrite_png)

