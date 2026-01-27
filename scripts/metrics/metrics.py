import sys, os, argparse, pygmt
sys.path.insert(0, '/home/581/da1339/AFIM/src/AFIM/src')
from sea_ice_toolbox import SeaIceToolbox, SeaIceToolboxManager
from pathlib import Path
import xarray as xr
import numpy as np
import pandas as pd

def main(sim_name, ispd_thresh, ice_type, BorC2T_type, dt0_str, dtN_str, overwrite_zarr, overwrite_png):
    print(f"ice_type passed is: {ice_type}")
    load_vars = ['aice','tarea','hi','uvel','vvel','strength','dvidtt','daidtt','dvidtd','daidtd']
    P_log     = Path(Path.home(), "logs", f"metrics_{sim_name}_ispd_thresh{ispd_thresh}.log")
    mgr       = SeaIceToolboxManager(P_log=P_log)
    tb        = mgr.get_toolbox(sim_name             = sim_name,
                                dt0_str              = dt0_str,
                                dtN_str              = dtN_str,
                                ice_type             = ice_type,
                                list_of_BorC2T       = BorC2T_type,
                                ice_speed_threshold  = ispd_thresh,
                                overwrite_zarr       = overwrite_zarr,
                                overwrite_saved_figs = overwrite_png)
    # Load classified ice masks and data
    tb.define_ice_mask_name(ice_type=ice_type)
    I_day   = tb.load_classified_ice(class_method="raw")[tb.mask_name]
    CICE_SO = tb.load_cice_zarr(slice_hem=True, variables=load_vars)
    if not ice_type=="SI":
        I_bin = tb.load_classified_ice(class_method="binary-days")[tb.mask_name]
        I_rol = tb.load_classified_ice(class_method="rolling-mean")[tb.mask_name]
    A = CICE_SO['tarea'].isel(time=0)
    # Apply the mask to the data
    I_daily = CICE_SO.where(I_day)
    if not ice_type=="SI":
        I_rolly = CICE_SO.where(I_rol)
        I_binly = CICE_SO.where(I_bin)
    # Create FI dictionaries for dy, rl, and bn
    I_dy = tb.metrics_data_dict(I_day, I_daily, A)
    if not ice_type=="SI":
        I_rl = tb.metrics_data_dict(I_rol, I_rolly, A)
        I_bn = tb.metrics_data_dict(I_bin, I_binly, A)
    # Compute metrics for each FI type
    tb.define_metrics_zarr(class_method = "raw")
    tb.compute_sea_ice_metrics(I_dy, tb.D_mets_zarr) 
    if not ice_type=="SI":
        tb.define_metrics_zarr(class_method = "rolling-mean")
        tb.compute_sea_ice_metrics(I_rl, tb.D_mets_zarr) 
        tb.define_metrics_zarr(class_method = "binary-days")
        tb.compute_sea_ice_metrics(I_bn, tb.D_mets_zarr) 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute FIA and FIP metrics, apply boolean mask, and plot spatial + temporal outputs.")
    parser.add_argument("--sim_name", type=str, required=True)
    parser.add_argument("--ispd_thresh", type=float, required=True)
    parser.add_argument("--ice_type", default="FI", help="either FI, PI, SI or MI")
    parser.add_argument("--BorC2T_type", default="Tc", help="must be Tc, Ta, Tb, Tx, B or BT")
    parser.add_argument("--start_date", default="1994-01-01", help="Start date (YYYY-MM-DD), which is then added to I_days as the first center-date")
    parser.add_argument("--end_date", default="1999-12-31", help="End date (YYYY-MM-DD), will stop processing when this end_date-I_days")
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

