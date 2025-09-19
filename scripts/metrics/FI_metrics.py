import sys, os, argparse, pygmt
sys.path.insert(0, '/home/581/da1339/AFIM/src/AFIM/src')
from sea_ice_toolbox import SeaIceToolbox, SeaIceToolboxManager
from pathlib import Path
import xarray as xr
import numpy as np
import pandas as pd

def main(sim_name, ispd_thresh, ice_type, dt0_str, dtN_str, smooth_FIA_days, overwrite_zarr, overwrite_png):
    print(f"ice_type passed is: {ice_type}")
    load_vars = ['aice','tarea','hi','strength','dvidtt','daidtt','dvidtd','daidtd']
    ispd_str  = f"{ispd_thresh:.1e}".replace("e-0", "e-")
    P_log     = Path(Path.home(), "logs", f"metrics_{sim_name}_ispd_thresh{ispd_thresh}.log")
    SI_mgr    = SeaIceToolboxManager(P_log=P_log)
    SI_tools  = SI_mgr.get_toolbox(sim_name             = sim_name,
                                   dt0_str              = dt0_str,
                                   dtN_str              = dtN_str,
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
    FI_types = [('FI_BT', FI_dy), ('FI_BT_roll', FI_rl), ('FI_BT_bin', FI_bn)]
    for prefix, FI_data in FI_types:
        P_mets_zarr = Path(SI_tools.D_ispd_thresh, f"{prefix}_mets.zarr")
        SI_tools.compute_sea_ice_metrics(FI_data, 
                                         ice_type="FI",
                                         P_mets_zarr=P_mets_zarr,
                                         ice_area_scale=SI_tools.FIC_scale)
    SI_tools.client.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute FIA and FIP metrics, apply boolean mask, and plot spatial + temporal outputs.")
    parser.add_argument("--sim_name", type=str, required=True)
    parser.add_argument("--ispd_thresh", type=float, required=True)
    parser.add_argument("--ice_type", nargs="+", default=["FI_BT"], help="must be of the format [TYPE_OF_ICE]_[TYPE_OF_GRID] (e.g. 'FI_BT')")
    parser.add_argument("--start_date", default="1994-01-01", help="Start date (YYYY-MM-DD), which is then added to FI_days as the first center-date")
    parser.add_argument("--end_date", default="1999-12-31", help="End date (YYYY-MM-DD), will stop processing when this end_date-FI_days")
    parser.add_argument("--overwrite_zarr", action="store_true")
    parser.add_argument("--overwrite_png", action="store_true")
    parser.add_argument("--smooth_FIA_days", type=int, default=0)
    parser.set_defaults(compute_boolean=False, overwrite_zarr=False, overwrite_png=False)
    args = parser.parse_args()

    # Flatten any space-separated or comma-separated strings into a proper list
    ice_type = args.ice_type
    if isinstance(ice_type, list):
        if len(ice_type) == 1:
            ice_type = ice_type[0].replace(",", " ").split()

    main(args.sim_name,
         args.ispd_thresh,
         ice_type,
         args.start_date,
         args.end_date,
         args.smooth_FIA_days,
         args.overwrite_zarr,
         args.overwrite_png)

