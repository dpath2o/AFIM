import sys, os, argparse, pygmt
sys.path.insert(0, '/home/581/da1339/AFIM/src/AFIM/src')
from sea_ice_toolbox import SeaIceToolbox
from pathlib         import Path
import xarray        as xr
import numpy         as np
import pandas        as pd

def main(sim_name, ispd_thresh, ice_type, dt0_str, dtN_str, log_file, smooth_FIA_days, overwrite_zarr, overwrite_png):
    print(f"ice_type passed is: {ice_type}")
    ispd_str = f"{ispd_thresh:.1e}".replace("e-0", "e-")
    SI_tools = SeaIceToolbox(sim_name       = sim_name,
                             dt0_str        = dt0_str,
                             dtN_str        = dtN_str,
                             P_log          = log_file,
                             ice_speed_threshold = ispd_thresh,
                             overwrite_zarr = overwrite_zarr,
                             overwrite_png  = overwrite_png)
    for itype in ice_type:
        SI_tools.sea_ice_metrics_wrapper(ice_type = itype)
    SI_tools.client.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute FIA and FIP metrics, apply boolean mask, and plot spatial + temporal outputs.")
    parser.add_argument("--sim_name"        , type=str, required=True)
    parser.add_argument("--ispd_thresh"     , type=float, required=True)
    parser.add_argument("--ice_type"         , nargs="+", default=["FI_B", "FI_Ta", "FI_Tx", "FI_BT"])
    parser.add_argument("--start_date"      , help="Start date (YYYY-MM-DD), which is then added to FI_days as the first center-date", default="1993-01-01")
    parser.add_argument("--end_date"        , help="End date (YYYY-MM-DD), will stop processing when this end_date-FI_days", default="1999-12-31")
    parser.add_argument("--log_file"        , help="Path to log file (default: use hard-coded JSON file in SeaIceToolbox)") 
    parser.add_argument("--overwrite_zarr"  , action="store_true")
    parser.add_argument("--overwrite_png"   , action="store_true")
    parser.add_argument("--smooth_FIA_days" , type=int, default=0)
    parser.set_defaults(compute_boolean=False, overwrite_zarr=False, overwrite_png=False)
    args     = parser.parse_args()
    ice_type = args.ice_type
    # Flatten any space-separated or comma-separated strings into a proper list
    if isinstance(ice_type, list):
        if len(ice_type) == 1:
            ice_type = ice_type[0].replace(",", " ").split()
    main(args.sim_name,
         args.ispd_thresh,
         ice_type,
         args.start_date,
         args.end_date,
         args.log_file,
         args.smooth_FIA_days,
         args.overwrite_zarr,
         args.overwrite_png)
