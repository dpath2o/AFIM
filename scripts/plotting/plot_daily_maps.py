import sys
import shutil
import pygmt
import argparse
import time
from pathlib import Path
import pandas as pd
import xarray as xr
mod_path = '/home/581/da1339/AFIM/src/AFIM/src'
if mod_path not in sys.path:
    sys.path.insert(0, mod_path)
from sea_ice_toolbox import SeaIceToolbox

def delete_gmt_session(D_gmt_session, age_threshold=120):
    """
    Delete whole directories that are older than the specified age (in seconds) from the given GMT session directory.
    
    Parameters:
    -----------
    D_gmt_session : str or Path
        The path to the directory containing GMT session files to be cleaned.
    age_threshold : int
        The age threshold in seconds. Directories older than this age will be deleted (default is 300 seconds, or 5 minutes).
    """
    current_time = time.time()
    directory = Path(D_gmt_session)
    # Check if the directory exists and is valid
    if not directory.exists() or not directory.is_dir():
        print(f"Directory {directory} does not exist or is not a valid directory.")
        return
    # Loop through directories in the session directory
    for file_path in directory.iterdir():
        if file_path.is_dir():  # Check if it's a directory
            file_mtime = file_path.stat().st_mtime
            # Check if the directory is older than the threshold
            if current_time - file_mtime > age_threshold:
                try:
                    shutil.rmtree(file_path)  # Delete the directory and its contents
                    print(f"Deleted: {file_path}")
                except Exception as e:
                    print(f"Failed to delete {file_path}: {e}")

def plot_monthly_variable_maps(sim_name, var_names, ispd_thresh=5.0e-4, dt0_str="1994-01-01", dtN_str="1999-12-31", script_log=None):
    if isinstance(var_names, str):
        var_names = [var_names]
    SI_tools      = SeaIceToolbox(sim_name             = sim_name,
                                  ice_speed_threshold  = ispd_thresh,
                                  dt0_str              = dt0_str,
                                  dtN_str              = dtN_str,
                                  P_log                = script_log,
                                  save_new_figs        = True,
                                  show_figs            = False,
                                  overwrite_saved_figs = True)
    FI_raw, CICE    = SI_tools.load_processed_cice(zarr_CICE=True)
    D_gmt_sess_base = Path(f"{Path.home()}/.gmt/sessions")
    
    # Loop through each time step
    for i in range(len(CICE['time'].values)):
        CICE_slc = CICE.isel(time=i, nj=SI_tools.hemisphere_dict['nj_slice'])
        dt       = pd.Timestamp(CICE.isel(time=i)['time'].values)
        dt_str   = f"{dt:%Y-%m-%d}"
        
        # Loop through each variable to plot
        for var_name in var_names:
            if var_name not in CICE_slc:
                SI_tools.logger.info(f"Skipping {var_name} as it is not found in {CICE_slc.keys()}")
                continue
            var_all = CICE_slc[var_name]
            SI_tools.logger.info(f"Plotting {var_name} for date {dt_str}")
            # Plot the variable
            SI_tools.pygmt_map_plot_one_var(var_all, var_name,
                                            sim_name       = sim_name,
                                            time_stamp     = dt_str,
                                            tit_str        = f"{sim_name} {dt_str}",
                                            plot_GI        = False,
                                            plot_iceshelves= False,
                                            plot_bathymetry= False,
                                            add_stat_annot = True,
                                            overwrite_fig  = False,
                                            show_fig       = False)
            
        # After plotting each var_name, delete old GMT session directories
        #delete_gmt_session(D_gmt_sess_base)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot monthly sea ice variable maps using SeaIcePlotter.")
    parser.add_argument("--sim_name"     , required=True, help="Simulation name (e.g. 'elps-min')")
    parser.add_argument("--var_names"    , required=True, nargs='+', help="One or more variable names to plot")
    parser.add_argument("--ispd_thresh"  , type=float, default=5.0e-4, help="Threshold for ispd plots")
    parser.add_argument("--start_date"   , help="Start date (YYYY-MM-DD), which is then added to FI_days as the first center-date", default="1994-01-01")
    parser.add_argument("--end_date"     , help="End date (YYYY-MM-DD), will stop processing when this end_date-FI_days", default="1999-12-31")
    parser.add_argument("--log_file"     , help="Path to log file to log the messages from SeaIceToolbox (default is ~/logs/plot_daily_maps/[SIM_NAME].log)") 
    args       = parser.parse_args()
    script_log = args.log_file if args.log_file is not None else Path(f"{Path.home()}/logs/plot_daily_maps/{args.sim_name}")
    plot_monthly_variable_maps(sim_name    = args.sim_name,
                               var_names   = args.var_names,
                               ispd_thresh = args.ispd_thresh,
                               dt0_str     = args.start_date,
                               dtN_str     = args.end_date,
                               script_log  = script_log)

