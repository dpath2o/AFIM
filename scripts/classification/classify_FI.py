import argparse
import xarray as xr
from pathlib  import Path
from datetime import datetime, timedelta
import sys
sys.path.insert(0, '/home/581/da1339/AFIM/src/AFIM/src')
from sea_ice_toolbox import SeaIceToolbox, SeaIceToolboxManager

def run_loop(sim_name,
             ispd_thresh          = None,
             ivec_type            = None,
             json_path            = None,
             start_date           = None,
             end_date             = None,
             log_file             = None,
             daily                = False,
             rolling              = False,
             roll_period          = None,
             overwrite_zarr       = False,
             delete_original_iceh = False):
    SI_tool_mgr = SeaIceToolboxManager(P_log=log_file)
    ispd_thresh = float(ispd_thresh)
    SI_tools    = SI_tool_mgr.get_toolbox(P_json              = json_path,
                                          dt0_str             = start_date,
                                          dtN_str             = end_date,
                                          sim_name            = sim_name,
                                          ice_speed_threshold = ispd_thresh,
                                          ice_vector_type     = ivec_type)
    SI_tools.daily_iceh_to_monthly_zarr(overwrite=overwrite_zarr, delete_original=delete_original_iceh)
    SI_tools.define_datetime_vars()
    for yr0_str,yrN_str in zip(SI_tools.yr0_strs,SI_tools.yrN_strs):
        yr_str = f"{yr0_str[:4]}"
        SI_tools.logger.info(f"looping over each month in year {yr_str}")
        FI_mo      = []
        FI_roll_mo = []
        FI_bin_mo  = []
        for mo0_str, moN_str in zip(SI_tools.mo0_strs, SI_tools.moN_strs):
            FI_raw,FI_bin,FI_roll = SI_tools.classify_fast_ice(dt0_str=mo0_str, dtN_str=moN_str, enable_rolling_output=True)
            FI_mo.append(FI_raw)
            FI_bin_mo.append(FI_bin)
            FI_roll_mo.append(FI_roll)
        SI_tools.logger.info("concatenating monthly datasets into yearly")
        FI_yr       = xr.concat(FI_mo, dim="time").chunk(SI_tools.CICE_dict["FI_chunks"])
        FI_bin_yr   = xr.concat(FI_bin_mo, dim="time").chunk(SI_tools.CICE_dict["FI_chunks"])
        FI_roll_yr  = xr.concat(FI_roll_mo, dim="time").chunk(SI_tools.CICE_dict["FI_chunks"])
        P_zarr      = Path(SI_tools.D_ispd_thresh, "FI_BT.zarr")
        P_zarr_bin  = Path(SI_tools.D_ispd_thresh, "FI_BT_bin.zarr")
        P_zarr_roll = Path(SI_tools.D_ispd_thresh, "FI_BT_roll.zarr")
        SI_tools.logger.info(f"writing zarr to disk {P_zarr}")
        FI_yr.to_zarr(P_zarr, group=yr_str, mode="w", consolidated=False)
        SI_tools.logger.info(f"writing zarr to disk {P_zarr_bin}")
        FI_bin_yr.to_zarr(P_zarr_bin, group=yr_str, mode="w", consolidated=False)
        SI_tools.logger.info(f"writing zarr to disk {P_zarr_roll}")
        FI_roll_yr.to_zarr(P_zarr_roll, group=yr_str, mode="w", consolidated=False)
    SI_tools.client.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run pack ice processor loop over time.")
    parser.add_argument("sim_name"              , help="Name of simulation")
    parser.add_argument("--ispd_thresh"         , help="ice speed threshold for masking fast ice (default: use hard-coded value in JSON file)")
    parser.add_argument("--ivec_type"           , nargs="+", choices=["B", "Ta", "Tx", "BT"],
                        help="List of ice speed masking types to use (choose one or more of: B, Ta, Tx). Example: --ivec_type B Ta")
    parser.add_argument("--rolling"             , action="store_true", 
                        help="Enable rolling temporal average (length [MEAN_PERIOD]) of sea ice concentration and sea ice speed prior to fast ice masking")
    parser.add_argument("--daily"               , action="store_true", help="Perform fast ice masking without ")
    parser.add_argument("--overwrite_zarr"      , action="store_true", help="overwrite zarrs")
    parser.add_argument("--delete_original_iceh", action="store_true", help="when creating monthly zarr files of original CICE files, if enabled this will then delete the original model ice history files")
    parser.add_argument("--mean_period"         , type=int , default=14 , help="averaging period in days over which to perform rolling calculation")
    parser.add_argument("--json_config"         , help="Path to JSON config file (default: use hard-coded JSON file in fast_ice_processor)")
    parser.add_argument("--log_file"            , help="Path to log file (default: use hard-coded JSON file in fast_ice_processor)")
    parser.add_argument("--start_date"          , help="Start date (YYYY-MM-DD), which is then added to FI_days as the first center-date", default="1993-01-01")
    parser.add_argument("--end_date"            , help="End date (YYYY-MM-DD), will stop processing when this end_date-FI_days", default="1999-12-31")
    args = parser.parse_args()
    run_loop(sim_name             = args.sim_name,
             ispd_thresh          = args.ispd_thresh,
             ivec_type            = args.ivec_type,
             json_path            = args.json_config,
             start_date           = args.start_date,
             end_date             = args.end_date,
             log_file             = args.log_file,
             daily                = args.daily,
             rolling              = args.rolling,
             roll_period          = args.mean_period,
             overwrite_zarr       = args.overwrite_zarr,
             delete_original_iceh = args.delete_original_iceh)
