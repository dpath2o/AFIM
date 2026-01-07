import argparse
import xarray as xr
from pathlib  import Path
from datetime import datetime, timedelta
import sys
sys.path.insert(0, '/home/581/da1339/AFIM/src/AFIM/src')
from sea_ice_toolbox import SeaIceToolbox, SeaIceToolboxManager

def run_loop(sim_name,
             ispd_thresh          = None,
             B2T_type             = None,
             json_path            = None,
             start_date           = None,
             end_date             = None,
             log_file             = None,
             netcdf_engine        = None,
             overwrite_zarr       = False,
             delete_original_iceh = False):
    SI_tool_mgr = SeaIceToolboxManager(P_log=log_file)
    ispd_thresh = float(ispd_thresh)
    SI_tools    = SI_tool_mgr.get_toolbox(P_json              = json_path,
                                          dt0_str             = start_date,
                                          dtN_str             = end_date,
                                          sim_name            = sim_name,
                                          list_of_B2T         = B2T_type,
                                          ice_speed_threshold = ispd_thresh)
    SI_tools.daily_iceh_to_monthly_zarr(netcdf_engine   = netcdf_engine,
                                        overwrite       = overwrite_zarr,
                                        delete_original = delete_original_iceh)
    SI_tools.define_datetime_vars()
    for yr0_str,yrN_str in zip(SI_tools.yr0_strs,SI_tools.yrN_strs):
        yr_str = f"{yr0_str[:4]}"
        SI_tools.logger.info(f"looping over each month in year {yr_str}")
        PI_mo      = []
        PI_roll_mo = []
        PI_bin_mo  = []
        for mo0_str, moN_str in zip(SI_tools.mo0_strs, SI_tools.moN_strs):
            PI_raw,PI_bin,PI_roll,_ = SI_tools.classify_pack_ice(dt0_str               = mo0_str,
                                                                 dtN_str               = moN_str, 
                                                                 enable_rolling_output = True)
            PI_mo.append(PI_raw)
            PI_bin_mo.append(PI_bin)
            PI_roll_mo.append(PI_roll)
        SI_tools.logger.info("concatenating monthly datasets into yearly")
        PI_yr       = xr.concat(PI_mo, dim="time").chunk(SI_tools.CICE_dict["PI_chunks"])
        PI_bin_yr   = xr.concat(PI_bin_mo, dim="time").chunk(SI_tools.CICE_dict["PI_chunks"])
        PI_roll_yr  = xr.concat(PI_roll_mo, dim="time").chunk(SI_tools.CICE_dict["PI_chunks"])
        F_prefix    = f"PI_{''.join(SI_tools.B2T_type)}"
        P_zarr      = Path(SI_tools.D_ispd_thresh, f"{F_prefix}.zarr")
        P_zarr_bin  = Path(SI_tools.D_ispd_thresh, f"{F_prefix}_bin.zarr")
        P_zarr_roll = Path(SI_tools.D_ispd_thresh, f"{F_prefix}_roll.zarr")
        SI_tools.logger.info(f"writing zarr to disk {P_zarr}")
        PI_yr.to_zarr(P_zarr, group=yr_str, mode="w", consolidated=False)
        SI_tools.logger.info(f"writing zarr to disk {P_zarr_bin}")
        PI_bin_yr.to_zarr(P_zarr_bin, group=yr_str, mode="w", consolidated=False)
        SI_tools.logger.info(f"writing zarr to disk {P_zarr_roll}")
        PI_roll_yr.to_zarr(P_zarr_roll, group=yr_str, mode="w", consolidated=False)
    SI_tools.client.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run pack-ice (PI) classification loop over time.")
    parser.add_argument("sim_name"              , help="Name of simulation")
    parser.add_argument("--ispd_thresh"         , help="Ice speed threshold (m/s) for PI classification (default: use value in JSON config)")
    parser.add_argument("--B2T_type"           , nargs="+", choices=["B", "Ta", "Tb", "Tx", "BT"], help="List of ice speed masking types to use (choose one or more of: B, Ta, Tb, Tx, BT). Example: --B2T_type B Ta")
    parser.add_argument("--overwrite_zarr"      , action="store_true", help="overwrite zarrs")
    parser.add_argument("--delete_original_iceh", action="store_true", help="when creating monthly zarr files of original CICE files, if enabled this will then delete the original model ice history files")
    parser.add_argument("--json_config"         , help="Path to JSON config file (default: use hard-coded JSON file in fast_ice_processor)")
    parser.add_argument("--log_file"            , help="Path to log file (default: use hard-coded JSON file in fast_ice_processor)")
    parser.add_argument("--start_date"          , help="Start date (YYYY-MM-DD), which is then added to PI_days as the first center-date", default="1993-01-01")
    parser.add_argument("--end_date"            , help="End date (YYYY-MM-DD), will stop processing when this end_date-PI_days", default="1999-12-31")
    parser.add_argument("--netcdf_engine"       , help="Engine used for xarray mfdataset method", default="netcdf4")
    args = parser.parse_args()
    run_loop(sim_name             = args.sim_name,
             ispd_thresh          = args.ispd_thresh,
             B2T_type             = args.B2T_type,
             json_path            = args.json_config,
             start_date           = args.start_date,
             end_date             = args.end_date,
             log_file             = args.log_file,
             netcdf_engine        = args.netcdf_engine,
             overwrite_zarr       = args.overwrite_zarr,
             delete_original_iceh = args.delete_original_iceh)
