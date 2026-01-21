import argparse
import xarray as xr
from pathlib  import Path
from datetime import datetime, timedelta
import sys
sys.path.insert(0, '/home/581/da1339/AFIM/src/AFIM/src')
from sea_ice_toolbox import SeaIceToolboxManager

def run_loop(sim_name,
             ispd_thresh          = None,
             BorC2T_type          = None,
             ice_type             = None,
             json_path            = None,
             start_date           = None,
             end_date             = None,
             log_file             = None,
             netcdf_engine        = None,
             overwrite_zarr       = False,
             delete_original_iceh = False):
    mgr = SeaIceToolboxManager(P_log=log_file)
    tb  = mgr.get_toolbox(P_json              = json_path,
                          dt0_str             = start_date,
                          dtN_str             = end_date,
                          sim_name            = sim_name,
                          list_of_BorC2T      = BorC2T_type,
                          ice_speed_threshold = float(ispd_thresh),
                          ice_type            = ice_type)
    tb.daily_iceh_to_monthly_zarr(netcdf_engine   = netcdf_engine,
                                  overwrite       = overwrite_zarr,
                                  delete_original = delete_original_iceh)
    tb.define_datetime_vars()
    tb.define_ice_class_name()
    tb.define_classification_dir()
    for yr0_str,yrN_str in zip(tb.yr0_strs,tb.yrN_strs):
        yr_str = f"{yr0_str[:4]}"
        tb.logger.info(f"looping over each month in year {yr_str}")
        I_mo, I_roll_mo, I_bin_mo = [], [], []
        I_yr, I_roll_yr, I_bin_yr = None, None, None
        for mo0_str, moN_str in zip(tb.mo0_strs, tb.moN_strs):
            I_raw,I_bin,I_roll,_ = tb.classify_ice(dt0_str               = mo0_str,
                                                dtN_str               = moN_str, 
                                                enable_rolling_output = True)
            tb.logger.info(f"concatenating monthly class methods:")
            tb.logger.info("    raw")
            I_mo.append(I_raw)
            if I_bin:
                tb.logger.info("    binary-days")
                I_bin_mo.append(I_bin)
            if I_roll:
                tb.logger.info("    rolling-mean")
                I_roll_mo.append(I_roll)
        tb.logger.info("concatenating monthly datasets into yearly")
        tb.logger.info("    raw")
        I_yr = xr.concat(I_mo, dim="time").chunk(tb.CICE_dict["FI_chunks"])
        if I_bin_mo:
            tb.logger.info("    binary-days")
            I_bin_yr = xr.concat(I_bin_mo, dim="time").chunk(tb.CICE_dict["FI_chunks"])
        if I_roll_mo:
            tb.logger.info("    rolling-mean")
            I_roll_yr = xr.concat(I_roll_mo, dim="time").chunk(tb.CICE_dict["FI_chunks"])
        P_zarr      = tb.define_classification_zarr(class_method = "raw")
        P_zarr_bin  = tb.define_classification_zarr(class_method = "binary-days")
        P_zarr_roll = tb.define_classification_zarr(class_method = "rolling-mean")
        toz         = dict(group=yr_str, mode="w", consolidated=False, zarr_format=2)
        tb.logger.info(f"writing {P_zarr}")
        I_yr.to_zarr(P_zarr, **toz)
        if I_bin_yr:
            tb.logger.info(f"writing {P_zarr_bin}")
            I_bin_yr.to_zarr(P_zarr_bin, **toz)
        if I_roll_yr:
            tb.logger.info(f"writing {P_zarr_roll}")
            I_roll_yr.to_zarr(P_zarr_roll, **toz)
    tb.client.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run pack ice processor loop over time.")
    parser.add_argument("sim_name"              , help="Name of simulation")
    parser.add_argument("--ispd_thresh"         , help="ice speed threshold for masking fast ice (default: use hard-coded value in JSON file)")
    parser.add_argument("--BorC2T_type"          , nargs="+", choices=["B", "Ta", "Tb", "Tc", "Tx"], help="List of ice speed masking types to use (choose one or more of: B, Ta, Tb, Tx). Example: --BorC2T_type B Ta")
    parser.add_argument("--ice_type"            , help="type of ice: 'FI', 'PI', 'SI', 'MI'")
    parser.add_argument("--overwrite_zarr"      , action="store_true", help="overwrite zarrs")
    parser.add_argument("--delete_original_iceh", action="store_true", help="when creating monthly zarr files of original CICE files, if enabled this will then delete the original model ice history files")
    parser.add_argument("--json_config"         , help="Path to JSON config file (default: use hard-coded JSON file in fast_ice_processor)")
    parser.add_argument("--log_file"            , help="Path to log file (default: use hard-coded JSON file in fast_ice_processor)")
    parser.add_argument("--start_date"          , help="Start date (YYYY-MM-DD), which is then added to FI_days as the first center-date", default="1993-01-01")
    parser.add_argument("--end_date"            , help="End date (YYYY-MM-DD), will stop processing when this end_date-FI_days", default="1999-12-31")
    parser.add_argument("--netcdf_engine"       , help="Engine used for xarray mfdataset method", default="netcdf4")
    args = parser.parse_args()
    run_loop(sim_name             = args.sim_name,
             ispd_thresh          = args.ispd_thresh,
             BorC2T_type          = args.BorC2T_type,
             ice_type             = args.ice_type,
             json_path            = args.json_config,
             start_date           = args.start_date,
             end_date             = args.end_date,
             log_file             = args.log_file,
             netcdf_engine        = args.netcdf_engine,
             overwrite_zarr       = args.overwrite_zarr,
             delete_original_iceh = args.delete_original_iceh)