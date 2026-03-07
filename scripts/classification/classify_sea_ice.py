#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path
import xarray as xr

# AFIM import path (as per your current scripts)
sys.path.insert(0, "/home/581/da1339/AFIM/src/AFIM/src")
from sea_ice_toolbox import SeaIceToolboxManager  # noqa: E402

DEFAULT_P_JSON = "/home/581/da1339/AFIM/src/AFIM/src/JSONs/sea_ice_config.json"

def _prep_for_zarr(ds: xr.Dataset) -> xr.Dataset:
    # 1) drop all non-dimension coordinates (this removes NLAT/NLON/etc if they are coords)
    ds = ds.reset_coords(drop=True)
    # 2) if any of these exist as data_vars, drop them too (belt + braces)
    drop_vars = ["NLAT","NLON","ELAT","ELON","ULAT","ULON","TLAT","TLON",
                 "ANGLE","ANGLET","tarea","uarea","HTE","HTN",
                 "dxt","dyt","dxu","dyu","dxe","dyn","dxn","dye"]
    ds = ds.drop_vars([v for v in drop_vars if v in ds], errors="ignore")
    # 3) clear conflicting encoding (esp. 'chunks') that can trigger the overlap error
    for v in ds.variables:
        ds[v].encoding.pop("chunks", None)
    return ds

def run_loop(sim_name,
             ispd_thresh          = None,
             BorC2T_type          = None,
             ice_type             = None,
             P_JSON               = None,
             start_date           = None,
             end_date             = None,
             log_file             = None,
             netcdf_engine        = None,
             overwrite_zarr       = False,
             delete_original_iceh = False,
             enable_rolling       = False):
    '''
    loop over years defined by start and stop date
    '''
    mgr = SeaIceToolboxManager(P_log=log_file)
    tb  = mgr.get_toolbox(P_json                       = P_JSON,
                          dt0_str                      = start_date,
                          dtN_str                      = end_date,
                          sim_name                     = sim_name,
                          list_of_BorC2T               = BorC2T_type,
                          ice_speed_threshold          = (float(ispd_thresh) if ispd_thresh is not None else None),
                          ice_type                     = ice_type,
                          overwrite                    = overwrite_zarr,
                          delete_original_cice_iceh_nc = delete_original_iceh)
    tb.daily_iceh_to_monthly_zarr(netcdf_engine=netcdf_engine)
    tb.define_datetime_vars()
    tb.define_ice_class_name()
    tb.define_classification_dir()
    for yr0_str, yrN_str in zip(tb.yr0_strs, tb.yrN_strs):
        yr_str = f"{yr0_str[:4]}"
        tb.logger.info(f"Looping over each month in year {yr_str}")
        I_mo, I_roll_mo, I_bin_mo = [], [], []
        I_yr, I_roll_yr, I_bin_yr = None, None, None
        for mo0_str, moN_str in zip(tb.mo0_strs, tb.moN_strs):
            # primary call for classifying sea ice
            I_raw, I_bin, I_roll, _ = tb.classify_ice(dt0_str               = mo0_str,
                                                      dtN_str               = moN_str,
                                                      enable_rolling_output = enable_rolling)
            tb.logger.info("Collecting monthly classifications")
            I_mo.append(I_raw)
            if I_bin is not None:
                I_bin_mo.append(I_bin)
            if I_roll is not None:
                I_roll_mo.append(I_roll)
        tb.logger.info("Concatenating monthly datasets into yearly")
        I_yr = xr.concat(I_mo, dim="time").chunk(tb.CICE_dict["FI_chunks"])
        if len(I_bin_mo) > 0:
            I_bin_yr = xr.concat(I_bin_mo, dim="time").chunk(tb.CICE_dict["FI_chunks"])
        if len(I_roll_mo) > 0:
            I_roll_yr = xr.concat(I_roll_mo, dim="time").chunk(tb.CICE_dict["FI_chunks"])
        P_zarr      = tb.define_classification_zarr(class_method="raw")
        P_zarr_bin  = tb.define_classification_zarr(class_method="binary-days")
        P_zarr_roll = tb.define_classification_zarr(class_method="rolling-mean")
        toz         = dict(group        = yr_str,
                           mode         = "w",
                           consolidated = True,
                           zarr_format  = 2,
                           safe_chunks  = True,
                           align_chunks = True)
        tb.logger.info(f"writing {P_zarr}")
        _prep_for_zarr(I_yr).to_zarr(P_zarr, **toz)
        if I_bin_yr:
            tb.logger.info(f"writing {P_zarr_bin}")
            _prep_for_zarr(I_bin_yr).to_zarr(P_zarr_bin, **toz)
        if I_roll_yr:
            tb.logger.info(f"writing {P_zarr_roll}")
            _prep_for_zarr(I_roll_yr).to_zarr(P_zarr_roll, **toz)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run ice classification loop over time.")
    parser.add_argument("sim_name", help="Name of simulation")
    parser.add_argument("--start_date", default="1993-01-01", help="Start date YYYY-MM-DD")
    parser.add_argument("--end_date", default="1999-12-31", help="End date YYYY-MM-DD")
    parser.add_argument("--P_JSON", "--P_json",
                        dest="P_JSON",
                        default=DEFAULT_P_JSON,
                        help=f"Path to AFIM config JSON (default: {DEFAULT_P_JSON})")
    parser.add_argument("--ispd_thresh",
                        default=None,
                        help="Ice speed threshold for masking fast ice (default: use value from JSON/toolbox defaults)")
    parser.add_argument("--BorC2T_type",
                        nargs="+",
                        choices=["B", "Ta", "Tb", "Tc", "Tx"],
                        help="One or more masking/regrid types: B, Ta, Tb, Tc, Tx (e.g. --BorC2T_type Tc)")
    parser.add_argument("--ice_type", default="FI", help="Ice type: FI, PI, SI, MI")
    parser.add_argument("--overwrite_zarr", action="store_true", help="Overwrite output Zarrs")
    parser.add_argument("--rolling_mean", action="store_true", help="turns on rolling-mean output in sea ice classification (defualt is 'binary-days' and 'daily' only)")
    parser.add_argument("--delete_original_iceh",
                        action="store_true",
                        help="After creating monthly iceh_*.zarr, delete original daily ice history files")
    parser.add_argument("--log_file", default=None, help="Path to log file")
    parser.add_argument("--netcdf_engine", default="netcdf4", help="Engine for xarray.open_mfdataset")
    args = parser.parse_args()
    # precedence: explicit --json_config overrides --P_JSON
    run_loop(sim_name             = args.sim_name,
             ispd_thresh          = args.ispd_thresh,
             BorC2T_type          = args.BorC2T_type,
             ice_type             = args.ice_type,
             P_JSON               = args.P_JSON,
             start_date           = args.start_date,
             end_date             = args.end_date,
             log_file             = args.log_file,
             netcdf_engine        = args.netcdf_engine,
             overwrite_zarr       = args.overwrite_zarr,
             delete_original_iceh = args.delete_original_iceh,
             enable_rolling       = args.rolling_mean)
