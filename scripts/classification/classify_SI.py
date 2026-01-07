import argparse
import xarray as xr
from pathlib  import Path
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

    # SI does not require a speed threshold, but we allow an optional value so the
    # output can still be written under a consistent ispd_thresh_* directory.
    if ispd_thresh is not None:
        try:
            ispd_thresh = float(ispd_thresh)
        except Exception as e:
            raise ValueError(f"--ispd_thresh must be float-like, got {ispd_thresh!r}") from e

    SI_tools = SI_tool_mgr.get_toolbox(P_json              = json_path,
                                       dt0_str             = start_date,
                                       dtN_str             = end_date,
                                       sim_name            = sim_name,
                                       list_of_B2T         = B2T_type,
                                       ice_speed_threshold = ispd_thresh)

    SI_tools.daily_iceh_to_monthly_zarr(netcdf_engine   = netcdf_engine,
                                        overwrite       = overwrite_zarr,
                                        delete_original = delete_original_iceh)

    SI_tools.define_datetime_vars()

    for yr0_str, yrN_str in zip(SI_tools.yr0_strs, SI_tools.yrN_strs):
        yr_str = f"{yr0_str[:4]}"
        SI_tools.logger.info(f"looping over each month in year {yr_str}")

        SI_mo = []
        for mo0_str, moN_str in zip(SI_tools.mo0_strs, SI_tools.moN_strs):
            SI_raw, _, _, _ = SI_tools.classify_sea_ice(dt0_str=mo0_str,
                                                        dtN_str=moN_str)
            SI_mo.append(SI_raw)

        SI_tools.logger.info("concatenating monthly datasets into yearly")
        chunks = SI_tools.CICE_dict.get("SI_chunks", SI_tools.CICE_dict.get("FI_chunks", None))
        SI_yr  = xr.concat(SI_mo, dim="time")
        if chunks is not None:
            SI_yr = SI_yr.chunk(chunks)

        F_prefix = f"SI"
        P_zarr   = Path(SI_tools.D_zarr, SI_tools.hemisphere_dict['abbreviation'], f"{F_prefix}.zarr")

        SI_tools.logger.info(f"writing zarr to disk {P_zarr}")
        SI_yr.to_zarr(P_zarr, group=yr_str, mode="w", consolidated=False)

    SI_tools.client.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run sea-ice (SI) classification loop over time.")
    parser.add_argument("sim_name", help="Name of simulation")
    parser.add_argument("--ispd_thresh", help="Optional: speed threshold only used for output directory naming (SI itself does not threshold speed)")
    parser.add_argument("--B2T_type", nargs="+", choices=["B", "Ta", "Tb", "Tx", "BT"], help="Optional: list of ice speed vector types for directory naming/metadata (SI itself does not threshold speed)")
    parser.add_argument("--overwrite_zarr", action="store_true", help="overwrite zarrs")
    parser.add_argument("--delete_original_iceh", action="store_true", help="after creating monthly iceh_*.zarr groups, delete original daily ice history files")
    parser.add_argument("--json_config", help="Path to JSON config file (default: use hard-coded JSON file in toolbox)")
    parser.add_argument("--log_file", help="Path to log file")
    parser.add_argument("--start_date", help="Start date (YYYY-MM-DD)", default="1993-01-01")
    parser.add_argument("--end_date", help="End date (YYYY-MM-DD)", default="1999-12-31")
    parser.add_argument("--netcdf_engine", help="Engine used for xarray mfdataset method", default="netcdf4")

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
