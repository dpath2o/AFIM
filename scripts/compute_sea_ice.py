import argparse
from datetime import datetime
import sys
sys.path.insert(0, '/home/581/da1339/AFIM/src/AFIM/src')
from sea_ice_processor import SeaIceProcessor

def run_loop(sim_name,
             json_path       = None,
             sea_ice         = False,
             pack_ice        = False,
             icon_thresh     = None,
             ispd_thresh     = None,
             roll_win        = None,
             extra_cice_vars = None,
             start_date      = None,
             end_date        = None,
             ow_zarrs        = False,
             hemisphere      = None,
             zarr_directory  = None):
    dt_start = datetime.strptime(start_date, "%Y-%m-%d")
    dt_end   = datetime.strptime(end_date, "%Y-%m-%d")
    if extra_cice_vars and isinstance(extra_cice_vars, str):
        extra_cice_vars = [v.strip() for v in extra_cice_vars.split(',')]
    icon_thresh = float(icon_thresh) if icon_thresh else None
    ispd_thresh = float(ispd_thresh) if ispd_thresh else None
    roll_win    = int(roll_win)      if roll_win   else None
    SI_proc     = SeaIceProcessor(sim_name,
                                  P_json                      = json_path,
                                  sea_ice                     = sea_ice,
                                  pack_ice                    = pack_ice,
                                  ice_concentration_threshold = icon_thresh,
                                  ice_speed_threshold         = ispd_thresh,
                                  extra_cice_vars             = extra_cice_vars,
                                  zarr_directory              = zarr_directory,
                                  hemisphere                  = hemisphere)
    DS          = SI_proc.process_window(dt0_str        = dt_start,
                                         dtN_str        = dt_end,
                                         rolling_window = roll_win,
                                         write_zarr     = True,
                                         ow_zarrs       = ow_zarrs)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run pack ice processor loop over time.")
    parser.add_argument("sim_name", help="Name of simulation")
    parser.add_argument("--json_config", help="Path to JSON config file (default: use hard-coded JSON file in fast_ice_processor)")
    parser.add_argument("--start_date", help="Start date (YYYY-MM-DD) to begin sea ice analysis", default="1993-01-01")
    parser.add_argument("--end_date", help="End date (YYYY-MM-DD), to end sea ice analysis", default="1999-12-31")
    parser.add_argument("--sea_ice", action="store_true", help="compute all sea ice")
    parser.add_argument("--pack_ice", action="store_true", help="compute pack ice, which is all mobile sea ice (i.e. not fast icea)")
    parser.add_argument("--icon_thresh", default=None, help="decimal number in range (0 to 1) for masking sea ice concentration *greater* than this value (default: None; which uses the default value in SeaIceProcessor)")
    parser.add_argument("--ispd_thresh", default=None, help="decimal number in range (1e-5 to 1e-2) for masking sea speeds *less* than this value (default: None; which uses the default value in SeaIceProcessor)")
    parser.add_argument("--roll_win", default=None, help="integer number in range (5 to 30) for computing rolling window average (default: None; which uses the default value in SeaIceProcessor)")
    parser.add_argument("--extra_cice_vars", default=None, help="Comma-separated list of valid CICE variables (caution! advanced user knowledge required) (default: None; which uses the default in SeaIceProcessor)")
    parser.add_argument("--hemisphere", default="south", help="Hemisphere to process (default: south)")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing Zarr files")
    parser.add_argument("--zarr_directory", default=None, help="directory location of where to store zarr files (default: none; which uses the default in SeaIceProcessor)")
    args = parser.parse_args()
    run_loop(sim_name        = args.sim_name,
             json_path       = args.json_config,
             sea_ice         = args.sea_ice,
             pack_ice        = args.pack_ice,
             icon_thresh     = args.icon_thresh,
             ispd_thresh     = args.ispd_thresh,
             roll_win        = args.roll_win,
             extra_cice_vars = args.extra_cice_vars,
             start_date      = args.start_date,
             end_date        = args.end_date,
             ow_zarrs        = args.overwrite,
             hemisphere      = args.hemisphere,
             zarr_directory  = args.zarr_directory)
