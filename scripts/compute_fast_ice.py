import argparse
from datetime import datetime, timedelta
import sys
sys.path.insert(0, '/home/581/da1339/AFIM/src/AFIM/src')
from sea_ice_processor import SeaIceProcessor

def run_loop(sim_name, json_path=None, sea_ice=False, pack_ice=False, start_date=None, end_date=None, ow_zarrs=False, hemisphere=None, zarr_directory=None):
    dt_start = datetime.strptime(start_date, "%Y-%m-%d")
    dt_end   = datetime.strptime(end_date, "%Y-%m-%d")
    SI_proc  = SeaIceProcessor(P_json                  = json_path,
                               pack_ice                = sea_ice,
                               pack_ice_minus_fast_ice = pack_ice
                               sim_name                = sim_name,
                               hemisphere              = hemisphere)
    SI       = SI_proc.process_window(dt0_str        = dt0_str,
                                      dtN_str        = dtN_str,
                                      write_zarr     = True,
                                      ow_zarrs       = ow_zarrs)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run pack ice processor loop over time.")
    parser.add_argument("sim_name", help="Name of simulation")
    parser.add_argument("--json_config", help="Path to JSON config file (default: use hard-coded JSON file in fast_ice_processor)")
    parser.add_argument("--sea_ice", action="store_true", help="compute all sea ice")
    parser.add_argument("--pack_ice", action="store_true", help="compute pack ice, which is all mobile sea ice (i.e. not fast ice)")
    parser.add_argument("--start_date", help="Start date (YYYY-MM-DD), which is then added to FI_days as the first center-date", default="1993-01-01")
    parser.add_argument("--end_date", help="End date (YYYY-MM-DD), will stop processing when this end_date-FI_days", default="1999-12-31")
    parser.add_argument("--hemisphere", default="south", help="Hemisphere to process (default: south)")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing Zarr files")
    parser.add_argument("--zarr_directory", default="FI", help="directory location (under the simulation directory) of where to store zarr files")

    args = parser.parse_args()
    run_loop(sim_name=args.sim_name,
             json_path=args.json_config,
             sea_ice=args.sea_ice,
             pack_ice=args.pack_ice,
             start_date=args.start_date,
             end_date=args.end_date,
             ow_zarrs=args.overwrite,
             hemisphere=args.hemisphere,
             zarr_directory=args.zarr_directory)
