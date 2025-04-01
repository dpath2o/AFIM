import argparse
from datetime import datetime, timedelta
from pathlib import Path
import sys
sys.path.insert(0, '/home/581/da1339/AFIM/src/AFIM/src')
from pack_ice_processor import PackIceProcessor  # Make sure this matches your actual module filename
import time

def run_loop(sim_name, json_path=None, start_date=None, end_date=None, ow_zarrs=False, hemisphere=None):
    print(f"🚀 Running {sim_name} from {start_date} to {end_date} with overwrite={ow_zarrs} and hemisphere={hemisphere}")
    dt_start = datetime.strptime(start_date, "%Y-%m-%d")
    dt_end   = datetime.strptime(end_date, "%Y-%m-%d")
    processor = PackIceProcessor(sim_name, dt0_str=dt_start, dtN_str=dt_end, mask_fast_ice=True, hemisphere=hemisphere)
    PI        = processor.process_window(save_zarr=True, ow_zarrs=ow_zarrs)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run pack ice processor loop over time.")
    parser.add_argument("sim_name", help="Name of simulation")
    parser.add_argument("--json_config", help="Path to JSON config file (default: use hard-coded JSON file in fast_ice_processor)")
    parser.add_argument("--start_date", help="Start date (YYYY-MM-DD), which is then added to FI_days as the first center-date", default="1993-01-01")
    parser.add_argument("--end_date", help="End date (YYYY-MM-DD), will stop processing when this end_date-FI_days", default="1999-12-31")
    parser.add_argument("--hemisphere", default="south", help="Hemisphere to process (default: south)")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing Zarr files")

    args = parser.parse_args()
    run_loop(sim_name=args.sim_name,
             json_path=args.json_config,
             start_date=args.start_date,
             end_date=args.end_date,
             ow_zarrs=args.overwrite,
             hemisphere=args.hemisphere)
