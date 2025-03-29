import argparse
from datetime import datetime, timedelta
from pathlib import Path
import sys
sys.path.insert(0, '/home/581/da1339/AFIM/src/AFIM/src')
from pack_ice_processor import PackIceProcessor  # Make sure this matches your actual module filename
import time

def run_loop(sim_name, json_path=None, start_date=None, end_date=None, ow_zarrs=False, hemisphere=None, PI_days=None, dry_run=False):
    dt_start = datetime.strptime(start_date, "%Y-%m-%d")
    dt_end   = datetime.strptime(end_date, "%Y-%m-%d")
    processor = PackIceProcessor(sim_name, json_path=json_path, roll_win=PI_days)
    current_date = dt_start+timedelta(days=FI_days)
    while current_date <= (dt_end-timedelta(days=FI_days)):
        print(f"\n → [{current_date.strftime('%Y-%m-%d')}] ", end="")
        if dry_run:
            print(f"(dry run) would process window with overwrite={ow_zarrs}")
        else:
            t0 = time.time()
            PI = processor.process_window(current_date, hemisphere=hemisphere, ow_zarrs=ow_zarrs)
            print(f"✓ Completed in {time.time() - t0:.1f} seconds\n")
        current_date += timedelta(days=FI_days)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run pack ice processor loop over time.")
    parser.add_argument("sim_name", help="Name of simulation")
    parser.add_argument("--json_config", help="Path to JSON config file (default: use hard-coded JSON file in fast_ice_processor)")
    parser.add_argument("--start_date", help="Start date (YYYY-MM-DD), which is then added to FI_days as the first center-date", default="1993-01-01")
    parser.add_argument("--end_date", help="End date (YYYY-MM-DD), will stop processing when this end_date-FI_days", default="1999-12-31")
    parser.add_argument("--hemisphere", default="south", help="Hemisphere to process (default: south)")
    parser.add_argument("--PI_days", default=15, help="number of days over which to compute pack ice (default: 15)", type=int)
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing Zarr files")
    parser.add_argument("--dry_run", action="store_true", help="Preview steps without running them")

    args = parser.parse_args()
    run_loop(sim_name=args.sim_name,
             json_path=args.json_config,
             start_date=args.start_date,
             end_date=args.end_date,
             ow_zarrs=args.overwrite,
             hemisphere=args.hemisphere,
             PI_days=args.PI_days,
             dry_run=args.dry_run)
