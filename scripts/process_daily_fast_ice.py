import argparse
from datetime import datetime, timedelta
import sys
sys.path.insert(0, '/home/581/da1339/AFIM/src/AFIM/src')
from sea_ice_processor import SeaIceProcessor

def run_loop(sim_name, ispd_thresh=None, ispd_type=None, json_path=None, start_date=None, end_date=None, log_file=None, rolling=True, roll_period=None):
    if ispd_type == "None":
        ispd_type = None
    dt_start = datetime.strptime(start_date, "%Y-%m-%d")
    dt_end   = datetime.strptime(end_date, "%Y-%m-%d")
    SI_proc  = SeaIceProcessor(P_json   = json_path,
                               P_log    = log_file,
                               dt0_str  = dt_start,
                               dtN_str  = dt_end,
                               sim_name = sim_name)
    SI_proc.process_daily_cice(ispd_thresh=ispd_thresh, ispd_type=ispd_type)
    if rolling:
        SI_proc.process_rolling_cice(mean_period=roll_period, ispd_thresh=ispd_thresh, ispd_type=ispd_type)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run pack ice processor loop over time.")
    parser.add_argument("sim_name"      , help="Name of simulation")
    parser.add_argument("--ispd_thresh" , help="ice speed threshold for masking fast ice (default: use hard-coded value in JSON file)")
    parser.add_argument("--ispd_type"   , default=None, help="name of the ice speed masking type: 'ispd_B', 'isp_Ta', 'ispd_Tx' (default=None; masks for each ice speed)")
    parser.add_argument("--no_rolling"  , action="store_true", help="Disable rolling averages")
    parser.add_argument("--mean_period" , type=int , default=14 , help="averaging period in days over which to perform rolling calculation")
    parser.add_argument("--json_config" , help="Path to JSON config file (default: use hard-coded JSON file in fast_ice_processor)")
    parser.add_argument("--log_file"    , help="Path to log file (default: use hard-coded JSON file in fast_ice_processor)")
    parser.add_argument("--start_date"  , help="Start date (YYYY-MM-DD), which is then added to FI_days as the first center-date", default="1993-01-01")
    parser.add_argument("--end_date"    , help="End date (YYYY-MM-DD), will stop processing when this end_date-FI_days", default="1999-12-31")
    args = parser.parse_args()
    run_loop(sim_name    = args.sim_name,
             ispd_thresh = args.ispd_thresh,
             ispd_type   = args.ispd_type,
             json_path   = args.json_config,
             start_date  = args.start_date,
             end_date    = args.end_date,
             log_file    = args.log_file,
             rolling     = not args.no_rolling,
             roll_period = args.mean_period)
