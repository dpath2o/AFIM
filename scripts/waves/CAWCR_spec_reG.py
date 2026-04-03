#!/usr/bin/env python3
import sys
from pathlib import Path
mod_path = "/home/581/da1339/AFIM/src/AFIM/src"
if mod_path not in sys.path:
    sys.path.insert(0, mod_path)
from sea_ice_toolbox import SeaIceToolboxManager

def main(year: int, month: int, overwrite_nc: bool = False, overwrite_weights: bool = False):
    D_cfg      = Path.home() / "AFIM" / "src" / "AFIM" / "src" / "JSONs"
    ice_type   = "FI"
    hemisphere = "south"
    sim_name   = "LD-waves"
    BorC2T     = "Tc"
    P_log      = Path.home() / "logs" / f"waves_setup_{year}{month:02d}.log"
    mgr        = SeaIceToolboxManager(P_log=P_log)
    tb         = mgr.get_toolbox(P_json         = D_cfg / "sea_ice_free-slip.json",
                                 sim_name       = sim_name,
                                 dt0_str        = f"{year:04d}-01-01",
                                 dtN_str        = f"{year:04d}-12-31",
                                 list_of_BorC2T = BorC2T,
                                 ice_type       = ice_type,
                                 hemisphere     = hemisphere)

    p_out      = tb.prepare_cawcr_wave_month(year, month,
                                             overwrite         = overwrite_nc,
                                             overwrite_weights = overwrite_weights,
                                             k                 = 5,
                                             power             = 2.5,
                                             radius_km         = 1000)
    print(p_out)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Prepare monthly CAWCR wave spectra remapped to the CICE grid.")
    parser.add_argument("year"    , type=int, help="4-digit year, e.g. 1993")
    parser.add_argument("month"   , type=int, help="Month number, 1-12")
    parser.add_argument("--ow_nc" , action="store_true", help="Overwrite existing monthly NetCDF output")
    parser.add_argument("--ow_wgt", action="store_true", help="Overwrite existing interpolation weights")
    args = parser.parse_args()
    if not (1 <= args.month <= 12):
        parser.error("month must be in the range 1..12")
    main(year=args.year, month=args.month, overwrite_nc=args.ow_nc, overwrite_weights=args.ow_wgt)
