#!/usr/bin/env python3
import sys
from pathlib import Path

mod_path = "/home/581/da1339/AFIM/src/AFIM/src"
if mod_path not in sys.path:
    sys.path.insert(0, mod_path)
from sea_ice_toolbox import SeaIceToolboxManager

def main(year: int, month: int):
    D_cfg      = Path.home() / "AFIM" / "src" / "AFIM" / "src" / "JSONs"
    ice_type   = "FI"
    hemisphere = "south"
    sim_name   = "waves"
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
    p_out      = tb.prepare_cawcr_wave_month(year, month, overwrite=False, overwrite_weights=False)
    print(p_out)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("Usage: prepare_cawcr_month.py YEAR MONTH")
    main(int(sys.argv[1]), int(sys.argv[2]))
