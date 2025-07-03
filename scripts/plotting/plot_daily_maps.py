#!/usr/bin/env python
import sys
import argparse
from pathlib import Path
import pandas as pd
import xarray as xr
mod_path = '/home/581/da1339/AFIM/src/AFIM/src'
if mod_path not in sys.path:
    sys.path.insert(0, mod_path)
from sea_ice_plotter import SeaIcePlotter

def plot_monthly_variable_maps(sim_name, ice_type, var_names, var_name_back=None, ispd_thresh=1.0e-3):
    if isinstance(var_names, str):
        var_names = [var_names]
    D_out   = Path(f"/g/data/gv90/da1339/afim_output/{sim_name}/zarr/{ice_type}")
    P_zarrs = sorted(D_out.glob(f"{ice_type}_*.zarr"))
    months  = [f.stem.split("_")[1] for f in P_zarrs]
    for yr_mo_str in months:
        P_zarr = D_out / f"{ice_type}_{yr_mo_str}.zarr"
        if not P_zarr.exists():
            print(f"\u26a0\ufe0f Missing dataset: {P_zarr}")
            continue
        ds              = xr.open_dataset(P_zarr, engine="zarr")
        dt_start        = pd.Timestamp(f"{yr_mo_str}-01")
        dt_end          = dt_start + pd.offsets.MonthEnd(0)
        yr_mo_start_str = dt_start.strftime("%Y-%m-%d")
        yr_mo_end_str   = dt_end.strftime("%Y-%m-%d")
        if ice_type=="SO":
            ice_type_plot = "SI"
        else:
            ice_type_plot = ice_type
        plotter = SeaIcePlotter(sim_name   = sim_name,
                                ice_type   = ice_type_plot,
                                plot_type  = 'regional',
                                hemisphere = 'south',
                                save_fig   = True,
                                show_fig   = False,
                                overwrite  = False )
        for var in var_names:
            if var not in ds:
                print(f"Skipping {var}: not found in {P_zarr}")
                continue
            print(f"Plotting {var} for {yr_mo_str}")
            extra_kwargs = {"ispd_thresh"   : ispd_thresh,
                            "var_name_back" : var_name_back,
                            "series"        : [0.0, 0.5] if var == "ispd" else None}
            if var == "aice" and ice_type=="FI":
                extra_kwargs = {"cmap": "viridis", "series": [0.9, 1], "cmap_reverse": True}
            elif var == "divu":
                extra_kwargs = {"cmap": "mag", "series": [-10, 10], "cmap_reverse": False}
            plotter.plot_map(ds            = ds,
                             var_name      = var,
                             var_name_back = var_name_back,
                             dt0_str       = yr_mo_start_str,
                             dtN_str       = yr_mo_end_str,
                             time_coord_name = "time",
                             lon_coord_name  = "lon",
                             lat_coord_name  = "lat",
                             **{k: v for k, v in extra_kwargs.items() if v is not None} )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot monthly sea ice variable maps using SeaIcePlotter.")
    parser.add_argument("--sim_name"     , required=True, help="Simulation name (e.g. 'baseline')")
    parser.add_argument("--ice_type"     , required=True, choices=["FI", "PI", "SO"], help="Ice type: FI, PI, or SO")
    parser.add_argument("--var_names"    , required=True, nargs='+', help="One or more variable names to plot")
    parser.add_argument("--var_name_back", default=None , help="Optional background variable name")
    parser.add_argument("--ispd_thresh"  , type=float, default=1.0e-3, help="Threshold for ispd plots")
    args = parser.parse_args()
    plot_monthly_variable_maps(sim_name=args.sim_name,
                               ice_type=args.ice_type,
                               var_names=args.var_names,
                               var_name_back=args.var_name_back,
                               ispd_thresh=args.ispd_thresh)
