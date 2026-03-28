#!/usr/bin/env python
import os
import sys
import argparse
from pathlib import Path
import xarray as xr
import numpy as np
import pandas as pd
import pygmt
sys.path.insert(0, "/home/581/da1339/AFIM/src/AFIM/src")
from sea_ice_toolbox import SeaIceToolboxManager, SeaIceToolbox

DEFAULT_P_JSON = "/home/581/da1339/AFIM/src/AFIM/src/JSONs/sea_ice_config.json"

def store_needs_compute(store, overwrite=False):
    store = Path(store)
    return overwrite or (not store.exists())

def main(P_JSON, sim_name, ispd_thresh, ice_type, BorC2T_type,
         dt0_str, dtN_str, rolling_mean, overwrite_zarr, overwrite_png,
         bin_win_days=None, bin_min_days=None, mean_period=None):
    D_dask    = os.environ.get("DASK_TEMPORARY_DIRECTORY", f"/scratch/gv90/{os.environ.get('USER', 'da1339')}/dask_tmp")
    load_vars = ["aice", "tarea", "hi", "uvel", "vvel", "strength",
                 "dvidtt", "daidtt", "dvidtd", "daidtd",
                 "KuxE", "KuxN", "KuyE", "KuyN",
                 "earea", "narea", "uarea"]
    P_log     = Path(Path.home(), "logs", f"metrics_{sim_name}_ispd_thresh{ispd_thresh}_BW{bin_win_days}_BM{bin_min_days}.log")
    mgr       = SeaIceToolboxManager(P_log     = P_log,
                                     n_workers = 1,
                                     n_threads = 4,
                                     mem_lim   = "48GB",
                                     process   = False,
                                     D_dask    = D_dask)
    tb        = mgr.get_toolbox(sim_name             = sim_name,
                                P_json               = P_JSON,
                                dt0_str              = dt0_str,
                                dtN_str              = dtN_str,
                                ice_type             = ice_type,
                                list_of_BorC2T       = BorC2T_type,
                                ice_speed_threshold  = ispd_thresh,
                                mean_period          = mean_period,
                                bin_win_days         = bin_win_days,
                                bin_min_days         = bin_min_days,
                                overwrite_zarr       = overwrite_zarr,
                                overwrite_saved_figs = overwrite_png)
    tb.define_toolbox_paths(ice_type     = ice_type,
                            BorC2T_type  = BorC2T_type,
                            ispd_thresh  = ispd_thresh,
                            bin_win_days = bin_win_days,
                            bin_min_days = bin_min_days,
                            mean_period  = mean_period)
    tb.define_ice_mask_name(ice_type=ice_type) # tb.mask_name
    # Resolve output stores by ice type
    if ice_type in {"FI", "PI"}:
        P_mets_raw = tb.P_zarrs[ice_type]["mets"]["raw"]
        P_mets_bin = tb.P_zarrs[ice_type]["mets"]["bin"]
        P_mets_roll = tb.P_zarrs[ice_type]["mets"]["roll"]
    elif ice_type in {"SI", "MIZ"}:
        P_mets_raw = tb.P_zarrs[ice_type]["mets"]
        P_mets_bin = None
        P_mets_roll = None
    else:
        raise ValueError(f"Unsupported ice_type: {ice_type}")
    tb.logger.info(f"write-to following zarr sources:\n\t{P_mets_raw}\n\t{P_mets_bin}\n\t{P_mets_roll}")
    # Decide what actually needs to be computed BEFORE loading large datasets
    do_raw  = store_needs_compute(P_mets_raw, overwrite=overwrite_zarr)
    do_bin  = (ice_type != "SI") and (P_mets_bin is not None) and store_needs_compute(P_mets_bin, overwrite=overwrite_zarr)
    do_roll = (ice_type != "SI") and rolling_mean and (P_mets_roll is not None) and store_needs_compute(P_mets_roll, overwrite=overwrite_zarr)
    if not any([do_raw, do_bin, do_roll]):
        tb.logger.info("All requested metric stores already exist; nothing to do.")
        return
    CICE_SO = tb.load_cice_zarr(slice_hem=True, variables=load_vars)
    A       = CICE_SO["tarea"]#.isel(time=0) if "time" in CICE_SO.dims else CICE_SO['tarea']
    I_day   = None
    I_bin   = None
    I_rol   = None
    if do_raw:
        I_day = tb.load_classified_ice(class_method="raw")[tb.mask_name]
        I_dly = CICE_SO.where(I_day)
        I_dy  = tb.metrics_data_dict(I_day, I_dly, A)
        tb.compute_sea_ice_metrics(I_dy, P_mets_raw, overwrite_zarr=overwrite_zarr)
    if do_bin:
        I_bin = tb.load_classified_ice(class_method = "binary-days",
                                       bin_win_days = bin_win_days,
                                       bin_min_days = bin_min_days,
                                       mean_period  = mean_period)[tb.mask_name]
        I_bly = CICE_SO.where(I_bin)
        I_bn  = tb.metrics_data_dict(I_bin, I_bly, A)
        tb.compute_sea_ice_metrics(I_bn, P_mets_bin, overwrite_zarr=overwrite_zarr)
    if do_roll:
        I_rol = tb.load_classified_ice(class_method = "rolling-mean",
                                       bin_win_days = bin_win_days,
                                       bin_min_days = bin_min_days,
                                       mean_period  = mean_period)[tb.mask_name]
        I_rly = CICE_SO.where(I_rol)
        I_rl  = tb.metrics_data_dict(I_rol, I_rly, A)
        tb.compute_sea_ice_metrics(I_rl, P_mets_roll, overwrite_zarr=overwrite_zarr)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute sea-ice metrics from classified masks and CICE Zarr.")
    parser.add_argument("--P_JSON", "--P_json",
                        dest="P_JSON",
                        default=DEFAULT_P_JSON,
                        help=f"Path to AFIM config JSON (default: {DEFAULT_P_JSON})")
    parser.add_argument("--sim_name", type=str, required=True)
    parser.add_argument("--ispd_thresh", type=float, required=True)
    parser.add_argument("--ice_type", default="FI", help="either FI, PI, SI or MIZ")
    parser.add_argument("--BorC2T_type", nargs="+", choices=["B", "Ta", "Tb", "Tc", "Tx"], default=["Tc"])
    parser.add_argument("--start_date", default="1994-01-01", help="Start date (YYYY-MM-DD)")
    parser.add_argument("--end_date", default="1999-12-31", help="End date (YYYY-MM-DD)")
    parser.add_argument("--bin_win_days", type=int, default=None, help="Binary-days window length")
    parser.add_argument("--bin_min_days", type=int, default=None, help="Binary-days minimum valid days")
    parser.add_argument("--mean_period", type=int, default=None, help="Rolling mean period")
    parser.add_argument("--rolling_mean", action="store_true",
                        help="Whether to compute rolling-mean metrics")
    parser.add_argument("--overwrite_zarr", action="store_true", help="Clobber existing metrics groups")
    parser.add_argument("--overwrite_png", action="store_true", help="Retained for compatibility")
    args = parser.parse_args()
    main(args.P_JSON,
         args.sim_name,
         args.ispd_thresh,
         args.ice_type,
         args.BorC2T_type,
         args.start_date,
         args.end_date,
         args.rolling_mean,
         args.overwrite_zarr,
         args.overwrite_png,
         bin_win_days=args.bin_win_days,
         bin_min_days=args.bin_min_days,
         mean_period=args.mean_period)
