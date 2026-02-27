#!/usr/bin/env python
import os
import sys
import argparse
from pathlib import Path
import xarray as xr
import numpy as np
import pandas as pd  # kept (present in current file)
import pygmt          # kept (present in current file)
# AFIM dev path (kept consistent with current codebase usage)
sys.path.insert(0, "/home/581/da1339/AFIM/src/AFIM/src")
from sea_ice_toolbox import SeaIceToolboxManager
DEFAULT_P_JSON = "/home/581/da1339/AFIM/src/AFIM/src/JSONs/sea_ice_config.json"

def safe_compute_sea_ice_metrics(tb, I_dict, P_mets_zarr):
    """
    Wrapper around tb.compute_sea_ice_metrics that handles edge cases where
    methods return None (e.g., no persistent fast-ice cells found).

    If the original method crashes with "'NoneType' object is not a mapping",
    this wrapper will catch it and write an empty metrics file instead.
    """
    try:
        return tb.compute_sea_ice_metrics(I_dict, P_mets_zarr)

    except TypeError as e:
        if "'NoneType' object is not a mapping" not in str(e):
            raise
        tb.logger.warning("compute_sea_ice_metrics encountered empty results "
                          "(likely no valid ice data for this classification). "
                          "Writing minimal metrics.\n"
                          f"Original error: {e}")
        DS_METS = xr.Dataset()
        # Add sim_config as attributes if present
        if hasattr(tb, "sim_config") and tb.sim_config:
            for k, v in tb.sim_config.items():
                DS_METS.attrs[k] = v
        ice_type = tb.ice_type.split("_")[0] if hasattr(tb, "ice_type") else "FI"
        placeholder_metrics = [f"{ice_type}A", f"{ice_type}V", f"{ice_type}T", f"{ice_type}P",
                               "persistence_stability_index",
                               "area_persistent_winter", "area_ever_winter",
                               "persistence_mean_distance", "persistence_max_distance",
                               "Bias", "RMSE", "MAE", "Corr", "SD_Model", "SD_Obs"]
        for metric in placeholder_metrics:
            DS_METS[metric] = xr.DataArray(np.nan, dims=())
        if P_mets_zarr:
            DS_METS.to_zarr(P_mets_zarr, mode="w", consolidated=True, zarr_format=2)
            tb.logger.info(f"Minimal (empty) metrics written to {P_mets_zarr}")
        return DS_METS

def main(P_JSON, sim_name, ispd_thresh, ice_type, BorC2T_type, dt0_str, dtN_str, compute_rolling_mean, overwrite_zarr, overwrite_png):
    print(f"ice_type passed is: {ice_type}")
    print(f"P_JSON: {P_JSON}")
    load_vars = ["aice", "tarea", "hi", "uvel", "vvel", "strength",
                 "dvidtt", "daidtt", "dvidtd", "daidtd",
                 "KuxE", "KuxN", "KuyE", "KuyN",
                 "earea", "narea", "uarea"]
    P_log = Path(Path.home(), "logs", f"metrics_{sim_name}_ispd_thresh{ispd_thresh}.log")
    mgr   = SeaIceToolboxManager(P_log=P_log)
    tb    = mgr.get_toolbox(sim_name             = sim_name,
                            P_json               = P_JSON, 
                            dt0_str              = dt0_str,
                            dtN_str              = dtN_str,
                            ice_type             = ice_type,
                            list_of_BorC2T       = BorC2T_type,
                            ice_speed_threshold  = ispd_thresh,
                            overwrite_zarr       = overwrite_zarr,
                            overwrite_saved_figs = overwrite_png)
    # Load classified ice masks and data
    tb.define_ice_mask_name(ice_type=ice_type)
    I_day   = tb.load_classified_ice(class_method="raw")[tb.mask_name]
    CICE_SO = tb.load_cice_zarr(slice_hem=True, variables=load_vars)
    if ice_type != "SI":
        I_bin = tb.load_classified_ice(class_method="binary-days")[tb.mask_name]
        if compute_rolling_mean:
            I_rol = tb.load_classified_ice(class_method="rolling-mean")[tb.mask_name]
    A = CICE_SO["tarea"].isel(time=0)
    # Apply mask(s)
    I_daily = CICE_SO.where(I_day)
    if ice_type != "SI":
        I_binly = CICE_SO.where(I_bin)
        if compute_rolling_mean:
            I_rolly = CICE_SO.where(I_rol)
    # Build dict(s)
    I_dy = tb.metrics_data_dict(I_day, I_daily, A)
    if ice_type != "SI":
        I_bn = tb.metrics_data_dict(I_bin, I_binly, A)
        if compute_rolling_mean:
            I_rl = tb.metrics_data_dict(I_rol, I_rolly, A)
    # Compute metrics
    tb.define_metrics_zarr(class_method="raw")
    safe_compute_sea_ice_metrics(tb, I_dy, tb.D_mets_zarr)
    if ice_type != "SI":
        tb.define_metrics_zarr(class_method="binary-days")
        safe_compute_sea_ice_metrics(tb, I_bn, tb.D_mets_zarr)
        if compute_rolling_mean:
            tb.define_metrics_zarr(class_method="rolling-mean")
            safe_compute_sea_ice_metrics(tb, I_rl, tb.D_mets_zarr)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute sea-ice metrics from classified masks and CICE Zarr.")
    parser.add_argument("--P_JSON", "--P_json",
                        dest="P_JSON",
                        default=DEFAULT_P_JSON,
                        help=f"Path to AFIM config JSON (default: {DEFAULT_P_JSON})")
    parser.add_argument("--sim_name", type=str, required=True)
    parser.add_argument("--ispd_thresh", type=float, required=True)
    parser.add_argument("--ice_type", default="FI", help="either FI, PI, SI or MI")
    parser.add_argument("--BorC2T_type", default="Tc", help="must be Tc, Ta, Tb, Tx, B or BT")
    parser.add_argument("--start_date", default="1994-01-01", help="Start date (YYYY-MM-DD)")
    parser.add_argument("--end_date", default="1999-12-31", help="End date (YYYY-MM-DD)")
    parser.add_argument("--compute_rolling_mean", action="store_true", help="whether to compute rolling-mean metrics in which case rolling-mean classification has to have been done")
    parser.add_argument("--overwrite_zarr", action="store_true", help="this will clobber existing metrics")
    parser.add_argument("--overwrite_png", action="store_true", help="pretty sure I can get rid of this ...")
    args = parser.parse_args()
    main(args.P_JSON,
         args.sim_name,
         args.ispd_thresh,
         args.ice_type,
         args.BorC2T_type,
         args.start_date,
         args.end_date,
         args.compute_rolling_mean,
         args.overwrite_zarr,
         args.overwrite_png)
