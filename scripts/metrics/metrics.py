import sys, os, argparse, pygmt
sys.path.insert(0, '/home/581/da1339/AFIM/src/AFIM/src')
from sea_ice_toolbox import SeaIceToolbox, SeaIceToolboxManager
from pathlib import Path
import xarray as xr
import numpy as np
import pandas as pd

def safe_compute_sea_ice_metrics(tb, I_dict, P_mets_zarr):
    """
    Wrapper around tb.compute_sea_ice_metrics that handles edge cases
    where methods return None (e.g., no persistent fast-ice cells found).
    
    If the original method crashes with "'NoneType' object is not a mapping",
    this wrapper will catch it and write an empty metrics file instead.
    """
    try:
        return tb.compute_sea_ice_metrics(I_dict, P_mets_zarr)
    except TypeError as e:
        if "'NoneType' object is not a mapping" in str(e):
            tb.logger.warning(
                f"compute_sea_ice_metrics encountered empty results "
                f"(likely no valid ice data for this classification). "
                f"Writing minimal metrics. Original error: {e}"
            )
            # Create minimal dataset with simulation config
            DS_METS = xr.Dataset()
            
            # Add sim_config as attributes
            if hasattr(tb, 'sim_config') and tb.sim_config:
                for k, v in tb.sim_config.items():
                    DS_METS.attrs[k] = v
            
            # Add placeholder NaN values for expected metrics
            ice_type = tb.ice_type.split('_')[0] if hasattr(tb, 'ice_type') else 'FI'
            placeholder_metrics = [
                f"{ice_type}A", f"{ice_type}V", f"{ice_type}T", f"{ice_type}P",
                "persistence_stability_index", "area_persistent_winter", "area_ever_winter",
                "persistence_mean_distance", "persistence_max_distance",
                "Bias", "RMSE", "MAE", "Corr", "SD_Model", "SD_Obs"
            ]
            for metric in placeholder_metrics:
                DS_METS[metric] = xr.DataArray(np.nan, dims=())
            
            # Write to zarr
            if P_mets_zarr:
                DS_METS.to_zarr(P_mets_zarr, mode="w", consolidated=True, zarr_format=2)
                tb.logger.info(f"Minimal (empty) metrics written to {P_mets_zarr}")
            
            return DS_METS
        else:
            # Re-raise if it's a different TypeError
            raise

def main(sim_name, ispd_thresh, ice_type, BorC2T_type, dt0_str, dtN_str, overwrite_zarr, overwrite_png):
    print(f"ice_type passed is: {ice_type}")
    load_vars = ['aice','tarea','hi','uvel','vvel','strength','dvidtt','daidtt','dvidtd','daidtd',
                 'KuxE','KuxN','KuyE','KuyN','earea','narea','uarea']
    P_log     = Path(Path.home(), "logs", f"metrics_{sim_name}_ispd_thresh{ispd_thresh}.log")
    mgr       = SeaIceToolboxManager(P_log=P_log)
    tb        = mgr.get_toolbox(sim_name             = sim_name,
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
    if not ice_type=="SI":
        I_bin = tb.load_classified_ice(class_method="binary-days")[tb.mask_name]
        I_rol = tb.load_classified_ice(class_method="rolling-mean")[tb.mask_name]
    A = CICE_SO['tarea'].isel(time=0)
    # Apply the mask to the data
    I_daily = CICE_SO.where(I_day)
    if not ice_type=="SI":
        I_rolly = CICE_SO.where(I_rol)
        I_binly = CICE_SO.where(I_bin)
    # Create FI dictionaries for dy, rl, and bn
    I_dy = tb.metrics_data_dict(I_day, I_daily, A)
    if not ice_type=="SI":
        I_rl = tb.metrics_data_dict(I_rol, I_rolly, A)
        I_bn = tb.metrics_data_dict(I_bin, I_binly, A)
    
    # Compute metrics for each FI type - USE SAFE WRAPPER
    tb.define_metrics_zarr(class_method = "raw")
    safe_compute_sea_ice_metrics(tb, I_dy, tb.D_mets_zarr)
    
    if not ice_type=="SI":
        tb.define_metrics_zarr(class_method = "rolling-mean")
        safe_compute_sea_ice_metrics(tb, I_rl, tb.D_mets_zarr)
        
        tb.define_metrics_zarr(class_method = "binary-days")
        safe_compute_sea_ice_metrics(tb, I_bn, tb.D_mets_zarr)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute FIA and FIP metrics, apply boolean mask, and plot spatial + temporal outputs.")
    parser.add_argument("--sim_name", type=str, required=True)
    parser.add_argument("--ispd_thresh", type=float, required=True)
    parser.add_argument("--ice_type", default="FI", help="either FI, PI, SI or MI")
    parser.add_argument("--BorC2T_type", default="Tc", help="must be Tc, Ta, Tb, Tx, B or BT")
    parser.add_argument("--start_date", default="1994-01-01", help="Start date (YYYY-MM-DD), which is then added to I_days as the first center-date")
    parser.add_argument("--end_date", default="1999-12-31", help="End date (YYYY-MM-DD), will stop processing when this end_date-I_days")
    parser.add_argument("--overwrite_zarr", action="store_true")
    parser.add_argument("--overwrite_png", action="store_true")
    parser.set_defaults(compute_boolean=False, overwrite_zarr=False, overwrite_png=False)
    args = parser.parse_args()

    main(args.sim_name,
         args.ispd_thresh,
         args.ice_type,
         args.BorC2T_type,
         args.start_date,
         args.end_date,
         args.overwrite_zarr,
         args.overwrite_png)