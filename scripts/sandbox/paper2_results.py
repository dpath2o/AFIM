import sys, os, pygmt, importlib
mod_path = '/home/581/da1339/AFIM/src/AFIM/src'
sys.path.insert(0, mod_path)
from sea_ice_toolbox         import SeaIceToolbox, SeaIceToolboxManager
import numpy                 as np
import pandas                as pd
import xarray                as xr
import xesmf                 as xe
import matplotlib.pyplot     as plt
from pathlib                 import Path

def maybe_get(ds, name, *, sim_name="", logger=None):
    """Return ds[name] if present; otherwise None (and optionally log)."""
    if isinstance(ds, dict):
        ok = name in ds
    else:
        ok = (name in ds.data_vars) or (name in ds.variables)
    if ok:
        return ds[name]
    if logger is not None:
        logger.info(f"{sim_name}: skipping '{name}' (not found in metrics dataset)")
    return None
    
def _has_any_data(da) -> bool:
    """True if DataArray has at least one non-NaN value (safe for dask)."""
    if da is None:
        return False
    if isinstance(da, xr.Dataset):
        return False
    try:
        ok = da.notnull().any()
        ok = ok.compute() if hasattr(ok, "compute") else ok
        return bool(ok.item())
    except Exception:
        # If any() can't be evaluated cheaply, fall back to "exists"
        return True

def prune_ts_dict(ts_dict: dict, primary_key: str, time_coord: str = "time") -> dict:
    """
    Keep only entries that:
      - are dict payloads,
      - contain primary_key,
      - have non-empty time dim (if present),
      - contain at least one non-NaN value.
    """
    out = {}
    for sim, payload in (ts_dict or {}).items():
        if not isinstance(payload, dict):
            continue
        if primary_key not in payload:
            continue
        da = payload[primary_key]
        if da is None:
            continue
        if hasattr(da, "dims") and (time_coord in da.dims):
            if da.sizes.get(time_coord, 0) == 0:
                continue
        if not _has_any_data(da):
            continue
        out[sim] = {primary_key: da}
    return out

def safe_timeseries(tb, ts_dict, *,
                    primary_key: str,
                    line_colors: dict[str, str] | None = None,
                    legend_labels: dict[str, str] | None = None,
                    **kwargs):
    """
    Prune ts_dict and skip plotting cleanly if nothing valid remains.
    Also supports per-key color/label overrides.
    """
    time_coord = kwargs.get("time_coord", "time")
    clean = prune_ts_dict(ts_dict, primary_key=primary_key, time_coord=time_coord)
    logger = getattr(tb, "logger", None)
    if not clean:
        msg = f"Skipping {primary_key}: no valid data across simulations after pruning."
        if logger is not None:
            logger.warning(msg)
        else:
            print(f"WARNING: {msg}")
        return None
    # Respect user-supplied ordering if provided, but keep only surviving keys.
    keys2plot_in = kwargs.get("keys2plot", None)
    if keys2plot_in is not None:
        keys2plot = [k for k in keys2plot_in if k in clean]
    else:
        keys2plot = list(clean.keys())
    if not keys2plot:
        msg = f"Skipping {primary_key}: keys2plot had no surviving series after pruning."
        if logger is not None:
            logger.warning(msg)
        else:
            print(f"WARNING: {msg}")
        return None
    kwargs["keys2plot"] = keys2plot
    # Filter overrides to plotted keys only (missing keys are fine).
    if line_colors:
        kwargs["line_colors"] = {k: v for k, v in line_colors.items() if k in keys2plot}
    if legend_labels:
        kwargs["legend_labels"] = {k: v for k, v in legend_labels.items() if k in keys2plot}
    return tb.pygmt_timeseries(clean, primary_key=primary_key, **kwargs)

#---------------------------------------------------------------------------------------------

ice_type   = 'FI' 
dt0, dtN   = '1993-01-01', '1995-12-31'
hemisphere = 'south'
sim_names  = ["elps-min", "LD-static-Cs1e-3", "LD-quad-Cq950", "LD-linear-CL0p25"] #"elps-min-GI", "elps-min-NS", "elps-min-NS-T24", "elps-min-NS-T7", "elps-min-NS-C", "elps-min-FS"]
#sim_names  = ["LD-static", "LD-static-Cs8e-4", "LD-static-Cs1e-3", "LD-quad", "LD-quad-Cq3p0", "LD-linear", "LD-linear-CL0p5", "LD-linear-CL0p1", "LD-sum"]
#sim_names  = ["LD-static", "LD-static-Cs1e-3", "LD-static-taub"]# "LD-linear-CL0p5", "LD-linear-CL0p1", "LD-sum"]
#sim_names  = ["LD-quad", "LD-quad-Cq3p0", "LD-quad-Cq15", "LD-quad-Cq45"]
BorC2T     = ['Tb', 'Tc', 'Tc', 'Tc']#, 'Tc', 'Tc']#, 'Tc', 'Tc']#, 'Tc', 'Tc', 'Tc', 'Tc']
clrs_dict  = {sim_names[0] : "#F72B02",
              sim_names[1] : "#009E73",
              sim_names[2] : "#D55E00",
              sim_names[3] : "#56B4E9",}
labs_dict  = {sim_names[0] : "elps-min (paper1)",
              sim_names[1] : "static @[C_s=10^{-3}@[",
              sim_names[2] : "quad @[C_q=950@[",
              sim_names[3] : "linear @[C_L=0.25@["}
comp_name  = 'LD-comp-GI'
FIA_dict, FIT_dict, FIS_dict, FIMVR_dict, FITVR_dict, FIMAR_dict, FITAR_dict, FIKuE_dict, FIKuN_dict, SIA_dict, SIT_dict = {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}
P_log = Path.home() / "logs" / f"LD_testing.log"
mgr   = SeaIceToolboxManager(P_log = P_log)
for i, sim_name in enumerate(sim_names):
    tb  = mgr.get_toolbox(sim_name       = sim_name, 
                          dt0_str        = dt0,
                          dtN_str        = dtN,
                          list_of_BorC2T = BorC2T[i],
                          ice_type       = ice_type,
                          hemisphere     = hemisphere)
    ds  = tb.load_cice_zarr(variables=['aice','hi','tarea'], slice_hem=True)
    if i == 0:
        NSIDC_ts           = tb.compute_NSIDC_metrics()
        SIA_dict['NSIDC']  = {'SIA': NSIDC_ts['SIA'].sel(time=slice(dt0, dtN))}
        FIA_dict['AF2020'] = {'FIA': tb.load_AF2020_FIA_summary()['FIA_obs']}
    SIA_dict[sim_name] = {'SIA': tb.compute_hemisphere_ice_area(ds['aice'], ds['tarea'],
                                                                ice_area_scale=tb.SIC_scale,
                                                                add_grounded_iceberg_area=False)}
    SIT_dict[sim_name] = {'SIT': tb.compute_hemisphere_ice_thickness(ds['aice'], ds['hi'], ds['tarea'])}
    FI_mets            = tb.load_computed_metrics(BorC2T_type=BorC2T[i]).sel(time=slice(dt0, dtN))
    FIA_dict[sim_name] = {'FIA': FI_mets['FIA']}
    FIT_dict[sim_name] = {'FIT': FI_mets['FIT']}
    FIS_dict[sim_name] = {'FIS': FI_mets['FIS']}
    # --- These two lines were swapped in your snippet (typo) ---
    FIMVR_dict[sim_name] = {'FIMVR': FI_mets['FIMVR']}   # <-- volume rate
    FITVR_dict[sim_name] = {'FITVR': FI_mets['FITVR']}
    # --- These were overwritten in your snippet; keep them separate ---
    FIMAR_dict[sim_name] = {'FIMAR': FI_mets['FIMAR'] / 1e9}
    FITAR_dict[sim_name] = {'FITAR': FI_mets['FITAR'] / 1e9}
    # --- Stress: choose keys that exist; skip if absent ---
    # Prefer magnitudes if you created them
    KuE = maybe_get(FI_mets, 'FIKuE_mag_abs_mean', sim_name=sim_name, logger=tb.logger)
    KuN = maybe_get(FI_mets, 'FIKuN_mag_abs_mean', sim_name=sim_name, logger=tb.logger)
    # else fall back to components if that’s what your metrics store
    if KuE is None:
        KuxE = maybe_get(FI_mets, 'FIKuxE', sim_name=sim_name, logger=tb.logger)
        KuyE = maybe_get(FI_mets, 'FIKuyE', sim_name=sim_name, logger=tb.logger)
        if (KuxE is not None) and (KuyE is not None):
            KuE = np.hypot(KuxE, KuyE)  # magnitude of the mean vector (OK for a quick look)
    if KuN is None:
        KuxN = maybe_get(FI_mets, 'FIKuxN', sim_name=sim_name, logger=tb.logger)
        KuyN = maybe_get(FI_mets, 'FIKuyN', sim_name=sim_name, logger=tb.logger)
        if (KuxN is not None) and (KuyN is not None):
            KuN = np.hypot(KuxN, KuyN)
    if KuE is not None:
        FIKuE_dict[sim_name] = {'FIKuE': KuE}
    if KuN is not None:
        FIKuN_dict[sim_name] = {'FIKuN': KuN}
    # These are scalar attrs/vars; guard them too if needed
    pmax = maybe_get(FI_mets, 'persistence_max_distance', sim_name=sim_name, logger=tb.logger)
    psi  = maybe_get(FI_mets, 'persistence_stability_index', sim_name=sim_name, logger=tb.logger)
    if pmax is not None:
        print(f"{sim_name} maximum persistence distance: {float(pmax.values):.2f} km")
    if psi is not None:
        print(f"{sim_name} persistence stability index: {float(psi.values):.4f}")

P_SIA_plot   = tb.D_graph / "timeseries" / "lateral_drag" / f"SIA_TS_{comp_name}_{hemisphere}_{dt0}_{dtN}.png"
P_SIT_plot   = tb.D_graph / "timeseries" / "lateral_drag" / f"SIT_TS_{comp_name}_{hemisphere}_{dt0}_{dtN}.png"
P_FIA_plot   = tb.D_graph / "timeseries" / "lateral_drag" / f"FIA_TS_{comp_name}_{hemisphere}_{dt0}_{dtN}.png"
P_FIT_plot   = tb.D_graph / "timeseries" / "lateral_drag" / f"FIT_TS_{comp_name}_{hemisphere}_{dt0}_{dtN}.png"
P_FIS_plot   = tb.D_graph / "timeseries" / "lateral_drag" / f"FIS_TS_{comp_name}_{hemisphere}_{dt0}_{dtN}.png"
P_FIMVR_plot = tb.D_graph / "timeseries" / "lateral_drag" / f"FIMVR_TS_{comp_name}_{hemisphere}_{dt0}_{dtN}.png"
P_FITVR_plot = tb.D_graph / "timeseries" / "lateral_drag" / f"FITVR_TS_{comp_name}_{hemisphere}_{dt0}_{dtN}.png"
P_FIMAR_plot = tb.D_graph / "timeseries" / "lateral_drag" / f"FIMAR_TS_{comp_name}_{hemisphere}_{dt0}_{dtN}.png"
P_FITAR_plot = tb.D_graph / "timeseries" / "lateral_drag" / f"FITAR_TS_{comp_name}_{hemisphere}_{dt0}_{dtN}.png"
P_FIKuE_plot = tb.D_graph / "timeseries" / "lateral_drag" / f"FIKuE_TS_{comp_name}_{hemisphere}_{dt0}_{dtN}.png"
P_FIKuN_plot = tb.D_graph / "timeseries" / "lateral_drag" / f"FIKuN_TS_{comp_name}_{hemisphere}_{dt0}_{dtN}.png"
# Make sure directory exists
for P in [P_SIA_plot, P_SIT_plot, P_FIA_plot, P_FIT_plot, P_FIS_plot,
          P_FIMVR_plot, P_FITVR_plot, P_FIMAR_plot, P_FITAR_plot,
          P_FIKuE_plot, P_FIKuN_plot]:
    P.parent.mkdir(parents=True, exist_ok=True)
# set across figures
legend_pos  = "JTL+jTL+o0.2c+w15c"
climatology = False
show_fig    = False
# plots
safe_timeseries(tb, SIA_dict,
                comp_name     = comp_name,
                primary_key   = "SIA",
                climatology   = climatology,
                ylabel        = "Sea Ice Area (@[1\\times10^6\\ \\mathrm{km}^2@[)",
                ylim          = [0, 20],
                P_png         = P_SIA_plot,
                show_fig      = show_fig,
                legend_pos    = legend_pos,
                line_colors   = clrs_dict,
                legend_labels = labs_dict)
safe_timeseries(tb, SIT_dict,
                comp_name     = comp_name,
                primary_key   = "SIT",
                climatology   = climatology,
                ylabel        = "Sea Ice Thickness (@[\\mathrm{m}@[)",
                ylim          = [0, 2],
                ytick_pri     = 1,
                ytick_sec     = 0.5,
                P_png         = P_SIT_plot,
                show_fig      = show_fig,
                legend_pos    = legend_pos,
                line_colors   = clrs_dict,
                legend_labels = labs_dict)
safe_timeseries(tb, FIA_dict,
                comp_name     = comp_name,
                primary_key   = "FIA",
                climatology   = climatology,
                ylabel        = "FIA (@[1\\times10^3\\ \\mathrm{km}^2@[)",
                ylim          = [0, 750],
                P_png         = P_FIA_plot,
                show_fig      = show_fig,
                legend_pos    = legend_pos,
                line_colors   = clrs_dict,
                legend_labels = labs_dict)
safe_timeseries(tb, FIT_dict,
                comp_name     = comp_name,
                primary_key   = "FIT",
                climatology   = climatology,
                ylabel        = "FI Thickness (@[\\mathrm{m}@[)",
                ylim          = [0, 3],
                ytick_pri     = 1,
                ytick_sec     = 0.5,
                P_png         = P_FIT_plot,
                show_fig      = show_fig,
                legend_pos    = legend_pos,
                line_colors   = clrs_dict,
                legend_labels = labs_dict,)
safe_timeseries(tb, FIS_dict,
                comp_name     = comp_name,
                primary_key   = "FIS",
                climatology   = climatology,
                ylabel        = "FI Area-wghtd Internal Ice Pressure (hPa)",
                ylim          = [0, 300],
                ytick_pri     = 50,
                ytick_sec     = 25,
                P_png         = P_FIS_plot,
                show_fig      = show_fig,
                legend_pos    = legend_pos,
                line_colors   = clrs_dict,
                legend_labels = labs_dict)
safe_timeseries(tb, FIMVR_dict,
                comp_name     = comp_name,
                primary_key   = "FIMVR",
                climatology   = climatology,
                ylabel        = "FI Volume Dyn-Growth Rate 10^6 km^3 / day",
                ylim          = [-2, 2],
                ytick_pri     = 0.5,
                ytick_sec     = 0.25,
                P_png         = P_FIMVR_plot,
                show_fig      = show_fig,
                legend_pos    = legend_pos,
                line_colors   = clrs_dict,
                legend_labels = labs_dict)
safe_timeseries(tb, FITVR_dict,
                comp_name     = comp_name,
                primary_key   = "FITVR",
                climatology   = climatology,
                ylabel        = "FI Volume Thermo-Growth Rate 10^6 km^3 / day",
                ylim          = [-10, 10],
                ytick_pri     = 1,
                ytick_sec     = 0.5,
                P_png         = P_FITVR_plot,
                show_fig      = show_fig,
                legend_pos    = legend_pos,
                line_colors   = clrs_dict,
                legend_labels = labs_dict)
safe_timeseries(tb, FIMAR_dict,
                comp_name     = comp_name,
                primary_key   = "FIMAR",
                climatology   = climatology,
                ylabel        = "FI Area Dyn-Growth Rate 10^6 km^2 / day",
                ylim          = [-5, .5],
                ytick_pri     = 0.5,
                ytick_sec     = 0.25,
                P_png         = P_FIMAR_plot,
                show_fig      = show_fig,
                legend_pos    = legend_pos,
                line_colors   = clrs_dict,
                legend_labels = labs_dict)
safe_timeseries(tb, FITAR_dict,
                comp_name     = comp_name,
                primary_key   = "FITAR",
                climatology   = climatology,
                ylabel        = "FI Area Thermo-Growth Rate 10^6 km^2 / day",
                ylim          = [-3, 5],
                ytick_pri     = 0.5,
                ytick_sec     = 0.25,
                P_png         = P_FITAR_plot,
                show_fig      = show_fig,
                legend_pos    = legend_pos,
                line_colors   = clrs_dict,
                legend_labels = labs_dict)
safe_timeseries(tb, FIKuE_dict,
                comp_name     = comp_name,
                primary_key   = "FIKuE",
                climatology   = climatology,
                ylabel        = "FI Area-wghtd Mean-Abs Estwd Lateral Stress (Pa)",
                ylim          = [0, 0.5],
                ytick_pri     = 0.1,
                ytick_sec     = 0.05,
                P_png         = P_FIKuE_plot,
                show_fig      = show_fig,
                legend_pos    = legend_pos,
                line_colors   = clrs_dict,
                legend_labels = labs_dict)
safe_timeseries(tb, FIKuN_dict,
                comp_name     = comp_name,
                primary_key   = "FIKuN",
                climatology   = climatology,
                ylabel        = "FI Area-wghtd Mean-Abs Nthwd Lateral Stress (Pa)",
                ylim          = [0, 0.5],
                ytick_pri     = 0.1,
                ytick_sec     = 0.05,
                P_png         = P_FIKuN_plot,
                show_fig      = show_fig,
                legend_pos    = legend_pos,
                line_colors   = clrs_dict,
                legend_labels = labs_dict)


