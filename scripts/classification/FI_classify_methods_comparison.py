#!/usr/bin/env python3
"""
Compare binary-day vs rolling-mean fast-ice classification for one window.

Usage (single window):
  python -u FI_classify_methods_comparison.py \
    --sim-name elps-min --dt0 2000-01-01 --dtN 2018-12-31 \
    --ispd 5e-4 --win 15 --min-days-span 3 \
    --out-dir /scratch/<proj>/<user>/afim/methods_comp/elps-min

Merge step (after multiple windows have written metrics_win*.csv to --out-dir):
  python -u FI_classify_methods_comparison.py --merge-only --out-dir <same>
"""

import os, sys, gc
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import xarray as xr
import dask
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- AFIM imports ------------------------------------------------------------
mod_path = str(Path.home() / "AFIM/src/AFIM/src")
if mod_path not in sys.path:
    sys.path.insert(0, mod_path)
from sea_ice_toolbox import SeaIceToolboxManager

# ---------------- CLI --------------------------------------------------------
ap = argparse.ArgumentParser()
ap.add_argument("--sim-name", default="elps-min")
ap.add_argument("--dt0", default="1994-01-01")
ap.add_argument("--dtN", default="2023-12-31")
ap.add_argument("--ispd", type=float, default=5e-4, help="ice_speed_threshold (m/s)")
ap.add_argument("--win", type=int, default=None, help="window length (days)")
ap.add_argument("--min-days-span", type=int, default=3,
               help="build min_days = win, win-1, ... (span elements)")
ap.add_argument("--out-dir", type=str,
               default=str(Path(os.environ.get("SCRATCH", "/tmp")) / "afim/methods_comp/elps-min"))
ap.add_argument("--merge-only", action="store_true",
               help="aggregate CSVs and write LaTeX + Taylor PNG")
args      = ap.parse_args()
OUTDIR    = Path(args.out_dir); OUTDIR.mkdir(parents=True, exist_ok=True)
HOME      = Path.home()
FINAL_TEX = HOME / "FI_classify_methods_comparison.tex"
FINAL_PNG = HOME / "graphical/AFIM/elps-min/TaylorDiagram_bin-day_v_roll-mean_methods.png"
FINAL_PNG.parent.mkdir(parents=True, exist_ok=True)

WINTER_MONTHS = [5, 6, 7, 8, 9, 10]

# ---------------- helpers (module-level; no nesting) -------------------------
def taylor_diagram(df: pd.DataFrame, out_png: Path):
    d = df.dropna(subset=["Corr","SD_Model","SD_Obs"]).copy()
    if d.empty:
        return
    d["sigma_norm"] = d["SD_Model"] / d["SD_Obs"]
    theta = np.arccos(np.clip(d["Corr"].to_numpy(), -1, 1))
    r = d["sigma_norm"].to_numpy()
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111, projection="polar")
    ax.set_theta_zero_location("E")
    ax.set_theta_direction(-1)
    ax.set_thetagrids([0, 30, 60, 90], labels=["1.0","0.866","0.5","0.0"])
    ax.set_rgrids(np.linspace(0, max(1.0, np.nanmax(r))*1.1, 5)[1:])
    ax.set_title("Taylor Diagram: Binary-day vs Rolling-mean (FIA vs AF2020)", pad=20)
    ax.plot([0], [1.0], marker="*", markersize=12)  # reference
    for th, rr, lbl in zip(theta, r, d["label"].to_numpy()):
        ax.plot(th, rr, marker="o", markersize=5)
        ax.text(th, rr, f" {lbl}", fontsize=7, ha="left", va="center")
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)

def make_latex_table(df: pd.DataFrame) -> str:
    df2 = df[["label","Bias","MAE","RMSE","SD_Model","FIPSI","FIPmax"]].copy()
    for c, fmt in [("Bias","{:.1f}"),("MAE","{:.1f}"),("RMSE","{:.1f}"),
                   ("SD_Model","{:.1f}"),("FIPSI","{:.3f}"),("FIPmax","{:.1f}")]:
        df2[c] = df2[c].map(lambda x: (fmt.format(x) if pd.notna(x) else ""))
    headers = [r"\textbf{Experiment}", r"\textbf{\gls{bias}}", r"\textbf{\gls{mae}}",
               r"\textbf{\gls{rmse}}", r"\textbf{STD model}" , r"\textbf{\gls{fipsi}}",
               r"\textbf{\gls{fipmax}}",]
    units   = ["", "", "", "", "", "", r"\si{\kilo\metre}"]
    align   = "lrrrrrr"
    body    = "\n".join(" & ".join(map(str, row.values)) + r" \\" for _, row in df2.iterrows())
    return rf"""\begin{table}[htbp]
\centering
\caption{{Binary-day and rolling-mean FIA skill and persistence geometry metrics for \gls{{elpsmin}} vs \glsentryshort{{af2020}}.}}
\label{{tab:FI_classify_methods_comparison}}
\begin{adjustbox}{{max width=\linewidth}}
\begin{tabular}{{{align}}}
\toprule
{ " & ".join(headers) } \\
{ " & ".join(units) } \\
\midrule
{body}
\bottomrule
\end{tabular}
\end{adjustbox}
\end{table}
"""

def write_partial_csv(rows: list[dict], out_csv: Path):
    pd.DataFrame(rows).to_csv(out_csv, index=False)

def read_all_partials(outdir: Path) -> pd.DataFrame:
    parts = sorted(outdir.glob("metrics_win*.csv"))
    if not parts:
        return pd.DataFrame()
    df = pd.concat([pd.read_csv(p) for p in parts], ignore_index=True)
    def sort_key(lbl: str):
        is_bin = "bindays" in lbl
        try:
            w = int(lbl.split("win")[1].split("-")[0])
        except Exception:
            w = 999
        return (0 if is_bin else 1, w, lbl)
    return df.sort_values(by="label", key=lambda s: s.map(sort_key)).reset_index(drop=True)

def _extract_number(obj, *keys):
    """Return a float from obj which may be a float, 0-D DataArray, or dict."""
    if isinstance(obj, (int, float, np.floating)):
        return float(obj)
    if isinstance(obj, xr.DataArray):
        v = obj.values
        if np.ndim(v) == 0:
            return float(v)
        raise TypeError(f"Expected 0-D DataArray, got shape {v.shape}")
    if isinstance(obj, dict):
        # Try provided keys first
        for k in keys:
            if k in obj:
                return float(_extract_number(obj[k]))
        # Fallback: first numeric-ish value
        for v in obj.values():
            try:
                return float(_extract_number(v))
            except Exception:
                pass
        raise KeyError(f"None of the keys {keys} found; and no numeric value in {obj}")
    raise TypeError(f"Unsupported type for metric extraction: {type(obj)}")

# ---------------- core compute ------------------------------------------------
def compute_for_window(win: int, min_days_span: int):
    csv_rows = []
    P_log  = Path.home() / "logs" / "FI_methods_comp.log"
    mgr    = SeaIceToolboxManager(P_log=P_log)
    SI     = mgr.get_toolbox(sim_name=args.sim_name, dt0_str=args.dt0, dtN_str=args.dtN, ice_speed_threshold=args.ispd)
    time_dim  = SI.CICE_dict["time_dim"]
    spat_dims = tuple(SI.CICE_dict["spatial_dims"])
    CICE_SO   = SI.load_cice_zarr(slice_hem=True, variables=["aice", "tarea"], dt0_str=args.dt0, dtN_str=args.dtN)
    A0        = CICE_SO["tarea"].isel(time=0)  # static 2-D area
    FIA_obs   = xr.open_dataset(SI.AF_FI_dict["P_AF2020_FIA"], engine="netcdf4")["AF2020"]
    # Build year list
    y0    = pd.to_datetime(args.dt0).year
    yN    = pd.to_datetime(args.dtN).year
    years = range(y0, yN + 1)
    # Min-days variants to evaluate
    min_days_list = list(range(win, max(win - (min_days_span - 1), 1) - 1, -1))
    # Accumulators
    FIA_roll_parts      : list[xr.DataArray] = []
    FIPSI_roll_parts    : list[float]        = []
    FIPDmax_roll_parts  : list[float]        = [] 
    FIA_bin_parts       : dict[int, list[xr.DataArray]] = {m: [] for m in min_days_list}
    FIPSI_bin_parts     : dict[int, list[float]]        = {m: [] for m in min_days_list}
    FIPDmax_bin_parts   : dict[int, list[float]]        = {m: [] for m in min_days_list}
    # ---------------- per-year loop: exactly one classify_fast_ice() per year ----------------
    for yr in years:
        y0_str    = f"{yr}-01-01"
        yN_str    = f"{yr}-12-31"
        CICE_SO_y = CICE_SO.sel(time=slice(y0_str, yN_str))
        aice_y    = CICE_SO_y["aice"]
        tarea_y   = CICE_SO_y["tarea"]
        for bin_min_day in min_days_list:
            _, FI_bin_y, FI_roll_y = SI.classify_fast_ice(dt0_str               = y0_str,
                                                          dtN_str               = yN_str,
                                                          bin_win_days          = win,
                                                          bin_min_days          = bin_min_day,  
                                                          enable_rolling_output = True,
                                                          roll_win_days         = win)
            bin_SO        = SI.slice_hemisphere(FI_bin_y)
            roll_SO       = SI.slice_hemisphere(FI_roll_y)
            bin_mask      = bin_SO["FI_mask"] > 0
            roll_mask     = roll_SO["FI_mask"] > 0
            FIA_bin_y     = SI.compute_hemisphere_ice_area(aice_y.where(bin_mask), A0.where(bin_mask), ice_area_scale=SI.FIC_scale)
            FIP_bin_y     = SI.compute_hemisphere_variable_aggregate(bin_mask)
            FIPSI_bin_y   = SI.persistence_stability_index(bin_mask, A0)
            FIPDmax_bin_y = SI.persistence_ice_distance_mean_max(FIP_bin_y)
            FIA_bin_parts[bin_min_day].append(FIA_bin_y)
            FIPSI_bin_parts[bin_min_day].append(_extract_number(FIPSI_bin_y, "FIPSI", "stability_index", "persistence_stability_index"))
            FIPDmax_bin_parts[bin_min_day].append(_extract_number(FIPDmax_bin_y, "persistence_max_distance", "max_distance"))
        FIA_roll_y     = SI.compute_hemisphere_ice_area(aice_y.where(roll_mask), A0.where(roll_mask), ice_area_scale=SI.FIC_scale)
        FIP_roll_y     = SI.compute_hemisphere_variable_aggregate(roll_mask)
        FIPSI_roll_y   = SI.persistence_stability_index(roll_mask, A0)
        FIPDmax_roll_y = SI.persistence_ice_distance_mean_max(FIP_roll_y)
        FIA_roll_parts.append(FIA_roll_y)
        FIPSI_roll_parts.append(_extract_number(FIPSI_roll_y, "FIPSI", "stability_index", "persistence_stability_index"))
        FIPDmax_roll_parts.append(_extract_number(FIPDmax_roll_y, "persistence_max_distance", "max_distance"))
    # Rolling
    if not FIA_roll_parts:
        raise RuntimeError("No yearly rolling FIA parts were produced.")
    FIA_roll     = xr.concat(FIA_roll_parts, dim=time_dim)
    FIPSI_roll   = float(np.mean(FIPSI_roll_parts))          # list[float] -> float
    FIPDmax_roll = float(np.mean(FIPDmax_roll_parts))        # list[float] -> float
    S_roll       = SI.compute_skill_statistics(FIA_roll, FIA_obs)
    roll_stats   = dict(label    = rf"\gls{{elpsmin}}-\gls{{rollmean}}-win{win}",
                        Bias     = S_roll["Bias"]    , MAE  = S_roll["MAE"] , RMSE   = S_roll["RMSE"],
                        SD_Model = S_roll["SD_Model"], Corr = S_roll["Corr"], SD_Obs = S_roll["SD_Obs"],
                        FIPSI    = FIPSI_roll        , FIPmax = FIPDmax_roll)
    csv_rows.append(roll_stats)
    for bin_min_day in min_days_list:
        FIA_bin     = xr.concat(FIA_bin_parts[bin_min_day], dim=time_dim)
        FIPSI_bin   = float(np.mean(FIPSI_bin_parts[bin_min_day]))
        FIPDmax_bin = float(np.mean(FIPDmax_bin_parts[bin_min_day]))
        S_bin       = SI.compute_skill_statistics(FIA_bin, FIA_obs)
        bin_stats   = dict(label     = rf"\gls{{elpsmin}}-\gls{{bindays}}-win{win}-min{bin_min_day}",
                            Bias     = S_bin["Bias"]    , MAE  = S_bin["MAE"] , RMSE   = S_bin["RMSE"],
                            SD_Model = S_bin["SD_Model"], Corr = S_bin["Corr"], SD_Obs = S_bin["SD_Obs"],
                            FIPSI    = FIPSI_bin        , FIPmax = FIPDmax_bin)
        csv_rows.append(bin_stats)
    write_partial_csv(csv_rows, OUTDIR / f"metrics_win{win:02d}.csv")

# ---------------- main --------------------------------------------------------
if __name__ == "__main__":
    # safer multiprocessing start for py3.11 + dask/distributed nanny
    try:
        import multiprocessing as mp
        if mp.get_start_method(allow_none=True) != "fork":
            mp.set_start_method("fork")
    except Exception:
        pass

    if args.merge_only:
        df = read_all_partials(OUTDIR)
        if df.empty:
            print(f"No partial CSVs found in {OUTDIR}", file=sys.stderr); sys.exit(1)
        FINAL_TEX.write_text(make_latex_table(df))
        taylor_diagram(df, FINAL_PNG)
        print(f"Wrote {FINAL_TEX} and {FINAL_PNG}")
    else:
        if args.win is None:
            print("No --win and not --merge-only; nothing to do.", file=sys.stderr)
            sys.exit(2)
        compute_for_window(args.win, args.min_days_span)
        print(f"Saved partial CSV to {OUTDIR}")
