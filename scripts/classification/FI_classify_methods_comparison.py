#!/usr/bin/env python3
import os, sys, json, math
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
from datetime import datetime
import matplotlib
matplotlib.use("Agg")  # headless on Gadi
import matplotlib.pyplot as plt

# --- AFIM imports ---
mod_path = str(Path.home() / "AFIM/src/AFIM/src")
if mod_path not in sys.path:
    sys.path.insert(0, mod_path)
from sea_ice_toolbox import SeaIceToolboxManager

# -------------------- CLI --------------------
import argparse
ap = argparse.ArgumentParser()
ap.add_argument("--sim-name", default="elps-min")
ap.add_argument("--dt0", default="1994-01-01")
ap.add_argument("--dtN", default="1999-12-31")
ap.add_argument("--ispd", type=float, default=5e-4, help="ice_speed_threshold")
ap.add_argument("--win", type=int, default=None, help="window length; if omitted, do merge-only or ALL")
ap.add_argument("--cnt-span", type=int, default=3, help="for bin-day: cnt=win,win-1,... span")
ap.add_argument("--out-dir", type=str, default=str(Path(os.environ.get("SCRATCH", "/tmp")) / "afim/methods_comp/elps-min"))
ap.add_argument("--merge-only", action="store_true", help="aggregate CSVs and write LaTeX + Taylor PNG")
args = ap.parse_args()

OUTDIR = Path(args.out_dir)
OUTDIR.mkdir(parents=True, exist_ok=True)
HOME = Path.home()
FINAL_TEX = HOME / "FI_classify_methods_comparison.tex"
FINAL_PNG = HOME / "graphical/AFIM/elps-min/TaylorDiagram_bin-day_v_roll-mean_methods.png"
FINAL_PNG.parent.mkdir(parents=True, exist_ok=True)

# -------------------- helpers --------------------
def iter_year_bounds(dt0_str: str, dtN_str: str):
    y0 = pd.to_datetime(dt0_str).year
    yN = pd.to_datetime(dtN_str).year
    for y in range(y0, yN + 1):
        yield f"{y}-01-01", f"{y}-12-31"

def ensure_float(x):
    try:
        return float(x)
    except Exception:
        return np.nan

def taylor_diagram(df: pd.DataFrame, out_png: Path):
    """
    Expect columns: label, Corr, SD_Model, SD_Obs.
    Plots normalized SD (sigma/sigma_obs) vs correlation (angle = arccos r).
    """
    d = df.dropna(subset=["Corr","SD_Model","SD_Obs"]).copy()
    if d.empty:
        return
    d["sigma_norm"] = d["SD_Model"] / d["SD_Obs"]
    # polar transform
    theta = np.arccos(np.clip(d["Corr"].values, -1, 1))
    r = d["sigma_norm"].values
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111, projection="polar")
    # Reference std = 1 circle and correlation spokes
    ax.set_theta_zero_location("E")
    ax.set_theta_direction(-1)  # 0° at (1,0), increasing clockwise
    ax.set_thetagrids([0, 30, 60, 90], labels=["1.0","0.866","0.5","0.0"])  # correlation labels
    # radial ticks for normalized std
    rticks = np.linspace(0, max(1.0, np.nanmax(r))*1.1, 5)[1:]
    ax.set_rgrids(rticks)
    ax.set_title("Taylor Diagram: Binary-day vs Rolling-mean (FIA vs AF2020)", pad=20)
    # reference point for perfect match
    ax.plot([0], [1.0], marker="*", markersize=12)
    # scatter points
    for (th, rr, lbl) in zip(theta, r, d["label"].values):
        ax.plot(th, rr, marker="o", markersize=5)
        ax.text(th, rr, f" {lbl}", fontsize=7, ha="left", va="center")
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)

def make_latex_table(df: pd.DataFrame) -> str:
    """
    Columns: label, Bias, MAE, RMSE, SD_Model, FIPSI, FIPmax
    """
    df2 = df[["label","Bias","MAE","RMSE","SD_Model","FIPSI","FIPmax"]].copy()
    # formatting
    for c, fmt in [("Bias","{:.1f}"),("MAE","{:.1f}"),("RMSE","{:.1f}"),("SD_Model","{:.1f}"),
                   ("FIPSI","{:.3f}"),("FIPmax","{:.1f}")]:
        df2[c] = df2[c].map(lambda x: (fmt.format(x) if pd.notna(x) else ""))
    # Build LaTeX
    headers = [r"\textbf{Experiment}",
               r"\textbf{\gls{bias}}",
               r"\textbf{\gls{mae}}",
               r"\textbf{\gls{rmse}}",
               r"\textbf{STD model}",
               r"\textbf{\gls{fipsi}}",
               r"\textbf{\gls{fipmax}}",]
    units = ["", "", "", "", "", "", r"\si{\kilo\metre}"]
    align = "lrrrrrr"
    body  = "\n".join(" & ".join(map(str, row.values)) + r" \\" for _, row in df2.iterrows())
    latex = rf"""\begin{table}[htbp]
\centering
\caption{{Binary-day and rolling-mean FIA skill and persistence geometry metrics (1994--2023) for \gls{{elpsmin}} vs \glsentryshort{{af2020}}.}}
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
    return latex

def write_partial_csv(rows: list[dict], out_csv: Path):
    df = pd.DataFrame(rows)
    df.to_csv(out_csv, index=False)

def read_all_partials(outdir: Path) -> pd.DataFrame:
    parts = sorted(outdir.glob("metrics_win*.csv"))
    if not parts:
        return pd.DataFrame()
    frames = [pd.read_csv(p) for p in parts]
    df = pd.concat(frames, ignore_index=True)
    # Nice ordering: bin-day first grouped by win asc, then rolling
    def sort_key(lbl: str):
        # \gls{elpsmin}-\gls{bindays}-win9-min8
        is_bin = "bindays" in lbl
        try:
            w = int(lbl.split("win")[1].split("-")[0])
        except Exception:
            w = 999
        return (0 if is_bin else 1, w, lbl)
    df = df.sort_values(by="label", key=lambda s: s.map(sort_key)).reset_index(drop=True)
    return df

# -------------------- core per-window compute --------------------
def compute_for_window(win: int, cnt_span: int):
    sim_name  = args.sim_name
    dt0_str   = args.dt0
    dtN_str   = args.dtN
    vars_load = ['aice','tarea']  # minimal
    P_log     = Path.home() / "logs" / "FI_methods_comp.log"
    SI_mgr    = SeaIceToolboxManager(P_log=P_log)
    SI        = SI_mgr.get_toolbox(sim_name=sim_name, dt0_str=dt0_str, dtN_str=dtN_str, ice_speed_threshold=args.ispd)
    # Load once
    CICE_all = SI.load_cice_zarr(slice_hem=False, variables=vars_load, dt0_str=dt0_str, dtN_str=dtN_str)
    CICE_SO  = SI.slice_hemisphere(CICE_all)
    A        = CICE_SO["tarea"].isel(time=0)
    # Obs FIA (AF2020)
    FIA_obs  = xr.open_dataset(SI.AF_FI_dict['P_AF2020_FIA'], engine="netcdf4")["AF2020"]
    yr0_strs = [f"{y}-01-01" for y in range(pd.to_datetime(dt0_str).year, pd.to_datetime(dtN_str).year + 1)]
    yrN_strs = [f"{y}-12-31" for y in range(pd.to_datetime(dt0_str).year, pd.to_datetime(dtN_str).year + 1)]
    rows = []
    # ---- cnt = win, win-1, win-2 (span) ----
    for cnt in range(win, max(win - (cnt_span - 1), 1) - 1, -1):
        label      = r"\gls{elpsmin}-\gls{bindays}-win%d-min%d" % (win, cnt)
        FI_bin_yr  = []
        FI_roll_yr = []
        for y0, yN in zip(yr0_strs, yrN_strs):
            _, FI_bin, FI_roll = SI.classify_fast_ice(dt0_str=y0,
                                                      dtN_str=yN,
                                                      bin_win_days=win,
                                                      bin_min_days=cnt,
                                                      roll_win_days=win,
                                                      enable_rolling_output=True)
            FI_bin_yr.append(FI_bin)
            FI_roll_yr.append(FI_roll)
        FI_bin     = xr.concat(FI_bin_yr, dim="time").chunk(SI.CICE_dict["FI_chunks"])
        FI_bin_SO  = SI.slice_hemisphere(FI_bin)
        FI_roll    = xr.concat(FI_roll_yr, dim="time").chunk(SI.CICE_dict["FI_chunks"])
        FI_roll_SO = SI.slice_hemisphere(FI_roll)
        aice_bin   = CICE_SO["aice"].where(FI_bin_SO)
        tarea_bin  = CICE_SO["tarea"].where(FI_bin_SO)
        FIA_bin    = SI.compute_hemisphere_ice_area(aice_bin, tarea_bin, ice_area_scale=SI.FIC_scale)
        FIP_bin    = SI.compute_hemisphere_variable_aggregate(aice_bin)
        aice_roll  = CICE_SO["aice"].where(FI_roll_SO)
        tarea_roll = CICE_SO["tarea"].where(FI_roll_SO)
        FIA_roll   = SI.compute_hemisphere_ice_area(aice_roll, tarea_roll, ice_area_scale=SI.FIC_scale)
        FIP_roll   = SI.compute_hemisphere_variable_aggregate(aice_roll)
        # Skills
        S_bin      = SI.compute_skill_statistics(FIA_bin, FIA_obs)
        Bias       = ensure_float(S_bin.get("Bias"))
        MAE        = ensure_float(S_bin.get("MAE"))
        RMSE       = ensure_float(S_bin.get("RMSE"))
        Corr       = ensure_float(S_bin.get("Corr"))
        SDm        = ensure_float(S_bin.get("SD_Model"))
        SDo        = ensure_float(S_bin.get("SD_Obs"))
        FIPSI_bin  = ensure_float(SI.persistence_stability_index(FIP_bin)["persistence_stability_index"])
        FIPmax_bin = ensure_float(SI.persistence_ice_distance_mean_max(FIP_bin)["persistence_max_distance"])
        rows.append(dict(label=label, Bias=Bias, MAE=MAE, RMSE=RMSE, SD_Model=SDm, FIPSI=FIPSI_bin, FIPmax=FIPmax_bin, Corr=Corr, SD_Obs=SDo))
        S_roll      = SI.compute_skill_statistics(FIA_roll, FIA_obs)
        Bias        = ensure_float(S_bin.get("Bias"))
        MAE         = ensure_float(S_bin.get("MAE"))
        RMSE        = ensure_float(S_bin.get("RMSE"))
        Corr        = ensure_float(S_bin.get("Corr"))
        SDm         = ensure_float(S_bin.get("SD_Model"))
        SDo         = ensure_float(S_bin.get("SD_Obs"))
        FIPSI_roll  = ensure_float(SI.persistence_stability_index(FIP_roll)["persistence_stability_index"])
        FIPmax_roll = ensure_float(SI.persistence_ice_distance_mean_max(FIP_roll)["persistence_max_distance"])
        rows.append(dict(label=label, Bias=Bias, MAE=MAE, RMSE=RMSE, SD_Model=SDm, FIPSI=FIPSI_roll, FIPmax=FIPmax_roll, Corr=Corr, SD_Obs=SDo))
    out_csv = OUTDIR / f"metrics_win{win:02d}.csv"
    write_partial_csv(rows, out_csv)

# -------------------- main --------------------
if args.merge_only:
    df = read_all_partials(OUTDIR)
    if df.empty:
        print(f"No partial CSVs found in {OUTDIR}", file=sys.stderr)
        sys.exit(1)
    # LaTeX
    latex = make_latex_table(df)
    FINAL_TEX.write_text(latex)
    print(f"Wrote LaTeX table to {FINAL_TEX}")
    # Taylor diagram
    taylor_diagram(df, FINAL_PNG)
    print(f"Wrote Taylor diagram to {FINAL_PNG}")
    sys.exit(0)
# array task: one window
if args.win is None:
    print("No --win and not --merge-only; nothing to do.", file=sys.stderr)
    sys.exit(2)
compute_for_window(args.win, args.cnt_span)
print(f"Saved partial CSV to {OUTDIR}")
