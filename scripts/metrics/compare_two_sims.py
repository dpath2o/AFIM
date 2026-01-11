#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare two CICE simulations (e.g., no-slip vs free-slip) for a given ice_type (default: SI).

Outputs:
- CSV of key daily time series (SIA, coastal speed, etc.)
- PNG plots for time series
- PNG maps of annual-mean differences (aice, speed, shear, sigP if present)

Notes:
- Designed for one-year experiments; uses lazy dask operations for time series.
- If you provide --nsidc_path, will also plot NSIDC SH SIA (native-grid aggregation; no regridding).
"""

import sys
from pathlib import Path
import argparse

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

# AFIM toolbox
sys.path.insert(0, str(Path.home() / "AFIM" / "src" / "AFIM" / "src"))
from sea_ice_toolbox import SeaIceToolboxManager  # noqa: E402

def _guess_spatial_dims(da: xr.DataArray):
    """Return likely horizontal dims for a DataArray (excluding time/cat dims)."""
    drop = {"time", "time_bounds", "nbnd", "bounds", "ncat", "NCAT"}
    return [d for d in da.dims if d not in drop]


def _coastal_band_from_tmask(tmask: xr.DataArray, band: int = 1) -> xr.DataArray:
    """
    Define coastal ocean cells: ocean cells within `band` grid-cells of land
    using 4-neighbour adjacency iterated `band` times.
    """
    # Ocean/land
    ocean = (tmask > 0)
    land = ~ocean

    # 4-neighbour land adjacency
    spatial = _guess_spatial_dims(tmask)
    if len(spatial) < 2:
        raise ValueError("tmask must be 2D (or 3D with time handled upstream).")

    ydim, xdim = spatial[-2], spatial[-1]
    nbr_land = (
        land.shift({xdim: 1}, fill_value=False)
        | land.shift({xdim: -1}, fill_value=False)
        | land.shift({ydim: 1}, fill_value=False)
        | land.shift({ydim: -1}, fill_value=False)
    )

    # Expand band iteratively (dilating land mask)
    dilated = nbr_land.copy()
    for _ in range(max(1, int(band)) - 1):
        dilated = (
            dilated
            | dilated.shift({xdim: 1}, fill_value=False)
            | dilated.shift({xdim: -1}, fill_value=False)
            | dilated.shift({ydim: 1}, fill_value=False)
            | dilated.shift({ydim: -1}, fill_value=False)
        )

    coastal = ocean & dilated
    coastal.name = "coastal_band"
    return coastal


def _area_weighted_mean(da: xr.DataArray, area: xr.DataArray, mask: xr.DataArray | None = None):
    spatial = _guess_spatial_dims(da)
    w = area
    if mask is not None:
        da = da.where(mask)
        w = w.where(mask)
    num = (da * w).sum(spatial, skipna=True)
    den = w.sum(spatial, skipna=True)
    return num / den


def _sum_over_mask(mask: xr.DataArray, area: xr.DataArray):
    spatial = _guess_spatial_dims(area)
    return (mask.astype("float32") * area).sum(spatial, skipna=True)


def _safe_get(ds: xr.Dataset, name: str):
    return ds[name] if name in ds.variables else None


def _make_timeseries_df(time, **series):
    df = pd.DataFrame({"time": pd.to_datetime(time)})
    for k, v in series.items():
        if v is None:
            continue
        df[k] = np.asarray(v)
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sim1", required=True)
    ap.add_argument("--sim2", required=True)
    ap.add_argument("--dt0", default="1993-01-01")
    ap.add_argument("--dtN", default="1993-12-31")
    ap.add_argument("--ice_type", default="SI", choices=["SI", "FI", "PI"])
    ap.add_argument("--hemisphere", default="SH", choices=["SH", "NH"])
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--coast_band", type=int, default=1, help="Coastal band width in grid cells (default 1)")
    ap.add_argument("--nsidc_path", default=None, help="Optional path to NSIDC daily SIC (zarr or netcdf)")
    ap.add_argument("--nsidc_var", default="siconc", help="NSIDC concentration var name (default siconc)")
    ap.add_argument("--nsidc_area_var", default="area", help="NSIDC cell-area var name (default area)")
    args = ap.parse_args()

    sim1, sim2 = args.sim1, args.sim2

    outdir = Path(args.outdir) if args.outdir else (Path.home() / "AFIM_archive" / "_comparisons" / f"{sim1}_VS_{sim2}" / args.ice_type)
    outdir.mkdir(parents=True, exist_ok=True)

    P_log = outdir / f"compare_{args.ice_type}_{sim1}_VS_{sim2}.log"
    mgr = SeaIceToolboxManager(P_log=P_log)

    # Keep separate toolboxes so metadata (paths, thresholds) remain intact
    tb1 = mgr.get_toolbox(sim_name=sim1, ice_type=args.ice_type, dt0_str=args.dt0, dtN_str=args.dtN, hemisphere=args.hemisphere)
    tb2 = mgr.get_toolbox(sim_name=sim2, ice_type=args.ice_type, dt0_str=args.dt0, dtN_str=args.dtN, hemisphere=args.hemisphere)

    mask_name = f"{args.ice_type}_mask" if args.ice_type in ("FI", "PI", "SI") else "SI_mask"
    # In your classifier, masks are named FI_mask / PI_mask / SI_mask
    mask_name = f"{args.ice_type}_mask".replace("si_mask", "SI_mask")  # safety
    mask_name = f"{args.ice_type}_mask".upper().replace("SI_MASK", "SI_mask").replace("FI_MASK","FI_mask").replace("PI_MASK","PI_mask")
    # Above can be messy; enforce explicitly:
    mask_name = {"SI": "SI_mask", "FI": "FI_mask", "PI": "PI_mask"}[args.ice_type]

    # Load classified daily masks
    msk1 = tb1.load_classified_ice(bin_days=False, roll_mean=False, variables=[mask_name])[mask_name]
    msk2 = tb2.load_classified_ice(bin_days=False, roll_mean=False, variables=[mask_name])[mask_name]

    # Load model variables needed for diagnostics
    vars_req = ["aice", "tarea", "tmask", "uvel", "vvel", "shear", "divu", "vort", "sigP", "sig1", "sig2", "daidtd", "daidtt", "dvidtd", "dvidtt"]
    ds1 = tb1.load_cice_zarr(variables=vars_req, slice_hem=True)
    ds2 = tb2.load_cice_zarr(variables=vars_req, slice_hem=True)

    # Align all on common time index
    ds1, ds2 = xr.align(ds1, ds2, join="inner")
    msk1, msk2 = xr.align(msk1, msk2, join="inner")

    # Shared time
    time = ds1["time"].values

    # Spatial area and masks
    area1 = ds1["tarea"]
    area2 = ds2["tarea"]
    # Use area from sim1 for all aggregated comparisons (should be identical grids)
    area = area1

    ice_union = (msk1 > 0) | (msk2 > 0)

    # Sea ice area (SIA) using the classified mask (icon_thresh already applied)
    sia1 = _sum_over_mask(msk1 > 0, area)
    sia2 = _sum_over_mask(msk2 > 0, area)
    # scale to million km^2 if toolbox provides SIC_scale, else assume 1e12 m^2
    scale = getattr(tb1, "SIC_scale", 1e12) or 1e12
    sia1_mkm2 = (sia1 / scale).compute()
    sia2_mkm2 = (sia2 / scale).compute()

    # Speed + coastal band diagnostics
    u1, v1 = ds1.get("uvel"), ds1.get("vvel")
    u2, v2 = ds2.get("uvel"), ds2.get("vvel")
    speed1 = xr.apply_ufunc(np.hypot, u1, v1) if (u1 is not None and v1 is not None) else None
    speed2 = xr.apply_ufunc(np.hypot, u2, v2) if (u2 is not None and v2 is not None) else None

    coastal = _coastal_band_from_tmask(ds1["tmask"], band=args.coast_band)

    # Means within (a) all-ice union and (b) coastal-ice union
    speed1_ice = _area_weighted_mean(speed1, area, mask=ice_union) if speed1 is not None else None
    speed2_ice = _area_weighted_mean(speed2, area, mask=ice_union) if speed2 is not None else None
    speed1_coast = _area_weighted_mean(speed1, area, mask=ice_union & coastal) if speed1 is not None else None
    speed2_coast = _area_weighted_mean(speed2, area, mask=ice_union & coastal) if speed2 is not None else None

    # Shear / divergence / pressure (area-weighted means)
    shear1 = _safe_get(ds1, "shear")
    shear2 = _safe_get(ds2, "shear")
    sigP1  = _safe_get(ds1, "sigP")
    sigP2  = _safe_get(ds2, "sigP")
    divu1  = _safe_get(ds1, "divu")
    divu2  = _safe_get(ds2, "divu")

    shear1_ice = _area_weighted_mean(shear1, area, mask=ice_union) if shear1 is not None else None
    shear2_ice = _area_weighted_mean(shear2, area, mask=ice_union) if shear2 is not None else None
    sigP1_ice  = _area_weighted_mean(sigP1,  area, mask=ice_union) if sigP1  is not None else None
    sigP2_ice  = _area_weighted_mean(sigP2,  area, mask=ice_union) if sigP2  is not None else None
    divu1_ice  = _area_weighted_mean(divu1,  area, mask=ice_union) if divu1  is not None else None
    divu2_ice  = _area_weighted_mean(divu2,  area, mask=ice_union) if divu2  is not None else None

    # Tendencies integrated over ice (useful for "started from zero ice" experiments)
    def _integrated_tendency(ds, name):
        if name not in ds:
            return None
        # tendency is per time step? we just area-integrate within union ice each day
        return (ds[name].where(ice_union) * area).sum(_guess_spatial_dims(area), skipna=True)

    daidtd1 = _integrated_tendency(ds1, "daidtd")
    daidtt1 = _integrated_tendency(ds1, "daidtt")
    daidtd2 = _integrated_tendency(ds2, "daidtd")
    daidtt2 = _integrated_tendency(ds2, "daidtt")

    # Collect to dataframe (compute only small 1D arrays)
    df = _make_timeseries_df(
        time,
        sia_mkm2_sim1=sia1_mkm2,
        sia_mkm2_sim2=sia2_mkm2,
        dsia_mkm2=(sia2_mkm2 - sia1_mkm2),
        speed_ice_sim1=(speed1_ice.compute() if speed1_ice is not None else None),
        speed_ice_sim2=(speed2_ice.compute() if speed2_ice is not None else None),
        speed_coast_sim1=(speed1_coast.compute() if speed1_coast is not None else None),
        speed_coast_sim2=(speed2_coast.compute() if speed2_coast is not None else None),
        shear_ice_sim1=(shear1_ice.compute() if shear1_ice is not None else None),
        shear_ice_sim2=(shear2_ice.compute() if shear2_ice is not None else None),
        sigP_ice_sim1=(sigP1_ice.compute() if sigP1_ice is not None else None),
        sigP_ice_sim2=(sigP2_ice.compute() if sigP2_ice is not None else None),
        divu_ice_sim1=(divu1_ice.compute() if divu1_ice is not None else None),
        divu_ice_sim2=(divu2_ice.compute() if divu2_ice is not None else None),
        daidtd_int_sim1=(daidtd1.compute() if daidtd1 is not None else None),
        daidtt_int_sim1=(daidtt1.compute() if daidtt1 is not None else None),
        daidtd_int_sim2=(daidtd2.compute() if daidtd2 is not None else None),
        daidtt_int_sim2=(daidtt2.compute() if daidtt2 is not None else None),
    )

    csv_path = outdir / "timeseries.csv"
    df.to_csv(csv_path, index=False)

    # Skill statistics (model-model) on SIA
    try:
        sia_skill = tb1.compute_skill_statistics(df["sia_mkm2_sim2"], df["sia_mkm2_sim1"])
    except Exception:
        sia_skill = {}

    # Optional NSIDC SIA
    nsidc_sia = None
    if args.nsidc_path:
        p = Path(args.nsidc_path).expanduser()
        if p.exists():
            if p.suffix in (".nc", ".nc4", ".cdf", ".netcdf"):
                dsN = xr.open_dataset(p)
            else:
                dsN = xr.open_zarr(p, consolidated=True)
            if args.nsidc_var in dsN and args.nsidc_area_var in dsN:
                conc = dsN[args.nsidc_var]
                areaN = dsN[args.nsidc_area_var]
                # try SH subset if a latitude coord exists; otherwise assume already SH
                if "lat" in conc.coords:
                    conc = conc.where(conc["lat"] < 0, drop=True)
                    areaN = areaN.where(areaN["lat"] < 0, drop=True)
                mskN = conc >= 0.15
                nsidc_sia = (_sum_over_mask(mskN, areaN) / scale).compute()
        # else silently skip (script still useful for model-model comparisons)

    # --- Plots ---
    # 1) SIA time series
    fig = plt.figure(figsize=(10, 4))
    ax = fig.add_subplot(111)
    ax.plot(df["time"], df["sia_mkm2_sim1"], label=sim1)
    ax.plot(df["time"], df["sia_mkm2_sim2"], label=sim2)
    if nsidc_sia is not None:
        ax.plot(df["time"][:len(nsidc_sia)], nsidc_sia, label="NSIDC")
    ax.set_ylabel("SIA (million km$^2$)")
    ax.set_title(f"{args.ice_type} SIA (SH) — {sim1} vs {sim2}")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "SIA_timeseries.png", dpi=200)
    plt.close(fig)

    # 2) Coastal speed time series
    if "speed_coast_sim1" in df.columns and df["speed_coast_sim1"].notna().any():
        fig = plt.figure(figsize=(10, 4))
        ax = fig.add_subplot(111)
        ax.plot(df["time"], df["speed_coast_sim1"], label=f"{sim1} coastal")
        ax.plot(df["time"], df["speed_coast_sim2"], label=f"{sim2} coastal")
        ax.set_ylabel("Area-weighted speed (m/s)")
        ax.set_title(f"Coastal-band ({args.coast_band} cells) mean speed within ice union")
        ax.legend()
        fig.tight_layout()
        fig.savefig(outdir / "coastal_speed_timeseries.png", dpi=200)
        plt.close(fig)

    # 3) Annual mean difference maps (free-slip - no-slip by convention: sim2 - sim1)

    def _xy_coords(A: xr.DataArray):
        spatial = _guess_spatial_dims(A)
        ydim, xdim = spatial[-2], spatial[-1]
        if xdim in A.coords:
            x = A.coords[xdim]
        else:
            x = xr.DataArray(np.arange(A.sizes[xdim]), dims=(xdim,))
        if ydim in A.coords:
            y = A.coords[ydim]
        else:
            y = xr.DataArray(np.arange(A.sizes[ydim]), dims=(ydim,))
        return x, y, xdim, ydim

    def _plot_mean_diff(varname, transform=None, mask_for_plot=None):
        if varname not in ds1 or varname not in ds2:
            return
        A = ds2[varname] - ds1[varname]
        if transform is not None:
            A = transform(A)
        if mask_for_plot is not None:
            A = A.where(mask_for_plot)
        Amean = A.mean("time").load()
        x, y, xdim, ydim = _xy_coords(Amean)
        fig = plt.figure(figsize=(6, 5))
        ax = fig.add_subplot(111)
        ax.pcolormesh(x, y, Amean, shading="auto")
        ax.set_title(f"Annual mean Δ{varname} ({sim2} − {sim1})")
        fig.tight_layout()
        fig.savefig(outdir / f"map_d{varname}_annual_mean.png", dpi=200)
        plt.close(fig)

    if speed1 is not None and speed2 is not None:
        _plot_mean_diff("aice", mask_for_plot=ice_union.mean("time") > 0)  # show where ice occurs at all
        # speed as derived: plot delta speed
        dspeed = (speed2 - speed1).where(ice_union)
        dspeed_mean = dspeed.mean("time").load()
        x, y, xdim, ydim = _xy_coords(dspeed_mean)
        fig = plt.figure(figsize=(6, 5))
        ax = fig.add_subplot(111)
        ax.pcolormesh(x, y, dspeed_mean, shading="auto")
        ax.set_title(f"Annual mean Δspeed ({sim2} − {sim1})")
        fig.tight_layout()
        fig.savefig(outdir / "map_dspeed_annual_mean.png", dpi=200)
        plt.close(fig)

    for vn in ["shear", "divu", "sigP"]:
        _plot_mean_diff(vn, mask_for_plot=ice_union.mean("time") > 0)

    # Write a small text summary
    summary_path = outdir / "summary.txt"
    with open(summary_path, "w") as f:
        f.write(f"Comparison: {sim1} vs {sim2}\n")
        f.write(f"ice_type: {args.ice_type} | hemisphere: {args.hemisphere} | dt0={args.dt0} dtN={args.dtN}\n\n")
        f.write("Model–model SIA skill (sim2 vs sim1):\n")
        for k, v in (sia_skill or {}).items():
            f.write(f"  {k}: {v}\n")
        f.write("\nOutputs:\n")
        f.write(f"  {csv_path}\n")
        f.write(f"  {outdir / 'SIA_timeseries.png'}\n")
        f.write(f"  {outdir / 'coastal_speed_timeseries.png'}\n")

    mgr.shutdown()


if __name__ == "__main__":
    main()