#!/usr/bin/env python3
"""
Create hemispheric PyGMT animations of CICE sea-ice fields (aice, hi)
for a given date range from daily mean files like:
  access-om3.cice.h.1day.mean.YYYY-MM-DD.nc

- Uses TLON/TLAT (curvilinear) -> regrids to regular lon/lat via xESMF.
- Plots each day using PyGMT (polar stereographic).
- Stitches PNG frames to MP4 via ffmpeg.

Example:
  python cice_pygmt_hemi_movie.py \
    --pattern "/scratch/ol01/da1339/access-om3/work/25km_jra_iaf_ld-dpath2o/dev-MC_25km_jra_iaf-95adec7f/access-om3.cice.h.1day.mean.%Y-%m-%d.nc" \
    --var aice --hemi sh \
    --date0 1993-01-01 --date1 1993-02-25 \
    --out aice_sh_19930101_19930225.mp4
"""

from __future__ import annotations

import argparse
import os
import subprocess
from pathlib import Path
from datetime import datetime, timedelta

import numpy as np
import xarray as xr

# xESMF is the simplest robust way to regrid curvilinear TLON/TLAT to regular lon/lat
import xesmf as xe

import pygmt


FILL_BIG = 1e20


def daterange(date0: str, date1: str):
    d0 = datetime.fromisoformat(date0)
    d1 = datetime.fromisoformat(date1)
    d = d0
    while d <= d1:
        yield d
        d += timedelta(days=1)


def find_existing_files(pattern: str, date0: str, date1: str):
    files = []
    dates = []
    for d in daterange(date0, date1):
        p = d.strftime(pattern)
        if os.path.exists(p):
            files.append(p)
            dates.append(d)
        else:
            print(f"[WARN] Missing: {p}")
    if not files:
        raise FileNotFoundError("No files found for the requested date range/pattern.")
    return files, dates


def _to_180(lon):
    """Convert lon to [-180, 180)"""
    lon = ((lon + 180.0) % 360.0) - 180.0
    return lon


def extract_aice_hi(ds: xr.Dataset, var: str) -> xr.DataArray:
    """
    Return 2D DataArray (nj, ni) for requested variable:
      - aice: prefer 'aice', else sum 'aicen' over category dim (nc/ncat)
      - hi  : prefer 'hi', else compute from 'vicen'/'vice' and aice

    Assumes each file has time length 1; handles both (time,nj,ni) and (nj,ni).
    """
    def _isel_time(da: xr.DataArray) -> xr.DataArray:
        return da.isel(time=0) if "time" in da.dims else da

    if var.lower() == "aice":
        if "aice" in ds.variables:
            da = _isel_time(ds["aice"])
        elif "aicen" in ds.variables:
            da = _isel_time(ds["aicen"])
            cat_dim = "nc" if "nc" in da.dims else ("ncat" if "ncat" in da.dims else None)
            if cat_dim is None:
                raise KeyError("Found aicen but could not identify category dim (nc/ncat).")
            da = da.sum(cat_dim)
        else:
            raise KeyError("Could not find 'aice' or 'aicen' in dataset.")
        da = da.where(np.isfinite(da) & (da < FILL_BIG))
        return da

    if var.lower() == "hi":
        if "hi" in ds.variables:
            da = _isel_time(ds["hi"])
            da = da.where(np.isfinite(da) & (da < FILL_BIG))
            return da

        # Fallback: compute hi = vice / aice (mean thickness over ice-covered area)
        # Try vice first, else sum vicen over categories
        if "aice" in ds.variables:
            aice = _isel_time(ds["aice"])
        elif "aicen" in ds.variables:
            aice = _isel_time(ds["aicen"])
            cat_dim = "nc" if "nc" in aice.dims else ("ncat" if "ncat" in aice.dims else None)
            if cat_dim is None:
                raise KeyError("Found aicen but could not identify category dim (nc/ncat).")
            aice = aice.sum(cat_dim)
        else:
            raise KeyError("Need aice/aicen to compute hi fallback, but neither found.")

        if "vice" in ds.variables:
            vice = _isel_time(ds["vice"])
        elif "vicen" in ds.variables:
            vice = _isel_time(ds["vicen"])
            cat_dim = "nc" if "nc" in vice.dims else ("ncat" if "ncat" in vice.dims else None)
            if cat_dim is None:
                raise KeyError("Found vicen but could not identify category dim (nc/ncat).")
            vice = vice.sum(cat_dim)
        else:
            raise KeyError("Could not find 'hi' and cannot fallback (need vice or vicen).")

        aice = aice.where(np.isfinite(aice) & (aice < FILL_BIG))
        vice = vice.where(np.isfinite(vice) & (vice < FILL_BIG))

        hi = vice / aice.where(aice > 0)
        return hi

    raise ValueError("var must be one of: aice, hi")


def build_regridder(example_file: str, hemi: str, res: float, weights_dir: Path, method: str):
    """
    Build (or reuse) an xESMF regridder from curvilinear TLON/TLAT to regular lon/lat.
    Separate weights per hemisphere because the output grid differs.
    """
    ds0 = xr.open_dataset(example_file)

    lon2 = ds0["TLON"].where(ds0["TLON"] < FILL_BIG)
    lat2 = ds0["TLAT"].where(ds0["TLAT"] < FILL_BIG)

    # Convert lon to [-180,180) for cleaner GMT region handling
    lon2 = xr.apply_ufunc(_to_180, lon2)

    grid_in = xr.Dataset(
        {
            "lon": lon2,
            "lat": lat2,
        }
    )

    if hemi.lower() == "nh":
        lat_out = np.arange(0.0, 90.0 + res, res)
        proj = "S0/90/12c"
        region = [-180, 180, 0, 90]
    elif hemi.lower() == "sh":
        lat_out = np.arange(-90.0, 0.0 + res, res)
        proj = "S0/-90/12c"
        region = [-180, 180, -90, 0]
    else:
        raise ValueError("hemi must be nh or sh")

    lon_out = np.arange(-180.0, 180.0 + res, res)

    grid_out = xr.Dataset(
        {
            "lon": (("lon",), lon_out),
            "lat": (("lat",), lat_out),
        }
    )

    weights_dir.mkdir(parents=True, exist_ok=True)
    weights_file = weights_dir / f"weights_{method}_{hemi.lower()}_{res:g}deg.nc"

    regridder = xe.Regridder(
        grid_in,
        grid_out,
        method=method,
        periodic=True,
        filename=str(weights_file),
        reuse_weights=weights_file.exists(),
    )

    ds0.close()
    return regridder, proj, region


def make_cpt(var: str):
    # Keep ranges fixed so the movie doesn't "breathe"
    if var == "aice":
        pygmt.makecpt(cmap="cividis", series=[0, 1, 0.05], continuous=True)
        label = "aice (fraction)"
    else:
        # Adjust if your hi routinely exceeds 5 m
        pygmt.makecpt(cmap="viridis", series=[0, 5, 0.25], continuous=True)
        label = "hi (m)"
    return label


def plot_frame(grid: xr.DataArray, hemi: str, proj: str, region, var: str, date_str: str, out_png: Path):
    fig = pygmt.Figure()

    # Basemap + coast
    fig.basemap(region=region, projection=proj, frame=["af", f"+t{var.upper()} {hemi.upper()}  {date_str}"])
    fig.coast(shorelines="0.5p,black", land="gray90", water="white")

    # Plot
    fig.grdimage(grid=grid, region=region, projection=proj, nan_transparent=True)

    # Colorbar
    if var == "aice":
        fig.colorbar(frame=['x+l"aice (fraction)"'])
    else:
        fig.colorbar(frame=['x+l"hi (m)"'])

    fig.savefig(str(out_png))
    return


def stitch_ffmpeg(frames_glob: str, out_mp4: str, fps: int):
    cmd = [
        "ffmpeg", "-y",
        "-framerate", str(fps),
        "-pattern_type", "glob",
        "-i", frames_glob,
        "-c:v", "libx264",
        "-pix_fmt", "yuv420p",
        "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2",
        out_mp4,
    ]
    print("[INFO] Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pattern", required=True,
                    help="strftime pattern for daily files, e.g. /path/access-om3.cice.h.1day.mean.%Y-%m-%d.nc")
    ap.add_argument("--var", required=True, choices=["aice", "hi"])
    ap.add_argument("--hemi", required=True, choices=["nh", "sh"])
    ap.add_argument("--date0", required=True, help="YYYY-MM-DD")
    ap.add_argument("--date1", required=True, help="YYYY-MM-DD")
    ap.add_argument("--res", type=float, default=0.25, help="Output lon/lat resolution in degrees (default 0.25)")
    ap.add_argument("--method", default="nearest_s2d",
                    choices=["nearest_s2d", "bilinear"],
                    help="xESMF regridding method (default nearest_s2d)")
    ap.add_argument("--fps", type=int, default=10, help="Frames per second for MP4 (default 10)")
    ap.add_argument("--out", required=True, help="Output MP4 filename")
    ap.add_argument("--workdir", default="movie_work", help="Working directory for frames/weights")
    ap.add_argument("--keep_frames", action="store_true", help="Do not delete PNG frames after stitching")
    args = ap.parse_args()

    files, dates = find_existing_files(args.pattern, args.date0, args.date1)

    workdir = Path(args.workdir) / f"{args.var}_{args.hemi}_{args.date0.replace('-','')}_{args.date1.replace('-','')}"
    frames_dir = workdir / "frames"
    weights_dir = workdir / "weights"
    frames_dir.mkdir(parents=True, exist_ok=True)

    regridder, proj, region = build_regridder(
        example_file=files[0],
        hemi=args.hemi,
        res=args.res,
        weights_dir=weights_dir,
        method=args.method,
    )

    make_cpt(args.var)  # defines current CPT in GMT session

    # Loop frames
    for i, (f, d) in enumerate(zip(files, dates), start=1):
        date_str = d.strftime("%Y-%m-%d")
        out_png = frames_dir / f"frame_{i:04d}.png"

        ds = xr.open_dataset(f)
        da = extract_aice_hi(ds, args.var)

        # Ensure TLON/TLAT-compatible dims; da should be (nj,ni)
        # Attach lon/lat for regridder and apply
        da_rg = regridder(da)

        # Mask hemisphere strictly (optional; mostly redundant due to output grid)
        if args.hemi == "nh":
            da_rg = da_rg.where(da_rg["lat"] >= 0)
        else:
            da_rg = da_rg.where(da_rg["lat"] <= 0)

        plot_frame(
            grid=da_rg,
            hemi=args.hemi,
            proj=proj,
            region=region,
            var=args.var,
            date_str=date_str,
            out_png=out_png,
        )

        ds.close()
        print(f"[INFO] Wrote {out_png.name}")

    # Stitch to MP4
    out_mp4 = str(Path(args.out).resolve())
    stitch_ffmpeg(frames_glob=str(frames_dir / "frame_*.png"), out_mp4=out_mp4, fps=args.fps)

    if not args.keep_frames:
        # Clean up frames (keep weights)
        for p in frames_dir.glob("frame_*.png"):
            p.unlink()
        try:
            frames_dir.rmdir()
        except OSError:
            pass

    print("[DONE] ->", out_mp4)


if __name__ == "__main__":
    main()
