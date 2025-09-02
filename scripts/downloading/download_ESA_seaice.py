#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mirror ESA-CCI Sea Ice Thickness data from CEDA over HTTP.

Matches the filename patterns you used previously, with robust downloading:
- HEAD exists check, then GET with Range resume to .part, atomic rename on success
- Retries with backoff, polite rate limiting
- Calendar-valid days only (handles leap years)
- Sensor-specific year and month bounds (CryoSat-2 vs Envisat; NH/SH L3C grid size)
- Mirrors into: DEST/<LEVEL>/<sensor>/<hem>/<filename.nc>

Defaults:
  base_url = https://dap.ceda.ac.uk/neodc/esacci/sea_ice/data/sea_ice_thickness
  dest     = /g/data/gv90/da1339/SeaIce/ESA_CCI
"""

from __future__ import annotations
import os, sys, time, math, argparse, itertools, datetime as dt
from pathlib import Path
from typing import Iterator, Optional, Tuple
from urllib.parse import urljoin, quote
import urllib.request
import urllib.error

# ---------- Defaults (can be overridden by CLI/env) ----------
BASE_URL_DEF = "https://dap.ceda.ac.uk/neodc/esacci/sea_ice/data/sea_ice_thickness"
DEST_DEF     = "/g/data/gv90/da1339/SeaIce/ESA/CCI"
LEVELS_DEF   = ["L2P", "L3C"]
SENSORS_DEF  = ["envisat", "cryosat2"]
HEMS_DEF     = ["nh", "sh"]
TIMEOUT_DEF  = 30
RETRIES_DEF  = 3
RATE_SEC_DEF = 0.5
# -------------------------------------------------------------

def log(msg: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)

def is_leap(y: int) -> bool:
    return (y % 4 == 0 and y % 100 != 0) or (y % 400 == 0)

def days_in_month(y: int, m: int) -> int:
    if m in (1,3,5,7,8,10,12): return 31
    if m in (4,6,9,11): return 30
    # February
    return 29 if is_leap(y) else 28

def month_iter(y: int, m0: int, mN: int) -> Iterator[int]:
    return range(m0, mN+1)

def date_iter(y: int, m0: int, mN: int) -> Iterator[Tuple[int,int,int]]:
    for m in month_iter(y, m0, mN):
        for d in range(1, days_in_month(y, m)+1):
            yield y, m, d

def sensor_year_bounds(level: str, sensor: str) -> Tuple[int,int]:
    # Matches your script: CS2: 2010–2017, Envisat: 2002–2012
    if sensor.lower() == "cryosat2":
        return 2010, 2017
    # envisat
    return 2002, 2012

def month_bounds(level: str, sensor: str, y: int) -> Tuple[int,int]:
    """
    Matches your script’s month bounds:
      - CS2: 2010: start Nov; 2017: end Apr; otherwise Jan–Dec
      - Envisat: 2002: start Oct; 2012: end Mar; otherwise Jan–Dec
    """
    s = sensor.lower()
    if s == "cryosat2":
        m0, mN = 1, 12
        if y == 2010: m0 = 11
        if y == 2017: mN = 4
        return m0, mN
    # Envisat
    m0, mN = 1, 12
    if y == 2002: m0 = 10
    if y == 2012: mN = 3
    return m0, mN

def build_filename(level: str, sensor: str, hem: str, y: int, m: int, d: Optional[int]) -> str:
    """
    Return the exact filename as per your patterns.
    """
    L = level.upper()
    S = sensor.upper()
    H = hem.upper()

    if L == "L2P":
        # daily
        date_str = f"{y:04d}{m:02d}{d:02d}"
        if sensor.lower() == "cryosat2":
            return f"ESACCI-SEAICE-L2P-SITHICK-SIRAL_{S}-{H}-{date_str}-fv2.0.nc"
        else:
            return f"ESACCI-SEAICE-L2P-SITHICK-RA2_{S}-{H}-{date_str}-fv2.0.nc"

    # L3C monthly (grid resolution differs NH/SH)
    date_str = f"{y:04d}{m:02d}"
    if sensor.lower() == "cryosat2":
        grid = "25KMEASE2" if hem.lower()=="nh" else "50KMEASE2"
        return f"ESACCI-SEAICE-L3C-SITHICK-SIRAL_{S}-{H}{grid}-{date_str}-fv2.0.nc"
    else:
        grid = "25KMEASE2" if hem.lower()=="nh" else "50KMEASE2"
        return f"ESACCI-SEAICE-L3C-SITHICK-RA2_{S}-{H}{grid}-{date_str}-fv2.0.nc"

def build_remote_path(base_url: str, level: str, sensor: str, hem: str, y: int, m: int, d: Optional[int]) -> str:
    """
    Remote directory layout:
      {base}/{LEVEL}/{sensor}/v2.0/{HEM}/{YYYY}/{MM}/filename.nc
    """
    L = level.upper(); H = hem.upper()
    subdir = f"{L}/{sensor}/v2.0/{H}/{y:04d}/{m:02d}/"
    fname  = build_filename(level, sensor, hem, y, m, d)
    # Ensure proper URL joining & encoding
    return urljoin(base_url.rstrip("/") + "/", quote(subdir + fname))

def head_request(url: str, timeout: int) -> Optional[int]:
    """
    Return Content-Length (int) if file exists (HTTP 200) else None.
    """
    req = urllib.request.Request(url, method="HEAD")
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            if 200 <= resp.status < 300:
                length = resp.headers.get("Content-Length")
                return int(length) if length is not None else None
            return None
    except urllib.error.HTTPError as e:
        # 404 etc.
        return None
    except Exception:
        return None

def get_with_resume(url: str, dest: Path, expected_size: Optional[int], timeout: int, retries: int, rate_sec: float) -> bool:
    """
    Download with HTTP Range resume to dest.part, then rename atomically.
    Returns True on success, False on failure.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_suffix(dest.suffix + ".part")
    # If final exists with expected size, skip
    if dest.exists() and (expected_size is None or dest.stat().st_size == expected_size):
        log(f"[SKIP] Exists: {dest}")
        return True

    resume_from = part.stat().st_size if part.exists() else 0
    headers = {}
    if resume_from > 0:
        headers["Range"] = f"bytes={resume_from}-"

    for attempt in range(1, retries+1):
        try:
            req = urllib.request.Request(url, headers=headers, method="GET")
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                # 206 if resume, 200 otherwise; stream to file in chunks
                mode = "ab" if resume_from > 0 else "wb"
                with open(part, mode) as f:
                    while True:
                        chunk = resp.read(1024 * 256)
                        if not chunk:
                            break
                        f.write(chunk)
            # Done: verify size if we know it
            if expected_size is not None:
                if part.stat().st_size != expected_size:
                    raise IOError(f"Size mismatch: got {part.stat().st_size}, expected {expected_size}")
            part.replace(dest)
            log(f"[OK ] {dest.name}")
            time.sleep(rate_sec)
            return True
        except Exception as e:
            log(f"[WARN] GET failed (attempt {attempt}/{retries}) for {dest.name}: {e}")
            time.sleep(rate_sec + 1.0 * attempt)
            # On retry, refresh resume offset (if partial written)
            resume_from = part.stat().st_size if part.exists() else 0
            headers = {"Range": f"bytes={resume_from}-"} if resume_from > 0 else {}

    # failure: clean partial to avoid confusion
    try:
        part.unlink(missing_ok=True)
    except Exception:
        pass
    log(f"[FAIL] {dest.name}")
    return False

def generate_targets(levels, sensors, hems, y_min: Optional[int], y_max: Optional[int]) -> Iterator[Tuple[str,str,str,int,int,Optional[int]]]:
    """
    Yields tuples: (level, sensor, hem, year, month, day_or_None)
    """
    for level in levels:
        for sensor in sensors:
            sy0, syN = sensor_year_bounds(level, sensor)
            if y_min is not None: sy0 = max(sy0, y_min)
            if y_max is not None: syN = min(syN, y_max)
            if sy0 > syN:
                continue
            for hem in hems:
                for y in range(sy0, syN+1):
                    m0, mN = month_bounds(level, sensor, y)
                    if level.upper() == "L2P":
                        for (yy, mm, dd) in date_iter(y, m0, mN):
                            yield level, sensor, hem, yy, mm, dd
                    else:
                        for mm in month_iter(y, m0, mN):
                            yield level, sensor, hem, y, mm, None

def path_local(dest_root: Path, level: str, sensor: str, hem: str, y: int, m: int, d: Optional[int]) -> Path:
    """
    Local path: DEST/<LEVEL>/<sensor>/<hem>/<filename.nc>
    """
    fname = build_filename(level, sensor, hem, y, m, d)
    return dest_root / level.upper() / sensor.lower() / hem.lower() / fname

def main():
    ap = argparse.ArgumentParser(description="Mirror ESA-CCI SIT from CEDA over HTTP with resume.")
    ap.add_argument("--base-url", default=os.environ.get("BASE_URL", BASE_URL_DEF))
    ap.add_argument("--dest",     default=os.environ.get("DLOC",     DEST_DEF))
    ap.add_argument("--levels",   default=os.environ.get("LEVELS", ",".join(LEVELS_DEF)))
    ap.add_argument("--sensors",  default=os.environ.get("SENSORS", ",".join(SENSORS_DEF)))
    ap.add_argument("--hems",     default=os.environ.get("HEMS",    ",".join(HEMS_DEF)))
    ap.add_argument("--year-min", type=int, default=int(os.environ.get("YEAR_MIN")) if os.environ.get("YEAR_MIN") else None)
    ap.add_argument("--year-max", type=int, default=int(os.environ.get("YEAR_MAX")) if os.environ.get("YEAR_MAX") else None)
    ap.add_argument("--timeout",  type=int, default=int(os.environ.get("TIMEOUT", TIMEOUT_DEF)))
    ap.add_argument("--retries",  type=int, default=int(os.environ.get("RETRIES", RETRIES_DEF)))
    ap.add_argument("--rate-sec", type=float, default=float(os.environ.get("RATE_SEC", RATE_SEC_DEF)))
    ap.add_argument("--dry-run",  action="store_true")
    args = ap.parse_args()

    base_url   = args.base_url.rstrip("/")
    dest_root  = Path(args.dest)
    dest_root.mkdir(parents=True, exist_ok=True)

    levels  = [s.strip() for s in args.levels.split(",")  if s.strip()]
    sensors = [s.strip() for s in args.sensors.split(",") if s.strip()]
    hems    = [s.strip() for s in args.hems.split(",")    if s.strip()]

    log(f"Base URL: {base_url}")
    log(f"Dest    : {dest_root}")
    log(f"Levels  : {levels}")
    log(f"Sensors : {sensors}")
    log(f"Hemis   : {hems}")
    log(f"Years   : {args.year_min}..{args.year_max} (None = default bounds)")
    log(f"Rate    : {args.rate_sec}s   Timeout: {args.timeout}s   Retries: {args.retries}")

    n_try = n_ok = n_skip = n_missing = 0
    for level, sensor, hem, y, m, d in generate_targets(levels, sensors, hems, args.year_min, args.year_max):
        url  = build_remote_path(base_url, level, sensor, hem, y, m, d)
        lpth = path_local(dest_root, level, sensor, hem, y, m, d)
        # Skip if exists already
        if lpth.exists():
            n_skip += 1
            continue
        # HEAD to confirm presence and size.  Some servers may reject HEAD
        # requests (e.g. with 403/405).  In that case, fall back to attempting
        # the GET anyway so that we do not incorrectly mark the file as
        # missing.  If the subsequent GET fails we will count it as missing.
        clen = head_request(url, timeout=args.timeout)
        if clen is None:
            # n_missing += 1
            # continue
            log(f"[WARN] HEAD failed for {url}; attempting GET anyway")
        n_try += 1
        if args.dry_run:
            log(f"[DRY] would fetch {url} -> {lpth}")
            time.sleep(args.rate_sec)
            continue
        # ok = get_with_resume(url, lpth, clen, timeout=args.timeout, retries=args.retries, rate_sec=args.rate_sec)
        # if ok: n_ok += 1
        ok = get_with_resume(url, lpth, clen, timeout=args.timeout, retries=args.retries, rate_sec=args.rate_sec)
        if ok:
            n_ok += 1
        else:
            n_missing += 1
    log(f"Summary: tried={n_try}, ok={n_ok}, skipped={n_skip}, missing={n_missing}")
    if n_ok == 0 and n_try > 0:
        sys.exit(2)

if __name__ == "__main__":
    main()
