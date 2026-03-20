#!/usr/bin/env python3
# download_nsidc_g02202.py
from __future__ import annotations

import argparse
import concurrent.futures as cf
import html
import re
import sys
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urljoin
from urllib.request import Request, urlopen

USER_AGENT = "Mozilla/5.0 (compatible; Gadi-NSIDC-downloader/1.0)"
HEMI_TO_GRID = {"north": "psn25", "south": "pss25"}


def fetch_text(url: str, timeout: int = 60) -> str:
    req = Request(url, headers={"User-Agent": USER_AGENT})
    with urlopen(req, timeout=timeout) as resp:
        return resp.read().decode("utf-8", "ignore")


def list_links(url: str) -> list[str]:
    text = fetch_text(url)
    hrefs = re.findall(r'href="([^"]+)"', text, flags=re.I)
    out: list[str] = []
    for href in hrefs:
        href = html.unescape(href)
        if href in ("../", "./") or href.startswith("?"):
            continue
        out.append(href)
    return sorted(set(out))


def list_nc(url: str) -> list[str]:
    return [name for name in list_links(url) if name.endswith(".nc")]


def choose_yearly_daily_aggregate(files: list[str], grid: str, year: int) -> str | None:
    rx = re.compile(rf"sic_{grid}_{year}0101-{year}\d{{4}}_v\d{{2}}r\d{{2}}\.nc$")
    cands = sorted(f for f in files if rx.fullmatch(f))
    return cands[-1] if cands else None


def choose_period_monthly_aggregate(files: list[str], grid: str) -> str | None:
    rx = re.compile(rf"sic_{grid}_\d{{6}}-\d{{6}}_v\d{{2}}r\d{{2}}\.nc$")
    cands = sorted(f for f in files if rx.fullmatch(f))
    return cands[-1] if cands else None


def choose_daily_individual(files: list[str], grid: str, year: int) -> list[str]:
    rx = re.compile(rf"sic_{grid}_{year}\d{{4}}_[A-Za-z0-9]+_v\d{{2}}r\d{{2}}\.nc$")
    return sorted(f for f in files if rx.fullmatch(f))


def choose_monthly_individual(files: list[str], grid: str, start_ym: int, end_ym: int) -> list[str]:
    rx = re.compile(rf"sic_{grid}_(\d{{6}})_[A-Za-z0-9]+_v\d{{2}}r\d{{2}}\.nc$")
    out: list[str] = []
    for f in files:
        m = rx.fullmatch(f)
        if not m:
            continue
        ym = int(m.group(1))
        if start_ym <= ym <= end_ym:
            out.append(f)
    return sorted(out)


def build_jobs(
    base_url: str,
    version: str,
    dest_root: Path,
    hemis: list[str],
    start_year: int,
    end_year: int,
    daily_mode: str,
    monthly_mode: str,
    include_ancillary: bool,
) -> list[tuple[str, Path]]:
    version_root = f"{base_url.rstrip('/')}/{version}"
    jobs: list[tuple[str, Path]] = []

    if include_ancillary:
        anc_url = f"{version_root}/ancillary/"
        anc_files = list_nc(anc_url)
        for hemi in hemis:
            grid = HEMI_TO_GRID[hemi]
            pats = [
                re.compile(rf"G02202-ancillary-{grid}-v\d{{2}}r\d{{2}}\.nc$"),
                re.compile(rf"G02202-ancillary-{grid}-daily-invalid-ice-v\d{{2}}r\d{{2}}\.nc$"),
            ]
            for name in anc_files:
                if any(p.fullmatch(name) for p in pats):
                    jobs.append((urljoin(anc_url, name), dest_root / version / hemi / "ancillary" / name))

    for hemi in hemis:
        grid = HEMI_TO_GRID[hemi]

        if daily_mode == "aggregate":
            agg_url = f"{version_root}/{hemi}/aggregate/"
            agg_files = list_nc(agg_url)
            for year in range(start_year, end_year + 1):
                name = choose_yearly_daily_aggregate(agg_files, grid, year)
                if name is None:
                    print(f"WARNING: no yearly daily aggregate found for {hemi} {year}", file=sys.stderr)
                    continue
                jobs.append((urljoin(agg_url, name), dest_root / version / hemi / "aggregate" / name))

        elif daily_mode == "individual":
            for year in range(start_year, end_year + 1):
                daily_url = f"{version_root}/{hemi}/daily/{year}/"
                try:
                    files = list_nc(daily_url)
                except HTTPError as exc:
                    print(f"WARNING: cannot list {daily_url}: {exc}", file=sys.stderr)
                    continue
                for name in choose_daily_individual(files, grid, year):
                    jobs.append((urljoin(daily_url, name), dest_root / version / hemi / "daily" / name))

        if monthly_mode == "aggregate":
            agg_url = f"{version_root}/{hemi}/aggregate/"
            agg_files = list_nc(agg_url)
            name = choose_period_monthly_aggregate(agg_files, grid)
            if name is None:
                print(f"WARNING: no period monthly aggregate found for {hemi}", file=sys.stderr)
            else:
                jobs.append((urljoin(agg_url, name), dest_root / version / hemi / "aggregate" / name))

        elif monthly_mode == "individual":
            mon_url = f"{version_root}/{hemi}/monthly/"
            files = list_nc(mon_url)
            start_ym = start_year * 100 + 1
            end_ym = end_year * 100 + 12
            for name in choose_monthly_individual(files, grid, start_ym, end_ym):
                jobs.append((urljoin(mon_url, name), dest_root / version / hemi / "monthly" / name))

    dedup: dict[str, tuple[str, Path]] = {}
    for url, path in jobs:
        dedup[str(path)] = (url, path)
    return list(dedup.values())


def download_one(job: tuple[str, Path], min_bytes: int, retries: int) -> tuple[str, str, str, str]:
    url, dest = job
    dest.parent.mkdir(parents=True, exist_ok=True)

    if dest.exists() and dest.stat().st_size >= min_bytes:
        return ("skip", url, str(dest), str(dest.stat().st_size))

    tmp = dest.with_suffix(dest.suffix + ".part")
    last_err: Exception | None = None

    for attempt in range(1, retries + 1):
        try:
            req = Request(url, headers={"User-Agent": USER_AGENT})
            with urlopen(req, timeout=180) as src, open(tmp, "wb") as out:
                while True:
                    chunk = src.read(1024 * 1024)
                    if not chunk:
                        break
                    out.write(chunk)

            size = tmp.stat().st_size
            if size < min_bytes:
                raise IOError(f"download too small: {size} bytes")

            tmp.replace(dest)
            return ("download", url, str(dest), str(size))

        except Exception as exc:
            last_err = exc
            if tmp.exists():
                try:
                    tmp.unlink()
                except OSError:
                    pass
            time.sleep(min(30, 2 * attempt))

    return ("error", url, str(dest), repr(last_err))


def main() -> int:
    p = argparse.ArgumentParser(description="Discover and download available NSIDC G02202 files.")
    p.add_argument("--base-url", default="https://noaadata.apps.nsidc.org/NOAA")
    p.add_argument("--version", default="G02202_V6")
    p.add_argument("--dest-root", required=True, type=Path)
    p.add_argument("--start-year", type=int, required=True)
    p.add_argument("--end-year", type=int, required=True)
    p.add_argument("--hemis", nargs="+", choices=["north", "south"], default=["north", "south"])
    p.add_argument("--daily", choices=["aggregate", "individual", "none"], default="aggregate")
    p.add_argument("--monthly", choices=["aggregate", "individual", "none"], default="aggregate")
    p.add_argument("--ancillary", action="store_true")
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--retries", type=int, default=4)
    p.add_argument("--min-bytes", type=int, default=10000)
    p.add_argument("--manifest-file", type=Path, default=None)
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args()

    jobs = build_jobs(
        base_url=args.base_url,
        version=args.version,
        dest_root=args.dest_root,
        hemis=args.hemis,
        start_year=args.start_year,
        end_year=args.end_year,
        daily_mode=args.daily,
        monthly_mode=args.monthly,
        include_ancillary=args.ancillary,
    )

    jobs = sorted(jobs, key=lambda x: str(x[1]))
    print(f"Planned files: {len(jobs)}")

    if args.manifest_file is not None:
        args.manifest_file.parent.mkdir(parents=True, exist_ok=True)
        with open(args.manifest_file, "w", encoding="utf-8") as f:
            for url, path in jobs:
                f.write(f"{url}\t{path}\n")
        print(f"Wrote manifest: {args.manifest_file}")

    if args.dry_run:
        for url, path in jobs:
            print(f"{url}\t{path}")
        return 0

    ok = 0
    skipped = 0
    failed = 0

    with cf.ThreadPoolExecutor(max_workers=args.workers) as ex:
        futures = [ex.submit(download_one, job, args.min_bytes, args.retries) for job in jobs]
        for fut in cf.as_completed(futures):
            status, url, path, info = fut.result()
            if status == "download":
                ok += 1
                print(f"DOWNLOADED\t{info}\t{path}")
            elif status == "skip":
                skipped += 1
                print(f"SKIPPED\t{info}\t{path}")
            else:
                failed += 1
                print(f"ERROR\t{path}\t{info}", file=sys.stderr)

    print(f"Summary: downloaded={ok}, skipped={skipped}, failed={failed}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
