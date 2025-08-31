#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Mirror AWI ESA-CCI tree via FTP (plain FTP, anonymous), verbose + resumable.

Mirrors everything UNDER REMOTE_BASE into LOCAL_ROOT preserving the relative
structure. Example:

  REMOTE_BASE = /sea_ice/projects/cci/crdp/v4p0
  REMOTE_ROOT = /sea_ice/projects/cci/crdp/v4p0

Remote: /sea_ice/projects/cci/crdp/v4p0/l2p_release/sh/cryosat2/2017/file.nc
Local :  ~/seaice/                               l2p_release/sh/cryosat2/2017/file.nc

You can also start deeper (REMOTE_ROOT inside REMOTE_BASE); the relative path
from REMOTE_BASE is still preserved.

Environment/flags:
  --remote-root REMOTE_ROOT        (env REMOTE_ROOT; default = REMOTE_BASE)
  --remote-base REMOTE_BASE        (env REMOTE_BASE; default = /sea_ice/projects/cci/crdp/v4p0)
  --local-root  LOCAL_ROOT         (env LOCAL_ROOT; default = ~/seaice)
  --year-min    YEAR_MIN           (env YEAR_MIN)
  --year-max    YEAR_MAX           (env YEAR_MAX)
  --include-glob '*.nc'            (repeatable; env INCLUDE_GLOBS comma-separated)
  --exclude-glob '*quicklook*'     (repeatable; env EXCLUDE_GLOBS comma-separated)
  --dry-run                        (env DRY_RUN=1)
  --quiet                          (env QUIET=1)
  --user USER  --pass PASS         (env FTP_USER / FTP_PASS)
  PASSIVE_MODE=1|0 (env)           start passive; set 0 to try active first
"""

from __future__ import annotations
import os, re, sys, time, fnmatch, socket, argparse, traceback
import ftplib
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import Optional, List

# ---------- defaults ----------
FTP_HOST     = "ftp.awi.de"
REMOTE_BASE_DEF = "/sea_ice/projects/cci/crdp/v4p0"
# If REMOTE_ROOT not provided, we default it to REMOTE_BASE
LOCAL_DEF    = os.path.expanduser("~/seaice/AWI")

MAX_RETRIES  = 5
TIMEOUT_SEC  = 180
PASSIVE_MODE = os.environ.get("PASSIVE_MODE", "1").lower() not in ("0", "false", "no")
# ------------------------------

@dataclass
class FtpEntry:
    name: str
    type: str                 # 'file' | 'dir'
    size: Optional[int] = None
    mtime: Optional[int] = None

def log(msg: str, quiet: bool = False) -> None:
    if quiet: return
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)

def connect(host: str, user: str, pwd: str, passive: bool, timeout: int, quiet: bool) -> ftplib.FTP:
    ftp = ftplib.FTP(host, timeout=timeout)  # plain FTP (anonymous allows no TLS)
    if not user: user = "anonymous"
    if not pwd:  pwd  = "anonymous@"
    ftp.login(user, pwd)
    ftp.set_pasv(passive)
    try: ftp.sock.settimeout(timeout)
    except Exception: pass
    log(f"Connected to {host} (TLS=False), user={user}", quiet)
    return ftp

def list_dir_abs(ftp: ftplib.FTP, path: str, timeout: int, quiet: bool) -> List[FtpEntry]:
    """Absolute LIST with passive→active one-time fallback, parsed minimally."""
    def _do_list() -> list[str]:
        lines: list[str] = []
        ftp.retrlines(f"LIST {path}", lines.append)
        return lines
    def _parse(lines: list[str]) -> List[FtpEntry]:
        out: List[FtpEntry] = []
        for line in lines:
            parts = line.split(maxsplit=8)
            if not parts: continue
            name = parts[-1]
            if name in (".", ".."): continue
            ftype = "dir" if parts[0].startswith("d") else "file"
            size = int(parts[4]) if (ftype == "file" and len(parts) >= 5 and parts[4].isdigit()) else None
            out.append(FtpEntry(name=name, type=ftype, size=size))
        return out
    for attempt in (1, 2):
        try:
            try: ftp.sock.settimeout(timeout)
            except Exception: pass
            return _parse(_do_list())
        except (ftplib.error_temp, ftplib.error_reply, ftplib.error_proto, socket.timeout, OSError) as e:
            log(f"[WARN] LIST {path} failed (attempt {attempt}/2): {e}", quiet)
            if attempt == 1:
                ftp.set_pasv(not PASSIVE_MODE)
                mode = "passive" if not PASSIVE_MODE else "active"
                log(f"Retrying LIST {path} in {mode} mode", quiet)
                time.sleep(1.5)
            else:
                raise

def want_file(remote_rel: str, include_globs: list[str], exclude_globs: list[str],
              year_min: Optional[int], year_max: Optional[int]) -> bool:
    base = os.path.basename(remote_rel.rstrip("/"))
    if include_globs and not any(fnmatch.fnmatch(base, g) for g in include_globs): return False
    if any(fnmatch.fnmatch(base, g) for g in exclude_globs): return False
    if year_min is None and year_max is None: return True
    year = None
    m = re.search(r"(\d{4})(\d{2})?(\d{2})?", base) or re.search(r"/(20\d{2})(/|$)", remote_rel)
    if m: year = int(m.group(1))
    if year is None: return True
    if year_min is not None and year < year_min: return False
    if year_max is not None and year > year_max: return False
    return True

def local_matches(local_path: Path, rsize: Optional[int]) -> bool:
    if not local_path.exists(): return False
    st = local_path.stat()
    return (rsize is None) or (st.st_size == rsize)

def download_file(ftp: ftplib.FTP, rdir: str, ent: FtpEntry, ldir: Path,
                  dry_run: bool, quiet: bool) -> None:
    """Absolute RETR (no cwd), resume to .part then rename."""
    ldir.mkdir(parents=True, exist_ok=True)
    lpath = ldir / ent.name
    ppath = ldir / (ent.name + ".part")
    rpath = f"{rdir}/{ent.name}"

    if local_matches(lpath, ent.size):
        log(f"SKIP (up-to-date): {lpath}", quiet); return
    if dry_run:
        log(f"[DRY] {rpath} → {lpath}", quiet); return

    resume_from = ppath.stat().st_size if ppath.exists() else 0
    if resume_from: log(f"RESUME {ppath.name} from {resume_from} bytes", quiet)
    bytes_done = resume_from; last = time.time()

    def _cb(chunk: bytes):
        nonlocal bytes_done, last
        with open(ppath, "ab") as f:  # safe for HPC NFS to reopen
            f.write(chunk)
        bytes_done += len(chunk)
        now = time.time()
        if not quiet and (now - last) > 1.5:
            if ent.size:
                pct = 100.0 * bytes_done / ent.size
                log(f"DOWN {ent.name}: {bytes_done}/{ent.size} ({pct:4.1f}%)", quiet)
            else:
                log(f"DOWN {ent.name}: {bytes_done} bytes", quiet)
            last = now

    if resume_from:
        try: ftp.sendcmd(f"REST {resume_from}")
        except Exception: pass

    for attempt in range(1, MAX_RETRIES + 1):
        try:
            ftp.retrbinary(f"RETR {rpath}", _cb, blocksize=1024 * 128)
            break
        except (ftplib.error_temp, ftplib.error_reply, ftplib.error_proto, socket.timeout, OSError) as e:
            log(f"[WARN] RETR fail {attempt}/{MAX_RETRIES} for {ent.name}: {e}", quiet)
            time.sleep(1.5 * attempt)
            if attempt == MAX_RETRIES: raise

    Path(ppath).replace(lpath)
    log(f"DONE {lpath}", quiet)

def local_path_for(rdir_abs: str, remote_base_abs: str, local_root: Path) -> Path:
    """Map absolute remote path to local path relative to remote_base_abs."""
    r = PurePosixPath(rdir_abs.rstrip("/"))
    b = PurePosixPath(remote_base_abs.rstrip("/"))
    try:
        rel = r.relative_to(b)           # e.g., l2p_release/sh/cryosat2/2017
    except ValueError:
        rel = PurePosixPath(r.name)      # fallback: just the last component
    return local_root / Path(str(rel))

def walk_and_mirror(ftp: ftplib.FTP, remote_root: str, remote_base: str, local_root: Path,
                    include_globs: list[str], exclude_globs: list[str],
                    year_min: Optional[int], year_max: Optional[int],
                    dry_run: bool, quiet: bool) -> None:
    remote_root = remote_root.rstrip("/")
    remote_base = (remote_base or remote_root).rstrip("/")

    start_local = local_path_for(remote_root, remote_base, local_root)
    start_local.mkdir(parents=True, exist_ok=True)
    log(f"Mirror start: {remote_root} (base={remote_base}) → {start_local}", quiet)

    stack: list[tuple[str, Path]] = [(remote_root, start_local)]
    while stack:
        rdir, ldir = stack.pop()
        try:
            entries = list_dir_abs(ftp, rdir, TIMEOUT_SEC, quiet)
        except Exception as e:
            log(f"[WARN] cannot list {rdir}: {e}", quiet)
            continue

        # queue subdirs
        for ent in entries:
            if ent.type == "dir":
                rsub = f"{rdir}/{ent.name}"
                lsub = local_path_for(rsub, remote_base, local_root)
                stack.append((rsub, lsub))

        # files
        for ent in entries:
            if ent.type != "file":
                continue
            rrel = f"{rdir}/{ent.name}"
            if not want_file(rrel, include_globs, exclude_globs, year_min, year_max):
                log(f"SKIP (filter): {rrel}", quiet)
                continue
            try:
                download_file(ftp, rdir, ent, ldir, dry_run, quiet)
            except Exception as e:
                log(f"[ERROR] {rrel}: {e}", quiet)

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Mirror AWI ESA-CCI (plain FTP, anonymous).")
    p.add_argument("--remote-base", default=os.environ.get("REMOTE_BASE", REMOTE_BASE_DEF))
    p.add_argument("--remote-root", default=os.environ.get("REMOTE_ROOT"))  # default to base in main()
    p.add_argument("--local-root",  default=os.environ.get("LOCAL_ROOT",  LOCAL_DEF))
    p.add_argument("--year-min", type=int, default=int(os.environ.get("YEAR_MIN")) if os.environ.get("YEAR_MIN") else None)
    p.add_argument("--year-max", type=int, default=int(os.environ.get("YEAR_MAX")) if os.environ.get("YEAR_MAX") else None)
    p.add_argument("--include-glob", action="append", default=(os.environ.get("INCLUDE_GLOBS","*.nc").split(",")))
    p.add_argument("--exclude-glob", action="append", default=(os.environ.get("EXCLUDE_GLOBS","").split(",") if os.environ.get("EXCLUDE_GLOBS") else []))
    p.add_argument("--dry-run", action="store_true", default=os.environ.get("DRY_RUN","").lower() in ("1","true","yes"))
    p.add_argument("--quiet",   action="store_true", default=os.environ.get("QUIET","").lower()   in ("1","true","yes"))
    p.add_argument("--user", default=os.environ.get("FTP_USER","anonymous"))
    p.add_argument("--pass", dest="pwd", default=os.environ.get("FTP_PASS","anonymous@"))
    return p.parse_args()

def flatten_globs(globs: List[List[str] | str]) -> list[str]:
    out: list[str] = []
    for item in globs or []:
        if isinstance(item, list):
            out.extend([g for g in item if g])
        elif item:
            out.append(item)
    # remove empties / whitespace
    return [g.strip() for g in out if g and g.strip()]

def main():
    ns = parse_args()
    try:
        remote_base = ns.remote_base.rstrip("/")
        remote_root = (ns.remote_root or ns.remote_base).rstrip("/")
        local_root  = Path(ns.local_root)
        local_root.mkdir(parents=True, exist_ok=True)

        include_globs = flatten_globs(ns.include_glob) or ["*.nc"]
        exclude_globs = flatten_globs(ns.exclude_glob)

        ftp = connect(FTP_HOST, ns.user, ns.pwd, PASSIVE_MODE, TIMEOUT_SEC, ns.quiet)

        # Early absolute LIST sanity check
        tmp: list[str] = []
        ftp.retrlines(f"LIST {remote_root}", tmp.append)
        log(f"Listing OK for {remote_root} ({len(tmp)} entries)", ns.quiet)

        walk_and_mirror(
            ftp,
            remote_root,
            remote_base,
            local_root,
            include_globs=include_globs,
            exclude_globs=exclude_globs,
            year_min=ns.year_min,
            year_max=ns.year_max,
            dry_run=ns.dry_run,
            quiet=ns.quiet,
        )

        log("All done.", ns.quiet)
        try: ftp.quit()
        except Exception: ftp.close()

    except SystemExit:
        raise
    except Exception as e:
        # Print full traceback to help debug if anything goes wrong
        print("".join(traceback.format_exception(e)), file=sys.stderr, flush=True)
        sys.exit(1)

if __name__ == "__main__":
    main()
