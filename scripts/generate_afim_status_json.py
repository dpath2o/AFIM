#!/usr/bin/env python3
import os, re, json
from tqdm import tqdm
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd

# === Paths ===
BASE_DIR      = Path("/g/data/gv90/da1339/afim_output")
GRAPHICAL_DIR = Path("/g/data/gv90/da1339/GRAPHICAL/AFIM")
OUTPUT_JSON   = BASE_DIR / "AFIM_archive_status.json"
OUTPUT_HTML   = Path("docs/AFIM_archive_status.html")

# === Globals to share between functions ===
total_files_all = 0
total_size_all = 0.0

# === Helpers ===
def get_model_date_range(sim_dir):
    daily_dir = sim_dir / "history" / "daily"
    if not daily_dir.exists():
        return None, None
    date_pattern = re.compile(r"iceh\.(\d{4}-\d{2}-\d{2})\.nc")
    dates = []
    for f in daily_dir.glob("iceh.*.nc"):
        m = date_pattern.match(f.name)
        if m:
            dates.append(m.group(1))
    return (min(dates), max(dates)) if dates else (None, None)

def get_file_count_and_size(sim_dir):
    total_files = 0
    total_size = 0
    for root, dirs, files in os.walk(sim_dir):
        total_files += len(files)
        for f in files:
            try:
                total_size += os.path.getsize(os.path.join(root, f))
            except OSError:
                continue
    size_gb = round(total_size / 1e9, 2)
    return total_files, size_gb

def scan_simulation(sim_dir):
    all_rows = []
    total_files, size_gb = get_file_count_and_size(sim_dir)
    start, end = get_model_date_range(sim_dir)

    zarr_dir = sim_dir / "zarr"
    if not zarr_dir.exists():
        all_rows.append({
            "sim_name": sim_dir.name,
            "threshold": "—",
            "start": start,
            "end": end,
            "file_count": total_files,
            "dir_size_gb": size_gb,
            "daily": False,
            "rolling": False,
            "FI": {},
            "PI": {},
            "SO": False,
            "FIA": False,
            "PIA": False,
            "SIA": False,
            "FIP": False,
            "FI_BOOL": {}
        })
        return all_rows

    for thresh_dir in zarr_dir.glob("ispd_thresh_*"):
        thresh_str = thresh_dir.name.replace("ispd_thresh_", "")
        result = {
            "sim_name": sim_dir.name,
            "threshold": thresh_str,
            "start": start,
            "end": end,
            "file_count": total_files,
            "dir_size_gb": size_gb,
            "daily": any(thresh_dir.glob("cice_daily_*.zarr")),
            "rolling": any(thresh_dir.glob("cice_rolling_*.zarr")),
            "FI": {g: False for g in ["FI_B", "FI_Ta", "FI_Tx", "FI_BT"]},
            "PI": {g: False for g in ["PI_B", "PI_Ta", "PI_Tx", "PI_BT"]},
            "SO": False,
            "FIA": False,
            "PIA": False,
            "SIA": False,
            "FIP": False,
            "FI_BOOL": {g: False for g in ["FI_B", "FI_Ta", "FI_Tx", "FI_BT"]}
        }

        for zarr_file in thresh_dir.glob("cice_*.zarr"):
            zpath = zarr_file
            for group in ["FI_B", "FI_Ta", "FI_Tx", "FI_BT"]:
                if (zpath / group).exists():
                    result["FI"][group] = True
            for group in ["PI_B", "PI_Ta", "PI_Tx", "PI_BT"]:
                if (zpath / group).exists():
                    result["PI"][group] = True
            if (zpath / "SO").exists():
                result["SO"] = True

        metric_dir = thresh_dir / "ice_metrics"
        if metric_dir.exists():
            for mfile in metric_dir.glob("*.zarr"):
                name = mfile.name
                if name.startswith("FIA"):
                    result["FIA"] = True
                if name.startswith("PIA"):
                    result["PIA"] = True
                if name.startswith("SIA"):
                    result["SIA"] = True
                if name.startswith("FIP"):
                    result["FIP"] = True
                for group in ["FI_B", "FI_Ta", "FI_Tx", "FI_BT"]:
                    if group in name and "_bool" in name:
                        result["FI_BOOL"][group] = True

        all_rows.append(result)

    return all_rows

def build_afim_archive_status(max_workers=8):
    global total_files_all, total_size_all
    sim_dirs = [d for d in BASE_DIR.iterdir() if d.is_dir()]
    status = []

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_sim = {executor.submit(scan_simulation, sim_dir): sim_dir.name for sim_dir in sim_dirs}
        for future in tqdm(as_completed(future_to_sim), total=len(future_to_sim), desc="Scanning simulations"):
            sim_name = future_to_sim[future]
            try:
                result = future.result()
                status.extend(result)
                for row in result:
                    total_files_all += row.get("file_count", 0)
                    total_size_all += row.get("dir_size_gb", 0.0)
                print(f"✅ Finished: {sim_name}")
            except Exception as exc:
                print(f"❌ Failed: {sim_name} — {exc}")

    with open(OUTPUT_JSON, "w") as f:
        json.dump(status, f, indent=2)
    print(f"\n✅ Archive status written to: {OUTPUT_JSON}")
    print(f"📦 Total files: {total_files_all:,}, Total size: {total_size_all:.2f} GB")
    return status

def render_status_html(status):
    global total_files_all, total_size_all
    rows = []
    for d in status:
        base_row = {
            "Simulation": d["sim_name"],
            "Start": d.get("start", "—"),
            "End": d.get("end", "—"),
            "Files": d.get("file_count", "—"),
            "Size (GB)": d.get("dir_size_gb", "—"),
            "Threshold": d.get("threshold", "—"),
            "daily": "🟢" if d.get("daily") else "🔴",
            "rolling": "🟢" if d.get("rolling") else "🔴",
            "FIA": "🟢" if d.get("FIA") else "⚪",
            "PIA": "🟢" if d.get("PIA") else "⚪",
            "SIA": "🟢" if d.get("SIA") else "⚪",
            "FIP": "🟢" if d.get("FIP") else "⚪"
        }
        for fi in ["FI_B", "FI_Ta", "FI_Tx", "FI_BT"]:
            base_row[fi] = "🟢" if d.get("FI", {}).get(fi) else "🔴"
        for pi in ["PI_B", "PI_Ta", "PI_Tx", "PI_BT"]:
            base_row[pi] = "🟢" if d.get("PI", {}).get(pi) else "🔴"
        for fi_bool in ["FI_B", "FI_Ta", "FI_Tx", "FI_BT"]:
            base_row[fi_bool + "_bool"] = "🟢" if d.get("FI_BOOL", {}).get(fi_bool) else "⚪"
        base_row["SO"] = "🟢" if d.get("SO") else "🔴"
        rows.append(base_row)

    df = pd.DataFrame(rows)
    df = df.sort_values(["Simulation", "Threshold"])
    df[["Simulation", "Start", "End", "Files", "Size (GB)"]] = df[["Simulation", "Start", "End", "Files", "Size (GB)"]].mask(df["Simulation"].duplicated(), '')
    df_html = df.to_html(index=False, escape=False)

    GRAPHICAL_DIR.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_HTML, "w") as f:
        f.write("""
        <!DOCTYPE html>
        <html lang=\"en\">
        <head>
        <meta charset=\"UTF-8\">
        <title>AFIM Archive Status</title>
        <style>
        body { font-family: sans-serif; }
        table { border-collapse: collapse; width: 100%; font-size: 12px; }
        th, td { border: 1px solid #ccc; padding: 4px 8px; text-align: center; }
        th { background-color: #f0f0f0; }
        </style>
        </head>
        <body>
        <h2>AFIM Archive Status</h2>
        """)
        f.write(df_html)
        f.write(f"<p><strong>📦 Total files:</strong> {total_files_all:,}, <strong>Total size:</strong> {total_size_all:.2f} GB</p>")
        f.write("</body></html>")

    print(f"✅ HTML status table written to: {OUTPUT_HTML}")

if __name__ == "__main__":
    status = build_afim_archive_status(max_workers=8)
    render_status_html(status)
