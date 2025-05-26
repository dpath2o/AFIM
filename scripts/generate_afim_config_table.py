#!/usr/bin/env python3

import os
import re
import json
from pathlib import Path
from datetime import datetime
import pandas as pd
from glob import glob

# === Configuration ===
BASE_DIR = Path("/g/data/gv90/da1339/afim_output")
OUTPUT_HTML = Path("docs/afim_simulation_configs.html")
ICE_IN_FILE = "ice_in"
HISTORY_DIR = "history/daily"

# === Parameters to Extract from ice_in ===
KEY_PARAMS = [
    "version_name", "grid_file", "kmt_file", "use_bathymetry",
    "grid_ice", "grid_atm", "grid_ocn",
    "ktherm", "kdyn", "ndte", "revised_evp",
    "kstrength", "Pstar", "Cstar", "Cf", "Ktens",
    "e_yieldcurve", "e_plasticpot", "seabed_stress",
    "atm_data_dir", "ocn_data_dir"
]


# === Helpers ===
def parse_ice_in(path):
    """Parse selected parameters from an ice_in namelist file."""
    config = {}
    if not path.exists():
        return config

    pattern = re.compile(r"^\s*(\w+)\s*=\s*(.*?)(!.*)?$")

    with path.open("r") as f:
        for line in f:
            match = pattern.match(line.strip())
            if match:
                key, val, _ = match.groups()
                val = val.strip().strip("'")
                if key in KEY_PARAMS:
                    config[key] = val
    return config


def get_date_range(sim_dir):
    """Get start and end dates from iceh files."""
    hist_path = sim_dir / HISTORY_DIR
    if not hist_path.exists():
        return None, None

    files = sorted(hist_path.glob("iceh.*.nc"))
    if not files:
        return None, None

    def extract_date(f):
        try:
            return datetime.strptime(f.name, "iceh.%Y-%m-%d.nc")
        except ValueError:
            return None

    dates = list(filter(None, map(extract_date, files)))
    if not dates:
        return None, None

    return dates[0].strftime("%Y-%m-%d"), dates[-1].strftime("%Y-%m-%d")


def build_simulation_config_table():
    records = []

    for sim_dir in sorted(BASE_DIR.iterdir()):
        if not sim_dir.is_dir():
            continue

        sim_name = sim_dir.name
        ice_in_path = sim_dir / ICE_IN_FILE
        config = parse_ice_in(ice_in_path)
        start_date, end_date = get_date_range(sim_dir)

        row = {"Simulation": sim_name, "Start": start_date, "End": end_date}
        for key in KEY_PARAMS:
            row[key] = config.get(key, "—")

        records.append(row)

    return pd.DataFrame(records)


def render_html_table(df):
    df = df.sort_values("Simulation")
    df_html = df.to_html(index=False, escape=False)

    OUTPUT_HTML.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_HTML, "w") as f:
        f.write("""
<!DOCTYPE html>
<html lang=\"en\">
<head>
  <meta charset=\"UTF-8\">
  <title>AFIM Simulation Configuration Table</title>
  <style>
    body { font-family: sans-serif; padding: 2rem; background: #f9f9f9; }
    table { border-collapse: collapse; width: 100%; font-size: 12px; }
    th, td { border: 1px solid #ccc; padding: 4px 8px; text-align: left; }
    th { background: #eee; }
  </style>
</head>
<body>
<h2>AFIM Simulation Configuration Table</h2>
<p>Extracted from <code>ice_in</code> files and model history outputs.</p>
""")
        f.write(df_html)
        f.write("\n</body></html>")
    print(f"✅ HTML table written to: {OUTPUT_HTML}")


if __name__ == "__main__":
    df = build_simulation_config_table()
    render_html_table(df)
