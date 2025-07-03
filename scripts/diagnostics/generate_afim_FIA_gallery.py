#!/usr/bin/env python3

from pathlib import Path
import shutil
import re
from collections import defaultdict
from datetime import datetime

# === CONFIGURATION ===
source_dir = Path("/g/data/gv90/da1339/GRAPHICAL/AFIM/timeseries/")
dest_dir = Path("/home/581/da1339/AFIM/src/AFIM/docs/figures/timeseries/")
html_output = Path("/home/581/da1339/AFIM/src/AFIM/docs/timeseries_gallery.html")
graphical_root = Path("/home/581/da1339/graphical/AFIM/")

# === SETUP ===
dest_dir.mkdir(parents=True, exist_ok=True)

# === COPY FILES ===
for file in source_dir.glob("FIA_*.png"):
    shutil.copy2(file, dest_dir / file.name)

# === PARSE FILENAMES ===
images_by_sim = defaultdict(lambda: defaultdict(list))
pattern = re.compile(r"FIA_([^-_]+)-?([^-_]+)?_([0-9.e+-]+)_smoothed_1993-1999\\.png")

for file in dest_dir.glob("FIA_*.png"):
    match = pattern.match(file.name)
    if match:
        sim1, sim2, thresh = match.groups()
        sim_name = sim1 if sim2 is None else f"{sim1}-{sim2}"
        images_by_sim[sim_name][thresh].append(file.name)

# === COPY ispd-thresh_vs_FIA-min-max.png IMAGES PER SIM ===
for sim_name in images_by_sim.keys():
    source_png = graphical_root / sim_name / "ispd-thresh_vs_FIA-min-max.png"
    if source_png.exists():
        dest_png_name = f"{sim_name}_ispd-thresh_vs_FIA-min-max.png"
        shutil.copy2(source_png, dest_dir / dest_png_name)
        images_by_sim[sim_name]["ispd-thresh-summary"] = [dest_png_name]

# === HTML HEADER ===
html = f"""<!DOCTYPE html>
<html lang=\"en\">
<head>
  <meta charset=\"UTF-8\">
  <title>AFIM Fast Ice Area Timeseries (1993–1999)</title>
  <style>
    body {{
        font-family: sans-serif;
        background: #f9f9f9;
        margin: 0 auto;
        padding: 1rem 2rem;
        max-width: 1400px;
    }}
    h1 {{
        text-align: center;
        color: #333;
    }}
    details {{
        margin: 1rem 0;
        border: 1px solid #ccc;
        padding: 0.5rem;
        background: white;
        box-shadow: 1px 1px 4px rgba(0,0,0,0.1);
    }}
    summary {{
        font-weight: bold;
        cursor: pointer;
    }}
    .plot {{
        text-align: center;
        margin-bottom: 1rem;
    }}
    .plot img {{
        max-width: 100%;
        height: auto;
    }}
    .plot p {{
        font-size: 14px;
        color: #555;
        margin-top: 0.3em;
    }}
    footer {{
        text-align: center;
        font-size: 13px;
        margin-top: 2rem;
        color: #888;
    }}
  </style>
</head>
<body>

<h1>AFIM Fast Ice Area Timeseries (1993–1999)</h1>
<p style=\"text-align:center; font-size: 1rem;\">All simulations — 15-day smoothed FIA from CICE model output</p>
"""

# === HTML BODY ===
for sim_name in sorted(images_by_sim.keys()):
    html += f"<details>\n<summary>{sim_name}</summary>\n"

    # Insert ispd-thresh summary plot at top of SIM_NAME section
    if "ispd-thresh-summary" in images_by_sim[sim_name]:
        filename = images_by_sim[sim_name]["ispd-thresh-summary"][0]
        html += f"""<div class=\"plot\">
  <img src=\"figures/timeseries/{filename}\" alt=\"{sim_name} - FIA vs Ice Speed Threshold\">
  <p>{sim_name} — FIA min/max vs Ice Speed Threshold</p>
</div>\n"""

    # Only plot numeric threshold sections
    numeric_thresh = [t for t in images_by_sim[sim_name].keys() if re.match(r'^[0-9.eE+-]+$', t)]
    numeric_thresh = sorted(numeric_thresh, key=lambda x: float(x))

    for thresh in numeric_thresh:
        html += f"<details style='margin-left: 1.5rem;'>\n<summary>Threshold: {thresh}</summary>\n"
        for filename in sorted(images_by_sim[sim_name][thresh]):
            label = f"{sim_name} — {thresh}"
            html += f"""<div class=\"plot\">
  <img src=\"figures/timeseries/{filename}\" alt=\"{label}\">
  <p>{label}</p>
</div>\n"""
        html += "</details>\n"
    html += "</details>\n"

# === HTML FOOTER ===
html += f"""
<footer>
  Last updated: {datetime.now().strftime("%B %d, %Y")} — AFIM archive visualisation by <a href=\"https://github.com/dpath2o/AFIM\">AFIM</a>
</footer>
</body>
</html>
"""

# === WRITE HTML FILE ===
html_output.parent.mkdir(parents=True, exist_ok=True)
html_output.write_text(html)

print(f"✅ HTML written to {html_output}")
