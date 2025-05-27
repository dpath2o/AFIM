import os
import shutil
from pathlib import Path
from collections import defaultdict
import re
from datetime import datetime

# === CONFIGURATION ===
source_dir = Path("/g/data/gv90/da1339/GRAPHICAL/AFIM/timeseries/")
dest_dir = Path("/home/581/da1339/AFIM/src/AFIM/docs/figures/timeseries/")
html_output = Path("/home/581/da1339/AFIM/src/AFIM/docs/timeseries_gallery.html")

# === SETUP ===
dest_dir.mkdir(parents=True, exist_ok=True)

# === COPY FILES ===
for file in source_dir.glob("FIA_*.png"):
    shutil.copy2(file, dest_dir / file.name)

# === PARSE FILENAMES ===
images_by_sim = defaultdict(lambda: defaultdict(list))
pattern = re.compile(r"FIA_(.+)_(1\.0e-3|5\.0e-4)_(smoothed|raw)_1993-1999\.png")

for file in dest_dir.glob("FIA_*.png"):
    match = pattern.match(file.name)
    if match:
        sim, thresh, smooth = match.groups()
        images_by_sim[sim][thresh].append((smooth, file.name))

# === HTML HEADER ===
html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
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
    .image-columns {{
        display: flex;
        justify-content: space-between;
        flex-wrap: wrap;
        margin-top: 1rem;
    }}
    .column {{
        width: 48%;
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
<p style="text-align:center; font-size: 1rem;">All simulations — 15-day smoothed FIA from CICE model output</p>
"""

# === HTML BODY ===
for sim_name in sorted(images_by_sim.keys()):
    html += f"<details>\n<summary>{sim_name}</summary>\n<div class='image-columns'>\n"

    for thresh in ["1.0e-3", "5.0e-4"]:
        html += "<div class='column'>\n"
        for smooth, filename in sorted(images_by_sim[sim_name].get(thresh, [])):
            label = f"{sim_name} — {thresh} — {smooth}"
            html += f"""<div class="plot">
  <img src="figures/timeseries/{filename}" alt="{label}">
  <p>{label}</p>
</div>\n"""
        html += "</div>\n"

    html += "</div>\n</details>\n"

# === HTML FOOTER ===
html += f"""
<footer>
  Last updated: {datetime.now().strftime("%B %d, %Y")} — AFIM archive visualisation by <a href="https://github.com/dpath2o/AFIM">AFIM</a>
</footer>
</body>
</html>
"""

# === WRITE HTML FILE ===
html_output.parent.mkdir(parents=True, exist_ok=True)
html_output.write_text(html)

print(f"✅ HTML written to {html_output}")
