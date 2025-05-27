#!/usr/bin/env python3

from pathlib import Path
from collections import defaultdict
import os

GRAPHICAL_DIR = Path("/g/data/gv90/da1339/GRAPHICAL/AFIM")
DOCS_DIR = Path("docs/figures/fip")
HTML_OUT_DIR = Path("docs")
REGIONS = ["AS", "Aus", "BS", "DML", "EIO", "VOL", "WIO", "WS"]
FIP_GROUPS = ["FI_BT_bool", "FI_Ta_bool", "FI_BT_roll", "FI_Ta_roll"]
ISPDS = ["1.0e-3", "5.0e-4"]
FILENAME = "FIP_1993-1999.png"

def find_all_fip_images():
    fip_index = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
    for sim_dir in GRAPHICAL_DIR.iterdir():
        if not sim_dir.is_dir():
            continue
        for region in REGIONS:
            region_dir = sim_dir / region / "FIP"
            if not region_dir.exists():
                continue
            for ispd in ISPDS:
                ispd_dir = region_dir / f"ispd_thresh_{ispd}"
                if not ispd_dir.exists():
                    continue
                for group in FIP_GROUPS:
                    fip_path = ispd_dir / group / FILENAME
                    if fip_path.exists():
                        fip_index[sim_dir.name][region][group][ispd] = fip_path
    return fip_index

def generate_tree_html(fip_index):
    html_path = HTML_OUT_DIR / "fip_tree_gallery.html"
    HTML_OUT_DIR.mkdir(parents=True, exist_ok=True)
    with open(html_path, "w") as f:
        f.write("""
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>AFIM FIP Tree Gallery</title>
  <style>
    body { font-family: sans-serif; padding: 2rem; background: #f9f9f9; }
    summary { cursor: pointer; font-weight: bold; }
    details { margin-left: 1rem; margin-top: 0.5rem; }
    img { max-height: 300px; display: block; margin: 1rem 0; }
    .img-group { display: flex; gap: 20px; flex-wrap: wrap; }
    .img-block { background: white; border: 1px solid #ccc; padding: 10px; }
  </style>
</head>
<body>
<h1>AFIM Fast Ice Persistence Gallery (Tree View)</h1>
""")
        for sim, region_dict in fip_index.items():
            f.write(f"<details>\n<summary>{sim}</summary>\n")
            for region, group_dict in region_dict.items():
                f.write(f"<details>\n<summary>{region}</summary>\n")
                for group, ispd_dict in group_dict.items():
                    f.write(f"<details>\n<summary>{group}</summary>\n<div class='img-group'>\n")
                    for ispd, path in ispd_dict.items():
                        rel_path = os.path.relpath(DOCS_DIR / sim / region / f"ispd_thresh_{ispd}" / group / FILENAME, HTML_OUT_DIR)
                        f.write(f"<div class='img-block'><img src='{rel_path}' alt='{sim} {region} {group} {ispd}'><br><b>{ispd}</b></div>\n")
                    f.write("</div>\n</details>\n")
                f.write("</details>\n")
            f.write("</details>\n")
        f.write("</body></html>")

if __name__ == "__main__":
    fip_index = find_all_fip_images()
    generate_tree_html(fip_index)
    print(f"✅ Generated HTML tree at: {HTML_OUT_DIR/'fip_tree_gallery.html'}")

