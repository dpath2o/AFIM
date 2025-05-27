#!/usr/bin/env python3

import os
from pathlib import Path
from collections import defaultdict

GRAPHICAL_DIR = Path("/g/data/gv90/da1339/GRAPHICAL/AFIM")
DOCS_DIR      = Path("docs/figures/fip")
HTML_OUT_DIR  = Path("docs")

REGIONS       = ["AS", "Aus", "BS", "DML", "EIO", "VOL", "WIO", "WS"]
FIP_GROUPS    = ["FI_BT_bool", "FI_Ta_bool", "FI_BT_roll", "FI_Ta_roll"]
ISPDS         = ["1.0e-3", "5.0e-4"]
FILENAME      = "FIP_1993-1999.png"


def find_all_fip_images():
    """Search the GRAPHICAL directory for all FIP image paths."""
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
                        fip_index[sim_dir.name][region][ispd][group] = fip_path

    return fip_index


def copy_pngs_to_docs(fip_index):
    """Copy all FIP PNGs to docs/figures/fip/[sim]/[region]/..."""
    for sim, region_dict in fip_index.items():
        for region, ispd_dict in region_dict.items():
            for ispd, group_dict in ispd_dict.items():
                for group, src_path in group_dict.items():
                    dest = DOCS_DIR / sim / region / f"ispd_thresh_{ispd}" / group
                    dest.mkdir(parents=True, exist_ok=True)
                    target = dest / FILENAME
                    if not target.exists():
                        os.system(f"cp {src_path} {target}")


def generate_intra_model_html(fip_index):
    """Each page: one sim, one region, one ispd_thresh -> multiple groups"""
    html_path = HTML_OUT_DIR / "fip_intra_model_gallery.html"
    with open(html_path, "w") as f:
        f.write("""
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>AFIM Intra-model FIP Gallery</title>
  <style>
    body { font-family: sans-serif; padding: 2rem; background: #f9f9f9; }
    h2 { margin-top: 2rem; color: #333; }
    .gallery { display: flex; flex-wrap: wrap; gap: 10px; }
    .gallery div { background: white; border: 1px solid #ccc; padding: 10px; }
    .gallery img { max-height: 300px; }
  </style>
</head>
<body>
<h1>AFIM Intra-model FIP Gallery</h1>
""")
        for sim, region_dict in fip_index.items():
            for region, ispd_dict in region_dict.items():
                for ispd, group_dict in ispd_dict.items():
                    f.write(f"<h2>{sim} — {region} — ispd_thresh {ispd}</h2>\n")
                    f.write("<div class='gallery'>\n")
                    for group, path in group_dict.items():
                        rel_path = os.path.relpath(DOCS_DIR / sim / region / f"ispd_thresh_{ispd}" / group / FILENAME, HTML_OUT_DIR)
                        f.write(f"  <div><img src='{rel_path}' alt='{group}'><br>{group}</div>\n")
                    f.write("</div>\n")

        f.write("</body></html>")


def generate_inter_model_html(fip_index):
    """Each page: one region, one ispd_thresh, one group -> multiple sims"""
    html_path = HTML_OUT_DIR / "fip_inter_model_gallery.html"
    with open(html_path, "w") as f:
        f.write("""
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>AFIM Inter-model FIP Gallery</title>
  <style>
    body { font-family: sans-serif; padding: 2rem; background: #f9f9f9; }
    h2 { margin-top: 2rem; color: #333; }
    .gallery { display: flex; flex-wrap: wrap; gap: 10px; }
    .gallery div { background: white; border: 1px solid #ccc; padding: 10px; }
    .gallery img { max-height: 300px; }
  </style>
</head>
<body>
<h1>AFIM Inter-model FIP Gallery</h1>
""")
        for region in REGIONS:
            for ispd in ISPDS:
                for group in FIP_GROUPS:
                    f.write(f"<h2>{region} — ispd_thresh {ispd} — {group}</h2>\n")
                    f.write("<div class='gallery'>\n")
                    for sim in fip_index:
                        path = fip_index[sim].get(region, {}).get(ispd, {}).get(group)
                        if path:
                            rel_path = os.path.relpath(DOCS_DIR / sim / region / f"ispd_thresh_{ispd}" / group / FILENAME, HTML_OUT_DIR)
                            f.write(f"  <div><img src='{rel_path}' alt='{sim}'><br>{sim}</div>\n")
                    f.write("</div>\n")

        f.write("</body></html>")

if __name__ == "__main__":
    fip_index = find_all_fip_images()
    print(f"Found {len(fip_index)} simulations with FIP plots.")
    for sim in fip_index:
        print(f" - {sim}: regions={list(fip_index[sim].keys())}")
    copy_pngs_to_docs(fip_index)
    generate_intra_model_html(fip_index)
    generate_inter_model_html(fip_index)
    print("✅ FIP galleries generated and copied to docs/")

