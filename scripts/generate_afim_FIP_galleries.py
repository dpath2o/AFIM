import os
from pathlib import Path
from collections import defaultdict

GRAPHICAL_DIR = Path("/g/data/gv90/da1339/GRAPHICAL/AFIM")
DOCS_DIR = Path("docs/figures/fip")
HTML_OUT_DIR = Path("docs")
FIP_GROUPS = ["FI_BT_bool", "FI_BT_roll"]
FILENAME = "FIP_1993-1999.png"

def find_all_fip_images():
    fip_index = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
    for sim_dir in GRAPHICAL_DIR.iterdir():
        if not sim_dir.is_dir():
            continue
        for region_dir in sim_dir.iterdir():
            region = region_dir.name
            fip_root = region_dir / "FIP"
            if not fip_root.exists():
                continue
            for ispd_dir in fip_root.iterdir():
                if not ispd_dir.is_dir() or not ispd_dir.name.startswith("ispd_thresh_"):
                    continue
                ispd = ispd_dir.name.replace("ispd_thresh_", "")
                for group in FIP_GROUPS:
                    fip_path = ispd_dir / group / FILENAME
                    if fip_path.exists():
                        fip_index[sim_dir.name][region][ispd][group] = fip_path
    return fip_index

def copy_pngs_to_docs(fip_index):
    for sim, region_dict in fip_index.items():
        for region, ispd_dict in region_dict.items():
            for ispd, group_dict in ispd_dict.items():
                for group, src_path in group_dict.items():
                    dest = DOCS_DIR / sim / region / f"ispd_thresh_{ispd}" / group
                    dest.mkdir(parents=True, exist_ok=True)
                    target = dest / FILENAME
                    if not target.exists():
                        os.system(f"cp {src_path} {target}")

def rel_path_for(sim, region, ispd, group):
    return f"figures/fip/{sim}/{region}/ispd_thresh_{ispd}/{group}/{FILENAME}"

def wrap_details(summary, content, level=0):
    indent = "  " * level
    return f"{indent}<details>\n{indent}<summary>{summary}</summary>\n{content}\n{indent}</details>\n"

def wrap_image_block(src, label):
    return f"<div class='img-block'><img src='{src}' alt='{label}'><br><b>{label}</b></div>"

def write_html(path, title, content):
    html = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="UTF-8">
<title>{title}</title>
<style>
body {{ font-family: sans-serif; padding: 2rem; background: #f9f9f9; }}
summary {{ cursor: pointer; font-weight: bold; }}
details {{ margin-left: 1rem; margin-top: 0.5rem; }}
img {{ max-height: 300px; display: block; margin: 1rem 0; }}
.img-group {{ display: flex; gap: 20px; flex-wrap: wrap; }}
.img-block {{ background: white; border: 1px solid #ccc; padding: 10px; }}
</style></head><body>
<h1>{title}</h1>\n{content}</body></html>"""
    path.write_text(html)

def generate_tree_html(fip_index):
    content = ""
    for sim, region_dict in sorted(fip_index.items()):
        sim_block = ""
        for region, ispd_dict in sorted(region_dict.items()):
            region_block = ""
            for ispd, group_dict in sorted(ispd_dict.items()):
                img_group = "<div class='img-group'>\n"
                for group in FIP_GROUPS:
                    if group in group_dict:
                        rel = rel_path_for(sim, region, ispd, group)
                        img_group += wrap_image_block(rel, group)
                img_group += "</div>"
                region_block += wrap_details(f"ispd_thresh {ispd}", img_group, 3)
            sim_block += wrap_details(region, region_block, 2)
        content += wrap_details(sim, sim_block, 1)
    write_html(HTML_OUT_DIR / "fip_tree_gallery.html", "AFIM Fast Ice Persistence Tree View", content)

def generate_intra_model_html(fip_index):
    content = ""
    for sim, region_dict in sorted(fip_index.items()):
        sim_block = ""
        for region, ispd_dict in sorted(region_dict.items()):
            img_group = "<div class='img-group'>\n"
            for ispd, group_dict in sorted(ispd_dict.items()):
                for group in FIP_GROUPS:
                    if group in group_dict:
                        rel = rel_path_for(sim, region, ispd, group)
                        label = f"{group} — ispd_thresh {ispd}"
                        img_group += wrap_image_block(rel, label)
            img_group += "</div>"
            sim_block += wrap_details(region, img_group, 2)
        content += wrap_details(sim, sim_block, 1)
    write_html(HTML_OUT_DIR / "fip_intra_model_gallery.html", "AFIM Intra-model FIP Gallery", content)

def generate_inter_model_html(fip_index):
    content = ""
    regions = sorted({r for sim in fip_index.values() for r in sim})
    ispds = sorted({i for sim in fip_index.values() for r in sim.values() for i in r})
    for region in regions:
        region_block = ""
        for ispd in ispds:
            img_group = "<div class='img-group'>\n"
            for sim in sorted(fip_index.keys()):
                group_dict = fip_index[sim].get(region, {}).get(ispd, {})
                for group in FIP_GROUPS:
                    if group in group_dict:
                        rel = rel_path_for(sim, region, ispd, group)
                        img_group += wrap_image_block(rel, f"{sim} — {group}")
            img_group += "</div>"
            region_block += wrap_details(f"ispd_thresh {ispd}", img_group, 2)
        content += wrap_details(region, region_block, 1)
    write_html(HTML_OUT_DIR / "fip_inter_model_gallery.html", "AFIM Inter-model FIP Gallery", content)

fip_index = find_all_fip_images()
copy_pngs_to_docs(fip_index)
generate_tree_html(fip_index)
generate_intra_model_html(fip_index)
generate_inter_model_html(fip_index)
print("✅ All FIP galleries generated and copied to docs/")


