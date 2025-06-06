#!/usr/bin/env python3

import re
import csv
import datetime
from pathlib import Path

# --- Paths ---
ARCHIVE_DIR = Path.home() / "AFIM_archive"
OUT_FILE    = ARCHIVE_DIR / "ice_diag_summary.csv"
HTML_PATH   = Path.home() / "AFIM/src/AFIM/docs/ice_diag_summary.html"

# --- Parameters to extract ---
PARAM_KEYS = [
    "dt", "ndtd", "ndte", "kdyn", "revised_evp", "e_yieldcurve", "e_plasticpot",
    "Ktens", "kstrength", "Pstar", "Cstar", "Cf", "visc_method", "kmt_file"
]

# --- Regex patterns ---
PATTERNS = {
    key: re.compile(rf"{key}\s*=\s*(.+?)\s*:") if key != "kmt_file"
    else re.compile(rf"{key}\s*=\s*(.+)$")
    for key in PARAM_KEYS
}

# --- DOI Links ---
DOI_LINKS = {
    "tsamados": "https://doi.org/10.1029/2012JC007990",
    "hibler": "https://doi.org/10.1175/1520-0485(1979)009<0815:ADTSIM>2.0.CO;2",
    "hunke": "https://doi.org/10.1006/jcph.2001.6710",
    "konig": "https://doi.org/10.1175/2009JPO4105.1",
    "rothrock": "https://doi.org/10.1029/JC080i033p04514",
    "kimmritz": "https://doi.org/10.1016/j.ocemod.2017.05.006",
    "bouillon": "https://doi.org/10.1016/j.ocemod.2013.05.013"
}

# --- Parse a single diagnostics file ---
def parse_diag_file(filepath):
    result = {key: "" for key in PARAM_KEYS}
    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        for i, line in enumerate(f):
            if i > 500:
                break
            for key, pattern in PATTERNS.items():
                if result[key] == "":
                    match = pattern.search(line)
                    if match:
                        value = match.group(1).strip()
                        if key == "kmt_file":
                            value = Path(value).name
                        result[key] = value
    return result

# --- Augment row with custom logic ---
def classify_simulation(sim_name, values):
    # rheology
    if sim_name.startswith("eap"):
        rheology = "anisotropic"
        citation = f"<a href='{DOI_LINKS['tsamados']}'>Tsamados et al. (2013)</a>"
    else:
        rheology = "isotropic"
        citation = ""

    # momentum term and additional citation
    momentum, extra_citation = "", ""
    if sim_name.startswith(("Cstar-", "Pstar-")):
        momentum = "internal ice strength"
        extra_citation = f"<a href='{DOI_LINKS['hibler']}'>Hibler (1979)</a>"
    elif sim_name.startswith("elps-"):
        momentum = "viscosities"
        extra_citation = f"<a href='{DOI_LINKS['hunke']}'>Hunke (2001)</a>"
    elif sim_name.startswith("ktens-"):
        momentum = "viscosities"
        extra_citation = f"<a href='{DOI_LINKS['konig']}'>König-Beatty & Holland (2010)</a>"
    elif sim_name.startswith("Roth-"):
        momentum = "internal ice strength"
        extra_citation = f"<a href='{DOI_LINKS['rothrock']}'>Rothrock (1975)</a>"
    elif sim_name == "visc-meth":
        momentum = "viscosities"
        extra_citation = f"<a href='{DOI_LINKS['kimmritz']}'>Kimmritz et al. (2016)</a>"
    elif sim_name == "re-evp-off":
        momentum = ""
        extra_citation = f"<a href='{DOI_LINKS['bouillon']}'>Bouillon et al. (2013)</a>"

    # override visc_method
    if sim_name == "visc-meth":
        values["visc_method"] = "average strength"
    else:
        values["visc_method"] = "average zeta"

    # post-process replacements
    if values["kdyn"] == "1":
        values["kdyn"] = "EVP"
    elif values["kdyn"] == "2":
        values["kdyn"] = "EAP"

    if values["kstrength"] == "0":
        values["kstrength"] = "Hibler"
    elif values["kstrength"] == "1":
        values["kstrength"] = "Rothrock"

    # Combine citations
    final_citation = citation or extra_citation
    return rheology, momentum, final_citation

# --- Main execution ---
def main():
    headers = ["rheology", "momentum term", "sim_name"] + PARAM_KEYS + ["relevant citations"]
    output_lines = [",".join(headers)]
    html_rows = []

    for diag_file in sorted(ARCHIVE_DIR.glob("*/ice_diag.d")):
        sim_name = diag_file.parent.name
        values = parse_diag_file(diag_file)
        rheo, mom, cite = classify_simulation(sim_name, values)
        row = [rheo, mom, sim_name] + [values[key] for key in PARAM_KEYS] + [cite]
        output_lines.append(",".join(row))
        html_rows.append(row)

    # Ensure CSV output directory exists
    OUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    OUT_FILE.write_text("\n".join(output_lines))
    print(f"✅ CSV summary written to: {OUT_FILE}")

    # --- HTML generation ---
    now_str = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    html = []
    html.append("<!DOCTYPE html><html><head><meta charset='UTF-8'>")
    html.append("<title>CICE Diagnostic Parameters</title><style>")
    html.append("""table {border-collapse: collapse; width: 100%; font-family: sans-serif;}
th {border: 1px solid #ccc; padding: 6px 10px; text-align: center; background-color: #f2f2f2; font-size: 1em;}
td {border: 1px solid #ccc; padding: 6px 10px; text-align: center; font-size: 0.85em;}""")
    html.append("</style></head><body>")
    html.append("<h1>CICE Diagnostic Parameters per Simulation</h1>")
    html.append("<p>Parsed from each simulation's <code>ice_diag.d</code> file.</p>")
    html.append(f"<p style='font-size: small;'>Last updated: {now_str}</p>")
    html.append("<table><thead><tr>" + "".join(f"<th>{h}</th>" for h in headers) + "</tr></thead><tbody>")
    for row in html_rows:
        html.append("<tr>" + "".join(f"<td>{cell}</td>" for cell in row) + "</tr>")
    html.append("</tbody></table></body></html>")

    HTML_PATH.parent.mkdir(parents=True, exist_ok=True)
    HTML_PATH.write_text("\n".join(html))
    print(f"✅ HTML table written to: {HTML_PATH}")

if __name__ == "__main__":
    main()

