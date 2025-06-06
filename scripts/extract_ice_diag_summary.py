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
    "dt", "ndtd", "ndte", "kdyn", "revised_evp", "e_yieldcurve", "e_plasticpot", "Ktens", "kstrength", "Pstar", "Cstar", "Cf", "visc_method", "kmt_file"
]

# --- Regex patterns ---
PATTERNS = {
    key: re.compile(rf"{key}\s*=\s*(.+?)\s*:") if key != "kmt_file"
    else re.compile(rf"{key}\s*=\s*(.+)$")
    for key in PARAM_KEYS
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

# --- Main execution ---
def main():
    output_lines = []
    headers = ["sim_name"] + PARAM_KEYS
    output_lines.append(",".join(headers))

    for diag_file in sorted(ARCHIVE_DIR.glob("*/ice_diag.d")):
        sim_name = diag_file.parent.name
        values = parse_diag_file(diag_file)
        row = [sim_name] + [values[key] for key in PARAM_KEYS]
        output_lines.append(",".join(row))

    # Write CSV summary
    OUT_FILE.write_text("\n".join(output_lines))
    print(f"✅ CSV summary written to: {OUT_FILE}")

    # --- Generate HTML table from CSV ---
    with open(OUT_FILE, newline="") as f:
        reader = csv.reader(f)
        rows = list(reader)

    header = rows[0]
    data = sorted(rows[1:], key=lambda row: row[0].lower())

    html_lines = []
    html_lines.append("<!DOCTYPE html>")
    html_lines.append("<html><head><meta charset='UTF-8'>")
    html_lines.append("<title>CICE Diagnostic Parameters</title>")
    html_lines.append("<style>")
    html_lines.append(""" table {border-collapse: collapse;
                                 width: 100%;
                                 font-family: sans-serif;}
                             th {border: 1px solid #ccc;
                                 padding: 6px 10px;
                                 text-align: center;
                                 background-color: #f2f2f2;
                                 font-size: 1em;}
                             td {border: 1px solid #ccc;
                                 padding: 6px 10px;
                                 text-align: center;
                                 font-size: 0.85em;}""")
    html_lines.append("</style></head><body>")
    html_lines.append("<h1>CICE Diagnostic Parameters per Simulation</h1>")
    html_lines.append("<p>Parsed from each simulation's <code>ice_diag.d</code> file.</p>")

    # Optional timestamp
    now_str = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    html_lines.append(f"<p style='font-size: small;'>Last updated: {now_str}</p>")

    html_lines.append("<table>")
    html_lines.append("<thead><tr>" + "".join([f"<th>{h}</th>" for h in header]) + "</tr></thead>")
    html_lines.append("<tbody>")
    for row in data:
        html_lines.append("<tr>" + "".join([f"<td>{cell}</td>" for cell in row]) + "</tr>")
    html_lines.append("</tbody></table>")
    html_lines.append("</body></html>")

    # Ensure output directory exists
    HTML_PATH.parent.mkdir(parents=True, exist_ok=True)
    HTML_PATH.write_text("\n".join(html_lines))
    print(f"✅ HTML table written to: {HTML_PATH}")

if __name__ == "__main__":
    main()
