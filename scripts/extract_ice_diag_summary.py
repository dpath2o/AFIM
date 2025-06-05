#!/usr/bin/env python3

import re
import csv
from pathlib import Path

# Directory containing the simulations
ARCHIVE_DIR = Path.home() / "AFIM_archive"
OUT_FILE    = ARCHIVE_DIR / "ice_diag_summary.csv"

# Parameters you want to extract
PARAM_KEYS = [
    "kdyn", "revised_evp", "ndtd", "ndte", "e_yieldcurve", "e_plasticpot",
    "Ktens", "kstrength", "Pstar", "Cstar", "dt", "kmt_file"
]

# Regex pattern to match key and value
PATTERNS = {key: re.compile(rf"{key}\s*=\s*(.+?)\s*:") for key in PARAM_KEYS}

# Helper to parse a single file
def parse_diag_file(filepath):
    result = {key: "" for key in PARAM_KEYS}
    with open(filepath, "r", encoding="utf-8", errors="replace") as f:    
        for i, line in enumerate(f):
            if i > 500:  # Stop after first 500 lines
                break
            for key, pattern in PATTERNS.items():
                if key in result and result[key] == "":
                    match = pattern.search(line)
                    if match:
                        result[key] = match.group(1).strip()
    return result

# Main loop over subdirectories
def main():
    output_lines = []
    headers = ["sim_name"] + PARAM_KEYS
    output_lines.append(",".join(headers))

    for diag_file in sorted(ARCHIVE_DIR.glob("*/ice_diag.d")):
        sim_name = diag_file.parent.name
        values = parse_diag_file(diag_file)
        row = [sim_name] + [values[key] for key in PARAM_KEYS]
        output_lines.append(",".join(row))

    # Write output CSV
    OUT_FILE.write_text("\n".join(output_lines))
    print(f"✅ Summary written to {OUT_FILE}")

    # Paths
    csv_path = Path.home() / "AFIM_archive/ice_diag_summary.csv"
    html_path = Path.home() / "AFIM/src/AFIM/docs/ice_diag_summary.html"

    # Load CSV
    with open(csv_path, newline="") as f:
        reader = csv.reader(f)
        rows = list(reader)

    header = rows[0]
    data = rows[1:]

    # Build HTML table
    html_lines = []
    html_lines.append("<!DOCTYPE html>")
    html_lines.append("<html><head><meta charset='UTF-8'>")
    html_lines.append("<title>CICE Diagnostic Parameters</title>")
    html_lines.append("<style>")
    html_lines.append("""
      table {border-collapse: collapse; width: 100%; font-family: sans-serif;}
      th, td {border: 1px solid #ccc; padding: 6px 10px; text-align: center;}
      th {background-color: #f2f2f2;}
    """)
    html_lines.append("</style></head><body>")
    html_lines.append("<h1>CICE Diagnostic Parameters per Simulation</h1>")
    html_lines.append("<p>Parsed from each simulation's <code>ice_diag.d</code> file.</p>")

    html_lines.append("<table>")
    html_lines.append("<thead><tr>" + "".join([f"<th>{h}</th>" for h in header]) + "</tr></thead>")
    html_lines.append("<tbody>")
    for row in data:
        html_lines.append("<tr>" + "".join([f"<td>{cell}</td>" for cell in row]) + "</tr>")
    html_lines.append("</tbody></table>")
    html_lines.append("</body></html>")

    # Write to file
    html_path.write_text("\n".join(html_lines))
    print(f"✅ HTML table written to: {html_path}")

if __name__ == "__main__":
    main()
