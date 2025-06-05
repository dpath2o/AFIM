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

    # === CONFIG ===
    csv_path = Path.home() / "AFIM_archive/ice_diag_summary.csv"
    md_path = Path.home() / "AFIM/src/AFIM/docs/ice_diag_summary.md"
    # === LOAD CSV AND CONVERT TO MARKDOWN TABLE ===
    with open(csv_path, newline="") as f:
        reader = csv.reader(f)
        rows = list(reader)
    header = rows[0]
    data = rows[1:]

    # Construct markdown table
    md_lines = []
    md_lines.append("| " + " | ".join(header) + " |")
    md_lines.append("|" + "|".join(["---"] * len(header)) + "|")
    for row in data:
        md_lines.append("| " + " | ".join(row) + " |")
    # Add optional title and note
    md_header = """# CICE Diagnostic Parameters per Simulation
    This table summarizes key dynamics-related parameters from each simulation's `ice_diag.d` file (parsed from first 500 lines).
    """
    md_output = md_header + "\n".join(md_lines) + "\n"
    # === WRITE TO MARKDOWN FILE ===
    md_path.write_text(md_output)
    print(f"✅ Markdown table written to {md_path}")

if __name__ == "__main__":
    main()
