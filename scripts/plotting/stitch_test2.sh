#!/bin/bash
set -euo pipefail

regions=(AS Aus BS DML EIO VOL WIO WS)
period="2000-2011"
var="FIP"
base="$HOME/graphical/AFIM"
outbase="$HOME/graphical/AFIM/_PANELS/${var}"
mkdir -p "$outbase"

sims=("elps-min" "LD-static-Cs1e-3" "LD-static-Cs7p5e-4" "LD-static-Cs5e-4" "LD-linear-CL0p25")

for region in "${regions[@]}"; do
    outdir="${outbase}/${region}"
    mkdir -p "$outdir"
    args=()
    for sim in "${sims[@]}"; do
        file="${base}/${sim}/${region}/${var}/${var}_${sim}_${region}_${period}.png"
        if [[ ! -f "$file" ]]; then
            echo "Missing: $file" >&2
            continue
        fi
        args+=("(" "$file" -background white -fill black -gravity north -splice 0x60 -pointsize 28 -annotate +0+10 "$sim" ")" )
    done
    if [[ ${#args[@]} -eq 0 ]]; then
        echo "No files found for ${region}" >&2
        continue
    fi
    magick "${args[@]}" -tile 5x1 -geometry +20+20 -background white montage "${outdir}/${var}_${region}_stitched.png"
    echo "Wrote ${outdir}/${var}_${region}_stitched.png"
done
