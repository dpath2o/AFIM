 #!/bin/bash
set -euo pipefail

period="2000-2011"
base="$HOME/graphical/AFIM"

declare -a regions=("AS" "Aus" "BS" "DML" "EIO" "VOL" "WIO" "WS")
declare -a sims=("elps-min" "LD-static-Cs1e-3" "LD-static-Cs7p5e-4" "LD-static-Cs5e-4" "LD-linear-CL0p25")
declare -a vars=("FIP" "FIHI" "FIST" "FIMAR" "FIMVR")

for var in $"${vars[@]}"; do
    for region in "${regions[@]}"; do
	      outdir="$HOME/graphical/AFIM/_PANELS/${var}/${region}"
	      tmpdir=$(mktemp -d)
	      mkdir -p "$outdir"
	      cleanup() { rm -rf "$tmpdir"; }
	      trap cleanup EXIT
	      for sim in "${sims[@]}"; do
	          infile="${base}/${sim}/${region}/${var}/${var}_${sim}_${region}_${period}.png"
	          outfile="${tmpdir}/${sim}.png"
	          [[ -f "$infile" ]] || { echo "Missing: $infile" >&2; exit 1; }
	          magick "$infile" -background white -fill black -gravity north -splice 0x60 -pointsize 50 -annotate +0+10 "$sim" "$outfile"
	      done
	      montage "${tmpdir}/elps-min.png" "${tmpdir}/LD-static-Cs1e-3.png"	"${tmpdir}/LD-static-Cs7p5e-4.png" "${tmpdir}/LD-static-Cs5e-4.png" "${tmpdir}/LD-linear-CL0p25.png" \
                -tile 3x2 \
		            -geometry +20+20 \
		            -background white \
		            "${outdir}/${var}_${region}_stitched.png"
    done
done
