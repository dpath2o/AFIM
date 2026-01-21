#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob
BASE="/g/data/gv90/da1339/GRAPHICAL/AFIM"
OUTBASE="${BASE}/animations"
VAR_NAME="fihi"
YR="1993"
sim_names=("elps-min" 
           "LD-GI-Cs-1e-2"
           "LD-Cs5e-4-mld60"
           "LD-Cs5e-4-mld10"
           "LD-quad-Cq-1p5"
           "LD-quad-cap-7p5e-2"
           "LD-sum-static-quad"
           "LD-sum-quad-cap"
           "LD-lin-C_L-1p1")
regions=("AS" "Aus" "BS" "DML" "EIO" "VOL" "WIO" "WS")
FPS=6
CRF=20
PRESET="medium"
for sim in "${sim_names[@]}"; do
  for reg in "${regions[@]}"; do
    in_dir="${BASE}/${sim}/${reg}/${VAR_NAME}"
    out_dir="${OUTBASE}/${sim}/${reg}/${VAR_NAME}"
    mkdir -p "${out_dir}"
    # Gather frames (brace expansion for months)
    frames=( "${in_dir}/${YR}"-{09,10}*.png )
    if (( ${#frames[@]} == 0 )); then
      echo "[skip] no frames: ${in_dir}/${YR}-{09,10}*.png"
      continue
    fi
    # Sort frames (lexicographic works for YYYY-MM-DD*)
    IFS=$'\n' frames_sorted=( $(printf '%s\n' "${frames[@]}" | sort) )
    unset IFS
    first_base="$(basename "${frames_sorted[0]}")"
    last_base="$(basename "${frames_sorted[-1]}")"
    # Extract YYYY-MM-DD from filename start (assumes it starts with date)
    first_date="${first_base:0:10}"
    last_date="${last_base:0:10}"
    out_mp4="${out_dir}/${first_date}-${last_date}_fice_${sim}_${reg}.mp4"
    echo "[make] ${out_mp4}"
    echo "       frames: ${#frames_sorted[@]}  (${first_date} .. ${last_date})"
    # Use glob input (simpler than writing a concat list). Requires glob pattern:
    # We run ffmpeg from inside in_dir so its glob finds files.
    (cd "${in_dir}"
      ffmpeg -hide_banner -loglevel warning -y \
        -framerate "${FPS}" \
        -pattern_type glob -i "${YR}-{09,10}*.png" \
        -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2,format=yuv420p" \
        -c:v libx264 -preset "${PRESET}" -crf "${CRF}" \
        "${out_mp4}")
  done
done

# #!/usr/bin/env bash
# set -Eeuo pipefail

# trap 'echo "ERROR: line ${LINENO}: ${BASH_COMMAND}" >&2' ERR

# SIM="elps-min"
# VAR="uvel"
# REGIONS=(DML WIO EIO Aus VOL AS BS WS)

# # Set MONTH="YYYY-MM" for a month; or MONTH="" for all frames
# MONTH="2001-09"

# SRC_ROOT="/g/data/gv90/da1339/GRAPHICAL/AFIM/${SIM}"
# DST_ROOT="${HOME}/graphical/animations/${SIM}/${VAR}"

# FPS=6
# CRF=18
# PRESET="medium"
# KEEP_TMP=0   # set to 1 to keep _frames_tmp for debugging

# command -v ffmpeg >/dev/null 2>&1 || {
#   echo "ERROR: ffmpeg not found in PATH. On Gadi: module load ffmpeg" >&2
#   exit 1
# }

# mkdir -p "${DST_ROOT}"

# for REGION in "${REGIONS[@]}"; do
#   SRC_DIR="${SRC_ROOT}/${REGION}/${VAR}"
#   OUT_DIR="${DST_ROOT}/${REGION}"
#   mkdir -p "${OUT_DIR}"

#   if [[ ! -d "${SRC_DIR}" ]]; then
#     echo "[SKIP] Missing source dir: ${SRC_DIR}"
#     continue
#   fi

#   if [[ -n "${MONTH}" ]]; then
#     PATTERN="${SRC_DIR}/${MONTH}-*_${SIM}_${REGION}_${VAR}.png"
#     TAG="${MONTH}"
#   else
#     PATTERN="${SRC_DIR}/*.png"
#     TAG="ALL"
#   fi

#   shopt -s nullglob
#   FRAMES=( ${PATTERN} )
#   shopt -u nullglob

#   if (( ${#FRAMES[@]} == 0 )); then
#     echo "[SKIP] No frames for ${REGION} (pattern: ${PATTERN})"
#     continue
#   fi

#   # Deterministic ordering
#   IFS=$'\n' FRAMES_SORTED=( $(printf '%s\n' "${FRAMES[@]}" | sort) )
#   unset IFS

#   OUT_MP4="${OUT_DIR}/${SIM}_${REGION}_${VAR}_${TAG}.mp4"
#   TMP_DIR="${OUT_DIR}/_frames_tmp"
#   rm -rf "${TMP_DIR}"
#   mkdir -p "${TMP_DIR}"

#   echo "[MAKE] ${REGION}: ${#FRAMES_SORTED[@]} frames -> ${OUT_MP4}"

#   i=0
#   for f in "${FRAMES_SORTED[@]}"; do
#     printf -v fname "frame_%05d.png" "${i}"
#     if ln -s "${f}" "${TMP_DIR}/${fname}" 2>/dev/null; then
#       :
#     else
#       cp -p "${f}" "${TMP_DIR}/${fname}"
#     fi
#     i=$((i+1))   # <-- FIX: avoids set -e early exit
#   done

#   ffmpeg -y \
#     -framerate "${FPS}" \
#     -start_number 0 \
#     -i "${TMP_DIR}/frame_%05d.png" \
#     -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2,format=yuv420p" \
#     -c:v libx264 \
#     -crf "${CRF}" \
#     -preset "${PRESET}" \
#     -movflags +faststart \
#     "${OUT_MP4}"

#   if [[ ! -s "${OUT_MP4}" ]]; then
#     echo "ERROR: ffmpeg did not produce a valid output: ${OUT_MP4}" >&2
#     exit 1
#   fi

#   if (( KEEP_TMP == 0 )); then
#     rm -rf "${TMP_DIR}"
#   else
#     echo "[INFO] Keeping temp frames at: ${TMP_DIR}"
#   fi
# done

# echo "Done. Animations under: ${DST_ROOT}"
