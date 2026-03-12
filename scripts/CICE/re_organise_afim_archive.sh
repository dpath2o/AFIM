#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

ARCHIVE_ROOT="${HOME}/AFIM_archive"
BIN_WIN_DAYS="11"
BIN_MIN_DAYS="09"
MEAN_PERIOD="15"
EXECUTE=0

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Options:
  --archive-root PATH   Root archive directory (default: ~/AFIM_archive)
  --bin-win NN          Binary window days, zero-padded or integer (default: 08)
  --bin-min NN          Binary minimum days (default: 05)
  --mean-period NN      Rolling mean period days (default: 05)
  --execute             Actually move files (default is dry-run)
  -h, --help            Show this help

Moves old stores like:
  SH/ispd_thresh_5.0e-4/FI_Tb.zarr
  SH/ispd_thresh_5.0e-4/FI_Tb_mets.zarr
  SH/ispd_thresh_5.0e-4/FI_Tb_bin.zarr
  SH/ispd_thresh_5.0e-4/FI_Tb_bin_mets.zarr
  SH/ispd_thresh_5.0e-4/FI_Tb_roll.zarr
  SH/ispd_thresh_5.0e-4/FI_Tb_roll_mets.zarr

to:
  SH/ispd_thresh_5.0e-4/FI/Tb/raw.zarr
  SH/ispd_thresh_5.0e-4/FI/Tb/mets.zarr
  SH/ispd_thresh_5.0e-4/FI/Tb/bin-win-08_bin-min-05/raw.zarr
  SH/ispd_thresh_5.0e-4/FI/Tb/bin-win-08_bin-min-05/mets.zarr
  SH/ispd_thresh_5.0e-4/FI/Tb/roll-days-05/raw.zarr
  SH/ispd_thresh_5.0e-4/FI/Tb/roll-days-05/mets.zarr
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --archive-root) ARCHIVE_ROOT="$2"; shift 2 ;;
        --bin-win)      BIN_WIN_DAYS="$(printf "%02d" "$2")"; shift 2 ;;
        --bin-min)      BIN_MIN_DAYS="$(printf "%02d" "$2")"; shift 2 ;;
        --mean-period)  MEAN_PERIOD="$(printf "%02d" "$2")"; shift 2 ;;
        --execute)      EXECUTE=1; shift ;;
        -h|--help)      usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
    esac
done

move_store() {
    local src="$1"
    local dst="$2"
    [[ -e "$src" ]] || return 0
    mkdir -p "$(dirname "$dst")"
    if [[ -e "$dst" ]]; then
        echo "SKIP (destination exists): $dst"
        return 0
    fi
    if [[ $EXECUTE -eq 1 ]]; then
        echo "MOVE: $src -> $dst"
        mv "$src" "$dst"
    else
        echo "DRYRUN: mv '$src' '$dst'"
    fi
}

for sim_dir in "${ARCHIVE_ROOT}"/*; do
    [[ -d "${sim_dir}/zarr/SH" ]] || continue
    sh_root="${sim_dir}/zarr/SH"
    echo "=== simulation: $(basename "$sim_dir") ==="
    # Optional SI / MIZ legacy layout
    move_store "${sh_root}/SI.zarr"      "${sh_root}/SI/raw.zarr"
    move_store "${sh_root}/SI_mets.zarr" "${sh_root}/SI/mets.zarr"
    move_store "${sh_root}/MIZ.zarr"      "${sh_root}/MIZ/raw.zarr"
    move_store "${sh_root}/MIZ_mets.zarr" "${sh_root}/MIZ/mets.zarr"
    for ispd_dir in "${sh_root}"/ispd_thresh_*; do
        [[ -d "$ispd_dir" ]] || continue
        for src in "${ispd_dir}"/*.zarr; do
            base="$(basename "$src")"
            if [[ "$base" =~ ^(FI|PI)_([A-Za-z0-9]+)(_bin|_roll)?(_mets)?\.zarr$ ]]; then
                ice_type="${BASH_REMATCH[1]}"
                reg_type="${BASH_REMATCH[2]}"
                method_suffix="${BASH_REMATCH[3]}"
                mets_suffix="${BASH_REMATCH[4]}"
                case "$method_suffix" in
                    "")
                        dst_dir="${ispd_dir}/${ice_type}/${reg_type}/raw"
                        ;;
                    "_bin")
                        dst_dir="${ispd_dir}/${ice_type}/${reg_type}/bin-win-${BIN_WIN_DAYS}_bin-min-${BIN_MIN_DAYS}"
                        ;;
                    "_roll")
                        dst_dir="${ispd_dir}/${ice_type}/${reg_type}/roll-days-${MEAN_PERIOD}"
                        ;;
                    *)
                        echo "SKIP (unhandled method): $src"
                        continue
                        ;;
                esac
                if [[ -n "$mets_suffix" ]]; then
                    dst="${dst_dir}/mets.zarr"
                else
                    dst="${dst_dir}/data.zarr"
                fi
                move_store "$src" "$dst"
            fi
        done
    done
done

if [[ $EXECUTE -eq 1 ]]; then
    find "${ARCHIVE_ROOT}" -depth -type d -empty -delete
else
    echo
    echo "Dry-run only. Re-run with --execute to perform the moves."
fi
