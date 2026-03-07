#!/bin/bash

# --- Defaults ---
SIM_NAME=""
ISPD_THRESH="5e-4"
ICE_TYPE="FI"
BORC2T_TYPE="Tc"
START_DATE="1993-01-01"
END_DATE="1999-12-31"
ROLL_MEAN=false
OVERWRITE_ZARR=false
OVERWRITE_PNG=false
DRY_RUN=false
PBS_SCRIPT="metrics.pbs"
P_JSON="/home/581/da1339/AFIM/src/AFIM/src/JSONs/sea_ice_config.json"

print_help() {
  echo ""
  echo "Usage: $0 -s SIM_NAME -t ISPD_THRESH [-i ICE_TYPE] [-g BORC2T_TYPE] [-S START_DATE] [-E END_DATE] [-J P_JSON] [-z] [-p] [-n]"
  echo ""
  echo "Options:"
  echo "  -s SIM_NAME        Simulation name (required)"
  echo "  -t ISPD_THRESH     Ice speed threshold (e.g. 5e-4) (required)"
  echo "  -i ICE_TYPE        Ice type: FI, PI, SI, MI (default: FI)"
  echo "  -g BORC2T_TYPE     Regrid method: Tc, Ta, Tb, Tx, B, BT (default: Tc)"
  echo "  -S START_DATE      Start date YYYY-MM-DD (default: 1993-01-01)"
  echo "  -E END_DATE        End date YYYY-MM-DD (default: 1999-12-31)"
  echo "  -J P_JSON          Path to AFIM config JSON (default: $P_JSON)"
  echo "  -r                 enable rolling-mean metrics (default: false)"
  echo "  -z                 Overwrite existing Zarr output"
  echo "  -p                 Overwrite existing PNG output"
  echo "  -n                 Dry run: print the qsub command without submitting"
  echo "  -h                 Show this help message and exit"
  echo ""
}

# --- Parse arguments ---
while getopts "s:t:i:g:S:E:J:rzpnh" opt; do
  case ${opt} in
    s) SIM_NAME="$OPTARG" ;;
    t) ISPD_THRESH="$OPTARG" ;;
    i) ICE_TYPE="$OPTARG" ;;
    g) BORC2T_TYPE="$OPTARG" ;;
    S) START_DATE="$OPTARG" ;;
    E) END_DATE="$OPTARG" ;;
    J) P_JSON="$OPTARG" ;;
    r) ROLL_MEAN=true ;;
    z) OVERWRITE_ZARR=true ;;
    p) OVERWRITE_PNG=true ;;
    n) DRY_RUN=true ;;
    h) print_help; exit 0 ;;
    \?) echo "Invalid option: -$OPTARG"; print_help; exit 1 ;;
    :)  echo "Option -$OPTARG requires an argument."; print_help; exit 1 ;;
  esac
done

# --- Validate required inputs ---
if [[ -z "$SIM_NAME" || -z "$ISPD_THRESH" ]]; then
  echo "-s SIM_NAME and -t ISPD_THRESH are required"
  print_help
  exit 1
fi

if [[ -z "$START_DATE" || -z "$END_DATE" ]]; then
  echo "Must provide both -S START_DATE and -E END_DATE." >&2
  exit 1
fi

# Force START_DATE to first of month; END_DATE to last day of month
START_DATE=$(date -d "$START_DATE" +%Y-%m-01)
END_DATE=$(date -d "$(date -d "$END_DATE +1 month" +%Y-%m-01) -1 day" +%F)

if [[ "$START_DATE" > "$END_DATE" ]]; then
  echo "START_DATE must be before END_DATE." >&2
  exit 1
fi

# --- Create temporary ENVFILE ---
TMP_ENV=$(mktemp /g/data/gv90/da1339/tmp/envfile_XXXXXX.sh)

cat > "$TMP_ENV" << EOF
export sim_name="${SIM_NAME}"
export ispd_thresh="${ISPD_THRESH}"
export ice_type="${ICE_TYPE}"
export borc2t_type="${BORC2T_TYPE}"
export start_date="${START_DATE}"
export end_date="${END_DATE}"
export rolling_mean=${ROLL_MEAN}
export overwrite_zarr=${OVERWRITE_ZARR}
export overwrite_png=${OVERWRITE_PNG}
export P_JSON="${P_JSON}"
EOF

JOB_NAME="${SIM_NAME}_${ICE_TYPE}_${BORC2T_TYPE}_${ISPD_THRESH}"
QSUB_CMD="qsub -N ${JOB_NAME} -v ENVFILE=$TMP_ENV $PBS_SCRIPT"

echo "$QSUB_CMD"
if [[ "$DRY_RUN" = true ]]; then
  echo "[DRY RUN] Would submit PBS job with ENVFILE: $TMP_ENV"
else
  eval "$QSUB_CMD"
fi
