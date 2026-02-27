#!/bin/bash

# --- Defaults ---
ICE_TYPE="FI"
ISPD_THRESH="5e-4"
BorC2T_TYPE_LIST=()
START_DATE="1993-01-01"
END_DATE=""
YEARLY=false
DRY_RUN=false
NC_ENG="netcdf4"
OVERWRITE_ZARR=false
DELETE_ORIGINAL_ICEH=false

P_JSON="/home/581/da1339/AFIM/src/AFIM/src/JSONs/sea_ice_config.json"

# --- PBS Script ---
PBS_SCRIPT=classify_sea_ice.pbs

print_help() {
  echo ""
  echo "Usage:"
  echo "  $0 -s SIM_NAME -t ISPD_THRESH [-i ICE_TYPE] [-g BorC2T_TYPE ...] [-S START_DATE] [-E END_DATE] [-x NC_ENG] [-J P_JSON] [-z] [-y] [-k] [-n]"
  echo ""
  echo "Options:"
  echo "  -s SIM_NAME        Simulation name (required)"
  echo "  -i ICE_TYPE        Ice type: FI, PI, SI, MI (default: FI)"
  echo "  -t ISPD_THRESH     Ice speed threshold (e.g. 1e-3) (required)"
  echo "  -g BorC2T_TYPE     B-grid->T-grid type(s): B, Ta, Tb, Tc, Tx (repeatable; default Tb)"
  echo "  -S START_DATE      Start date YYYY-MM-DD (forced to first of month; default 1993-01-01)"
  echo "  -E END_DATE        End date YYYY-MM-DD (forced to last of month) (required)"
  echo "  -x NC_ENG          NetCDF engine for xarray (default: netcdf4)"
  echo "  -J P_JSON          Path to AFIM config JSON (default: $P_JSON)"
  echo "  -z                 Overwrite zarrs"
  echo "  -y                 Process yearly loop (else monthly)"
  echo "  -k                 Delete original daily ice history files after monthly zarr creation"
  echo "  -n                 Dry run (print qsub without submitting)"
  echo "  -h                 Help"
  echo ""
}

# --- Parse options ---
while getopts "s:i:t:g:S:E:x:J:zkynh" opt; do
  case ${opt} in
    s) SIM_NAME="$OPTARG" ;;
    i) ICE_TYPE="$OPTARG" ;;
    t) ISPD_THRESH="$OPTARG" ;;
    g) BorC2T_TYPE_LIST+=("$OPTARG") ;;
    S) START_DATE="$OPTARG" ;;
    E) END_DATE="$OPTARG" ;;
    x) NC_ENG="$OPTARG" ;;
    J) P_JSON="$OPTARG" ;;
    z) OVERWRITE_ZARR=true ;;
    k) DELETE_ORIGINAL_ICEH=true ;;
    y) YEARLY=true ;;
    n) DRY_RUN=true ;;
    h) print_help; exit 0 ;;
    \?) echo "Invalid option: -$OPTARG" >&2; print_help; exit 1 ;;
    :)  echo "Option -$OPTARG requires an argument." >&2; print_help; exit 1 ;;
  esac
done

# --- Check required arguments ---
if [[ -z "$SIM_NAME" ]]; then
  echo "❌ Error: -s SIM_NAME is required." >&2
  print_help
  exit 1
fi

if [[ -z "$START_DATE" || -z "$END_DATE" ]]; then
  echo "❌ Must provide both -S START_DATE and -E END_DATE." >&2
  exit 1
fi

# --- Default BorC2T types if none specified ---
if [ ${#BorC2T_TYPE_LIST[@]} -eq 0 ]; then
  BorC2T_TYPE_LIST=("Tb")
fi

# --- Normalize dates ---
START_DATE=$(date -d "$START_DATE" +%Y-%m-01)
END_DATE=$(date -d "$(date -d "$END_DATE +1 month" +%Y-%m-01) -1 day" +%F)

if [[ "$START_DATE" > "$END_DATE" ]]; then
  echo "❌ START_DATE must be before END_DATE." >&2
  exit 1
fi

current_period=$(date -d "$START_DATE" +%Y-%m-01)
end_period=$(date -d "$END_DATE" +%Y-%m-01)

while [[ "$current_period" < "$end_period" || "$current_period" == "$end_period" ]]; do
  if [ "$YEARLY" = true ]; then
    DT0=$(date -d "$current_period" +%Y-01-01)
    DTN=$(date -d "$DT0 +1 year -1 day" +%F)
    JOB_SUFFIX=$(date -d "$DT0" +%Y)
    next_period=$(date -d "$current_period +1 year" +%Y-%m-01)
  else
    DT0="$current_period"
    DTN=$(date -d "$DT0 +1 month -1 day" +%F)
    JOB_SUFFIX=$(date -d "$DT0" +%Y_%m)
    next_period=$(date -d "$current_period +1 month" +%Y-%m-01)
  fi

  YEAR=$(date -d "$DT0" +%Y)
  MONTH=$(date -d "$DT0" +%m)

  VAR_PASS="SIM_NAME=${SIM_NAME},MONTH=${MONTH},YEAR=${YEAR},START_DATE=${DT0},END_DATE=${DTN},ICE_TYPE=${ICE_TYPE},ISPD_THRESH=${ISPD_THRESH},NC_ENG=${NC_ENG},P_JSON=${P_JSON}"

  [ "${#BorC2T_TYPE_LIST[@]}" -gt 0 ] && VAR_PASS+=",BorC2T_TYPE=${BorC2T_TYPE_LIST[*]}"
  [ "${OVERWRITE_ZARR}" = true ] && VAR_PASS+=",OVERWRITE_ZARR=true"
  [ "${DELETE_ORIGINAL_ICEH}" = true ] && VAR_PASS+=",DELETE_ORIGINAL_ICEH=true"

  JOB_NAME="${SIM_NAME}_${ICE_TYPE}_${ISPD_THRESH}_${JOB_SUFFIX}"
  QSUB_CMD="qsub -N ${JOB_NAME} -v \"${VAR_PASS}\" ${PBS_SCRIPT}"

  if [ "$DRY_RUN" = true ]; then
    echo " [DRY RUN] Would submit: $QSUB_CMD"
  else
    eval $QSUB_CMD
  fi

  current_period=$next_period
done
