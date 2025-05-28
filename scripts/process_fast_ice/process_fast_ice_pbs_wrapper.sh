#!/bin/bash

# --- Defaults ---
ROLLING=false
DAILY=false
ISPD_TYPE_LIST=()
START_DATE="1993-01-01"
END_DATE="1993-01-31"
DRY_RUN=false
OVERWRITE_ZARR=false
DELETE_ORIGINAL_ICEH=false

# --- PBS Script ---
PBS_SCRIPT=process_fast_ice.pbs

print_help() {
    echo ""
    echo "Usage:"
    echo "  $0 -s SIM_NAME -t ISPD_THRESH [-i ISPD_TYPE ...] [-r] [-d]"
    echo ""
    echo "Options:"
    echo "  -s SIM_NAME      Simulation name (required)"
    echo "  -t ISPD_THRESH   Ice speed threshold (e.g. 1e-3) (required)"
    echo "  -i ISPD_TYPE     Ice speed type(s): ispd_B, ispd_Ta, ispd_Tx (repeatable)"
    echo "  -S START_DATE    Start date (YYYY-MM-DD, forced to first of month)"
    echo "  -E END_DATE      End date (YYYY-MM-DD, forced to last of month)"
    echo "  -r               Enable rolling fast ice masking"
    echo "  -d               Enable daily fast ice masking"
    echo "  -z               overwrite zarrs"
    echo "  -k               if enabled then when creating new iceh_*.zarr monthly files, this will then delete the original daily ice history files for a particular month"
    echo "  -n               Dry run: print PBS command without submitting"
    echo "  -h               Show this help message and exit"
    echo ""
    echo "Examples:"
    echo "  $0 -s gi-mid -t 1e-3 -i ispd_Ta -d"
    echo "  $0 -s gi-mid -t 1e-3 -i ispd_B -i ispd_Ta -r -d"
    echo "  $0 -s gi-mid -t 1e-3 -r -d     # runs all ispd types by default"
    echo ""
}

# --- Parse options ---
while getopts "s:t:i:rdzkhS:E:n" opt; do
    case ${opt} in
        s) SIM_NAME="$OPTARG" ;;
        t) ISPD_THRESH="$OPTARG" ;;
        i) ISPD_TYPE_LIST+=("$OPTARG") ;;
        r) ROLLING=true ;;
        d) DAILY=true ;;
        z) OVERWRITE_ZARR=true ;;
        k) DELETE_ORIGINAL_ICEH=true ;;
        S) START_DATE="$OPTARG" ;;
        E) END_DATE="$OPTARG" ;;
        n) DRY_RUN=true ;;  # <-- Dry run flag
        h) print_help; exit 0 ;;
        \?) echo "Invalid option: -$OPTARG" >&2; print_help; exit 1 ;;
        :) echo "Option -$OPTARG requires an argument." >&2; print_help; exit 1 ;;
    esac
done

# --- Check required arguments ---
if [[ -z "$SIM_NAME" || -z "$ISPD_THRESH" ]]; then
    echo "❌ Error: -s SIM_NAME and -t ISPD_THRESH are required." >&2
    print_help
    exit 1
fi

# --- Default ISPD types if none specified ---
if [ ${#ISPD_TYPE_LIST[@]} -eq 0 ]; then
    ISPD_TYPE_LIST=("ispd_B" "ispd_Ta" "ispd_Tx" "ispd_BT")
fi

# Validate start and end date
if [[ -z "$START_DATE" || -z "$END_DATE" ]]; then
    echo "❌ Must provide both -S START_DATE and -E END_DATE." >&2
    exit 1
fi
# Force START_DATE to first of the month
START_DATE=$(date -d "$START_DATE" +%Y-%m-01)
# Force END_DATE to last day of its month
END_DATE=$(date -d "$(date -d "$END_DATE +1 month" +%Y-%m-01) -1 day" +%F)
# Ensure start is before end
if [[ "$START_DATE" > "$END_DATE" ]]; then
    echo "❌ START_DATE must be before END_DATE." >&2
    exit 1
fi

current_month=$(date -d "$START_DATE" +%Y-%m-01)
while [[ "$current_month" < "$END_DATE" || "$current_month" == "$END_DATE" ]]; do
    DT0="$current_month"
    DTN=$(date -d "$DT0 +1 month -1 day" +%F)
    YEAR=$(date -d "$DT0" +%Y)
    MONTH=$(date -d "$DT0" +%m)
    VAR_PASS="SIM_NAME=${SIM_NAME},MONTH=${MONTH},YEAR=${YEAR},START_DATE=${DT0},END_DATE=${DTN},ISPD_THRESH=${ISPD_THRESH}"
    [ "${#ISPD_TYPE_LIST[@]}" -gt 0 ] && VAR_PASS+=",ISPD_TYPE=${ISPD_TYPE_LIST[*]}"
    [ "${ROLLING}" = true ] && VAR_PASS+=",ROLLING=true"
    [ "${DAILY}" = true ] && VAR_PASS+=",DAILY=true"
    [ "${OVERWRITE_ZARR}" = true ] && VAR_PASS+=",OVERWRITE_ZARR=true"
    [ "${DELETE_ORIGINAL_ICEH}" = true ] && VAR_PASS+=",DELETE_ORIGINAL_ICEH=true"
    JOB_NAME="fi_${SIM_NAME}_${YEAR}_${MONTH}"
    QSUB_CMD="qsub -N ${JOB_NAME} -v \"${VAR_PASS}\" ${PBS_SCRIPT}"
    if [ "$DRY_RUN" = true ]; then
        echo "🧪 [DRY RUN] Would submit: $QSUB_CMD"
    else
        eval $QSUB_CMD
    fi
    current_month=$(date -d "$current_month +1 month" +%Y-%m-01)
done
