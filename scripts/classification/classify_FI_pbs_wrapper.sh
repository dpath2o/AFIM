#!/bin/bash

# --- Defaults ---
ROLLING=false
DAILY=false
IVEC_TYPE_LIST=()
START_DATE="1993-01-01"
END_DATE=""
YEARLY=false
DRY_RUN=false
OVERWRITE_ZARR=false
DELETE_ORIGINAL_ICEH=false

# --- PBS Script ---
PBS_SCRIPT=classify_FI.pbs

print_help() {
    echo ""
    echo "Usage:"
    echo "  $0 -s SIM_NAME -t ISPD_THRESH [-i IVEC_TYPE ...] [-r] [-d]"
    echo ""
    echo "Options:"
    echo "  -s SIM_NAME      Simulation name (required)"
    echo "  -t ISPD_THRESH   Ice speed threshold (e.g. 1e-3) (required)"
    echo "  -i IVEC_TYPE     Ice speed type(s): B, Ta, Tx, BT (repeatable)"
    echo "  -S START_DATE    Start date (YYYY-MM-DD, forced to first of month; default 1993-01-01)"
    echo "  -E END_DATE      End date (YYYY-MM-DD, forced to last of month)"
    echo "  -r               Enable rolling fast ice masking"
    echo "  -d               Enable daily fast ice masking"
    echo "  -z               overwrite zarrs"
    echo "  -y               process yearly loop"
    echo "  -k               if enabled then when creating new iceh_*.zarr monthly files, this will then delete the original daily ice history files for a particular month"
    echo "  -n               Dry run: print PBS command without submitting"
    echo "  -h               Show this help message and exit"
    echo ""
    echo "Examples:"
    echo "  $0 -s gi-mid -t 1e-3 -i Ta -d"
    echo "  $0 -s gi-mid -t 1e-3 -i B -i Ta -r -d"
    echo "  $0 -s gi-mid -t 1e-3 -r -d     # runs all ispd types by default"
    echo ""
}

# --- Parse options ---
while getopts "s:t:i:rdzkhS:E:ny" opt; do
    case ${opt} in
        s) SIM_NAME="$OPTARG" ;;
        t) ISPD_THRESH="$OPTARG" ;;
        i) IVEC_TYPE_LIST+=("$OPTARG") ;;
        r) ROLLING=true ;;
        d) DAILY=true ;;
        z) OVERWRITE_ZARR=true ;;
        k) DELETE_ORIGINAL_ICEH=true ;;
        S) START_DATE="$OPTARG" ;;
        E) END_DATE="$OPTARG" ;;
        y) YEARLY=true ;;
        n) DRY_RUN=true ;;  # <-- Dry run flag
        h) print_help; exit 0 ;;
        \?) echo "Invalid option: -$OPTARG" >&2; print_help; exit 1 ;;
        :) echo "Option -$OPTARG requires an argument." >&2; print_help; exit 1 ;;
    esac
done

# --- Check required arguments ---
if [[ -z "$SIM_NAME" ]]; then
    echo "❌ Error: -s SIM_NAME is required." >&2
    print_help
    exit 1
fi

if { [[ "$ROLLING" == true ]] || [[ "$DAILY" == true ]]; } && [[ -z "$ISPD_THRESH" ]]; then
    echo "❌ Error: -t ISPD_THRESH is required when using -r (rolling) or -d (daily) modes." >&2
    print_help
    exit 1
fi


# --- Default ISPD types if none specified ---
if [ ${#ISPD_TYPE_LIST[@]} -eq 0 ]; then
    IVEC_TYPE_LIST=("B" "Ta" "Tx" "BT")
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
    VAR_PASS="SIM_NAME=${SIM_NAME},MONTH=${MONTH},YEAR=${YEAR},START_DATE=${DT0},END_DATE=${DTN},ISPD_THRESH=${ISPD_THRESH}"
    [ "${#IVEC_TYPE_LIST[@]}" -gt 0 ] && VAR_PASS+=",ISPD_TYPE=${IVEC_TYPE_LIST[*]}"
    [ "${ROLLING}" = true ] && VAR_PASS+=",ROLLING=true"
    [ "${DAILY}" = true ] && VAR_PASS+=",DAILY=true"
    [ "${OVERWRITE_ZARR}" = true ] && VAR_PASS+=",OVERWRITE_ZARR=true"
    [ "${DELETE_ORIGINAL_ICEH}" = true ] && VAR_PASS+=",DELETE_ORIGINAL_ICEH=true"
    JOB_NAME="fi_${SIM_NAME}_${JOB_SUFFIX}_${ISPD_THRESH}"
    QSUB_CMD="qsub -N ${JOB_NAME} -v \"${VAR_PASS}\" ${PBS_SCRIPT}"
    if [ "$DRY_RUN" = true ]; then
        echo "🧪 [DRY RUN] Would submit: $QSUB_CMD"
    else
        eval $QSUB_CMD
    fi
    current_period=$next_period
done
