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

# --- PBS Script ---
PBS_SCRIPT=classify_sea_ice.pbs

print_help() {
    echo ""
    echo "Usage:"
    echo "  $0 -s SIM_NAME -t ISPD_THRESH [-i BorC2T_TYPE ...] [-r] [-d]"
    echo ""
    echo "Options:"
    echo "  -s SIM_NAME      Simulation name (required)"
    echo "  -i ICE_TYPE      ice classification type 'FI', 'PI', 'SI', 'MI'"
    echo "  -t ISPD_THRESH   Ice speed threshold (e.g. 1e-3) (required)"
    echo "  -g BorC2T_TYPE   B-grid regridding to T-grid method(s): B (no regridding), Ta, Tb, Tc, Tx (some combination)"
    echo "  -S START_DATE    Start date (YYYY-MM-DD, forced to first of month; default 1993-01-01)"
    echo "  -E END_DATE      End date (YYYY-MM-DD, forced to last of month)"
    echo "  -x NC_ENG        NetCDF engine used by xarray mfdataset method; default 'netcdf4'"
    echo "  -z               overwrite zarrs"
    echo "  -y               process yearly loop"
    echo "  -k               if enabled then when creating new iceh_*.zarr monthly files, this will then delete the original daily ice history files for a particular month"
    echo "  -n               Dry run: print PBS command without submitting"
    echo "  -h               Show this help message and exit"
    echo ""
    echo "Examples:"
    echo "  $0 -s gi-mid -t 1e-3 -i Ta"
    echo "  $0 -s gi-mid -t 1e-3 -i Tb -i Tx"
    echo "  $0 -s gi-mid # runs defaults"
    echo ""
}

# --- Parse options ---
while getopts "s:i:t:g:zkhS:E:ny" opt; do
    case ${opt} in
        s) SIM_NAME="$OPTARG" ;;
        i) ICE_TYPE="$OPTARG" ;;
        t) ISPD_THRESH="$OPTARG" ;;
        g) BorC2T_TYPE_LIST+=("$OPTARG") ;;
        z) OVERWRITE_ZARR=true ;;
        k) DELETE_ORIGINAL_ICEH=true ;;
        S) START_DATE="$OPTARG" ;;
        E) END_DATE="$OPTARG" ;;
        x) NC_ENG="$OPTARG" ;;
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

# --- Default ISPD types if none specified ---
if [ ${#BorC2T_TYPE_LIST[@]} -eq 0 ]; then
    BorC2T_TYPE_LIST=("Tb") #"B" "Ta" "Tx" "BT"
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
    VAR_PASS="SIM_NAME=${SIM_NAME},MONTH=${MONTH},YEAR=${YEAR},START_DATE=${DT0},END_DATE=${DTN},ICE_TYPE=${ICE_TYPE},ISPD_THRESH=${ISPD_THRESH},NC_ENG=${NC_ENG}"
    [ "${#BorC2T_TYPE_LIST[@]}" -gt 0 ] && VAR_PASS+=",BorC2T_TYPE=${BorC2T_TYPE_LIST[*]}"
    [ "${OVERWRITE_ZARR}" = true ] && VAR_PASS+=",OVERWRITE_ZARR=true"
    [ "${DELETE_ORIGINAL_ICEH}" = true ] && VAR_PASS+=",DELETE_ORIGINAL_ICEH=true"
    JOB_NAME="${SIM_NAME}_${ICE_TYPE}_${ISPD_THRESH}_${JOB_SUFFIX}"
    QSUB_CMD="qsub -N ${JOB_NAME} -v \"${VAR_PASS}\" ${PBS_SCRIPT}"
    if [ "$DRY_RUN" = true ]; then
        echo "🧪 [DRY RUN] Would submit: $QSUB_CMD"
    else
        eval $QSUB_CMD
    fi
    current_period=$next_period
done
