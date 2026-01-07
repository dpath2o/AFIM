#!/bin/bash

# --- Defaults ---
B2T_TYPE_LIST=()
START_DATE="1993-01-01"
END_DATE=""
YEARLY=false
DRY_RUN=false
NC_ENG="netcdf4"
OVERWRITE_ZARR=false
DELETE_ORIGINAL_ICEH=false

# --- PBS Script ---
PBS_SCRIPT=classify_SI.pbs

print_help() {
    echo ""
    echo "Usage:"
    echo "  $0 -s SIM_NAME [-t ISPD_THRESH] [-i B2T_TYPE ...] [-S START_DATE -E END_DATE] [-y] [-z] [-k] [-n]"
    echo ""
    echo "Options:"
    echo "  -s SIM_NAME      Simulation name (required)"
    echo "  -t ISPD_THRESH   Optional: speed threshold only used for output directory naming (SI does not threshold speed)"
    echo "  -i B2T_TYPE      Optional: ice-speed vector type(s) used for naming/metadata: B, Ta, Tb, Tx, BT"
    echo "  -S START_DATE    Start date (YYYY-MM-DD, forced to first of month; default 1993-01-01)"
    echo "  -E END_DATE      End date (YYYY-MM-DD, forced to last of month)"
    echo "  -x NC_ENG         NetCDF engine used by xarray mfdataset method; default 'netcdf4'"
    echo "  -z               overwrite zarrs"
    echo "  -y               process yearly loop (submit one job per year)"
    echo "  -k               after creating monthly iceh_*.zarr groups, delete original daily ice history files"
    echo "  -n               Dry run: print PBS command without submitting"
    echo "  -h               Show this help message and exit"
    echo ""
}

# --- Parse options ---
while getopts "s:t:i:rdzkhS:E:nyx:" opt; do
    case ${opt} in
        s) SIM_NAME="$OPTARG" ;;
        t) ISPD_THRESH="$OPTARG" ;;
        i) B2T_TYPE_LIST+=("$OPTARG") ;;
        z) OVERWRITE_ZARR=true ;;
        k) DELETE_ORIGINAL_ICEH=true ;;
        S) START_DATE="$OPTARG" ;;
        E) END_DATE="$OPTARG" ;;
        x) NC_ENG="$OPTARG" ;;
        y) YEARLY=true ;;
        n) DRY_RUN=true ;;
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

# --- Default vector types if none specified ---
if [ ${#B2T_TYPE_LIST[@]} -eq 0 ]; then
    B2T_TYPE_LIST=("Tb")
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

    VAR_PASS="SIM_NAME=${SIM_NAME},MONTH=${MONTH},YEAR=${YEAR},START_DATE=${DT0},END_DATE=${DTN},NC_ENG=${NC_ENG}"
    [ -n "${ISPD_THRESH}" ] && VAR_PASS+=",ISPD_THRESH=${ISPD_THRESH}"
    [ "${#B2T_TYPE_LIST[@]}" -gt 0 ] && VAR_PASS+=",B2T_TYPE=${B2T_TYPE_LIST[*]}"
    [ "${OVERWRITE_ZARR}" = true ] && VAR_PASS+=",OVERWRITE_ZARR=true"
    [ "${DELETE_ORIGINAL_ICEH}" = true ] && VAR_PASS+=",DELETE_ORIGINAL_ICEH=true"

    JOB_NAME="si_${SIM_NAME}_${JOB_SUFFIX}"
    [ -n "${ISPD_THRESH}" ] && JOB_NAME+="_${ISPD_THRESH}"

    QSUB_CMD="qsub -N ${JOB_NAME} -v "${VAR_PASS}" ${PBS_SCRIPT}"

    if [ "$DRY_RUN" = true ]; then
        echo "🧪 [DRY RUN] Would submit: $QSUB_CMD"
    else
        eval $QSUB_CMD
    fi

    current_period=$next_period
done
