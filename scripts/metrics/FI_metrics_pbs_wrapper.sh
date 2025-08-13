#!/bin/bash
# --- Defaults ---
SIM_NAME=""
ISPD_THRESH="5e-4"
ICE_TYPE_LIST=("FI_B" "FI_Ta" "FI_Tx" "FI_BT")
START_DATE="1993-01-01"
END_DATE="1999-12-31"
OVERWRITE_ZARR=false
OVERWRITE_PNG=false
SMOOTH_FIA_DAYS=15
DRY_RUN=false
PBS_SCRIPT="FI_metrics.pbs"

# --- Help message ---
print_help() {
    echo ""
    echo "Usage: $0 -s SIM_NAME -t ISPD_THRESH [-i ICE_TYPE ...] [-b] [-z] [-p] [-d DAYS] [-n]"
    echo ""
    echo "Options:"
    echo "  -s SIM_NAME        Simulation name (required)"
    echo "  -t ISPD_THRESH     Ice speed threshold (e.g. 5e-4) (required)"
    echo "  -i ICE_TYPE        Ice type(s): FI_B, FI_Ta, FI_Tx, FI_BT (repeatable or comma-separated)"
    echo "  -S START_DATE      Start date (YYYY-MM-DD, default 1993-01-01)"
    echo "  -E END_DATE        End date (YYYY-MM-DD, default 1999-12-31)"
    echo "  -z                 Overwrite existing Zarr output"
    echo "  -p                 Overwrite existing PNG output"
    echo "  -d DAYS            Days to smooth FIA time series (default: 15)"
    echo "  -n                 Dry run: print the qsub command without submitting"
    echo "  -h                 Show this help message and exit"
    echo ""
    echo "Examples:"
    echo "  $0 -s gi-mid -t 1e-3 -i FI_Ta -b -z -p"
    echo "  $0 -s gi-mid -t 1e-3 -i FI_B,FI_Tx -d 7"
    echo ""
}

# --- Parse arguments ---
while getopts "s:t:i:S:E:bzpd:nh" opt; do
    case ${opt} in
        s) SIM_NAME="$OPTARG" ;;
        t) ISPD_THRESH="$OPTARG" ;;
        i)
            IFS=',' read -ra ICE_TYPE_LIST <<< "$OPTARG"
            ;;
        S) START_DATE="$OPTARG" ;;
        E) END_DATE="$OPTARG" ;;
        z) OVERWRITE_ZARR=true ;;
        p) OVERWRITE_PNG=true ;;
        d) SMOOTH_FIA_DAYS="$OPTARG" ;;
        n) DRY_RUN=true ;;
        h) print_help; exit 0 ;;
        \?) echo "Invalid option: -$OPTARG"; print_help; exit 1 ;;
        :) echo "Option -$OPTARG requires an argument."; print_help; exit 1 ;;
    esac
done

# --- Validate required inputs ---
if [[ -z "$SIM_NAME" || -z "$ISPD_THRESH" ]]; then
    echo "-s SIM_NAME and -t ISPD_THRESH are required"
    print_help
    exit 1
fi
# Validate start and end date
if [[ -z "$START_DATE" || -z "$END_DATE" ]]; then
    echo "Must provide both -S START_DATE and -E END_DATE." >&2
    exit 1
fi
# Force START_DATE to first of the month
START_DATE=$(date -d "$START_DATE" +%Y-%m-01)
# Force END_DATE to last day of its month
END_DATE=$(date -d "$(date -d "$END_DATE +1 month" +%Y-%m-01) -1 day" +%F)
# Ensure start is before end
if [[ "$START_DATE" > "$END_DATE" ]]; then
    echo "START_DATE must be before END_DATE." >&2
    exit 1
fi

# --- Create temporary ENVFILE ---
TMP_ENV=$(mktemp /g/data/gv90/da1339/tmp/envfile_XXXXXX.sh)
cat <<EOF > "$TMP_ENV"
export sim_name="$SIM_NAME"
export ispd_thresh="$ISPD_THRESH"
export ice_type=("${ICE_TYPE_LIST[@]}")
export start_date=$START_DATE
export end_date=$END_DATE
export overwrite_zarr=$OVERWRITE_ZARR
export overwrite_png=$OVERWRITE_PNG
export smooth_FIA_days=$SMOOTH_FIA_DAYS
EOF
>> "$TMP_ENV"

JOB_NAME="fi_mets_${SIM_NAME}_${ISPD_THRESH}"
QSUB_CMD="qsub -N ${JOB_NAME} -v ENVFILE=$TMP_ENV $PBS_SCRIPT"
echo "🔧 $QSUB_CMD"
if [[ "$DRY_RUN" = true ]]; then
    echo "[DRY RUN] Would submit PBS job with ENVFILE: $TMP_ENV"
else
    eval $QSUB_CMD
fi
