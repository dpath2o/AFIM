#!/bin/bash
# --- Defaults ---
SIM_NAME=""
ISPD_THRESH="1e-3"
ICE_TYPE_LIST=("FI_B" "FI_Ta" "FI_Tx" "FI_BT")
COMPUTE_BOOLEAN=false
OVERWRITE_ZARR=false
OVERWRITE_PNG=false
SMOOTH_FIA_DAYS=15
DRY_RUN=false
PBS_SCRIPT="compute_fast_ice_metrics.pbs"

# --- Help message ---
print_help() {
    echo ""
    echo "Usage: $0 -s SIM_NAME -t ISPD_THRESH [-i ICE_TYPE ...] [-b] [-z] [-p] [-d DAYS] [-n]"
    echo ""
    echo "Options:"
    echo "  -s SIM_NAME        Simulation name (required)"
    echo "  -t ISPD_THRESH     Ice speed threshold (e.g. 1e-3) (required)"
    echo "  -i ICE_TYPE        Ice type(s): FI_B, FI_Ta, FI_Tx, FI_BT (repeatable or comma-separated)"
    echo "  -b                 Enable boolean fast ice computation"
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
while getopts "s:t:i:bzpd:nh" opt; do
    case ${opt} in
        s) SIM_NAME="$OPTARG" ;;
        t) ISPD_THRESH="$OPTARG" ;;
        i)
            IFS=',' read -ra ICE_TYPE_LIST <<< "$OPTARG"
            ;;
        b) COMPUTE_BOOLEAN=true ;;
        z) OVERWRITE_ZARR=true ;;
        p) OVERWRITE_PNG=true ;;
        d) SMOOTH_FIA_DAYS="$OPTARG" ;;
        n) DRY_RUN=true ;;
        h) print_help; exit 0 ;;
        \?) echo "❌ Invalid option: -$OPTARG"; print_help; exit 1 ;;
        :) echo "❌ Option -$OPTARG requires an argument."; print_help; exit 1 ;;
    esac
done

# --- Validate required inputs ---
if [[ -z "$SIM_NAME" || -z "$ISPD_THRESH" ]]; then
    echo "❌ -s SIM_NAME and -t ISPD_THRESH are required"
    print_help
    exit 1
fi

# --- Create temporary ENVFILE ---
TMP_ENV=$(mktemp /g/data/gv90/da1339/tmp/envfile_XXXXXX.sh)
cat <<EOF > "$TMP_ENV"
export sim_name="$SIM_NAME"
export ispd_thresh="$ISPD_THRESH"
export ice_type=("${ICE_TYPE_LIST[@]}")
export compute_boolean=$COMPUTE_BOOLEAN
export overwrite_zarr=$OVERWRITE_ZARR
export overwrite_png=$OVERWRITE_PNG
export smooth_FIA_days=$SMOOTH_FIA_DAYS
EOF
>> "$TMP_ENV"

JOB_NAME="fi_mets_${SIM_NAME}_${ISPD_THRESH}"
QSUB_CMD="qsub -N ${JOB_NAME} -v ENVFILE=$TMP_ENV $PBS_SCRIPT"
echo "🔧 $QSUB_CMD"
if [[ "$DRY_RUN" = true ]]; then
    echo "🧪 [DRY RUN] Would submit PBS job with ENVFILE: $TMP_ENV"
else
    eval $QSUB_CMD
fi
