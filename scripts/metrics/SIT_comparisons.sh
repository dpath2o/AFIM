#!/bin/bash
SIM_NAME=""
START_YEAR=2002
END_YEAR=2012
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --sim_name) SIM_NAME="$2"; shift ;;
        --start_year) START_YEAR="$2"; shift ;;
        --end_year) END_YEAR="$2"; shift ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done
if [ -z "$SIM_NAME" ]; then
    echo "Error: --sim_name is required"
    exit 1
fi
echo "Submitting SIT computation for $SIM_NAME from $START_YEAR to $END_YEAR"
for YEAR in $(seq $START_YEAR $END_YEAR); do
    qsub -v SIM_NAME=$SIM_NAME,YEAR=$YEAR <<EOF
#!/bin/bash
#PBS -l walltime=03:00:00
#PBS -l mem=64GB
#PBS -l ncpus=10
#PBS -j oe
#PBS -q normalbw
#PBS -P gv90
#PBS -l storage=gdata/xp65+gdata/gv90
#PBS -N SIT_${SIM_NAME}_${YEAR}
module use /g/data/xp65/public/modules
module load conda/analysis3-25.05
python3 ~/AFIM/src/AFIM/scripts/metrics/SIT_comparisons.py --sim_name ${SIM_NAME} --year ${YEAR}
EOF

done
