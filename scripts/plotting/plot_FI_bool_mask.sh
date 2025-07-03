#PBS -P gv90
#PBS -q normalbw
#PBS -l walltime=24:00:00
#PBS -l mem=64GB
#PBS -l ncpus=8
#PBS -l jobfs=10GB
#PBS -l storage=gdata/gv90+scratch/gv90+gdata/xp65+gdata/jk72
#PBS -N plot_FI_mask
#PBS -o plot_FI_mask.o
#PBS -e plot_FI_mask.e
#PBS -l wd
module purge
module use /g/data/xp65/public/modules
module load conda/analysis3-25.05
python /home/581/da1339/AFIM/src/AFIM/scripts/plot_FI_bool_mask.py