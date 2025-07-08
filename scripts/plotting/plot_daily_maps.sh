#!/bin/bash
# List of simulations
simulations=("elps-min" "elps-ext" "elps-max" "elps-mid" "ktens-min" "ktens-max" "ktens-ext" "Pstar-max" "Cstar-min" "gi-mid" "gi-nil" "gi-nil-def" "ndte-max" "ndte-max-re-off") 
# Directory for plotting scripts
cd ~/AFIM/src/AFIM/scripts/plotting
# Loop through each simulation and submit a job
for sim_name in "${simulations[@]}"; do
    # Submit a job for each simulation
    qsub -N "plot_${sim_name}" <<EOF
    #!/bin/bash
    #PBS -P jk72
    #PBS -q normalbw
    #PBS -N ${sim_name}-plot-daily-maps
    #PBS -l walltime=24:00:00
    #PBS -l mem=128GB
    #PBS -l ncpus=28
    #PBS -l storage=gdata/gv90+gdata/xp65
    #PBS -o plot_daily_maps_${sim_name}.out
    #PBS -e plot_daily_maps_${sim_name}.err
    # Load modules
    module use /g/data/xp65/public/modules
    module load conda/analysis3-25.05
    source activate
    # Change to the directory containing the plotting script
    cd ~/AFIM/src/AFIM/scripts/plotting
    # Run the Python script for the simulation
    python plot_daily_maps.py --sim_name ${sim_name} --var_names aice hi divu shear strength dvidtt
EOF
done
