#!/bin/bash
#PBS -P gv90
#PBS -q normalbw
#PBS -l ncpus=28
#PBS -l mem=64GB
#PBS -l walltime=24:00:00
#PBS -l jobfs=10GB
#PBS -l storage=gdata/gv90+gdata/xp65+gdata/jk72
#PBS -N AF2020db
#PBS -l wd

module purge
module use /g/data/xp65/public/modules
module load conda/analysis3-25.05
python3 - <<EOF
import sys, os, pygmt, importlib
mod_path = '/home/581/da1339/AFIM/src/AFIM/src'
sys.path.insert(0, mod_path)
from sea_ice_processor          import SeaIceProcessor
from sea_ice_plotter            import SeaIcePlotter
from grounded_iceberg_processor import GroundedIcebergProcessor
from datetime                   import timedelta, date, datetime
from pathlib                    import Path
from dask.distributed           import Client, LocalCluster
import numpy                    as np
import pandas                   as pd
import xarray                   as xr

sim_name        = "elps-min"
dt0_str         = "2000-01-01"
dtN_str         = "2018-12-31"
ispd_thresh     = 5.0e-4
ice_type        = "FI_BT"
ispd_str        = f"{ispd_thresh:.1e}".replace("e-0", "e-")
smooth_FIA_days = 15
overwrite_zarr  = False
overwrite_png   = True
SI_proc         = SeaIceProcessor(sim_name = sim_name,
                                  dt0_str  = dt0_str,
                                  dtN_str  = dtN_str)
FI_obs = SI_proc.filter_AF2020_FI_by_date(dt0_str, dtN_str)
EOF
