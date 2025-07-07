#!/usr/bin/env python3
import sys, os, pygmt, importlib
mod_path = '/home/581/da1339/AFIM/src/AFIM/src'
sys.path.insert(0, mod_path)
from sea_ice_toolbox import SeaIceToolbox
SI_tools = SeaIceToolbox(sim_name="elps-min")
FI_obs = SI_tools.load_AF2020db(overwrite_zarr=True)
