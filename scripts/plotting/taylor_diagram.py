#!/usr/bin/env python3
import os, sys
import numpy             as np
import pandas            as pd
import xarray            as xr
import matplotlib.pyplot as plt

def main():
    mod_path = '/home/581/da1339/AFIM/src/AFIM/src'
    sys.path.insert(0, mod_path)
    from sea_ice_toolbox import SeaIceToolbox
    dt0_str = "1993-01-01"
    dtN_str = "1999-12-31"
    out_dir = "/g/data/gv90/da1339/GRAPHICAL/AFIM/taylor_diagram"
    sim_name = "gi-nil"
    os.makedirs(out_dir, exist_ok=True)
    SI_cice = SeaIceToolbox(sim_name, dt0_str, dtN_str, True, True, True, True)
    _, CICE = SI_cice.load_processed_cice(zarr_CICE=True)
    CICE    = CICE.isel(nj=SI_cice.hemisphere_dict['nj_slice'])
    SI_aom2 = SeaIceToolbox("AOM2-ERA5", dt0_str, dtN_str, True, True, True, True)
    _, AOM2 = SI_aom2.load_processed_cice(zarr_CICE=True)
    AOM2    = AOM2.isel(nj=SI_aom2.hemisphere_dict['nj_slice'])
    OSI_SAF = xr.open_mfdataset("/home/581/da1339/seaice/OSI_SAF/ispd_reG_SH*")
    da_obs  = OSI_SAF.isel(nj=slice(0, 540)).ice_speed * 1000
    da_cice = CICE.isel(nj=slice(0, 540)).ispd_BT
    da_aom2 = AOM2.isel(nj=slice(0, 540)).ispd_BT
    da_obs, da_cice, da_aom2 = SI_cice.align_and_subset(da_obs, da_cice, da_aom2)
    stats   = {"CICE vs OSI-SAF": SI_cice.compute_taylor_stats(da_cice, da_obs),
               "AOM2 vs OSI-SAF": SI_cice.compute_taylor_stats(da_aom2, da_obs),
               "CICE vs AOM2"   : SI_cice.compute_taylor_stats(da_cice, da_aom2)}
    SI_cice.plot_taylor(stats, os.path.join(out_dir, f"{sim_name}_taylor_diagram_seaice_speed.png"))

if __name__ == "__main__":
    from multiprocessing import freeze_support
    freeze_support()
    main()
