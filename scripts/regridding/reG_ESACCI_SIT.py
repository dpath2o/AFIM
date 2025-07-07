#!/usr/bin/env python3

def sort_time(ds):
    return ds.sortby("time")

def main():
    import sys, glob
    import xarray as xr
    from pathlib import Path
    # === Setup AFIM environment ===
    mod_path = '/home/581/da1339/AFIM/src/AFIM/src'
    sys.path.insert(0, mod_path)
    from sea_ice_toolbox import SeaIceToolbox
    # === Load model grid from SeaIceToolbox ===
    D_ESACCI = "/g/data/gv90/da1339/SeaIce/ESA_CCI/L2P/envisat/sh"
    SI_tools = SeaIceToolbox(sim_name='elps-min', dt0_str="2002-01-01", dtN_str="2012-12-31", P_log="/g/data/gv90/da1339/logs/ESACCI_SIT_reG.log")
    CICE_all = SI_tools.load_iceh_zarr()
    CICE_SO  = CICE_all.isel(nj=SI_tools.hemisphere_dict['nj_slice'])
    # === Load ESA CCI L2P ===
    files   = sorted(glob.glob(f"{D_ESACCI}/ESACCI-SEAICE-L2P-SITHICK-RA2_ENVISAT-SH-*.nc"))
    ESA_CCI = xr.open_mfdataset(files, combine="nested", concat_dim="time", preprocess=sort_time, parallel=False)
    # === Regrid and Save ===
    ds_out = SI_tools.bin_esa_thickness_to_cice_grid(ESA_CCI, CICE_SO)
    output_path = Path(f"{D_ESACCI}/reG/ESA_CCI_SIT_regridded.zarr")
    ds_out.to_zarr(output_path, mode="w")
    print(f"✅ Regridded ESA CCI saved to: {output_path}")

if __name__ == '__main__':
    from multiprocessing import freeze_support
    freeze_support()
    main()
