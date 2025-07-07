def main():
    import sys, os, time
    from pathlib import Path
    import xarray as xr
    import numpy as np
    import xesmf as xe
    from dask.diagnostics import ProgressBar
    mod_path = '/home/581/da1339/AFIM/src/AFIM/src'
    sys.path.insert(0, mod_path)
    from sea_ice_toolbox import SeaIceToolbox
    P_wghts  = "/g/data/gv90/da1339/grids/weights/CMEMS-ORAS-0p083_to_AOM2-025_tgrid_bil.nc"
    P_reG    = "/g/data/gv90/da1339/SeaIce/CMEMS/0p083/daily/reG/CMEMS-ORAS-SI.zarr"
    reG_vars = ['siconc', 'sithick', 'usi', 'vsi']
    cnk_sz   = 15
    G_CMEMS  = xr.open_dataset("/g/data/gv90/da1339/grids/GLORYS/CMEMS_0p083_grid.nc")
    CMEMS_SI = xr.open_mfdataset("/g/data/gv90/da1339/SeaIce/CMEMS/0p083/daily/*_CMEMS*_org.nc", combine='by_coords', chunks={'time': 1})
    SI_tools = SeaIceToolbox(sim_name="elps-min")
    SI_tools.load_bgrid()
    G_src = {"lon": G_CMEMS["longitude"], "lat": G_CMEMS["latitude"]}
    G_dst = {"lon": SI_tools.G_t["lon"], "lat": SI_tools.G_t["lat"]}
    reG = xe.Regridder(G_src, G_dst,
                       method="bilinear",
                       periodic=True,
                       reuse_weights=True,
                       weights=P_wghts)
    n_T = CMEMS_SI.sizes["time"]
    n_bat = n_T // cnk_sz + (1 if n_T % cnk_sz else 0)
    SI_tools.logger.info(f"📦 Processing {n_T} timesteps in {n_bat} batches of size {cnk_sz}")
    for i in range(n_bat):
        t_loop_start = time.time()
        t0 = i * cnk_sz
        t1 = min(t0 + cnk_sz, n_T)
        times = CMEMS_SI.isel(time=slice(t0, t1)).time.values
        SI_tools.logger.info(f"\n🔁 Batch {i+1}/{n_bat}: {times[0]} to {times[-1]}")
        regridded_vars = {}
        for var in reG_vars:
            t_var_start = time.time()
            SI_tools.logger.info(f"  ↪ Regridding {var}...")
            sliced = CMEMS_SI[var].isel(time=slice(t0, t1)).load()
            regridded = reG(sliced)
            regridded["time"] = times
            regridded_vars[var] = regridded
            SI_tools.logger.info(f"{time.time() - t_var_start:.2f}s")
        usi, vsi = regridded_vars["usi"], regridded_vars["vsi"]
        ispd = (usi**2 + vsi**2)**0.5
        ispd["time"] = times
        ispd.attrs = {"long_name": "ice speed", "units": "m/s"}
        regridded_vars["ispd"] = ispd
        ds_out = xr.Dataset(regridded_vars).chunk({"time": cnk_sz})
        SI_tools.logger.info("  💾 Writing to Zarr... ")
        t_write_start = time.time()
        with ProgressBar():
            if i == 0:
                ds_out.to_zarr(P_reG, mode='w', consolidated=False)
            else:
                ds_out.to_zarr(P_reG, mode='a', append_dim='time', consolidated=False)
        SI_tools.logger.info(f"  💾 Write complete in {time.time() - t_write_start:.2f}s ✅")
    SI_tools.logger.info("🎉 All batches completed successfully.")

if __name__ == '__main__':
    from multiprocessing import freeze_support
    freeze_support()  # safe on all platforms, even if not frozen
    main()  