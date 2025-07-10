import sys, os, pygmt, argparse
import numpy  as np
import pandas as pd
import xarray as xr
from datetime import datetime
from pathlib  import Path
from glob     import glob

def main(year=None):
    mod_path = '/home/581/da1339/AFIM/src/AFIM/src'
    sys.path.insert(0, mod_path)
    from sea_ice_toolbox import SeaIceToolbox
    # --- Setup ---
    sim_name = "elps-min"
    dt0_str = f"{year}-01-01"
    dtN_str = f"{year}-12-31"
    print(f"Initializing SeaIceToolbox for simulation {sim_name}")
    SI_tools = SeaIceToolbox(sim_name            = sim_name,
                             dt0_str             = dt0_str,
                             dtN_str             = dtN_str,
                             P_log               = Path(Path.home(),"logs",f"reG_OSISAF_ispd_{year}.log"),
                            overwrite_zarr       = True,
                            save_new_figs        = False,
                            show_figs            = False,
                            overwrite_saved_figs = False)
    # --- Load OSI-SAF and compute observed ice speed in m/s ---
    SI_tools.logger.info("Loading OSI SAF data...")
    D_OSI_SAF_raw    = Path("/g/data/gv90/da1339/SeaIce/OSI_SAF/ice_drift_455m/sh/")
    F_OSI_SAF_glob   = f"ice_drift_sh_ease2-750_cdr-v1p0_24h-{year}*"
    P_OSI_SAF_raw    = sorted(glob(str(Path(D_OSI_SAF_raw, F_OSI_SAF_glob))))
    OSI_SAF_raw      = xr.open_mfdataset(P_OSI_SAF_raw, combine='by_coords')    
    dt               = (OSI_SAF_raw["time_bnds"].isel(nv=1) - OSI_SAF_raw["time_bnds"].isel(nv=0)) / np.timedelta64(1, 's')
    disp_mag         = np.sqrt(OSI_SAF_raw['dX']**2 + OSI_SAF_raw['dY']**2)  # km
    ispd_osisaf      = (disp_mag * 1000) / dt.broadcast_like(disp_mag)  # m/s
    ispd_osisaf.name = "ispd"
    ispd_osisaf.attrs["units"] = "m s-1"
    ispd_osisaf.attrs["long_name"] = "Sea ice drift speed"
    # --- Load CICE ice speed ---
    SI_tools.logger.info("Loading processed CICE data...")
    _, CICE       = SI_tools.load_processed_cice(zarr_CICE=True)
    CICE_isp      = SI_tools.compute_ice_magnitude_from_ice_components_on_Bgrid(CICE, ivec_type="BT")
    ispd_obs_algn = ispd_osisaf.sel(time=slice(dt0_str, dtN_str))
    ispd_sim_algn = CICE_isp['ispd_BT'].sel(time=slice(dt0_str, dtN_str))
    # --- Loop over time and regrid ---
    obs_list, sim_list, diff_list = [], [], []
    for i in range(len(ispd_obs_algn.time)):
        dt = pd.Timestamp(ispd_obs_algn.time[i].values)
        dt_str = f"{dt:%Y-%m-%d}"
        SI_tools.logger.info(f"Processing {dt_str}")
        # --- OBS ---
        obs_vals = ispd_obs_algn.isel(time=i).values.flatten()
        obs_lat  = ispd_obs_algn.lat.values.flatten()
        obs_lon  = ispd_obs_algn.lon.values.flatten()
        obs_df   = pd.DataFrame({"longitude": obs_lon, "latitude": obs_lat, "z": obs_vals})
        obs_reG  = pygmt.nearneighbor(data=obs_df, region=[0, 360, -90, -50], spacing="0.25/0.25", search_radius="100k")
        obs_grid = xr.DataArray(obs_reG.values, coords={"lat": obs_reG.lat.values, "lon": obs_reG.lon.values, "time": dt}, dims=("lat", "lon"))
        obs_list.append(obs_grid)
        # --- SIM ---
        sim_vals = ispd_sim_algn.isel(time=i).values.flatten()
        sim_lat  = ispd_sim_algn.TLAT.values.flatten()
        sim_lon  = ispd_sim_algn.TLON.values.flatten()
        sim_df   = pd.DataFrame({"longitude": sim_lon, "latitude": sim_lat, "z": sim_vals})
        sim_reG  = pygmt.nearneighbor(data=sim_df, region=[0, 360, -90, -50], spacing="0.25/0.25", search_radius="50k")
        sim_grid = xr.DataArray(sim_reG.values, coords={"lat": sim_reG.lat.values, "lon": sim_reG.lon.values, "time": dt}, dims=("lat", "lon"))
        sim_list.append(sim_grid)
        # --- DIFF ---
        diff_grid = obs_grid - sim_grid
        diff_list.append(diff_grid)
    # --- Create Dataset ---
    print("Creating final dataset...")
    obs_da  = xr.concat(obs_list, dim="time")
    sim_da  = xr.concat(sim_list, dim="time")
    diff_da = xr.concat(diff_list, dim="time")
    ds      = xr.Dataset({"ispd_obs"  : obs_da,
                          "ispd_sim"  : sim_da,
                          "ispd_diff" : diff_da})
    ds.attrs["description"] = "Daily regridded sea ice speed from OSI-SAF and model, 1994–1999"
    ds.attrs["units"]       = "m/s"
    P_reG                   = Path(SI_tools.D_sim, f"ispd_{sim_name}_and_OSI-SAF_{year}.nc")
    SI_tools.logger.info(f"Saving to {P_reG}")
    ds.to_netcdf(P_reG, unlimited_dims="time")

if __name__ == '__main__':
    from multiprocessing import freeze_support
    freeze_support()
    parser = argparse.ArgumentParser()
    parser.add_argument("--year", type=int, required=True, help="Year to process (e.g. 1994)")
    args = parser.parse_args()
    main(year=args.year)
