import sys, os
import numpy              as np
import pandas             as pd
import xarray             as xr
import matplotlib.pyplot  as plt
from pathlib              import Path
import warnings
warnings.filterwarnings("ignore", message="Sending large graph of size", category=UserWarning, module="distributed.client")

def convert_stats_dict_to_xarray(stats_dict):
    sim_names = list(stats_dict.keys())
    scalar_vars = {}
    for sim in sim_names:
        for metric_group, metrics in stats_dict[sim].items():
            if isinstance(metrics, xr.DataArray):  # skip time series here
                continue
            for metric_name, value in metrics.items():
                var_key = f"{metric_group}_{metric_name}".replace(" ", "_").replace("-", "_")
                if var_key not in scalar_vars:
                    scalar_vars[var_key] = []
                scalar_vars[var_key].append(value)
    scalar_da_dict = {key: xr.DataArray(data=np.array(vals),
                                        dims=["simulation"],
                                        coords={"simulation": sim_names},
                                        name=key) for key, vals in scalar_vars.items()}
    time_series_da = {}
    for var_name in ['SIA', 'SIV']:
        ts_list = []
        for sim in sim_names:
            da = stats_dict[sim].get(var_name)
            if da is None:
                raise ValueError(f"Missing {var_name} for {sim}")
            da = da.expand_dims({"simulation": [sim]})
            da = da.chunk({"time": 180})  # uniform rechunking
            ts_list.append(da)
        time_series_da[var_name] = xr.concat(ts_list, dim="simulation")
    ds = xr.Dataset({**scalar_da_dict, **time_series_da})
    return ds

def main():
    mod_path = '/home/581/da1339/AFIM/src/AFIM/src'
    sys.path.insert(0, mod_path)
    from sea_ice_toolbox import SeaIceToolbox
    P_log                         = Path(Path.home(),"logs","big_stats_run.log")
    SI_tools_initial              = SeaIceToolbox(sim_name = 'elps-min', P_log=P_log) #dummy to load tools
    SI_growth_days                = (110,200)
    SI_retreat_days1              = (330,360)
    SI_retreat_days2              = (1,60)
    PI_stats                      = {}
    FI_stats                      = {}
    sim_names                     = ["gi-max","gi-mid","gi-nil","gi-nil-def","ndte-max","ndte-min","re-evp-off",'ndte-max-re-off'
                                    "Cstar-max","Cstar-min","Pstar-max","Pstar-min","elps-ext","elps-min","elps-mid","elps-max",
                                    "ktens-ext","ktens-max","ktens-min","ktens-nil","ry93","AOM2-ERA5"]
    G_CMEMS                       = xr.open_dataset("/g/data/gv90/da1339/grids/GLORYS/CMEMS_0p25_grid.nc").rename_dims({'x':'longitude','y':'latitude'})
    CMEMS_SI                      = xr.open_mfdataset("/g/data/gv90/da1339/SeaIce/CMEMS/0p25/daily/199*_CMEMS_org.nc")
    G_CMEMS_SO                    = G_CMEMS.isel(latitude=slice(0,340))
    CMEMS_SI_SO                   = CMEMS_SI.isel(latitude=slice(0,340))
    PI_stats['CMEMS-ORAS']        = {}
    PI_stats['ESA_CCI']           = {}
    PI_stats['NSIDC']             = {}
    PI_stats['CMEMS-ORAS']['SIA'] = SI_tools_initial.compute_ice_area(CMEMS_SI_SO['siconc_oras'], G_CMEMS_SO['area'], ice_area_scale=1e6, spatial_dim_names=("latitude","longitude"))
    PI_stats['CMEMS-ORAS']['SIV'] = SI_tools_initial.compute_ice_volume(CMEMS_SI_SO['siconc_oras'], CMEMS_SI_SO['sithick_oras'], G_CMEMS_SO['area'], ice_volume_scale=1e6, spatial_dim_names=("latitude","longitude"))
    PI_stats['CMEMS-ORAS']['SIT'] = SI_tools_initial.compute_ice_thickness(CMEMS_SI_SO['sithick_oras'], CMEMS_SI_SO['siconc_oras'], G_CMEMS_SO['area'], spatial_dim_names=("latitude","longitude"))
    SIT_comps                     = xr.open_mfdataset("/g/data/gv90/da1339/afim_output/elps-min/SIT_elps-min_*.nc")
    PI_stats['ESA_CCI']['SIT']    = SIT_comps['SIT_obs']
    NSIDC                         = SI_tools_initial.compute_NSIDC_metrics()
    PI_stats['NSIDC']['SIA']      = NSIDC['SIA']
    for sim in sim_names:
        print(f"\n\nprocessing {sim}")
        PI_stats[sim]               = {}
        FI_stats[sim]               = {}
        SI_tools                    = SeaIceToolbox(sim_name = sim, 
                                                    dt0_str  = "1994-01-01",
                                                    dtN_str  = "1999-12-31",
                                                    client   = SI_tools_initial.client,
                                                    P_log    = P_log)
        FId, CICE_org               = SI_tools.load_processed_cice( zarr_CICE = True )
        CICE_SO                     = CICE_org.isel(nj=SI_tools.hemisphere_dict['nj_slice'])
        PI_stats[sim]['SIA']        = SI_tools.compute_ice_area(CICE_SO['aice'], CICE_SO['tarea'], ice_area_scale=SI_tools.SIC_scale)
        PI_stats[sim]['SIV']        = SI_tools.compute_ice_volume(CICE_SO['aice'], CICE_SO['hi'], CICE_SO['tarea'])
        PI_stats[sim]['SIT']        = SI_tools.compute_ice_thickness(CICE_SO['hi'], CICE_SO['aice'], CICE_SO['tarea'])
        PI_stats[sim]['SIA_season'] = SI_tools.compute_seasonal_statistics(PI_stats[sim]['SIA'], 
                                                                        growth_range        = SI_growth_days,
                                                                        retreat_early_range = SI_retreat_days1,
                                                                        retreat_late_range1 = SI_retreat_days1,
                                                                        retreat_late_range2 = SI_retreat_days2)
        PI_stats[sim]['SIV_season'] = SI_tools.compute_seasonal_statistics(PI_stats[sim]['SIV'], 
                                                                        growth_range        = SI_growth_days,
                                                                        retreat_early_range = SI_retreat_days1,
                                                                        retreat_late_range1 = SI_retreat_days1,
                                                                        retreat_late_range2 = SI_retreat_days2)
        PI_stats[sim]['SIA_skills'] = SI_tools.compute_skill_statistics(PI_stats[sim]['SIA'], PI_stats['NSIDC']['SIA'])
        FI_bool                     = SI_tools.boolean_fast_ice( FId['FI_mask'] )
        aice_bool                   = CICE_SO['aice'].where(FI_bool)
        hi_bool                     = CICE_SO['hi'].where(FI_bool)
        tarea_bool                  = CICE_SO['tarea'].where(FI_bool)
        P_bool                      = CICE_SO['strength'].where(FI_bool)
        divu_bool                   = CICE_SO['divu'].where(FI_bool)
        dvidtt_bool                 = CICE_SO['dvidtt'].where(FI_bool)
        FI_stats[sim]['FIA']        = SI_tools.compute_ice_area(FI_bool, tarea_bool)
        FI_stats[sim]['FIV']        = SI_tools.compute_ice_volume(FI_bool, hi_bool, tarea_bool)
        FI_stats[sim]['FIT']        = SI_tools.compute_ice_thickness(hi_bool, FI_bool, tarea_bool)
        FI_stats[sim]['FIA_season'] = SI_tools.compute_seasonal_statistics(FI_stats[sim]['FIA'], 
                                                                        growth_range        = SI_growth_days,
                                                                        retreat_early_range = SI_retreat_days1,
                                                                        retreat_late_range1 = SI_retreat_days1,
                                                                        retreat_late_range2 = SI_retreat_days2)
        FI_stats[sim]['FIV_season'] = SI_tools.compute_seasonal_statistics(FI_stats[sim]['FIV'], 
                                                                        growth_range        = SI_growth_days,
                                                                        retreat_early_range = SI_retreat_days1,
                                                                        retreat_late_range1 = SI_retreat_days1,
                                                                        retreat_late_range2 = SI_retreat_days2)
        AF_clim                     = SI_tools.load_AF2020_FIA_summary(start="1994-01-01", end="1999-12-31")
        obs_fia                     = SI_tools.AF2020_clim_to_model_time(FI_stats[sim]['FIA'] ,
                                                                        AF_clim["FIA_clim"].sel(region="circumpolar") )
        FI_stats[sim]['FIA_stats']  = SI_tools.compute_skill_statistics(FI_stats[sim]['FIA'] , obs_fia )
    PI_DS = convert_stats_dict_to_xarray(PI_stats)
    PI_DS.to_zarr("/g/data/gv90/da1339/afim_output/pack_ice.zarr", mode="w")
    FI_DS = convert_stats_dict_to_xarray(FI_stats)
    FI_DS.to_zarr("/g/data/gv90/da1339/afim_output/fast_ice.zarr", mode="w")

if __name__ == "__main__":
    from multiprocessing import freeze_support
    freeze_support()
    main()