import os
import afim
import time
import cmocean 
import numpy             as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.path   as mpath
import cartopy.crs       as ccrs
import cartopy.feature   as cfeature
cice_anal = afim.cice_analysis('/home/581/da1339/AFIM/src/python/afim_cice_analysis.json')
c0        = time.process_time()
for run_inc,model_run_name in enumerate(cice_anal.afim_runs.keys()):
    cice_anal.build_data_path(cice_anal.afim_runs[model_run_name]['run_dir'], 'daily', mf_switch=True)
    cice_anal.load_run_and_time_convert(mf_switch=True, monthly=False)
    print(f'main: finished opening {cice_anal.P_afim}, {(time.process_time()-c0):.3f}')
    cice_anal.CICE['icespd']        = np.sqrt(cice_anal.CICE['uvel'] ** 2 + cice_anal.CICE['vvel'] ** 2)
    condition                       = cice_anal.CICE['aice'] < cice_anal.aice_thresh
    cice_anal.CICE['icespd']        = cice_anal.CICE['icespd'].where(~condition, 99999)
    cice_anal.CICE['zero_vel_days'] = (cice_anal.CICE['icespd'] <= cice_anal.FI_thresh).groupby('time.month').sum(dim='time')
    for reg_inc,region_key in enumerate(cice_anal.regions_info.keys()):
        for t_inc in cice_anal.CICE['zero_vel_days'].month.values:
            if region_key=='circumpolar':
                proj = ccrs.SouthPolarStereo()
                reg  = [-180, 180, -90, -50]
                hemisphere = 'both'
            elif region_key=="prince olav coast":
                proj = ccrs.Stereographic(central_longitude=42.5, central_latitude=-68)
                reg  = [60, 75, -69, -65]
                hemisphere = 'sh'
            elif region_key=="mawson coast":
                proj = ccrs.Stereographic(central_longitude=67.5, central_latitude=-67)
                reg  = [35, 50, -71, -65]
                hemisphere = 'sh'
            elif region_key=="sabrina coast":
                proj = ccrs.Stereographic(central_longitude=119.5, central_latitude=-65)
                reg  = [112, 127, -68, -63]
                hemisphere = 'sh'
            cice_anal.construct_date_string_from_model(t_inc, 'zero_vel_days')
            cice_anal.construct_graphical_output_directory(D_base          = '',
                                                           model_name      = cice_anal.afim_runs[model_run_name]['long_name'],
                                                           coast_name      = cice_anal.regions_info[region_key]['coast_name'],
                                                           mean_length_str = 'monthly',
                                                           var_name        = 'zero_vel_days')
            cice_anal.construct_graphical_output_filename(dt_str     = cice_anal.fig_dt_str,
                                                          model_name = cice_anal.afim_runs[model_run_name]['long_name'],
                                                          var_name   = 'zero_vel_days',
                                                          hemisphere = hemisphere,
                                                          img_type   = '')
            cice_anal.P_fig = os.path.join(cice_anal.D_fig,'alt',cice_anal.F_fig)
            if not os.path.exists(os.path.join(cice_anal.D_fig,'alt')): os.makedirs(os.path.join(cice_anal.D_fig,'alt'))
            print(f"main: plotting {model_run_name}, variable: zero_vel_days, time increment: {t_inc}, processing time (s): {(time.process_time()-c0):.3f}")
            cice_anal.find_nearest_coordinates(cice_anal.CICE['ULAT'].values, cice_anal.CICE['ULON'].values, cice_anal.regions_info[region_key]['gi_locations'])
            gi_locs = cice_anal.GI_locations
            da      = cice_anal.CICE['zero_vel_days'].sel(month=t_inc)
            fig, ax = plt.subplots(figsize=(15,9),subplot_kw={'projection': proj}) #  
            ax.set_extent(reg, ccrs.PlateCarree())
            print(f'main: calling pcolormesh, processing time (s): {(time.process_time()-c0):.3f}')
            cf_da = ax.pcolormesh(da.ULON.values, da.ULAT.values, da.values, cmap=cmocean.cm.amp, transform=ccrs.PlateCarree())
            if not hemisphere=='both':
                print(f'main: calling scatter, processing time (s): {(time.process_time()-c0):.3f}')
                for gi_loc in gi_locs:
                    ax.scatter(x=gi_loc[0], y=gi_loc[1], color="black", s=1, alpha=0.5, transform=ccrs.PlateCarree())
            cbar = plt.colorbar(cf_da, ax=ax, orientation='vertical')
            cbar.set_label('near-zero velocity days', rotation=90, labelpad=15)
            print(f'main: calling coastlines, processing time (s): {(time.process_time()-c0):.3f}')
            ax.coastlines(resolution='10m')
            ax.add_feature(cfeature.LAND, edgecolor='black', facecolor='#666666')
            var_desc = cice_anal.plt_stat_dict['var_descs'][0]
            tit_str = f"{var_desc} {cice_anal.afim_runs[model_run_name]['long_name']} {cice_anal.fig_dt_str}"
            ax.set_title(tit_str)
            print(f'main: calling savefig, processing time (s): {(time.process_time()-c0):.3f}')
            plt.savefig(cice_anal.P_fig)
            print(f'main: completed {model_run_name}, variable: zero_vel_days, time increment: {t_inc}, processing time (s): {(time.process_time()-c0):.3f}\n\n')
            plt.close(fig)