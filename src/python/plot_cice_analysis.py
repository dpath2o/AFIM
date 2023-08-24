import afim
import time
import numpy  as np
from datetime import date

cice_anal = afim.cice_analysis('/home/581/da1339/src/python/afim_cice_analysis.json')

c0 = time.process_time()
for run_inc,model_run_name in enumerate(cice_anal.afim_runs.keys()):
    cice_anal.build_data_path(cice_anal.afim_runs[model_run_name]['run_dir'], 'monthly', mf_switch=True)
    cice_anal.load_run_and_time_convert(mf_switch=True, monthly=True)
    print(f'main: finished opening {cice_anal.P_afim}, {(time.process_time()-c0):.3f}')
    for t_inc,dt_obj in enumerate(cice_anal.CICE.time.values):
        for var_inc,var_key in enumerate(cice_anal.plt_mo_dict['var_keys']):
            for region_key in cice_anal.regions_info.keys():
                if region_key=='circumpolar':
                    hemisphere = 'both'
                else:
                    hemisphere = 'sh'
                if not cice_anal.plt_mo_dict['var_vmins'][var_inc]:
                    vmin = cice_anal.CICE[var_key].isel(time=t_inc).min().values
                else:
                    vmin = cice_anal.plt_mo_dict['var_vmins'][var_inc]
                if not cice_anal.plt_mo_dict['var_vmaxs'][var_inc]:
                    vmax = cice_anal.CICE[var_key].isel(time=t_inc).max().values
                else:
                    vmax = cice_anal.plt_mo_dict['var_vmaxs'][var_inc]
                if cice_anal.CICE[var_key].isel(time=t_inc).min().values==cice_anal.CICE[var_key].isel(time=t_inc).max().values:
                    print('vmin = vmax, hence nothing to plot, skipping')
                    continue
                cice_anal.construct_date_string_from_model(t_inc,var_key)
                cice_anal.construct_date_string_from_model(t_inc,var_key)
                tit_str = f"{cice_anal.plt_mo_dict['var_descs'][var_inc]} {cice_anal.fig_dt_str}"
                cice_anal.plot_cice_map(cice_anal.CICE.isel(time=t_inc), t_inc,
                                        which_grid      = cice_anal.plt_mo_dict['var_grids'][var_inc],
                                        var_name        = var_key,
                                        cbar_lab        = cice_anal.plt_mo_dict['var_labs'][var_inc],
                                        cbar_units      = cice_anal.plt_mo_dict['var_units'][var_inc],
                                        model_name      = cice_anal.afim_runs[model_run_name]['long_name'],
                                        vmin            = vmin,
                                        vmax            = vmax,
                                        tit_str         = tit_str,
                                        coast_name      = cice_anal.regions_info[region_key]['coast_name'],
                                        text_str        = cice_anal.regions_info[region_key]['text_string'],
                                        text_loc        = cice_anal.regions_info[region_key]['text_location'],
                                        region          = cice_anal.regions_info[region_key]['region'],
                                        projection      = cice_anal.regions_info[region_key]['projection'],
                                        mean_length_str = 'monthly',
                                        hemisphere      = hemisphere,
                                        spacing         = cice_anal.spacing,
                                        search_radius   = cice_anal.search_radius,
                                        cmap            = cice_anal.plt_mo_dict['var_cmaps'][var_inc],
                                        overwrite       = False)
    #### Now work on geographical statistical threshold figures
    cice_anal.build_data_path(cice_anal.afim_runs[model_run_name]['run_dir'], 'daily', mf_switch=True)
    cice_anal.load_run_and_time_convert(mf_switch=True, monthly=False)
    print(f'main: finished opening {cice_anal.P_afim}, {(time.process_time()-c0):.3f}')
    cice_anal.CICE['icespd']           = np.sqrt(cice_anal.CICE['uvel'] ** 2 + cice_anal.CICE['vvel'] ** 2)
    cice_anal.CICE['zero_vel_days']    = (cice_anal.CICE['icespd'] <= cice_anal.FI_thresh).groupby('time.month').sum(dim='time')
    cice_anal.CICE['frazil_grow_days'] = (cice_anal.CICE['frazil'] >= cice_anal.frazil_thresh).groupby('time.month').sum(dim='time')
    plt_vars = ['zero_vel_days','frazil_grow_days']
    for var_inc,var_key in enumerate(cice_anal.plt_stat_dict['var_keys']):
        for reg_inc,region_key in enumerate(cice_anal.regions_info.keys()):
            if region_key=='circumpolar':
                hemisphere = 'both'
            else:
                hemisphere       = 'sh'
                glon_str         = f"{cice_anal.plt_stat_dict['var_grids'][var_inc-1].upper()}LON"
                glat_str         = f"{cice_anal.plt_stat_dict['var_grids'][var_inc-1].upper()}LAT"
                cice_anal.find_nearest_coordinates(cice_anal.CICE[glat_str].values,
                                                   cice_anal.CICE[glon_str].values,
                                                   cice_anal.regions_info[region_key]['gi_locations'])
            for t_inc in cice_anal.CICE[var_key].month.values:
                if not cice_anal.plt_stat_dict['var_vmins'][var_inc-1]:
                    vmin = cice_anal.CICE[var_key].isel(month=t_inc-1).min().values
                else:
                    vmin = cice_anal.plt_stat_dict['var_vmins'][var_inc-1]
                if not cice_anal.plt_stat_dict['var_vmaxs'][var_inc-1]:
                    vmax = cice_anal.CICE[var_key].isel(month=t_inc-1).max().values
                else:
                    vmax = cice_anal.plt_stat_dict['var_vmaxs'][var_inc-1]
                if cice_anal.CICE[var_key].isel(month=t_inc-1).min().values==cice_anal.CICE[var_key].isel(month=t_inc-1).max().values:
                    print('vmin = vmax, hence nothing to plot, skipping')
                    continue
                cice_anal.construct_date_string_from_model(t_inc,var_key)
                tit_str = f"{cice_anal.plt_stat_dict['var_descs'][var_inc]}, {date(1,t_inc,1).strftime('%b')}"
                cice_anal.plot_cice_map(cice_anal.CICE.isel(month=t_inc-1), t_inc,
                                        which_grid      = cice_anal.plt_stat_dict['var_grids'][var_inc],
                                        var_name        = var_key,
                                        cbar_lab        = cice_anal.plt_stat_dict['var_labs'][var_inc],
                                        cbar_units      = cice_anal.plt_stat_dict['var_units'][var_inc],
                                        model_name      = cice_anal.afim_runs[model_run_name]['long_name'],
                                        vmin            = vmin,
                                        vmax            = vmax,
                                        tit_str         = tit_str,
                                        coast_name      = cice_anal.regions_info[region_key]['coast_name'],
                                        text_str        = cice_anal.regions_info[region_key]['text_string'],
                                        text_loc        = cice_anal.regions_info[region_key]['text_location'],
                                        region          = cice_anal.regions_info[region_key]['region'],
                                        projection      = cice_anal.regions_info[region_key]['projection'],
                                        gi_locs         = cice_anal.GI_locations,
                                        mean_length_str = 'monthly',
                                        hemisphere      = hemisphere,
                                        spacing         = cice_anal.spacing,
                                        search_radius   = cice_anal.search_radius,
                                        cmap            = cice_anal.plt_stat_dict['var_cmaps'][var_inc],
                                        overwrite       = False)
