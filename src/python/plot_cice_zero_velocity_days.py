import afim
import time
import numpy  as np
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
    #cice_anal.CICE['frazil_grow_days'] = (cice_anal.CICE['frazil'] >= cice_anal.frazil_thresh).groupby('time.month').sum(dim='time')
    plt_vars = ['zero_vel_days']#,'frazil_grow_days']
    for var_inc,var_key in enumerate(cice_anal.plt_stat_dict['var_keys']):
        for reg_inc,region_key in enumerate(cice_anal.regions_info.keys()):
            for t_inc in cice_anal.CICE[var_key].month.values:
                print(f'main: plotting {model_run_name}, variable: {var_key}, time increment: {t_inc}, processing time (s): {(time.process_time()-c0):.3f}')
                cice_anal.plot_cice_map(model_run_name, t_inc,
                                        region_key      = region_key,
                                        var_name        = var_key,
                                        var_inc         = var_inc,
                                        nodata          = '0',
                                        mo_or_ti        = 'month',
                                        overwrite       = False)
                print(f'main: completed {model_run_name}, variable: {var_key}, time increment: {t_inc}, processing time (s): {(time.process_time()-c0):.3f}\n\n')