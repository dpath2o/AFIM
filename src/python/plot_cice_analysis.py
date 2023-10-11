import afim
import time
import numpy  as np
from datetime import date
cice_anal = afim.cice_analysis('/home/581/da1339/AFIM/src/python/afim_cice_analysis.json')
c0        = time.process_time()
for run_inc,model_run_name in enumerate(cice_anal.afim_runs.keys()):
    cice_anal.build_data_path(cice_anal.afim_runs[model_run_name]['run_dir'], 'monthly', mf_switch=True)
    cice_anal.load_run_and_time_convert(mf_switch=True, monthly=True)
    print(f'main: finished opening {cice_anal.P_afim}, {(time.process_time()-c0):.3f}')
    for var_inc,var_key in enumerate(cice_anal.plt_mo_dict['var_keys']):
        for region_key in cice_anal.regions_info.keys():
            for t_inc,dt_obj in enumerate(cice_anal.CICE.time.values):
                cice_anal.plot_cice_map(model_run_name, t_inc,
                                        region_key      = region_key,
                                        var_name        = var_key,
                                        var_inc         = var_inc,
                                        nodata          = '0',
                                        mo_or_ti        = 'time',
                                        overwrite       = False)
