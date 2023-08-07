import os
import afim
import numpy as np
#####################################################################################
model_name      = 'CICE6'
grid_res_str    = '0p25'
fig_type        = 'png'
ani_type        = 'gif'
overwrite       = False
D_plt_base      = '/g/data/jk72/da1339/GRAPHICAL'
D_afim_data     = '/g/data/jk72/da1339/afim_output'
D_ani_base      = '/g/data/jk72/da1339/GRAPHICAL/animations'
model_run_date  = ['20230806',
                   '20230807']
mean_length_str = ['monthly',
                   'monthly']
atm_frcg_name   = ['jra55','jra55']
ocn_frcg_name   = ['aom2','aom2']
intervals       = ['-delay 100','-delay 110']
D_afim_out      = ['/g/data/jk72/da1339/afim_output/20230807_jra55_aom2_3yr_bathy/history/monthly/',
                   '/g/data/jk72/da1339/afim_output/20230807_jra55_aom2_3yr_FI/history/monthly/']
var_names       = ['aice_m',
                   'frazil_m',
                   'hi_m',
                   'uvel_m',
                   'congel_m',
                   'snowfrac_m',
                   'Tsfc_m',
                   'iage_m']
cmaps           = ['cmocean/ice',
                   'cmocean/haline',
                   'cmocean/dense',
                   'cmocean/speed',
                   'cmocean/amp',
                   'cmocean/amp',
                   'cmocean/thermal',
                   'cmocean/amp']
cbar_lab        = ['sea ice concentration',
                   'frazil',
                   'ice volume',
                   'ice speed',
                   'congelation ice growth',
                   'snow fraction',
                   'snow/ice surface temperature',
                   'age of ice']
cbar_units      = ['1/100',
                   'cm/day',
                   'cm',
                   'm/s',
                   'cm',
                   'cm',
                   'C',
                   'years']
vmin            = [0,
                   0,
                   0,
                   -.5,
                   0,
                   0,
                   -30,
                   0]
vmax            = [1,
                   10,
                   10,
                   .5,
                   5,
                   1,
                   1,
                   5]
cnt_data        = 0
for D_ in D_afim_out:
      F_l = np.sort([f for f in os.listdir(D_) if f.endswith(".nc")])
      cnt_vars = 0
      for V_ in var_names:
            for F_ in F_l:
                  afim.plot_cice_map(F_name         = os.path.join(D_,F_),
                                     D_plt_base     = D_plt_base,
                                     model_name     = model_name,
                                     img_type       = fig_type,
                                     atm_frcg_name  = atm_frcg_name[cnt_data],
                                     ocn_frcg_name  = ocn_frcg_name[cnt_data],
                                     model_run_date = model_run_date[cnt_data],
                                     var_name       = V_,
                                     cmap           = cmaps[cnt_vars],
                                     cbar_lab       = cbar_lab[cnt_vars],
                                     cbar_units     = cbar_units[cnt_vars],
                                     vmin           = vmin[cnt_vars],
                                     vmax           = vmax[cnt_vars],
                                     overwrite      = overwrite)
            D_seq_plt = afim.directory_name_from_model_parameters(directory_base  = D_plt_base,
                                                                  model_name      = model_name,
                                                                  model_run_date  = model_run_date[cnt_data],
                                                                  mean_length_str = mean_length_str[cnt_data],
                                                                  atm_frcg_name   = atm_frcg_name[cnt_data],
                                                                  ocn_frcg_name   = ocn_frcg_name[cnt_data],
                                                                  var_name        = V_)
            afim.animate_sequenced_figures(D_seq_fig      = D_seq_plt,
                                          D_ani_base      = D_ani_base,
                                          interval        = intervals[cnt_data],
                                          img_in_type     = fig_type,
                                          img_out_type    = ani_type,
                                          var_name        = V_,
                                          mean_length_str = mean_length_str[cnt_data],
                                          model_name      = model_name,
                                          model_run_date  = model_run_date[cnt_data],
                                          atm_frcg_name   = atm_frcg_name[cnt_data],
                                          ocn_frcg_name   = ocn_frcg_name[cnt_data],
                                          grid_res_str    = grid_res_str,
                                          overwrite       = overwrite)
            cnt_vars+=1
      cnt_data+=1