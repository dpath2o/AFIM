import afim
import numpy  as np
from datetime import date

dt0            = '2005-01-01'
dtN            = '2007-12-31'
proc_period    = 'daily' #or 'monthly'
P_nsdic        = '/g/data/jk72/da1339/SeaIce/nsdic/G02202_V4/south/seaice_conc_monthly_sh_2005-2015.nc'
P_cices        = [f'/g/data/jk72/da1339/afim_output/20230802_jra55_aom2_10yr/history/{proc_period}/*.nc',
                  f'/g/data/jk72/da1339/afim_output/20230807_jra55_aom2_3yr_bathy/history/{proc_period}/*.nc',
                  f'/g/data/jk72/da1339/afim_output/20230807_jra55_aom2_3yr_FI/history/{proc_period}/*.nc',
                  f'/g/data/jk72/da1339/afim_output/20230813_jra55_aom2_3yr_GI10m/history/{proc_period}/*.nc',
                  f'/g/data/jk72/da1339/afim_output/20230814_jra55_aom2_3yr_GIneg0p2m/history/{proc_period}/*.nc']
ln_olav        = [ 49.1, 46.9, 46.3, 45.9, 45.1, 44.7, 44.3, 43.9, 43.1, 42.7, 42.1, 41.5, 41.0, 40.2, 39.7]
lt_olav        = [-66.5,-66.9,-67.1,-67.3,-67.5,-67.5,-67.6,-67.7,-67.9,-67.7,-67.9,-68.0,-68.3,-68.3,-68.5]
ln_maws        = [ 62.3, 62.5, 62.9, 63.1, 63.5, 64.0, 65.0, 65.2, 65.5, 66.5, 66.7, 67.0, 67.7, 68.1, 70.0, 70.0, 70.3, 70.3, 70.3]
lt_maws        = [-66.8,-66.8,-66.8,-66.9,-67.0,-67.2,-67.0,-67.0,-67.0,-67.0,-67.0,-67.1,-67.5,-67.5,-67.5,-67.7,-67.2,-67.5,-67.7]
FI_thresh      = 0.001
spacing        = '15m'
search_radius  = '15m'
cmap           = 'cmocean/amp'
cice_labels    = ['CICE6-baseline','CICE6-bathy','CICE6-FI-noGI','CICE6-GI-10m','CICE6-GI-neg0p2m']
regions        = ['prince olav coast', 'mawson coast', 'sabrina coast']
projections    = {'prince olav coast' : 'S45/-90/25c',
                  'mawson coast'      : 'S67/-90/25c',
                  'sabrina coast'     : 'S119/-90/25c'}
regions        = {'prince olav coast' : (35,50,-71,-65),
                  'mawson coast'      : (60,75,-69,-65),
                  'sabrina coast'     : (112,127,-69,-65)}
coast_names    = {'prince olav coast' : 'olav',
                  'mawson coast'      : 'mawson',
                  'sabrina coast'     : 'sabrina'}
text_strings   = {'prince olav coast' : 'Lutzow-Holm Bay',
                  'mawson coast'      : 'Cape Darnley',
                  'sabrina coast'     : 'Sabrina Coast'}
text_locations = {'prince olav coast' : [40,-70.5],
                  'mawson coast'      : [68,-68.5],
                  'sabrina coast'     : [119,-68.5]}

for i,model_run_name in enumerate(cice_labels):
    CICE_SH               = afim.cice_load_and_mask(P_cices[i], dt0=dt0, dtN=dtN, mf_switch=True, which_grid='u', daily_or_monthly='daily')
    speed                 = np.sqrt(CICE_SH['uvel'] ** 2 + CICE_SH['vvel'] ** 2)
    vel0_days             = speed <= FI_thresh
    vel0_days_monthly_sum = vel0_days.groupby('time.month').sum(dim='time')
    for region in regions:
        if region=='prince olav coast':
            gi_lons, gi_lats = afim.find_nearest_coordinates(CICE_SH['TLAT'].values, CICE_SH['TLON'].values, lt_olav, ln_olav)
        elif region=='mawson coast':
            gi_lons, gi_lats = afim.find_nearest_coordinates(CICE_SH['TLAT'].values, CICE_SH['TLON'].values, lt_maws, ln_maws)
        # iceberg_indices = get_iceberg_grid_indices(lon_maws, lat_maws, CICE_SH['TLON'].values, CICE_SH['TLAT'].values)
        # GI_mask = create_iceberg_mask(iceberg_indices, CICE_SH['TLON'].values.shape)
        for month in vel0_days_monthly_sum.month.values:
            afim.plot_zero_ice_speed_day_count(vel0_days_monthly_sum.sel(month=month), #GI_mask,
                                               proj           = projections[region],
                                               reg            = regions[region],
                                               cmap           = cmap,
                                               spac           = spacing,
                                               sr             = search_radius,
                                               cst_name       = coast_names[region],
                                               text_str       = text_strings[region],
                                               text_loc       = text_locations[region],
                                               month_label    = date(1,month,1).strftime('%b'), 
                                               model_run_name = model_run_name, 
                                               gi_lons        = gi_lons,
                                               gi_lats        = gi_lats)