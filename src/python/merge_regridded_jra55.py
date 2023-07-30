import os
import numpy  as np
import xarray as xr
from datetime import datetime

########################################################

reG_meth  = 'bilinear'
reG_perd  = True
reG_reuse = True
vars_map  = {'airtmp':'tas',
             'spchmd':'huss',
             'wndewd':'uas',
             'wndnwd':'vas',
             'ttlpcp':'prsn',#hmmmmm, what about prsa?
             'dlwsfc':'rlds',
             'glbrad':'rsds'}
D_reG     = '/scratch/jk72/da1339/cice-dirs/input/AFIM/forcing/0p25/JRA55/regridded_raw'
D_frcg    = '/scratch/jk72/da1339/cice-dirs/input/AFIM/forcing/0p25/JRA55/8XDAILY/'
F_G_CICE  = '/g/data/jk72/da1339/grids/aom2_cice_0p25_deg.nc'
YR0       = 2005
YRN       = 2019
G_CICE    = xr.open_dataset(F_G_CICE)

########################################################

for yr in np.arange(YR0,YRN):
    airtmp = xr.open_dataset('{D:s}/jra55do_v1p5_{var:s}_{yr:s}.nc'.format(D=D_reG,var=vars_map['airtmp'],yr=str(yr)))
    dlwsfc = xr.open_dataset('{D:s}/jra55do_v1p5_{var:s}_{yr:s}.nc'.format(D=D_reG,var=vars_map['dlwsfc'],yr=str(yr)))
    glbrad = xr.open_dataset('{D:s}/jra55do_v1p5_{var:s}_{yr:s}.nc'.format(D=D_reG,var=vars_map['glbrad'],yr=str(yr)))
    spchmd = xr.open_dataset('{D:s}/jra55do_v1p5_{var:s}_{yr:s}.nc'.format(D=D_reG,var=vars_map['spchmd'],yr=str(yr)))
    ttlpcp = xr.open_dataset('{D:s}/jra55do_v1p5_{var:s}_{yr:s}.nc'.format(D=D_reG,var=vars_map['ttlpcp'],yr=str(yr)))
    wndewd = xr.open_dataset('{D:s}/jra55do_v1p5_{var:s}_{yr:s}.nc'.format(D=D_reG,var=vars_map['wndewd'],yr=str(yr)))
    wndnwd = xr.open_dataset('{D:s}/jra55do_v1p5_{var:s}_{yr:s}.nc'.format(D=D_reG,var=vars_map['wndnwd'],yr=str(yr)))
    dicdat = {"airtmp" : (['time','nj','ni'],airtmp[vars_map['airtmp']].values,
                          {'long_name' :"2 metre temperature",
                           'units'     :"Kelvin",
                           '_FillValue':-2e8}),
              "dlwsfc" : (['time','nj','ni'],dlwsfc[vars_map['dlwsfc']].values,
                          {'long_name':"Mean surface downward long-wave radiation flux",
                           'units'    :"W m**-2",
                           '_FillValue':-2e8}),
              "glbrad" : (['time','nj','ni'],glbrad[vars_map['glbrad']].values,
                          {'long_name':"Mean surface downward short-wave radiation flux",
                           'units'    :"W m**-2",
                           '_FillValue':-2e8}),
              "spchmd" : (['time','nj','ni'],spchmd[vars_map['spchmd']].values,
                          {'long_name':"specific humidity",
                           'units'    :"kg/kg",
                           '_FillValue':-2e8}),
              "ttlpcp" : (['time','nj','ni'],ttlpcp[vars_map['ttlpcp']].values,
                          {'long_name':"Mean total precipitation rate",
                           'units'    :"kg m**-2 s**-1",
                           '_FillValue':-2e8}),
              "wndewd" : (['time','nj','ni'],wndewd[vars_map['wndewd']].values,
                          {'long_name':"10 metre meridional wind component",
                           'units'    :"m s**-1",
                           '_FillValue':-2e8}),
              "wndnwd" : (['time','nj','ni'],wndnwd[vars_map['wndnwd']].values,
                          {'long_name':"10 metre zonal wind component",
                           'units'    :"m s**-1",
                           '_FillValue':-2e8}) }
    coords = {"LON"  : (["nj","ni"],G_CICE.lon.values,{'units':'degrees_east'}),
              "LAT"  : (["nj","ni"],G_CICE.lat.values,{'units':'degrees_north'}),
              "time" : (["time"],airtmp.time.values)}#,'calendar':'gregorian','axis':'T','long_name':'time','standard_name':'time'})}
    attrs = {'creation_date': datetime.now().strftime('%Y-%m-%d %H'),
             'conventions'  : "CCSM data model domain description -- for CICE6 standalone 'JRA55' atmosphere option",
             'title'        : "re-gridded JRA55 for CICE6 standalone atmospheric forcing",
             'source'       : "JRA55-do 1.5.0, doi:10.22033/ESGF/input4MIPs.15017, ",
             'comment'      : "source files found on gadi, /g/data/qv56/inputs",
             'author'       : 'Daniel Patrick Atwater',
             'email'        : 'daniel.atwater@utas.edu.au'}
    enc_dict = {'shuffle':True,'zlib':True,'complevel':5}
    JRA55    = xr.Dataset(data_vars=dicdat,coords=coords,attrs=attrs)
    F_write  = f'{D_frcg}/JRA55_03hr_forcing_tx1_{str(yr):s}.nc'
    if not os.path.exists(F_write):
        write_job = JRA55.to_netcdf(F_write, unlimited_dims=['time'], compute=False, encoding={'glbrad':enc_dict})
        with ProgressBar():
            print(f"Writing to {F_write}")
            write_job.compute()