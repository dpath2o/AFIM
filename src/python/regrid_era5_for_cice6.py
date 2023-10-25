import os
import afim
import logging
import xarray as xr
import numpy  as np
import xesmf  as xe
from datetime         import datetime
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
file_handler = logging.FileHandler('/home/581/da1339/logs/regrid_era5_for_cice6.log')
formatter    = logging.Formatter('%(asctime)s : %(levelname)s : %(name)s : %(message)s')
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)
yr_str        = '2005'
G_origin      = 'aom2'
G_res         = '0p25'
regrid_method = 'nearest_s2d'
regrid_prdcty = False
F_weights     = "/g/data/jk72/da1339/grids/weights/map_ERA5_{origin:s}_{res:s}_{method:s}.nc".format(origin=G_origin,res=G_res,method=regrid_method)
reuse_weights = False
G_CICE        = xr.open_dataset('/g/data/ik11/inputs/access-om2/input_20200530/cice_025deg/grid.nc')
G_CICE['lat'] = (['ny','nx'],G_CICE.tlat.values*(180/np.pi),{'units':'degrees_north'})
G_CICE['lon'] = (['ny','nx'],G_CICE.tlon.values*(180/np.pi),{'units':'degrees_east'})
logger.debug('Loading and regridding 2t')
ERA5_oG       = xr.open_mfdataset( os.path.join('/g/data/rt52/era5/single-levels/reanalysis','2t'     ,yr_str,'*0101-*.nc') , chunks={'time': 10, 'ny': 100, 'nx': 100} )
LON,LAT       = np.meshgrid(ERA5_oG.longitude.values,ERA5_oG.latitude.values)
G_ERA5        = xr.Dataset(data_vars={'lon':(('ny','nx'),LON),
                                      'lat':(('ny','nx'),LAT)})
rg            = xe.Regridder(G_ERA5, G_CICE, method=regrid_method,  periodic=regrid_prdcty, filename=F_weights, reuse_weights=reuse_weights)
t2m           = rg(ERA5_oG)
t2m_long_name = ERA5_oG['t2m'].attrs['long_name']
t2m_units     = ERA5_oG['t2m'].attrs['units']
logger.debug('Loading and regridding msdwlwrf')
ERA5_oG       = xr.open_mfdataset( os.path.join('/g/data/rt52/era5/single-levels/reanalysis','msdwlwrf',yr_str,'*0101-*.nc') , chunks={'time': 10, 'ny': 100, 'nx': 100} )
lw            = rg(ERA5_oG)
lw_long_name = ERA5_oG['msdwlwrf'].attrs['long_name']
lw_units     = ERA5_oG['msdwlwrf'].attrs['units']
logger.debug('Loading and regridding msdwswrf')
ERA5_oG       = xr.open_mfdataset( os.path.join('/g/data/rt52/era5/single-levels/reanalysis','msdwswrf',yr_str,'*0101-*.nc') , chunks={'time': 10, 'ny': 100, 'nx': 100} )
sw            = rg(ERA5_oG)
sw_long_name = ERA5_oG['msdwswrf'].attrs['long_name']
sw_units     = ERA5_oG['msdwswrf'].attrs['units']
logger.debug('Loading and regridding mtpr')
ERA5_oG       = xr.open_mfdataset( os.path.join('/g/data/rt52/era5/single-levels/reanalysis','mtpr'    ,yr_str,'*0101-*.nc') , chunks={'time': 10, 'ny': 100, 'nx': 100} )
mtpr          = rg(ERA5_oG)
mtpr_long_name = ERA5_oG['mtpr'].attrs['long_name']
mtpr_units     = ERA5_oG['mtpr'].attrs['units']
logger.debug('Loading and regridding 10u')
ERA5_oG       = xr.open_mfdataset( os.path.join('/g/data/rt52/era5/single-levels/reanalysis','10u'     ,yr_str,'*0101-*.nc') , chunks={'time': 10, 'ny': 100, 'nx': 100} )
u10           = rg(ERA5_oG)
u10_long_name = ERA5_oG['u10'].attrs['long_name']
u10_units     = ERA5_oG['u10'].attrs['units']
logger.debug('Loading and regridding 10v')
ERA5_oG       = xr.open_mfdataset( os.path.join('/g/data/rt52/era5/single-levels/reanalysis','10v'     ,yr_str,'*0101-*.nc') , chunks={'time': 10, 'ny': 100, 'nx': 100} )
v10           = rg(ERA5_oG)
v10_long_name = ERA5_oG['v10'].attrs['long_name']
v10_units     = ERA5_oG['v10'].attrs['units']
logger.debug('Loading and regridding 2d')
ERA5_oG       = xr.open_mfdataset( os.path.join('/g/data/rt52/era5/single-levels/reanalysis','2d'      ,yr_str,'*0101-*.nc') , chunks={'time': 10, 'ny': 100, 'nx': 100} )
d2m           = rg(ERA5_oG)
logger.debug('Loading and regridding sp')
ERA5_oG       = xr.open_mfdataset( os.path.join('/g/data/rt52/era5/single-levels/reanalysis','sp'      ,yr_str,'*0101-*.nc') , chunks={'time': 10, 'ny': 100, 'nx': 100} )
sp            = rg(ERA5_oG)
logger.debug('Computing qsat')
qsat          = afim.compute_sfc_qsat(d2m.d2m,sp.sp)
qsat_long_name = "specific humidity"
qsat_units     = "kg/kg"
logger.debug("Data array sizes (GB):")
logger.debug("\t airtmp: ", t2m.t2m.astype(np.single).nbytes / (1024**3))
logger.debug("\t dlwsfc: ", lw.msdwlwrf.astype(np.single).nbytes / (1024**3))
logger.debug("\t glbrad: ", sw.msdwswrf.astype(np.single).nbytes / (1024**3))
logger.debug("\t spchmd: ", qsat.astype(np.single).nbytes / (1024**3))
logger.debug("\t ttlpcp: ", mtpr.mtpr.astype(np.single).nbytes / (1024**3))
logger.debug("\t wndewd: ", u10.u10.astype(np.single).nbytes / (1024**3))
logger.debug("\t wndnwd: ", v10.v10.astype(np.single).nbytes / (1024**3))
d_vars = {"airtmp" : (['time','ny','nx'], t2m.t2m.astype(np.single).data,
                      {'long_name' : t2m_long_name,
                       'units'     : t2m_units}),
          "dlwsfc" : (['time','ny','nx'], lw.msdwlwrf.astype(np.single).data,
                      {'long_name' : lw_long_name,
                       'units'     : lw_units}),
          "glbrad" : (['time','ny','nx'], sw.msdwswrf.astype(np.single).data,
                      {'long_name' : sw_long_name,
                       'units'     : sw_units}),
          "spchmd" : (['time','ny','nx'], qsat.astype(np.single).data,
                      {'long_name' : qsat_long_name,
                       'units'     : qsat_units}),
          "ttlpcp" : (['time','ny','nx'], mtpr.mtpr.astype(np.single).data,
                      {'long_name' : mtpr_long_name,
                       'units'     : mtpr_units}),
          "wndewd" : (['time','ny','nx'], u10.u10.astype(np.single).data,
                      {'long_name' : u10_long_name,
                       'units'     : u10_units}),
          "wndnwd" : (['time','ny','nx'], v10.v10.astype(np.single).data,
                      {'long_name' : v10_long_name,
                       'units'     : v10_units}) }
coords = {"LON"  : (["ny","nx"], G_CICE.lon.data, {'units':'degrees_east'}),
          "LAT"  : (["ny","nx"], G_CICE.lat.data, {'units':'degrees_north'}),
          "time" : (["time"], t2m.time.data)}
attrs = {'creation_date': datetime.now().strftime('%Y-%m-%d %H'),
         'conventions'  : "CCSM data model domain description -- for CICE6 standalone 'JRA55' atmosphere option",
         'title'        : "re-gridded ERA5 for CICE6 standalone atmosphere forcing",
         'source'       : "ERA5, https://doi.org/10.1002/qj.3803, ",
         'comment'      : "source files found on gadi, /g/data/rt52/era5/single-levels/reanalysis",
         'note1'        : "ERA5 documentation, https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation",
         'note2'        : "regridding weight file, {:s}".format(F_weights),
         'note3'        : "re-gridded using ESMF_RegridGenWeights",
         'author'       : 'Daniel P Atwater',
         'email'        : 'daniel.atwater@utas.edu.au'}
ATM   = xr.Dataset(data_vars=d_vars,coords=coords,attrs=attrs)
ATM   = ATM.compute()
logger.debug('Writing out netcdf file ... which can take a very long time')
D_out = '/scratch/jk72/da1339/cice-dirs/input/AFIM/forcing/0p25/ERA5/24XDAILY'
F_out = "era5_for_cice6_{yr:s}.nc".format(yr=yr_str)
P_out = os.path.join(D_out,F_out)
logger.debug(f"Writing to {P_out}")
ATM.to_netcdf(P_out,unlimited_dims=['time'])