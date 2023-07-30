import os
import afim
import xarray as xr
import numpy  as np
import xesmf  as xe
from datetime         import datetime, timedelta
from dask.diagnostics import ProgressBar
#cice_prep = afim.cice_prep(os.path.join(os.path.expanduser('~'), 'src','python', 'afim_on_gadi.json'))
#cice_prep.regrid_era5_for_cice6()
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
print('Loading and regridding 2t')
ERA5_oG       = xr.open_mfdataset( os.path.join('/g/data/rt52/era5/single-levels/reanalysis','2t'     ,yr_str,'*.nc') , chunks={"time":1} )
LON,LAT       = np.meshgrid(ERA5_oG.longitude.values,ERA5_oG.latitude.values)
G_ERA5        = xr.Dataset(data_vars={'lon':(('ny','nx'),LON),'lat':(('ny','nx'),LAT)})
rg            = xe.Regridder(G_ERA5, G_CICE, method=regrid_method,  periodic=regrid_prdcty, filename=F_weights, reuse_weights=reuse_weights)
t2m           = rg(ERA5_oG)
print('Loading and regridding msdwlwrf')
ERA5_oG       = xr.open_mfdataset( os.path.join('/g/data/rt52/era5/single-levels/reanalysis','msdwlwrf',yr_str,'*.nc') , chunks={"time":1} )
lw            = rg(ERA5_oG)
print('Loading and regridding msdwswrf')
ERA5_oG       = xr.open_mfdataset( os.path.join('/g/data/rt52/era5/single-levels/reanalysis','msdwswrf',yr_str,'*.nc') , chunks={"time":1} )
sw            = rg(ERA5_oG)
print('Loading and regridding mtpr')
ERA5_oG       = xr.open_mfdataset( os.path.join('/g/data/rt52/era5/single-levels/reanalysis','mtpr'    ,yr_str,'*.nc') , chunks={"time":1} )
mtpr          = rg(ERA5_oG)
print('Loading and regridding 10u')
ERA5_oG       = xr.open_mfdataset( os.path.join('/g/data/rt52/era5/single-levels/reanalysis','10u'     ,yr_str,'*.nc') , chunks={"time":1} )
u10           = rg(ERA5_oG)
print('Loading and regridding 10v')
ERA5_oG       = xr.open_mfdataset( os.path.join('/g/data/rt52/era5/single-levels/reanalysis','10v'     ,yr_str,'*.nc') , chunks={"time":1} )
v10           = rg(ERA5_oG)
print('Loading and regridding 2d')
ERA5_oG       = xr.open_mfdataset( os.path.join('/g/data/rt52/era5/single-levels/reanalysis','2d'      ,yr_str,'*.nc') , chunks={"time":1} )
d2m           = rg(ERA5_oG)
print('Loading and regridding sp')
ERA5_oG       = xr.open_mfdataset( os.path.join('/g/data/rt52/era5/single-levels/reanalysis','sp'      ,yr_str,'*.nc') , chunks={"time":1} )
sp            = rg(ERA5_oG)
print('Computing qsat')
qsat          = afim.compute_sfc_qsat(d2m.d2m,sp.sp)
print("Data array sizes (GB):")
print("\t airtmp: ", t2m.t2m.astype(np.single).nbytes / (1024**3))
print("\t dlwsfc: ", lw.msdwlwrf.astype(np.single).nbytes / (1024**3))
print("\t glbrad: ", sw.msdwswrf.astype(np.single).nbytes / (1024**3))
print("\t spchmd: ", qsat.astype(np.single).nbytes / (1024**3))
print("\t ttlpcp: ", mtpr.mtpr.astype(np.single).nbytes / (1024**3))
print("\t wndewd: ", u10.u10.astype(np.single).nbytes / (1024**3))
print("\t wndnwd: ", v10.v10.astype(np.single).nbytes / (1024**3))
d_vars = {"airtmp" : (['time','ny','nx'],t2m.t2m.astype(np.single).data,
                      {'long_name' :"2 metre temperature",
                       'units'     :"Kelvin"}),
          "dlwsfc" : (['time','ny','nx'],lw.msdwlwrf.astype(np.single).data,
                      {'long_name':"Mean surface downward long-wave radiation flux",
                       'units'    :"W m**-2"}),
          "glbrad" : (['time','ny','nx'],sw.msdwswrf.astype(np.single).data,
                      {'long_name':"Mean surface downward short-wave radiation flux",
                       'units'    :"W m**-2"}),
          "spchmd" : (['time','ny','nx'],qsat.astype(np.single).data,
                      {'long_name':"specific humidity",
                       'units'    :"kg/kg"}),
          "ttlpcp" : (['time','ny','nx'],mtpr.mtpr.astype(np.single).data,
                      {'long_name':"Mean total precipitation rate",
                       'units'    :"kg m**-2 s**-1"}),
          "wndewd" : (['time','ny','nx'],u10.u10.astype(np.single).data,
                      {'long_name':"10 metre meridional wind component",
                       'units'    :"m s**-1"}),
          "wndnwd" : (['time','ny','nx'],v10.v10.astype(np.single).data,
                      {'long_name':"10 metre zonal wind component",
                       'units'    :"m s**-1"}) }
coords = {"LON"  : (["ny","nx"],G_CICE.lon.data,{'units':'degrees_east'}),
          "LAT"  : (["ny","nx"],G_CICE.lat.data,{'units':'degrees_north'}),
          "time" : (["time"],t2m.time.data)}
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
print('Writing out netcdf file ... which can take a very long time')
D_out = '/scratch/jk72/da1339/cice-dirs/input/AFIM/forcing/0p25/ERA5/8XDAILY'
F_out = "ERA5_{yr:s}_with_jra55do_var_names_reG_{meth:s}_{res:s}_{org:s}.nc".format(yr=yr_str,meth=regrid_method,res=G_res,org=G_origin)
P_out = os.path.join(D_out,F_out)
write = ATM.to_netcdf(P_out,unlimited_dims=['time'],compute=False)
print(f"Writing to {P_out}")
with ProgressBar(): write.compute()