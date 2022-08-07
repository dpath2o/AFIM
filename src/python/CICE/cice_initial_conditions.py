
import numpy  as np
import xarray as xr

tgrid        = xr.open_dataset('/Volumes/ioa02/model_input/grids/0p25/g0p25_cice_tgrid.nc')
ugrid        = xr.open_dataset('/Volumes/ioa02/model_input/grids/0p25/g0p25_cice_ugrid.nc')
nxny         = list(np.shape(tgrid.tlat))
ny           = nxny[0]
nx           = nxny[1]
ncat         = 5

uvel         = xr.open_dataset('/Volumes/ioa02/model_input/grids/0p25/regrid/SOSE/regrid_SIuvel.nc').SIuice.values
vvel         = xr.open_dataset('/Volumes/ioa02/model_input/grids/0p25/regrid/SOSE/regrid_SIvvel.nc').SIvice.values
scale_factor = np.zeros((ny,nx))
swvdr        = xr.open_dataset('/Volumes/ioa02/model_input/grids/0p25/regrid/MERRA/regrid_SWdirect.nc').SWGNT.values
swvdf        = xr.open_dataset('/Volumes/ioa02/model_input/grids/0p25/regrid/MERRA/regrid_SWdiffuse.nc').__xarray_dataarray_variable__.values
swidr        = xr.open_dataset('/Volumes/ioa02/model_input/grids/0p25/regrid/MERRA/regrid_LWdirect.nc').LWGNT.values
swidf        = xr.open_dataset('/Volumes/ioa02/model_input/grids/0p25/regrid/MERRA/regrid_LWdiffuse.nc').__xarray_dataarray_variable__.values
strocnxT     = xr.open_dataset('/Volumes/ioa02/model_input/grids/0p25/regrid/SOSE/regrid_taux.nc').oceTAUX.values
strocnyT     = xr.open_dataset('/Volumes/ioa02/model_input/grids/0p25/regrid/SOSE/regrid_tauy.nc').oceTAUY.values
stressp_1    = np.zeros((ny,nx))
stressp_2    = np.zeros((ny,nx))
stressp_3    = np.zeros((ny,nx))
stressp_4    = np.zeros((ny,nx))
stressm_1    = np.zeros((ny,nx))
stressm_2    = np.zeros((ny,nx))
stressm_3    = np.zeros((ny,nx))
stressm_4    = np.zeros((ny,nx))
stress12_1   = np.zeros((ny,nx))
stress12_2   = np.zeros((ny,nx))
stress12_3   = np.zeros((ny,nx))
stress12_4   = np.zeros((ny,nx))
iceumask     = xr.open_dataset('/Volumes/ioa02/model_input/grids/0p25/regrid/SOSE/regrid_aicen.nc').SIarea > 0
iceumask     = iceumask.values
sst          = xr.open_dataset('/Volumes/ioa02/model_input/grids/0p25/regrid/SOSE/regrid_sst.nc').THETA.values
frzmlt       = xr.open_dataset('/Volumes/ioa02/model_input/grids/0p25/regrid/SOSE/regrid_oceFreez.nc').oceFreez.values
frz_onset    = np.zeros((ny,nx))
fsnow        = np.zeros((ny,nx))
tlat         = tgrid.tlat.values
tlon         = tgrid.tlon.values
ulat         = ugrid.ulat.values
ulon         = ugrid.ulon.values

aicen        = np.zeros( (ncat-1,ny,nx) )
aicen        = np.insert(aicen, 1, xr.open_dataset('/Volumes/ioa02/model_input/grids/0p25/regrid/SOSE/regrid_aicen.nc').SIarea.values, axis=0)
vicen        = np.zeros( (ncat-1,ny,nx) )
vicen        = np.insert(vicen, 1, xr.open_dataset('/Volumes/ioa02/model_input/grids/0p25/regrid/SOSE/regrid_SIheff.nc').SIheff.values, axis=0)
vsnon        = np.zeros( (ncat-1,ny,nx) )
vsnon        = np.insert(vsnon, 1, xr.open_dataset('/Volumes/ioa02/model_input/grids/0p25/regrid/SOSE/regrid_SIhsnow.nc').SIhsnow.values, axis=0)
Tsfcn        = np.zeros( (ncat,ny,nx) )
iage         = np.zeros( (ncat,ny,nx) )
FY           = np.zeros( (ncat,ny,nx) )
alvl         = np.zeros( (ncat,ny,nx) )
vlvl         = np.zeros( (ncat,ny,nx) )
apnd         = np.zeros( (ncat,ny,nx) )
hpnd         = np.zeros( (ncat,ny,nx) )
dhs          = np.zeros( (ncat,ny,nx) )
ffrac        = np.zeros( (ncat,ny,nx) )
fbrn         = np.zeros( (ncat,ny,nx) )
first_ice    = np.zeros( (ncat,ny,nx) )
sice001      = np.zeros( (ncat,ny,nx) )
sice002      = np.zeros( (ncat,ny,nx) )
sice003      = np.zeros( (ncat,ny,nx) )
sice004      = np.zeros( (ncat,ny,nx) )
sice005      = np.zeros( (ncat,ny,nx) )
sice006      = np.zeros( (ncat,ny,nx) )
sice007      = np.zeros( (ncat,ny,nx) )
qice001      = np.zeros( (ncat,ny,nx) )
qice002      = np.zeros( (ncat,ny,nx) )
qice003      = np.zeros( (ncat,ny,nx) )
qice004      = np.zeros( (ncat,ny,nx) )
qice005      = np.zeros( (ncat,ny,nx) )
qice006      = np.zeros( (ncat,ny,nx) )
qice007      = np.zeros( (ncat,ny,nx) )
qsno001      = np.zeros( (ncat,ny,nx) )

# 
IC = ''
IC = xr.Dataset({'uvel'         : (['nj','ni'],uvel),
                 'vvel'         : (['nj','ni'],vvel),
                 'scale_factor' : (['nj','ni'],scale_factor),
                 'swvdr'        : (['nj','ni'],swvdr),
                 'swvdf'        : (['nj','ni'],swvdf),
                 'swidr'        : (['nj','ni'],swidr),
                 'swidf'        : (['nj','ni'],swidf),
                 'strocnxT'     : (['nj','ni'],strocnxT),
                 'strocnyT'     : (['nj','ni'],strocnyT),
                 'stressp_1'    : (['nj','ni'],stressp_1),
                 'stressp_2'    : (['nj','ni'],stressp_2),
                 'stressp_3'    : (['nj','ni'],stressp_3),
                 'stressp_4'    : (['nj','ni'],stressp_4),
                 'stressm_1'    : (['nj','ni'],stressm_1),
                 'stressm_2'    : (['nj','ni'],stressm_2),
                 'stressm_3'    : (['nj','ni'],stressm_3),
                 'stressm_4'    : (['nj','ni'],stressm_4),
                 'stress12_1'   : (['nj','ni'],stress12_1),
                 'stress12_2'   : (['nj','ni'],stress12_2),
                 'stress12_3'   : (['nj','ni'],stress12_3),
                 'stress12_4'   : (['nj','ni'],stress12_4),
                 'iceumask'     : (['nj','ni'],iceumask),
                 'sst'          : (['nj','ni'],sst),
                 'frzmlt'       : (['nj','ni'],frzmlt),
                 'frz_onset'    : (['nj','ni'],frz_onset),
                 'fsnow'        : (['nj','ni'],fsnow),
                 'tlat'         : (['nj','ni'],tlat),
                 'tlon'         : (['nj','ni'],tlon),
                 'ulat'         : (['nj','ni'],ulat),
                 'ulon'         : (['nj','ni'],ulon),
                 'aicen'        : (['ncat','nj','ni'],aicen),
                 'vicen'        : (['ncat','nj','ni'],vicen),
                 'vsnon'        : (['ncat','nj','ni'],vsnon),
                 'Tsfcn'        : (['ncat','nj','ni'],Tsfcn),
                 'iage'         : (['ncat','nj','ni'],iage),
                 'FY'           : (['ncat','nj','ni'],FY),
                 'alvl'         : (['ncat','nj','ni'],alvl),
                 'vlvl'         : (['ncat','nj','ni'],vlvl),
                 'apnd'         : (['ncat','nj','ni'],apnd),
                 'hpnd'         : (['ncat','nj','ni'],hpnd),
                 'dhs'          : (['ncat','nj','ni'],dhs),
                 'ffrac'        : (['ncat','nj','ni'],ffrac),
                 'fbrn'         : (['ncat','nj','ni'],fbrn),
                 'first_ice'    : (['ncat','nj','ni'],first_ice),
                 'sice001'      : (['ncat','nj','ni'],sice001),
                 'sice002'      : (['ncat','nj','ni'],sice002),
                 'sice003'      : (['ncat','nj','ni'],sice003),
                 'sice004'      : (['ncat','nj','ni'],sice004),
                 'sice005'      : (['ncat','nj','ni'],sice005),
                 'sice006'      : (['ncat','nj','ni'],sice006),
                 'sice007'      : (['ncat','nj','ni'],sice007),
                 'qice001'      : (['ncat','nj','ni'],qice001),
                 'qice002'      : (['ncat','nj','ni'],qice002),
                 'qice003'      : (['ncat','nj','ni'],qice003),
                 'qice004'      : (['ncat','nj','ni'],qice004),
                 'qice005'      : (['ncat','nj','ni'],qice005),
                 'qice006'      : (['ncat','nj','ni'],qice006),
                 'qice007'      : (['ncat','nj','ni'],qice007),
                 'qsno001'      : (['ncat','nj','ni'],qsno001)})
IC.to_netcdf('/Users/dpath2o/cice-dirs/input/CICE_data/ic/0p25/cice6_v6_2013-01.nc')
