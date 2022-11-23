import os
import pdb
import afim
import pygmt
import xarray            as xr
import numpy             as np
import pandas            as pd
import metpy.calc        as mpc
import matplotlib.pyplot as plt

##################################################################################
# Directories, files and switches
plot_gx1   = False
plot_acom2 = False
cice_prep  = afim.cice_prep('/Users/dpath2o/PHD/src/python/afim.json')
sttdt      = '2010-01-01'
months     = 12
years      = 1
stdts      = pd.date_range(sttdt, freq='MS', periods=months*years)
spdts      = pd.date_range(sttdt, freq='M', periods=months*years)
D_02       = cice_prep.D_02
D_SOSE     = cice_prep.D_SOSE
D_MERRA    = cice_prep.D_MERRA
G_type     = cice_prep.G_type
D_graph    = cice_prep.D_graph
Facout     = os.path.join(D_02,'model_input','ACCESS-OM2','output','025deg_jra55_iaf_omip2_cycle6','iceh.2010-01-daily.nc')
Ftgrd      = os.path.join(D_02,'model_input','grids',G_type,'g{:s}_cice_tgrid.nc'.format(G_type))
Fugrd      = os.path.join(D_02,'model_input','grids',G_type,'g{:s}_cice_ugrid.nc'.format(G_type))
D_reG      = os.path.join(D_02,'model_input','grids',G_type,'regrid')
D_ic       = os.path.join(cice_prep.D_IC,G_type)

##################################################################################
# Load grid files
Gu = xr.open_dataset(Fugrd)
Gt = xr.open_dataset(Ftgrd)

##################################################################################
# initial condition NON-tripolar grid provided from CICE6
if plot_gx1:
    Bgx1  = xr.open_dataset(cice_prep.F_gx1bath)
    lon   = np.ravel(Bgx1.TLON.values.T)
    lat   = np.ravel(Bgx1.TLAT.values.T)
    ICgx1 = xr.open_dataset(cice_prep.F_gx1ic)
    # plot these initial conditions and save the images in a reasonable place
    D_plt = os.path.join(D_graph,'cice','ic','gx1',sttdt)
    if not os.path.exists(D_plt): os.mkdir(D_plt)
    for i in ICgx1.keys():
        fig = pygmt.Figure()
        fig.basemap(region='g',projection='W180/15c') # Cyl_stere/30/-20/12c
        if np.ndim(ICgx1[i])==2:
            tit_str  = 'CICE{:d}, Init.Cond. {:s}, {:s}'.format(cice_prep.CICE_ver,i,sttdt)
            print(tit_str)
            dat      = np.ravel(ICgx1[i].values.T)
            DAT      = np.ndarray((len(dat),3))
            DAT[:,0] = lon
            DAT[:,1] = lat
            DAT[:,2] = dat
            bmean    = pygmt.blockmean( data=DAT, region="g", spacing="30m" )
            grd      = pygmt.xyz2grd(data=DAT,region='g',spacing=[1,1])
            pygmt.makecpt(cmap="batlow", series=[ICgx1[i].min().values, ICgx1[i].max().values])
            fig.grdimage(grid=grd, region="g", frame="a",transparency=25)
            fig.coast(region="g", land="black",frame=["af", f'WSne+t"{tit_str}"'])
            fig.colorbar(frame='af')
            fig.savefig(os.path.join(D_plt,'{:s}.png'.format(i)))
        elif np.ndim(ICgx1[i])==3:
            for j in range(np.shape(ICgx1[i])[0]):
                tit_str  = 'CICE{:d}, Init.Cond. {:s}, ncat{:d}, {:s}'.format(cice_prep.CICE_ver,i,j,sttdt)
                print(tit_str)
                dat      = np.ravel(ICgx1[i][j,:,:].values.T)
                DAT      = np.ndarray((len(dat),3))
                DAT[:,0] = lon
                DAT[:,1] = lat
                DAT[:,2] = dat
                bmean    = pygmt.blockmean( data=DAT, region="g", spacing="30m" )
                grd      = pygmt.xyz2grd(data=DAT,region='g',spacing=[1,1])
                pygmt.makecpt(cmap="batlow", series=[ICgx1[i].min().values, ICgx1[i].max().values])
                fig.grdimage(grid=grd, region="g", frame="a",transparency=25)
                fig.coast(region="g", land="black",frame=["af", f'WSne+t"{tit_str}"'])
                fig.colorbar(frame='af')
                fig.savefig(os.path.join(D_plt,'{:s}_ncat{:d}.png'.format(i,j)))

##################################################################################
if plot_acom2:
    ICaco2 = xr.open_dataset(Facout).isel(time=0)
    lon    = np.ravel(ICaco2.TLON.values.T)
    lat    = np.ravel(ICaco2.TLAT.values.T)
    # plot these initial conditions and save the images in a reasonable place
    D_plt = os.path.join(D_graph,'ACOM2','output','tx0p25',sttdt)
    if not os.path.exists(D_plt): os.mkdir(D_plt)
    for i in ICaco2.keys():
        fig = pygmt.Figure()
        fig.basemap(region='g',projection='W180/15c') # Cyl_stere/30/-20/12c
        if np.ndim(ICaco2[i])==2:
            tit_str  = 'CICE_5, AC-OM2 out {:s}, {:s}'.format(i,sttdt)
            print(tit_str)
            dat      = np.ravel(ICaco2[i].values.T)
            DAT      = np.ndarray((len(dat),3))
            DAT[:,0] = lon
            DAT[:,1] = lat
            DAT[:,2] = dat
            bmean    = pygmt.blockmean( data=DAT, region="g", spacing="25m" )
            grd      = pygmt.xyz2grd( data=DAT, region='g', spacing="25m" )
            pygmt.makecpt(cmap="batlow", series=[ICaco2[i].min().values, ICaco2[i].max().values])
            fig.grdimage(grid=grd, region="g", frame="a",transparency=25)
            fig.coast(region="g", land="black",frame=["af", f'WSne+t"{tit_str}"'])
            fig.colorbar(frame='af')
            fig.savefig(os.path.join(D_plt,'{:s}.png'.format(i)))
        elif np.ndim(ICaco2[i])==3:
            for j in range(np.shape(ICaco2[i])[0]):
                tit_str  = 'CICE_5, AC-OM2 out {:s}, ncat{:d}, {:s}'.format(i,j,sttdt)
                dat      = np.ravel(ICaco2[i][j,:,:].values.T)
                DAT      = np.ndarray((len(dat),3))
                DAT[:,0] = lon
                DAT[:,1] = lat
                DAT[:,2] = dat
                bmean    = pygmt.blockmean( data=DAT, region="g", spacing="25m" )
                grd      = pygmt.xyz2grd( data=DAT, region='g', spacing="25m" )
                pygmt.makecpt(cmap="batlow", series=[ICaco2[i].min().values, ICaco2[i].max().values])
                fig.grdimage(grid=grd, region="g", frame="a",transparency=25)
                fig.coast(region="g", land="black",frame=["af", f'WSne+t"{tit_str}"'])
                fig.colorbar(frame='af')
                fig.savefig(os.path.join(D_plt,'{:s}_ncat{:d}.png'.format(i,j)))

##################################################################################
# output shape
sttdt = '2013-01-01'
nxny = list(np.shape(Gt.tlat))
ny   = nxny[0]
nx   = nxny[1]
ncat = 5

##################################################################################
# initial condition external files
# ice velocity
# uvel, vvel
F_uvel = os.path.join(D_SOSE,'bsose_i122_2013to2017_1day_SeaIceUvel.nc')
Pin    = os.path.join(F_uvel)
Dout   = os.path.join(D_reG,'SOSE')
Fout   = 'g{Gtype:s}_{field:s}_{yrmo:s}'.format(Gtype=G_type,field='uvel',yrmo=stdts[0].strftime('%Y%M'))
Pout   = os.path.join(Dout,Fout)
afim.regrid_netcdf_with_cdo(Pin, Pout, Fugrd, cdo_options='-f nc -b F64') # remember to get the right grid
uvel   = xr.open_dataset(Pout).sel(time='2013-01-01')
lon    = np.ravel(uvel.lon.values.T)
lat    = np.ravel(uvel.lat.values.T)
tit_str  = "SOSE, uvel, {sttdt}"
D_plt    = os.path.join(D_graph,'SOSE','tx0p25',sttdt)
fig      = pygmt.Figure()
fig.basemap(region='g',projection='W180/15c') # Cyl_stere/30/-20/12c
dat      = np.ravel(uvel.SIuice.values.T)
DAT      = np.ndarray((len(dat),3))
DAT[:,0] = lon
DAT[:,1] = lat
DAT[:,2] = dat
bmean    = pygmt.blockmean( data=DAT, region="g", spacing="25m" )
grd      = pygmt.xyz2grd( data=DAT, region='g', spacing="25m" )
pygmt.makecpt(cmap="batlow", series=[uvel.SIuice.min().values, uvel.SIuice.max().values])
fig.grdimage(grid=grd, region="g", frame="a",transparency=25)
fig.coast(region="g", land="black",frame=["af", f'WSne+t"{tit_str}"'])
fig.colorbar(frame='af')
fig.savefig(os.path.join(D_plt,'{:s}.png'.format(i)))
F_vvel = os.path.join(D_SOSE,'bsose_i122_2013to2017_1day_SeaIceVvel.nc')
Pin    = os.path.join(F_vvel)
Dout   = os.path.join(D_reG,'SOSE')
Fout   = 'g{Gtype:s}_{field:s}_{yrmo:s}'.format(Gtype=G_type,field='vvel',yrmo=stdts[0].strftime('%Y%M'))
Pout   = os.path.join(Dout,Fout)
afim.regrid_netcdf_with_cdo(Pin, Pout, Fugrd, cdo_options='-f nc -b F64') # remember to get the right grid
vvel   = xr.open_dataset(Pout).SIvice.values
            
# scaling factor for shortwave radiation components
# scale_factor
F_sf         = ''
scale_factor = np.ones( (ny,nx) )*-9999.
# incoming shortwave radiation, visible (near IR), direct (diffuse), W/m^2
# swvdr, swvdf, swidr, swidf
F_sw  = os.path.join(D_MERRA,'rad','2013','01','MERRA2_400.tavg1_2d_rad_Nx.20130101.nc4')
Pin   = os.path.join(F_sw)
Dout  = os.path.join(D_reG,'SOSE')
Fout  = 'g{Gtype:s}_{field:s}_{yrmo:s}'.format(Gtype=G_type,field='merra_rad',yrmo=stdts[0].strftime('%Y%M'))
Pout  = os.path.join(Dout,Fout)
afim.regrid_netcdf_with_cdo(Pin, Pout, Ftgrd, cdo_options='-f nc -b F64') # remember to get the right grid
swvdr = xr.open_dataset(Pout).SWGNT.values
swvdf = np.ones( (ny,nx) )*-9999.
swidr = xr.open_dataset(Pout).LWGNT.values
swidf = np.ones( (ny,nx) )*-9999.
# ice–ocean stress, x(y)-dir. (T-cell), N/m^2
# strocnxT, strocnyT
F_strocnxT = os.path.join(D_SOSE,'bsose_i122_2013to2017_1day_oceTAUX.nc')
Pin        = os.path.join(F_strocnxT)
Dout       = os.path.join(D_reG,'SOSE')
Fout       = 'g{Gtype:s}_{field:s}_{yrmo:s}'.format(Gtype=G_type,field='taux',yrmo=stdts[0].strftime('%Y%M'))
Pout       = os.path.join(Dout,Fout)
afim.regrid_netcdf_with_cdo(Pin, Pout, Fugrd, cdo_options='-f nc -b F64') # remember to get the right grid
strocnxT   = xr.open_dataset(Pout).oceTAUX.values
F_strocnyT = os.path.join(D_SOSE,'bsose_i122_2013to2017_1day_oceTAUY.nc')
Pin        = os.path.join(F_strocnyT)
Dout       = os.path.join(D_reG,'SOSE')
Fout       = 'g{Gtype:s}_{field:s}_{yrmo:s}'.format(Gtype=G_type,field='tauy',yrmo=stdts[0].strftime('%Y%M'))
Pout       = os.path.join(Dout,Fout)
afim.regrid_netcdf_with_cdo(Pin, Pout, Fugrd, cdo_options='-f nc -b F64') # remember to get the right grid
strocnyT   = xr.open_dataset(Pout).oceTAUY.values
# internal ice stress, sigma_11 + sigma_22, N/m
# stressp_1-4
F_stressp = ''
stressp_1 = np.ones( (ny,nx) )*-9999.
stressp_2 = np.ones( (ny,nx) )*-9999.
stressp_3 = np.ones( (ny,nx) )*-9999.
stressp_4 = np.ones( (ny,nx) )*-9999.
# internal ice stress, sigma_11 - sigma_22, N/m
# stressm_1-4
F_stressm = ''
stressm_1 = np.ones( (ny,nx) )*-9999.
stressm_2 = np.ones( (ny,nx) )*-9999.
stressm_3 = np.ones( (ny,nx) )*-9999.
stressm_4 = np.ones( (ny,nx) )*-9999.
# internal ice stress, sigma_12, N/m
# stress12_1-4
F_stress12 = ''
stress12_1 = np.ones( (ny,nx) )*-9999.
stress12_2 = np.ones( (ny,nx) )*-9999.
stress12_3 = np.ones( (ny,nx) )*-9999.
stress12_4 = np.ones( (ny,nx) )*-9999.
#
# aicen
F_aicen = os.path.join(D_SOSE,'bsose_i122_2013to2017_1day_SeaIceArea.nc')
Pin     = os.path.join(F_aicen)
Dout    = os.path.join(D_reG,'SOSE')
Fout    = 'g{Gtype:s}_{field:s}_{yrmo:s}'.format(Gtype=G_type,field='taux',yrmo=stdts[0].strftime('%Y%M'))
Pout    = os.path.join(Dout,Fout)
afim.regrid_netcdf_with_cdo(Pin, Pout, Ftgrd, cdo_options='-f nc -b F64') # remember to get the right grid
aicen   = np.ones( (ncat-1,ny,nx) )*-9999.
aicen   = np.insert(aicen, 1, xr.open_dataset(Pout).SIarea.values, axis=0)
# ice extent mask (U-cell)
# iceumask
F_iceumask = ''
iceumask   = xr.open_dataset(Pout).SIarea > 0
# sea surface temperature, degree C
# sst
F_sst = os.path.join(D_SOSE,'bsose_i122_2013to2017_daily_Theta.nc')
# reezing/melting potential, W/m^2
# frzmlt
F_frzmlt = os.path.join(D_SOSE,'bsose_i122_2013to2017_daily_oceFreez.nc')
# day of year that freezing begins
# frz_onset
F_frz_onset = ''
# 




iceumask     = iceumask.values
sst          = xr.open_dataset('/Volumes/ioa02/model_input/grids/0p25/regrid/SOSE/regrid_sst.nc').THETA.values
frzmlt       = xr.open_dataset('/Volumes/ioa02/model_input/grids/0p25/regrid/SOSE/regrid_oceFreez.nc').oceFreez.values
frz_onset    = np.ones( (ny,nx) )*-9999.
fsnow        = np.ones( (ny,nx) )*-9999.
tlat         = tgrid.tlat.values
tlon         = tgrid.tlon.values
ulat         = ugrid.ulat.values
ulon         = ugrid.ulon.values


vicen        = np.ones( (ncat-1,ny,nx) )*-9999.
vicen        = np.insert(vicen, 1, xr.open_dataset('/Volumes/ioa02/model_input/grids/0p25/regrid/SOSE/regrid_SIheff.nc').SIheff.values, axis=0)
vsnon        = np.ones( (ncat-1,ny,nx) )*-9999.
vsnon        = np.insert(vsnon, 1, xr.open_dataset('/Volumes/ioa02/model_input/grids/0p25/regrid/SOSE/regrid_SIhsnow.nc').SIhsnow.values, axis=0)
Tsfcn        = np.ones( (ncat,ny,nx) )*-9999.
iage         = np.ones( (ncat,ny,nx) )*-9999.
FY           = np.ones( (ncat,ny,nx) )*-9999.
alvl         = np.ones( (ncat,ny,nx) )*-9999.
vlvl         = np.ones( (ncat,ny,nx) )*-9999.
apnd         = np.ones( (ncat,ny,nx) )*-9999.
hpnd         = np.ones( (ncat,ny,nx) )*-9999.
dhs          = np.ones( (ncat,ny,nx) )*-9999.
ffrac        = np.ones( (ncat,ny,nx) )*-9999.
fbrn         = np.ones( (ncat,ny,nx) )*-9999.
first_ice    = np.ones( (ncat,ny,nx) )*-9999.
sice001      = np.ones( (ncat,ny,nx) )*-9999.
sice002      = np.ones( (ncat,ny,nx) )*-9999.
sice003      = np.ones( (ncat,ny,nx) )*-9999.
sice004      = np.ones( (ncat,ny,nx) )*-9999.
sice005      = np.ones( (ncat,ny,nx) )*-9999.
sice006      = np.ones( (ncat,ny,nx) )*-9999.
sice007      = np.ones( (ncat,ny,nx) )*-9999.
qice001      = np.ones( (ncat,ny,nx) )*-9999.
qice002      = np.ones( (ncat,ny,nx) )*-9999.
qice003      = np.ones( (ncat,ny,nx) )*-9999.
qice004      = np.ones( (ncat,ny,nx) )*-9999.
qice005      = np.ones( (ncat,ny,nx) )*-9999.
qice006      = np.ones( (ncat,ny,nx) )*-9999.
qice007      = np.ones( (ncat,ny,nx) )*-9999.
qsno001      = np.ones( (ncat,ny,nx) )*-9999.

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
