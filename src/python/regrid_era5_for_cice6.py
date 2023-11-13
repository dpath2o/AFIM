#!/usr/bin/env python
# author: dpath2o, daniel.atwater@utas.edu.au, Sep.23
import os, sys, logging, calendar
import xarray         as xr
import numpy          as np
import xesmf          as xe
#import pandas         as pd
#import dask.array     as da
#import dask.dataframe as dd
from dask.diagnostics import ProgressBar
#from dask.distributed import Client, LocalCluster
from datetime         import datetime
# cluster = LocalCluster()
# client = Client(cluster)

#############################################################################
# HELPFUL FUNCTIONS
def get_last_day_of_month(year, month):
    return calendar.monthrange(year, month)[1]
####################################################
def compute_sfc_qsat(d2m, sp):
    """
    Computes specific humidity at 2-meters based on dewpoint and surface pressure.
    
    Parameters:
        d2m (float): Dewpoint temperature at 2 meters.
        sp (float): Surface pressure.
        
    Returns:
        float: Specific humidity at 2-meters.
    """
    Rdry = 287.0597
    Rvap = 461.5250
    a1   = 611.21
    a3   = 17.502
    a4   = 32.19
    T0   = 273.16
    E    = a1 * np.exp(a3 * (d2m-T0) / (d2m-a4) )
    return (Rdry/Rvap) * E / (sp - ( (1-Rdry/Rvap) * E) )

#############################################################################
# LOCAL VARIABLES / PARAMETERS
user_name       = "da1339"
project_name    = "jk72"
D_longterm      = os.path.join("/","g","data",project_name,user_name)
D_scratch       = os.path.join("/","scratch",project_name,user_name)
D_home          = os.path.join("/","home","581",user_name)
D_afim_input    = os.path.join(D_scratch,"afim_input","ERA5","0p25")
D_grids         = os.path.join(D_longterm,"grids")
F_log           = os.path.join(D_home,"logs","regrid_era5_for_cice6.log")
D_ERA5          = os.path.join("/","g","data","rt52","era5","single-levels","reanalysis")
ERA5_vars_D     = [ "2t"    , "msdwlwrf", "msdwswrf", "mtpr"  , "10u"   , "10v"   , "2d"     ]
ERA5_vars_short = [ "t2m"   , "msdwlwrf", "msdwswrf", "mtpr"  , "u10"   , "v10"   , "d2m"    ]
CICE6_var_names = [ "airtmp", "dlwsfc"  , "glbrad"  , "ttlpcp", "wndewd", "wndnwd", "spchmd" ]
F_G_CICE5       = os.path.join(D_grids,"CICE5_0p25_grid_with_t_deg.nc")
G_CICE5         = xr.open_dataset(F_G_CICE5)
F_G_CICE5_SCRIP = os.path.join(D_grids,"CICE5_SCRIP_0p25.nc")
G_G_ERA5        = os.path.join(D_grids,"ERA5_grid.nc")
F_G_ERA5_SCRIP  = os.path.join(D_grids,"ERA5_SCRIP_native.nc")
F_wgt           = os.path.join(D_grids,"weights","map_ERA5_to_CICE5_0p25_patch_extrap_neareststod.nc")
reG_meth        = "patch"
reG_extrap_meth = "neareststod"
yrs             = np.arange(2005,2019)
mos             = np.arange(1,13)
F_ERA5_fmt      = "{var}_era5_oper_sfc_{yr:04d}{mo:02d}01-{yr:04d}{mo:02d}{dyN:02d}.nc"
netcdf_enc      = {'time': dict(unlimited=True)}
# CHUNKING
# I did some testing with various chunking of the dimensions and the following was the fatest
# Any chunking that I attempted to do on latitude and longitude dimensions, either failed
# the regridding or really slowed things down. Likewise if I went to small on the time
# dimension that also really slowed things down.
chunk_in = {'time': 10}

#############################################################################
# LOGGING
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
file_handler = logging.FileHandler(F_log)
formatter    = logging.Formatter('%(asctime)s : %(levelname)s : %(name)s : %(message)s')
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)

#############################################################################
# GRID WEIGHT GENERATION:
# only needs to be un-commented if weights don't exist or changing the weights
# F_wgt = ["ESMF_RegridWeightGen", "-s", F_G_ERA5_SCRIP, "-d", F_G_CICE5_SCRIP, "-w", F_wgt, "-m", reG_meth, "--extrap_method", reG_extrap_meth, "-i"]
# subprocess.run(gen_F_wgt)

#############################################################################
logger.debug("starting loop")
logger.debug("looping order: year, month, variable")
for y in yrs:
    for m in mos:
        dyN = get_last_day_of_month(y,m)
        for k,var in enumerate(ERA5_vars_D):
            var_D   = var
            var_sh  = ERA5_vars_short[k]
            var_out = CICE6_var_names[k]
            D_out   = os.path.join(D_afim_input,reG_meth)
            if not os.path.exists(D_out): os.makedirs(D_out)
            F_out  = f"{var_out}_{y:04d}_{m:02d}.nc"
            P_out  = os.path.join(D_out,F_out)
            if os.path.exists(P_out):
                logger.debug(f"{P_out} exists ... skipping")
                continue
            F_ERA5 = F_ERA5_fmt.format(var=var,yr=y,mo=m,dyN=dyN)
            P_ERA5 = os.path.join(D_ERA5,var,"{y:04d}".format(y=y),F_ERA5)
            logger.debug(f"opening: {P_ERA5} with chunks {chunk_in}")
            ERA5      = xr.open_dataset(P_ERA5).chunk(chunk_in)
            long_name = ERA5[var_sh].attrs['long_name']
            units     = ERA5[var_sh].attrs['units']
            logger.debug(f"creating re-gridding object using XESMF")
            reG = xe.Regridder(ERA5, G_CICE5, reG_meth, weights=F_wgt)
            logger.debug(f"re-gridding {var} using XESMF re-gridder object")
            if var=="2d":
                d2m_reG = reG(ERA5[var_sh])
                F_ERA5  = F_ERA5_fmt.format(var="sp",yr=y,mo=m,dyN=dyN)
                P_ERA5  = os.path.join(D_ERA5,"sp","{:04d}".format(y),F_ERA5)
                logger.debug(f"opening: {P_ERA5} with chunks {chunk_in}")
                SP = xr.open_dataset(P_ERA5).chunk(chunk_in)
                logger.debug(f"re-gridding sp using XESMF re-gridder object")
                sp_reG = reG(SP["sp"])
                logger.debug(f"computing specific humidity from sp and d2m")
                qsat      = compute_sfc_qsat(d2m_reG,sp_reG)
                ds_reG    = xr.DataArray( qsat , dims=("time","ny","nx") )
                long_name = "specific humidity"
                units     = "kg/kg"
            else:
                ds_reG = reG(ERA5[var_sh])
            logger.debug(f"{var_out} byte size: {ds_reG.astype(np.single).nbytes / (1024**3):.2f}GB")
            logger.debug(f"creating dataset for saving to temporary file")
            var_dict = { var_out : (['time','ny','nx'], ds_reG.astype(np.single).data,
                                    {'long_name' : long_name,
                                     'units'     : units})}
            coords = {"LON"  : (["ny","nx"], G_CICE5["lon"].data, {'units':'degrees_east'}),
                      "LAT"  : (["ny","nx"], G_CICE5["lat"].data, {'units':'degrees_north'}),
                      "time" : (["time"], ERA5["time"].data)}
            attrs = {'creation_date': datetime.now().strftime('%Y-%m-%d %H'),
                     'conventions'  : "CCSM data model domain description",
                     'title'        : "re-gridded ERA5 for CICE6 standalone atmosphere forcing",
                     'source'       : "ERA5, https://doi.org/10.1002/qj.3803, ",
                     'comment'      : "source files found on gadi, /g/data/rt52/era5/single-levels/reanalysis",
                     'note1'        : "ERA5 documentation, https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation",
                     'note2'        : "regridding weight file, {:s}".format(F_wgt),
                     'note3'        : "re-gridded using XESMF",
                     'author'       : 'Daniel P Atwater',
                     'email'        : 'daniel.atwater@utas.edu.au'}
            ds_out = xr.Dataset(data_vars=var_dict,coords=coords,attrs=attrs).chunk(chunk_in)
            logger.debug('Wri')
            logger.debug(f"Writing to {P_out}")
            nc_write = ds_out.to_netcdf(P_out,unlimited_dims=['time'],compute=False)
            with ProgressBar():
                logger.debug(f"Writing to {P_out}")
                nc_write.compute()
