'''
This is a python module for functions and objects (classes) that are beneficial
to the Antarctica Fast Ice Modelling that I'm doing for my PhD.

author: dpath2o, daniel.atwater@utas.edu.au, May 2022 
'''

import os
import re
import json
import pdb
import numpy             as np
import xarray            as xr
import xesmf             as xe
import pandas            as pd
import metpy.calc        as mpc
from datetime            import datetime, timedelta
from dask.distributed    import Client
from dask.diagnostics    import ProgressBar

############################################################################
# globals
ncrcat = os.path.join('/','apps','nco','5.0.5','bin','ncrcat')

################################################################################################################################################
######################################################## NON-CLASS FUNCTIONS #################################################################
################################################################################################################################################
def compute_sfc_qsat(d2m, sp):
    '''
    compute specific humidity at 2-metres based on dewpoint and surface pressure
    '''
    Rdry = 287.0597
    Rvap = 461.5250
    a1   = 611.21
    a3   = 17.502
    a4   = 32.19
    T0   = 273.16
    E    = a1 * np.exp(a3 * (d2m-T0) / (d2m-a4) )
    return (Rdry/Rvap) * E / (sp - ( (1-Rdry/Rvap) * E) )

############################################################################
def compute_ocn_sfc_slope(eta, dxdy_array, direction='x', grid_scale_factor=100):
    '''
    grid_scale_factor assumes dxdy_array is in cm and eta is m
    '''
    if direction=='x':
        dhdx                    = eta / (dxdy_array/grid_scale_factor)
        dhdx.attrs['units']     = 'meters'
        dhdx.attrs['long_name'] = 'sea surface slope in x-direction'
        return dhdx
    if direction=='y':
        dhdy                    = eta / (dxdy_array/grid_scale_factor)
        dhdy.attrs['units']     = 'meters'
        dhdy.attrs['long_name'] = 'sea surface slope in y-direction'
        return dhdy
    
############################################################################
def compute_ocn_heat_flux_at_depth(rho_D,cp_D,D,F_net,dTdt_D,
                                   time_unit_to_seconds=3600):
    '''
    Returns the ocean heat flux at depth in W/m^2
    
    Requires:
        rho, rho_D, density at depth, in kg/m^2
        cp, cp_D, heat capacity at depth, in J/(kg*C)
        depth, D, in metres
        dTdt_D, time derivative of temperature at depth,
                assumes in C/hr, but if time derivative is different
                then the option to provide the correct unit of time in
                seconds is provided, (i.e. default is 3600 seconds)
        F_net, atmospheric heat flux at surface, W/m^2
        
    Unit analysis on assumption of time derivative:
        
        W/m^2 - (kg/m^3) * (J/(kg*C)) * m * (C/hr)
        
        reduces to:
        
        W/m^2 - J/(m^2 * hr)
        
        or
        
        (J * m^-2 * s^-1) - (J * m^-2 * 3600*s^-1)
    
    '''
    # cp_D is in J/(kg*C) and W is equal to J 
    return F_net - (rho_D*cp_D*D*dTdt_D)/time_unit_to_seconds

############################################################################
def compute_ocn_pressure_at_depth(D,latitude):
    '''
    Calculates pressure in dbars from depth in meters.
    
    REFERENCE:
    Saunders, P.M. 1981
    "Practical conversion of Pressure to Depth"
    Journal of Physical Oceanography, 11, 573-574
    
    '''
    x  = np.sin(abs(latitude)*np.pi/180)
    c1 = 5.92E-3+(x**2)*5.25E-3;
    return ((1-c1)-np.sqrt(((1-c1)**2)-(8.84E-6*D)))/4.42E-6;

############################################################################
def compute_ocn_secant_bulk_modulus(S,T,P):
    '''
     Secant Bulk Modulus (K) of Sea Water using Equation of state 1980. 
     UNESCO polynomial implementation.
     
     P in dbar
     T in Celsius
     S in PSU
     
     REFERENCES:
      Fofonoff, P. and Millard, R.C. Jr
      Unesco 1983. Algorithms for computation of fundamental properties of 
      seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
      Eqn.(15) p.18
      
      Millero, F.J. and  Poisson, A.
      International one-atmosphere equation of state of seawater.
      Deep-Sea Res. 1981. Vol28A(6) pp625-629.
    '''
    # pressure dbar to atmospheric
    P = P/10
    
    #--------------------------------------------------------------------
    # Pure water terms of the secant bulk modulus at atmos pressure.
    # UNESCO eqn 19 p 18
    h3 = -5.77905e-7
    h2 =  1.16092e-4
    h1 =  1.43713e-3
    h0 =  3.239908
    AW = h0 + (h1 + (h2 + h3*T)*T)*T

    k2 =  5.2787e-8
    k1 = -6.12293e-6
    k0 =  8.50935e-5
    BW = k0 + (k1 + k2*T)*T;

    e4 = -5.155288e-5
    e3 =  1.360477e-2
    e2 = -2.327105
    e1 =  1.484206e2
    e0 =  1.96522e4
    KW = e0 + (e1 + (e2 + (e3 + e4*T)*T)*T)*T
    
    #--------------------------------------------------------------------
    # SEA WATER TERMS OF SECANT BULK MODULUS AT ATMOS PRESSURE.
    j0 =  1.91075e-4
    i2 = -1.6078e-6
    i1 = -1.0981e-5
    i0 =  2.2838e-3

    SR = np.sqrt(S);
    A  = AW + (i0 + (i1 + i2*T)*T + j0*SR)*S
    
    m2 =  9.1697e-10
    m1 =  2.0816e-8
    m0 = -9.9348e-7
    
    #eqn.18
    B  = BW + (m0 + (m1 + m2*T)*T)*S
    
    f3 = -6.1670e-5
    f2 =  1.09987e-2
    f1 = -0.603459
    f0 = 54.6746
    g2 = -5.3009e-4
    g1 =  1.6483e-2
    g0 =  7.944e-2
    
    #eqn.16
    K0 = KW + (f0 + (f1 + (f2 + f3*T)*T)*T + (g0 + (g1 + g2*T)*T)*SR)*S

    #eqn.15
    return K0 + (A + B*P)*P

############################################################################
def compute_ocn_density_standard_mean_ocean_water(T):
    '''
    Denisty of Standard Mean Ocean Water (Pure Water) using EOS 1980. 
    Returns kg/m^3
    
    Requires input, T, temperature, in Celsius
    
    REFERENCES:
        Fofonoff, P. and Millard, R.C. Jr
        Unesco 1983. Algorithms for computation of fundamental properties of 
        seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
        UNESCO 1983 p17  Eqn(14)
         
        Millero, F.J & Poisson, A.
        International one-atmosphere equation of state for seawater.
        Deep-Sea Research Vol28A No.6. 1981 625-629.    Eqn (6)
    '''
    a0 =  9.99842594e2
    a1 =  6.793952e-2
    a2 = -9.095290e-3
    a3 =  1.001685e-4
    a4 = -1.120083e-6
    a5 =  6.536332e-9
    
    return a0 + (a1 + (a2 + (a3 + (a4 + a5*T)*T)*T)*T)*T

############################################################################
def compute_ocn_density_at_surface(S,T):
    '''
    Density of Sea Water at atmospheric pressure using UNESCO 1983 (EOS 1980) polynomial.
    Returns kg/m^3
    
    Requires, S (salinity) to be in practical salinity units (psu), and, T (temperature)
    to be in Celsius
    
    REFERENCES:
         Unesco 1983. Algorithms for computation of fundamental properties of 
         seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
         UNESCO 1983 p17
         
         Millero, F.J & Poisson, A.
         International one-atmosphere equation of state for seawater.
         Deep-Sea Research Vol28A No.6. 1981 625-629.
    '''
    # UNESCO 1983 eqn(13) p17.
    b0 =  8.24493e-1
    b1 = -4.0899e-3
    b2 =  7.6438e-5
    b3 = -8.2467e-7
    b4 =  5.3875e-9
    c0 = -5.72466e-3
    c1 =  1.0227e-4
    c2 = -1.6546e-6
    d0 =  4.8314e-4
    d_smow = compute_ocn_density_standard_mean_ocean_water(T)
    
    return d_smow + (b0 + (b1 + (b2 + (b3 + b4*T)*T)*T)*T)*S + (c0 + (c1 + c2*T)*T)*S*np.sqrt(S) + d0*(S**2)

############################################################################
def compute_ocn_density_at_depth(S,T,D,latitude):
    '''
    Density of Sea Water using UNESCO 1983 (EOS 80) polynomial.
    Returns kg/m^3
    
    Requires:
        S (salinity) to be in practical salinity units (psu)
        T (temperature) to be in Celsius
        D (depth) to be in metres
        latitude to be in degrees
    
    REFERENCES:
        Fofonoff, P. and Millard, R.C. Jr
        Unesco 1983. Algorithms for computation of fundamental properties of 
        seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
        UNESCO 1983 p17  Eqn(14)
         
        Millero, F.J & Poisson, A.
        International one-atmosphere equation of state for seawater.
        Deep-Sea Research Vol28A No.6. 1981 625-629.    Eqn (6)
    '''
    # compute pressure at depth and convert to bars
    P    = compute_ocn_pressure_at_depth(D,latitude)/10 
    rho0 = compute_ocn_density_at_surface(S,T)
    K    = compute_ocn_secant_bulk_modulus(S,T,P)
    
    return rho0/(1-P/K)

############################################################################
def compute_ocn_heat_capacity_at_depth(S,T,D,latitude):
    '''
    Heat Capacity of Sea Water using UNESCO 1983 polynomial.
    Returns [J kg^-1 C^-1]
    
    Requires:
        S (salinity) to be in practical salinity units (psu)
        T (temperature) to be in Celsius
        D (depth) to be in metres
        latitude to be in degrees
    
    REFERENCES:
        Fofonoff, P. and Millard, R.C. Jr
        Unesco 1983. Algorithms for computation of fundamental properties of 
        seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
        UNESCO 1983 p17  Eqn(14)
         
        Millero, F.J & Poisson, A.
        International one-atmosphere equation of state for seawater.
        Deep-Sea Research Vol28A No.6. 1981 625-629.    Eqn (6)
    '''
    # compute pressure at depth and convert to bars
    P    = compute_ocn_pressure_at_depth(D,latitude)/10
    #-----------------
    # eqn.26, p.32
    a0   = -7.64357
    a1   =   .1072763
    a2   = -1.38385e-3
    b0   =   .1770383
    b1   = -4.07718e-3
    b2   =  5.148e-5
    c0   =  4.2174e3
    c1   = -3.720283
    c2   =   .1412855
    c3   = -2.654387e-3
    c4   = 2.093236e-5
    A1   = c0 + c1*T + c2*T**2 + c3*T**3 + c4*T**4
    A2   = (a0 + a1*T + a2*T**2)*S
    A3   = (b0 + b1*T + b2*T**2)*S*np.sqrt(S)
    Cp0  = A1 + A2 + A3
    #-----------------
    # eqn.28, p.33
    d0   = -4.9592e-1
    d1   =  1.45747e-2
    d2   = -3.13885e-4
    d3   =  2.0357e-6
    d4   =  1.7168e-8
    e0   =  2.4931e-4
    e1   = -1.08645e-5
    e2   =  2.87533e-7
    e3   = -4.0027e-9
    e4   =  2.2956e-11
    f0   = -5.422e-8
    f1   =  2.6380e-9
    f2   = -6.5637e-11
    f3   =  6.136e-13
    B1   = (d0 + d1*T + d2*T**2 + d3*T**3 + d4*T**4)*P
    B2   = (e0 + e1*T + e2*T**2 + e3*T**3 + e4*T**4)*P**2
    B3   = (f0 + f1*T + f2*T**2 + f3*T**3)*P**3
    dCp0 =  B1 + B2 + B3
    #-----------------
    # eqn.29, p.34
    d0   =  4.9247e-3
    d1   = -1.28315e-4
    d2   =  9.802e-7
    d3   =  2.5941e-8
    d4   = -2.9179e-10
    e0   = -1.2331e-4
    e1   = -1.517e-6
    e2   =  3.122e-8
    f0   = -2.9558e-6
    f1   =  1.17054e-7
    f2   = -2.3905e-9
    f3   =  1.8448e-11
    g0   =  9.971e-8
    h0   =  5.540e-10
    h1   = -1.7682e-11
    h2   =  3.513e-13
    j1   = -1.4300e-12
    S_3  = S*np.sqrt(S)
    C1   = ((d0 + d1*T + d2*T**2 + d3*T**3 + d4*T**4)*S + (e0 + e1*T + e2*T**2)*S_3)*P
    C2   = ((f0 + f1*T + f2*T**2 + f3*T**3)*S + g0*S_3)*P**2
    C3   = ((h0 + h1*T + h2*T**2)*S + j1*T*S_3)*P**3
    dCp  = C1 + C2 + C3
    
    return (Cp0 + dCp0 + dCp)

############################################################################
def compute_ocn_density_at_depth_alternative_method(S,T,D,latitude):
    '''
    '''
    # depth to pressure
    P = compute_ocn_pressure_at_depth(D,latitude)
    # standar ocean water
    a0       =   999.842594
    a1       =     6.793953e-2
    a2       =    -9.095290e-3
    a3       =     1.001685e-4
    a4       =    -1.120083e-6
    a5       =     6.536332e-9
    rho_SMOW = a0 + a1*T + a2*T**2 + a3*T**3 + a4*T**4 +a5*T**5
    b0       =     8.2449e-1
    b1       =    -4.0899e-3
    b2       =     7.6438e-5
    b3       =    -8.2467e-7
    b4       =     5.3875e-9
    c0       =    -5.7246e-3
    c1       =     1.0227e-4
    c2       =    -1.6546e-6
    d0       =     4.8314e-4
    B1       = b0 + b1*T + b2*T**2 + b3*T**3 + b4*T**4
    C1       = c0 + c1*T + c2*T**2
    rho_p0   = rho_SMOW + B1*S + C1*S**1.5 + d0*S**2
    e0       = 19652.21
    e1       =   148.4206
    e2       =    -2.327105
    e3       =     1.360477e-2
    e4       =    -5.155288e-5
    f0       =    54.674600
    f1       =    -0.603459
    f2       =     1.099870e-2
    f3       =    -6.167e-5
    g0       =     7.9440e-2
    g1       =     1.6483e-2
    g2       =    -5.3009e-4
    G1       = g0 + g1*T + g2*T**2
    F1       = f0 + f1*T + f2*T**2 + f3*T**3
    Kw       = e0 + e1*T + e2*T**2 + e3*T**3 + e4*T**4
    K0       = Kw + F1*S + G1*S**1.5
    h0       =     3.23990
    h1       =     1.43713e-3
    h2       =     1.16092e-4
    h3       =    -5.77905e-7
    i0       =     2.28380e-3
    i1       =    -1.09810e-5
    i2       =    -1.60780e-6
    j0       =     1.91075e-4
    k0       =     8.50935e-5
    k1       =    -6.12293e-6
    k2       =     5.27870e-8
    m0       =    -9.9348e-7
    m1       =     2.0816e-8
    m2       =     9.1697e-10
    Bw       = k0 + k1*T + k2*T**2
    B2       = Bw + (m0 + m1*T + m2*T**2)*S
    Aw       = h0 + h1*T + h2*T**2 + h3*T**3
    A1       = Aw + (i0 + i1*T + i2*T**2)*S + j0*S**1.5
    K        = K0 + A1*P + B2*P**2
    return ( rho_p0 / (1 - ( P / K )) )

############################################################################
def compute_diffuse_sfc_em(em_sfc, em_toa):
    '''
    a simple function to compute the diffuse short/long wave EM at surface
    '''
    k_t = em_sfc / em_toa
    return 0.952 - 1.041 * np.exp(np.exp((2.3 - 4.702*k_t)))

############################################################################
def compute_u2_from_u10(u10):
    '''
    a simple function transfer/compute wind speed at 2-metres from reported wind speed at 10-metres
    '''
    return (u10 * 4.87) / np.log((67.8 * 10) - 5.42)

############################################################################
def compute_sfc_airrho(t2m, d2m, sp):
    '''
    compute atmospheric air density at the surface based on temperature, dewpoint and surface pressure using metpy calculation
    '''
    RH  = mpc.relative_humidity_from_dewpoint(t2m,d2m)
    r   = mpc.mixing_ratio_from_relative_humidity(sp, t2m, RH)
    return mpc.density(sp, t2m, r)

############################################################################
def compute_sfc_qsat(d2m, sp):
    '''
    compute specific humidity at 2-metres based on dewpoint and surface pressure
    '''
    Rdry = 287.0597
    Rvap = 461.5250
    a1   = 611.21
    a3   = 17.502
    a4   = 32.19
    T0   = 273.16
    E    = a1 * np.exp(a3 * (d2m-T0) / (d2m-a4) )
    return (Rdry/Rvap) * E / (sp - ( (1-Rdry/Rvap) * E) )

############################################################################
def ocean_temp_from_theta(p0, p):
    '''
    compute ocean temperature from potential temperature
    p0 :: reference pressure [db]
    p  :: pressure at level [db]

    '''
    dp = p0 - p

############################################################################
def xesmf_regrid_dataset(DS_src, DS_dst, F_DS_na,
                         method='bilinear',
                         periodic=True,
                         reuse_weights=True,
                         F_weights='',
                         multi_file=False):
    '''
    '''
    print("regridding file: ",F_DS_na)
    rg = xe.Regridder(DS_src, DS_dst, 
                      method=method, 
                      periodic=periodic, 
                      filename=F_weights, 
                      reuse_weights=reuse_weights)
    if multi_file:
        DS_na = xr.open_mfdataset(F_DS_na, parallel=True, chunks={'time':1})
    else:
        DS_na = xr.open_dataset(F_DS_na)
    DS_rg = rg(DS_na)
    return DS_rg

##############################################################################################################################################
############################################################### CLASSES ######################################################################
##############################################################################################################################################

# DIRECTORIES
class cice_prep:
    '''
    This class is written to prepare datasets for CICE6 input, either forcing files or initial conditions. It
    works from a supplied JSON file for defining the parameters required.

    External software dependencies:
    - CDO
    *** With this dependency, the user must edit this module's global variable "cdo" (just a few lines above these words),
        and provide the path location to the binary cdo executables.

    Internal software dependencies:
    *** Author recommends using package manager "Anaconda" to download, install and manage the following packages. Version
        dependencies with the following packages are not being rigorously accounted for by the author of this module and
        hence errors and crashes of internal functions and classes of this module may occur do to differences in versions
        in which the author wrote the package and versions that the current user is using.
    - xarray
    - metpy
    - pygmt
    - numpy
    - pandas
    - json
    
    '''
    def __init__(self,FJSON='',**kwargs):
        '''
        Below is a description of the required/expected fields of the JSON file and the contents of those fields.
        Initialising this class in essence gives the user access to those fields and sub-modules that rely
        heavily on the definitions.

        Given a JSON file define the following:

        CICE version        :: 'CICE_ver'
                            :: version number of CICE. Currently on CICE6 is supported
                            :: default is 6 
        Start date          :: 'start_date'
                            :: the date that is used to determine when to start looking for data
                            :: default is 2010-01-01
        Number of Months    :: 'n_months'
                            :: the number of dates from the start date with which to consider biulding
                            :: a period within a year. For instance if one was just considering June, July
                            :: and August, then they would have a start date of 2010-06-01 and set this
                            :: value to 3
                            :: defualt is 12
        Number of Years     :: much like the number of months, but this defines the duration
                            :: default is 2
        CICE grid resolution:: this is string that defines the grid resolution of CICE and is strictly 
                            :: used in naming of files and directories
                            :: default is "0p1"
        Interpolation method:: "regrid_interp_type"
        for re-gridding     :: this is a string that is passed to ESMF_RegridWeightGen for defining the
                            :: type of re-gridding method used for interpolation. See 
                            :: ESMF_RegridWeightGen for more details on methods of interpolation. This 
                            :: is used in the method cice_prep.esmf_generate_weights(), however, that
                            :: method does have an over-ride option 'interp_type' that when provided
                            :: will over-ride the setting in this JSON file.
                            :: default is "bilinear"
        Directories:
           The following inputs are useful directories
           base data        :: "D_data"
                            :: defines the base directory where data will be output
           BRAN data        :: "D_BRAN", location of native BRAN climatology data
           ERA5 data        :: "D_ERA5", location of native ERA5 climatology data
           Model Input      :: "D_modin", probably no longer necessary as it was used in previous 
                            :: version of this module
           Graphical        :: "D_graph", this is figures/plots/animations are stored
           ACCESS-OM2 output:: "D_access-om2_out", this defines which ACCESS-OM2 cycle data will be 
                            :: in the concatenation routines and then subsequently forcing the CICE
                            :: stand-alone model
        Files:
           The followin are specific files that are used by the module
           Bathymetry       :: "F_amo2bath", this file is no presently used in any of the methods,
                            :: but is highly relevant and might want to be used for plotting
           CICE grid        :: "F_G_CICE", this contains the spherical grid (lat/lon in decimal degrees) 
                            :: and is an AFIM-compatible-version of the grid that stand-alone CICE uses
                            :: when running. It is used here in the  generation of re-gridding weights and 
                            :: re-gridding of climatology datasets
           BRAN grid        :: "F_G_BRAN", this contains the spherical grid (lat/lon in decimal degrees) 
                            :: of the BRAN climatology and is used in the generation of re-gridding weights
                            :: and re-gridding of climatology datasets
           ERA5 grid        :: "F_G_ERA5", this contains the spherical grid (lat/lon in decimal degrees) 
                            :: of the ERA5 climatology and is used in the generation of re-gridding weights
                            :: and re-gridding of climatology datasets
           BRAN weights     :: "F_BRAN_weights", once created via the AFIM python module
                            :: cice_prep.esmf_generate_weights this is the name of the file used in the 
                            :: method cice_prep.regrid_climatology() for BRAN climatology
           ERA5 weights     :: "F_ERA5_weights", once created via the AFIM python module
                            :: cice_prep.esmf_generate_weights this is the name of the file used in the 
                            :: method cice_prep.regrid_climatology() for ERA5 climatology
           
        '''
        # read in the JSON file
        self.FJSON = FJSON
        with open(FJSON) as f:
            PARAMS = json.load(f)
        # miscellaneous
        self.CICE_ver           = PARAMS['CICE_ver']
        self.start_date         = PARAMS['start_date']
        self.n_months           = PARAMS['n_months']
        self.n_years            = PARAMS['n_years']
        self.ncrcat_opts        = PARAMS["ncrcat_opts"]
        # dataset specific
        self.acom2_model_run    = PARAMS["acom2_model_run"]
        self.acom2_start_date   = PARAMS['acom2_start_date']
        self.acom2_chunking     = PARAMS['acom2_chunking']
        self.ERA5_chunking      = PARAMS['ERA5_chunking']
        self.ERA5_start_date    = PARAMS['ERA5_start_date']
        self.ERA5_regen_weights = PARAMS['ERA5_regen_weights']
        self.ERA5_regrid_method = PARAMS['ERA5_regrid_method']
        self.ERA5_periodicity   = PARAMS['ERA5_periodicity']
        self.BRAN_chunking      = PARAMS['BRAN_chunking']
        self.BRAN_start_date    = PARAMS['BRAN_start_date']
        self.BRAN_regen_weights = PARAMS['BRAN_regen_weights']
        self.BRAN_regrid_method = PARAMS['BRAN_regrid_method']
        self.BRAN_periodicity   = PARAMS['BRAN_periodicity']
        # grid
        self.G_res        = PARAMS['G_res']
        # directories
        self.D01          = PARAMS['D01']
        self.D_data       = os.path.join(PARAMS['D01'],'data')
        self.D_grids      = os.path.join('/','g','data','jk72','da1339','grids')
        self.D_weights    = os.path.join(self.D_grids,'weights')
        self.D_BRAN       = PARAMS['D_BRAN']
        self.D_ERA5       = PARAMS['D_ERA5']
        self.D_CICE       = PARAMS['D_CICE']
        self.D_graph      = PARAMS['D_graph']
        self.D_modin      = PARAMS['D_modin']
        self.D_frcg       = PARAMS['D_frcg']
        self.D_acom2_out = PARAMS['D_access-om2_out']
        self.D_CICE_IC    = os.path.join(PARAMS['D_CICE'],'input','AFIM','ic',self.G_res)
        self.D_CICE_force = os.path.join(PARAMS['D_CICE'],'input','AFIM','forcing',self.G_res)
        self.D_CICE_grid  = os.path.join(PARAMS['D_CICE'],'input','AFIM','grid',self.G_res)
        # files
        self.F_gx1ic         = os.path.join(self.D_CICE_grid,'gx1','iced_gx1_v6.2005-01-01.nc')
        self.F_gx1bath       = os.path.join(self.D_CICE_grid,'gx1','global_gx1.bathy.nc')
        self.F_aom2bath      = PARAMS['F_aom2bath']
        self.F_G_CICE        = PARAMS['F_G_CICE']
        self.F_G_CICE_vars   = PARAMS['F_G_CICE_vars']
        self.F_G_CICE_original = PARAMS['F_G_CICE_original']
        self.F_G_BRAN        = PARAMS['F_G_BRAN']
        self.F_G_ERA5        = PARAMS['F_G_ERA5']
        self.F_ERA5_weights  = PARAMS['F_ERA5_weights']
        self.F_BRAN_weights  = PARAMS['F_BRAN_weights']
        self.F_ERA5_reG_form = PARAMS['F_ERA5_reG_form']
        self.F_ERA5_reG_form_tmp = PARAMS["F_ERA5_reG_form_tmp"]
        self.F_BRAN_reG_form = PARAMS['F_BRAN_reG_form']

    ############################################################################################
    def esmf_generate_weights(self, ds_name='ERA5', weight_file_prefix='access-om2_cice',
                              interp_type='', regrid_options='-p none'):
        '''
        '''
        if not interp_type:
            interp_type = self.regrid_interp_type
        if ds_name=='ERA5':
            P_wgt = os.path.join( self.D_weights, 'map_{from_name:s}_{to_name:s}_{grid_res:s}_{type:s}.nc'.format(from_name=ds_name,
                                                                                                                  to_name=weight_file_prefix,
                                                                                                                  grid_res=self.G_res,
                                                                                                                  type=interp_type))
            str_esmf_gen_weights = 'ESMF_RegridWeightGen -m {type:s} {opts:s} -s {src:s} -d {dst:s} -w {wgt:s}'.format(type=interp_type,
                                                                                                              opts=regrid_options,
                                                                                                              src=self.F_G_ERA5,
                                                                                                              dst=self.F_G_CICE,
                                                                                                              wgt=P_wgt)
        elif ds_name=='BRAN':
            P_wgt = os.path.join( self.D_weights, 'map_{from_name:s}_{to_name:s}_{grid_res:s}_{type:s}.nc'.format(from_name=ds_name,
                                                                                                                  to_name=weight_file_prefix,
                                                                                                                  grid_res=self.G_res,
                                                                                                                  type=interp_type))

            str_esmf_gen_weights = 'ESMF_RegridWeightGen -m {type:s} {opts:s} -s {src:s} -d {dst:s} -w {wgt:s}'.format(type=interp_type,
                                                                                                              opts=regrid_options,
                                                                                                              src=self.F_G_BRAN,
                                                                                                              dst=self.F_G_CICE,
                                                                                                              wgt=P_wgt)
        print("generating weights with the command:\n",str_esmf_gen_weights)
        os.system(str_esmf_gen_weights)
        
    ############################################################################################
    def time_series_day_month_start_or_end(self,start_date='',n_months='',n_years='',month_start=True):
        '''
        '''
        if not n_months:
            n_mos = self.n_months
        else:
            n_mos = n_months
        if not n_years:
            n_yrs = self.n_years
        else:
            n_yrs = n_years
        if not start_date:
            stt_dt = self.start_date
        else:
            stt_dt = start_date
        if month_start:
            freq = 'MS'
        else:
            freq = 'M'
        return pd.date_range(stt_dt, freq=freq, periods=n_mos*n_yrs)

    ############################################################################################
    def define_datetime_object(self,start_date='',year_offset=0):
        '''
        '''
        if not start_date: start_date = self.start_date
        return datetime.strptime(start_date,'%Y-%m-%d').date() + pd.DateOffset(years=year_offset)
    
    ############################################################################################
    def year_string_from_datetime_object(self,dt_obj):
        '''
        '''
        return dt_obj.strftime('%Y')
    
    ############################################################################################
    def define_era5_variable_path(self,var_name='',stt='',stp=''):
        '''
        '''
        if not stt:
            stt = self.time_series_day_month_start_or_end()
        if not stp:
            ptp = self.time_series_day_month_end()
        F = '{var:s}_era5_oper_sfc_{stt:s}-{stp:s}.nc'.format(var=var_name,stt=stt,stp=stp)
        return os.path.join(self.D_ERA5,var_name,self.yr_str,F)
    
    ############################################################################################
    def era5_load_and_regrid(self,era5_var='2t',yr_str='', D_base='',
                             regrid_method='', regrid_periodicity='',cnk_dict='',
                             generate_weights=False, weights_file_name=''):
        '''
        '''
        # optional inputs
        if not D_base:             D_ERA5             = self.D_ERA5
        if not regrid_method:      regrid_method      = self.ERA5_regrid_method
        if not regrid_periodicity: regrid_periodicity = self.ERA5_periodicity
        if not cnk_dict:           cnk_dict           = self.ERA5_chunking
        if not weights_file_name: 
            F_weights = self.F_ERA5_weights
        else:
            F_weights = weights_file_name
        if not yr_str:
            dt_obj = define_datetime_object(start_date=self.ERA5_start_date, year_offset=0)
            yr_str = year_string_from_datetime_object(dt_obj)
        # pull-in the grids
        G_ERA5 = self.era5_grid_prep()
        print("\nSource file looks like this: ",G_ERA5)
        G_CICE = self.cice_grid_prep()
        print("\nDestination file looks like this: ",G_CICE)
        if generate_weights or not os.path.exists(F_weights):
            print("User has requested to generate weights file *or* {:s} does not exist".format(F_weights))
            print("Creating weights can take some time ... sometimes it can take a long time!!")
            print("Intialising xesmf regridder")
            print("Using xesmf method {:s}".format(regrid_method))
            print(f"Periodicity is set to {regrid_periodicity}")
            print("Full path to weight file is: {:s}".format(F_weights))
            rg = xe.Regridder(G_ERA5, G_CICE, method=regrid_method,  periodic=regrid_periodicity, filename=F_weights, reuse_weights=False)
        else:
            print("REUSING EXISTING WEIGHT FILE: {:s}".format(F_weights))
            print("Using xesmf method {:s}".format(regrid_method))
            print(f"Periodicity is set to {regrid_periodicity}")
            rg = xe.Regridder(G_ERA5, G_CICE, method=regrid_method,  periodic=regrid_periodicity, filename=F_weights, reuse_weights=True)
        dat_n = self.era5_load(era5_var=era5_var, D_base=D_ERA5, yr_str=yr_str, cnk_dict=cnk_dict )
        print("regridding ERA5 ",dat_n)
        return rg(dat_n)
    
    ############################################################################################
    def era5_load(self,era5_var='2t',D_base='',yr_str='',cnk_dict=''):
        '''
        '''
        #optional inputs
        if not D_base:
            D_ERA5 = self.D_ERA5
        else:
            D_ERA5 = D_base
        if not cnk_dict: cnk_dict = self.ERA5_chunking
        if not yr_str:
            dt_obj = define_datetime_object(start_date=self.ERA5_start_date, year_offset=0)
            yr_str = year_string_from_datetime_object(dt_obj)
        P_dat = os.path.join(D_ERA5,era5_var,yr_str,'*.nc')
        print("loading ERA5: {:s}".format(P_dat))
        prt_str = "".join(str(key) + str(value) for key, value in cnk_dict.items())
        print("\t with chunking dictionary: {:s}".format(prt_str))
        return xr.open_mfdataset( P_dat , chunks=cnk_dict )
     
    ############################################################################################
    def bran_load_and_regrid(self,bran_var='temp',yr_str='', D_base='',
                             regrid_method='', regrid_periodicity='',cnk_dict='',
                             generate_weights=False, weights_file_name=''):
        '''
        '''
        # optional inputs
        if not D_base:             D_BRAN             = self.D_BRAN
        if not regrid_method:      regrid_method      = self.BRAN_regrid_method
        if not regrid_periodicity: regrid_periodicity = self.BRAN_periodicity
        if not cnk_dict:           cnk_dict           = self.BRAN_chunking
        if not weights_file_name: 
            F_weights = self.F_BRAN_weights
        else:
            F_weights = weights_file_name
        if not yr_str:
            dt_obj = define_datetime_object(start_date=self.BRAN_start_date, year_offset=0)
            yr_str = year_string_from_datetime_object(dt_obj)
        # pull-in the grids
        G_BRAN = self.bran_grid_prep()
        print("\nSource file looks like this: ",G_BRAN)
        G_CICE = self.cice_grid_prep()
        print("\nDestination file looks like this: ",G_CICE)
        if generate_weights or not os.path.exists(F_weights):
            print("User has requested to generate weights file *or* {:s} does not exist".format(F_weights))
            print("Creating weights can take some time ... sometimes it can take a long time!!")
            print("Intialising xesmf regridder")
            print("Using xesmf method {:s}".format(regrid_method))
            print(f"Periodicity is set to {regrid_periodicity}")
            print("Full path to weight file is: {:s}".format(F_weights))
            rg = xe.Regridder(G_BRAN, G_CICE, method=regrid_method,  periodic=regrid_periodicity, filename=F_weights, reuse_weights=False)
        else:
            print("REUSING EXISTING WEIGHT FILE: {:s}".format(F_weights))
            print("Using xesmf method {:s}".format(regrid_method))
            print(f"Periodicity is set to {regrid_periodicity}")
            rg = xe.Regridder(G_BRAN, G_CICE, method=regrid_method,  periodic=regrid_periodicity, filename=F_weights, reuse_weights=True)
        dat_n = self.bran_load(bran_var=bran_var, D_base=D_BRAN, yr_str=yr_str, cnk_dict=cnk_dict )
        print("regridding BRAN ",dat_n)
        return rg(dat_n)

    ############################################################################################
    def bran_load(self,bran_var='temp',D_base='',yr_str='',cnk_dict=''):
        '''
        '''
        #optional inputs
        if not D_base:
            D_BRAN = self.D_BRAN
        else:
            D_BRAN = D_base
        if not cnk_dict: cnk_dict = self.BRAN_chunking
        if not yr_str:
            dt_obj = define_datetime_object(start_date=self.BRAN_start_date, year_offset=0)
            yr_str = year_string_from_datetime_object(dt_obj)
        P_dat = os.path.join(D_BRAN,'ocean_{var:s}_{yr:s}*.nc'.format(var=bran_var,yr=yr_str))
        print("loading BRAN: {:s}".format(P_dat))
        prt_str = "".join(str(key) + str(value) for key, value in cnk_dict.items())
        print("\t with chunking dictionary: {:s}".format(prt_str))
        return xr.open_mfdataset( P_dat , chunks=cnk_dict )
        
    ############################################################################################
    def era5_grid_prep(self):
        '''
        '''
        return xr.open_dataset(self.F_G_ERA5)
    
    ############################################################################################
    def bran_grid_prep(self):
        '''
        '''
        G_BRAN        = xr.open_dataset(self.F_G_BRAN)
        LN,LT         = np.meshgrid(G_BRAN.xt_ocean.values,G_BRAN.yt_ocean.values)
        G_BRAN['lon'] = (['yt_ocean','xt_ocean'],LN,{'units':'degrees_east','_FillValue':-2e8})
        G_BRAN['lat'] = (['yt_ocean','xt_ocean'],LT,{'units':'degrees_north','_FillValue':-2e8})
        G_BRAN        = G_BRAN.drop(('xu_ocean','yu_ocean','hu','kmu','umask','tmask','st_edges_ocean','Time','st_ocean'))
        return G_BRAN

    ############################################################################################
    def cice_grid_prep(self):
        '''
        '''
        G_CICE        = xr.open_dataset(self.F_G_CICE_original)
        G_CICE['lat'] = (['ny','nx'],G_CICE.tlat.values*(180/np.pi),{'units':'degrees_north','_FillValue':-2e8})
        G_CICE['lon'] = (['ny','nx'],G_CICE.tlon.values*(180/np.pi),{'units':'degrees_east','_FillValue':-2e8})
        G_CICE        = G_CICE.drop(('ulat','ulon','tlat','tlon','clon_t','clat_t','clat_u','clon_u','angle','uarea'))
        G_CICE        = G_CICE.assign_coords(xt_ocean=G_CICE['lon'][0,:],yt_ocean=G_CICE['lat'][:,0])
        G_CICE        = G_CICE.set_index(nx='xt_ocean',ny='yt_ocean')
        return G_CICE

    ############################################################################################
    def regrid_era5_for_cice6(self,start_date='', regrid_method='', regrid_periodicity='',
                              cnk_dict='', generate_weights='', weights_file_name=''):
        '''
        '''
        if not start_date:         start_date         = self.ERA5_start_date
        if not regrid_method:      regrid_method      = self.ERA5_regrid_method
        if not regrid_periodicity: regrid_periodicity = self.ERA5_periodicity
        if not cnk_dict:           cnk_dict           = self.ERA5_chunking
        if not generate_weights:   generate_weights   = self.ERA5_regen_weights
        if not weights_file_name: 
            F_weights = self.F_ERA5_weights
        else:
            F_weights = weights_file_name
        for i in range(self.n_years):
            dt_obj = self.define_datetime_object(start_date=start_date, year_offset=i)
            yr_str = self.year_string_from_datetime_object(dt_obj)
            t2m   = self.era5_load_and_regrid(era5_var           = '2t', 
                                          yr_str             = yr_str,
                                          regrid_method      = regrid_method,
                                          regrid_periodicity = regrid_periodicity,
                                          cnk_dict           = cnk_dict,
                                          generate_weights   = generate_weights,
                                          weights_file_name  = F_weights)
            lw    = self.era5_load_and_regrid(era5_var            = 'msdwlwrf', 
                                          yr_str             = yr_str,
                                          regrid_method      = regrid_method,
                                          regrid_periodicity = regrid_periodicity,
                                          cnk_dict           = cnk_dict,
                                          generate_weights   = generate_weights,
                                          weights_file_name  = F_weights)
            sw    = self.era5_load_and_regrid(era5_var            = 'msdwswrf', 
                                          yr_str             = yr_str,
                                          regrid_method      = regrid_method,
                                          regrid_periodicity = regrid_periodicity,
                                          cnk_dict           = cnk_dict,
                                          generate_weights   = generate_weights,
                                          weights_file_name  = F_weights)
            mtpr  = self.era5_load_and_regrid(era5_var            = 'mtpr', 
                                          yr_str             = yr_str,
                                          regrid_method      = regrid_method,
                                          regrid_periodicity = regrid_periodicity,
                                          cnk_dict           = cnk_dict,
                                          generate_weights   = generate_weights,
                                          weights_file_name  = F_weights)
            u10   = self.era5_load_and_regrid(era5_var            = 'u10', 
                                          yr_str             = yr_str,
                                          regrid_method      = regrid_method,
                                          regrid_periodicity = regrid_periodicity,
                                          cnk_dict           = cnk_dict,
                                          generate_weights   = generate_weights,
                                          weights_file_name  = F_weights)
            v10   = self.era5_load_and_regrid(era5_var            = 'v10', 
                                          yr_str             = yr_str,
                                          regrid_method      = regrid_method,
                                          regrid_periodicity = regrid_periodicity,
                                          cnk_dict           = cnk_dict,
                                          generate_weights   = generate_weights,
                                          weights_file_name  = F_weights)
            d2m   = self.era5_load_and_regrid(era5_var            = '2d', 
                                          yr_str             = yr_str,
                                          regrid_method      = regrid_method,
                                          regrid_periodicity = regrid_periodicity,
                                          cnk_dict           = cnk_dict,
                                          generate_weights   = generate_weights,
                                          weights_file_name  = F_weights)
            sp    = self.era5_load_and_regrid(era5_var            = 'sp', 
                                          yr_str             = yr_str,
                                          regrid_method      = regrid_method,
                                          regrid_periodicity = regrid_periodicity,
                                          cnk_dict           = cnk_dict,
                                          generate_weights   = generate_weights,
                                          weights_file_name  = F_weights)
            qsat   = compute_sfc_qsat(d2m.d2m, sp.sp)
            print("Data array sizes (GB):")
            print("\t airtmp: ", t2m.t2m.astype(np.single).nbytes / (1024**3))
            print("\t dlwsfc: ", lw.msdwlwrf.astype(np.single).nbytes / (1024**3))
            print("\t glbrad: ", sw.msdwswrf.astype(np.single).nbytes / (1024**3))
            print("\t spchmd: ", qsat.astype(np.single).nbytes / (1024**3))
            print("\t ttlpcp: ", mtpr.mtpr.astype(np.single).nbytes / (1024**3))
            print("\t wndewd: ", u10.u10.astype(np.single).nbytes / (1024**3))
            print("\t wndnwd: ", v10.v10.astype(np.single).nbytes / (1024**3))
            d_vars = {"airtmp" : (['time','nj','ni'],t2m.t2m.astype(np.single).data,
                                      {'long_name' :"2 metre temperature",
                                       'units'     :"Kelvin",
                                       '_FillValue':-2e8}),
                          "dlwsfc" : (['time','nj','ni'],lw.msdwlwrf.astype(np.single).data,
                                      {'long_name':"Mean surface downward long-wave radiation flux",
                                       'units'    :"W m**-2",
                                       '_FillValue':-2e8}),
                          "glbrad" : (['time','nj','ni'],sw.msdwswrf.astype(np.single).data,
                                      {'long_name':"Mean surface downward short-wave radiation flux",
                                       'units'    :"W m**-2",
                                       '_FillValue':-2e8}),
                          "spchmd" : (['time','nj','ni'],qsat.astype(np.single).data,
                                      {'long_name':"specific humidity",
                                       'units'    :"kg/kg",
                                       '_FillValue':-2e8}),
                          "ttlpcp" : (['time','nj','ni'],mtpr.mtpr.astype(np.single).data,
                                      {'long_name':"Mean total precipitation rate",
                                       'units'    :"kg m**-2 s**-1",
                                       '_FillValue':-2e8}),
                          "wndewd" : (['time','nj','ni'],u10.u10.astype(np.single).data,
                                      {'long_name':"10 metre meridional wind component",
                                       'units'    :"m s**-1",
                                       '_FillValue':-2e8}),
                          "wndnwd" : (['time','nj','ni'],v10.v10.astype(np.single).data,
                                      {'long_name':"10 metre zonal wind component",
                                       'units'    :"m s**-1",
                                       '_FillValue':-2e8}) }
            coords = {"LON"  : (["nj","ni"],G_CICE.lon.data,{'units':'degrees_east'}),
                          "LAT"  : (["nj","ni"],G_CICE.lat.data,{'units':'degrees_north'}),
                          "time" : (["time"],t2m.time.data)}
            attrs = {'creation_date': datetime.now().strftime('%Y-%m-%d %H'),
                         'conventions'  : "CCSM data model domain description -- for CICE6 standalone 'JRA55' atmosphere option",
                         'title'        : "re-gridded ERA5 for CICE6 standalone ocean forcing",
                         'source'       : "ERA5, https://doi.org/10.1002/qj.3803, ",
                         'comment'      : "source files found on gadi, /g/data/rt52/era5/single-levels/reanalysis",
                         'note1'        : "ERA5 documentation, https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation",
                         'note2'        : "regridding weight file, /g/data/jk72/da1339/grids/weights/map_ERA5_access-om2_cice_0p1_bilinear.nc",
                         'note3'        : "re-gridded using ESMF_RegridGenWeights",
                         'author'       : 'Daniel P Atwater',
                         'email'        : 'daniel.atwater@utas.edu.au'}
            ATM = xr.Dataset(data_vars=d_vars,coords=coords,attrs=attrs)
            #write to monthly files
            cnt = 0
            #for i in wks_arr:
            mo0_dates = pd.date_range(start=dt_obj,freq='MS',periods=self.n_months)
            moN_dates = pd.date_range(start=dt_obj,freq='M',periods=self.n_months)
            enc_dict  = {'shuffle':True,'zlib':True,'complevel':5}
            for i in mo0_dates:
                dt0_str   = i.strftime('%Y-%m-%d %H:%M')
                dtN_str   = moN_dates[cnt].strftime('%Y-%m-%d %H:%M')
                yrmo0_str = i.strftime('%Y_%m')
                P_ATM_tmp = os.path.join(self.D_data,'ERA5','monthly',self.F_ERA5_reG_form_tmp.format(dt_str=yrmo0_str))
                if not os.path.exists(P_ATM_tmp):
                    ATM_tmp   = ATM.sel(time=slice(dt0_str,dtN_str))
                    write_job = ATM_tmp.to_netcdf(P_ATM_tmp,unlimited_dims=['time'],compute=False,encoding={'airtmp':enc_dict,
                                                                                                            'dlwsfc':enc_dict,
                                                                                                            'glbrad':enc_dict,
                                                                                                            'spchmd':enc_dict,
                                                                                                            'ttlpcp':enc_dict,
                                                                                                            'wndewd':enc_dict,
                                                                                                            'wndnwd':enc_dict})
                    with ProgressBar():
                        print(f"Writing to {P_ATM_tmp}")
                        write_job.compute()
                cnt+=1
                
    ############################################################################################
    def regrid_bran_for_cice6(self,start_date='', regrid_method='', regrid_periodicity='',
                              cnk_dict='', generate_weights='', weights_file_name=''):
        '''
        '''
        if not start_date:         start_date         = self.BRAN_start_date
        if not regrid_method:      regrid_method      = self.BRAN_regrid_method
        if not regrid_periodicity: regrid_periodicity = self.BRAN_periodicity
        if not cnk_dict:           cnk_dict           = self.BRAN_chunking
        if not generate_weights:   generate_weights   = self.BRAN_regen_weights
        if not weights_file_name: 
            F_weights = self.F_BRAN_weights
        else:
            F_weights = weights_file_name
        for i in range(self.n_years):
            dt_obj = self.define_datetime_object(start_date=start_date, year_offset=i)
            yr_str = self.year_string_from_datetime_object(dt_obj)
            temp   = self.bran_load_and_regrid(bran_var           = 'temp', 
                                          yr_str             = yr_str,
                                          regrid_method      = regrid_method,
                                          regrid_periodicity = regrid_periodicity,
                                          cnk_dict           = cnk_dict,
                                          generate_weights   = generate_weights,
                                          weights_file_name  = F_weights)
            salt   = self.bran_load_and_regrid(bran_var           = 'salt', 
                                          yr_str             = yr_str,
                                          regrid_method      = regrid_method,
                                          regrid_periodicity = regrid_periodicity,
                                          cnk_dict           = cnk_dict,
                                          generate_weights   = generate_weights,
                                          weights_file_name  = F_weights)
            mld    = self.bran_load_and_regrid(bran_var           = 'mld', 
                                          yr_str             = yr_str,
                                          regrid_method      = regrid_method,
                                          regrid_periodicity = regrid_periodicity,
                                          cnk_dict           = cnk_dict,
                                          generate_weights   = generate_weights,
                                          weights_file_name  = F_weights)
            uocn   = self.bran_load_and_regrid(bran_var           = 'u', 
                                          yr_str             = yr_str,
                                          regrid_method      = regrid_method,
                                          regrid_periodicity = regrid_periodicity,
                                          cnk_dict           = cnk_dict,
                                          generate_weights   = generate_weights,
                                          weights_file_name  = F_weights)
            vocn   = self.bran_load_and_regrid(bran_var           = 'v', 
                                          yr_str             = yr_str,
                                          regrid_method      = regrid_method,
                                          regrid_periodicity = regrid_periodicity,
                                          cnk_dict           = cnk_dict,
                                          generate_weights   = generate_weights,
                                          weights_file_name  = F_weights)
            eta    = self.bran_load_and_regrid(bran_var           = 'eta_t', 
                                          yr_str             = yr_str,
                                          regrid_method      = regrid_method,
                                          regrid_periodicity = regrid_periodicity,
                                          cnk_dict           = cnk_dict,
                                          generate_weights   = generate_weights,
                                          weights_file_name  = F_weights)
            sst     = temp.sel(st_ocean=0,method='nearest')
            sss     = salt.sel(st_ocean=0,method='nearest')
            u       = uocn.sel(st_ocean=0,method="nearest")
            v       = vocn.sel(st_ocean=0,method="nearest")
            dhdx    = compute_ocn_sfc_slope(eta, G_CICE, direction='x', grid_scale_factor=100)
            dhdy    = compute_ocn_sfc_slope(eta, G_CICE, direction='y', grid_scale_factor=100)
            print("computing temperature derivative over time")
            dTdt    = dTdt.differentiate("time")
            print("computing atmospheric surface heat flux")
            F_net   = self.compute_era5_net_atmospheric_surface_heat_flux(yr_str             = yr_str,
                                                                              return_daily_mean  = True,
                                                                              regrid_method      = regrid_method,
                                                                              regrid_periodicity = regrid_periodicity,
                                                                              cnk_dict           = self.ERA5_chunking,
                                                                              generate_weights   = generate_weights,
                                                                              weights_file_name  = F_weights)
            F_net = F_net.rename({'nx':'xt_ocean','ny':'yt_ocean'})
            print("computing ocean heat capacity at the mixed layer depth")
            cp_j   = compute_ocn_heat_capacity_at_depth(salt_j, temp_j, mld_j, mld.yt_ocean).astype(np.single)
            print("computing ocean density at the mixed layer depth")
            rho_j  = compute_ocn_density_at_depth(salt_j, temp_j, mld_j, mld.yt_ocean).astype(np.single)
            print("computing ocean heat flux at the mixed layer depth")
            qdp_j  = compute_ocn_heat_flux_at_depth(rho_j, cp_j, mld_j, Fnet_j, dTdt_j)
            qdp_j  = qdp_j.assign_coords(time=temp.time.isel(time=j).values).expand_dims('time').astype(np.single).to_dataset(name='qdp')
            ln    = sst.lon
            lt    = sst.lat
            LN,LT = np.meshgrid(ln,lt)
            time  = sst.Time
            OCN   = xr.Dataset({ "T"    : sst,
                                 "S"    : sss,
                                 "hblt" : mld,
                                 "u"    : u,
                                 "v"    : v,
                                 "dhdx" : dhdx,
                                 "dhdy" : dhdy,
                                 "qdp"  : qdp},
                               coords = { "LON"  : (["ny","nx"],LN),
                                          "LAT"  : (["ny","nx"],LT),
                                          "time" : time.values })
            OCN   = OCN.rename_dims({"ny":"nj","nx":"ni"})
            P_ATM = os.path.join(self.D_frcg,self.F_BRAN_reG_form.format(yr=yr_str))
            OCN.to_netcdf(P_ocn)

    ############################################################################################
    def acom2_compute_qdp(self,start_date='',acom2_model_run='', ncrcat_opts="", freq='1 daily',
                          regrid_method='', regrid_periodicity='', cnk_dict='', generate_weights='', weights_file_name=''):
        '''
        start_date='', regrid_method='', regrid_periodicity='',
                              cnk_dict='', generate_weights='', weights_file_name=''
                              
                if not start_date:         start_date         = self.BRAN_start_date
        if not regrid_method:      regrid_method      = self.BRAN_regrid_method
        if not regrid_periodicity: regrid_periodicity = self.BRAN_periodicity
        if not cnk_dict:           cnk_dict           = self.BRAN_chunking
        if not generate_weights:   generate_weights   = self.BRAN_regen_weights
        '''
        import cosima_cookbook as cc
        if not ncrcat_opts:        ncrcat_opts        = self.ncrcat_opts
        if not acom2_model_run:    acom2_model_run    = self.acom2_model_run
        if not start_date:         start_date         = self.acom2_start_date
        if not regrid_method:      regrid_method      = self.ERA5_regrid_method
        if not regrid_periodicity: regrid_periodicity = self.ERA5_periodicity
        if not cnk_dict:           cnk_dict           = self.acom2_chunking
        if not generate_weights:   generate_weights   = self.ERA5_regen_weights
        if not weights_file_name: 
            F_weights = self.F_ERA5_weights
        else:
            F_weights = weight_file_name
        G_ERA5 = self.era5_grid_prep()
        P_qdp  = os.path.join(self.D01,'qdp','tmp')
        for i in range(self.n_years):
            if i>0:
                dt0 = datetime.strptime(start_date,'%Y-%m-%d %H:%M') + timedelta(days=(i*365))
                dtN = datetime.strptime(start_date,'%Y-%m-%d %H:%M') + timedelta(days=(i*364*2))
            else:
                dt0 = datetime.strptime(start_date,'%Y-%m-%d %H:%M')
                dtN = datetime.strptime(start_date,'%Y-%m-%d %H:%M') + timedelta(days=(364))
            dt0_str = dt0.strftime('%Y-%m-%d %H:%M')
            dtN_str = dtN.strftime('%Y-%m-%d %H:%M')
            P_out  = os.path.join(self.D01,'qdp','ac-om2-2d-qdp-1-daily-mean-ym_{yr:s}.nc'.format(yr=str(dt0.year)))
            if os.path.exists(P_out):
                print('File exists: {:s}\n SKIPPING concatentation')
            else:
                print("extracting ac-om2 data from model run: ",acom2_model_run," from start dt: ",dt0_str," to stop dt: ",dtN_str)
                cnk     = {'time':3,'xt_ocean':100,'yt_ocean':100}
                cc_sess = cc.database.create_session()
                temp    = cc.querying.getvar(expt=acom2_model_run,
                                             variable='temp',
                                             session=cc_sess,
                                             frequency=freq,
                                             start_time=dt0_str,
                                             end_time=dtN_str).chunk(cnk)
                salt    = cc.querying.getvar(expt=acom2_model_run,
                                             variable='salt',
                                             session=cc_sess,
                                             frequency=freq,
                                             start_time=dt0_str,
                                             end_time=dtN_str).chunk(cnk)
                mld     = cc.querying.getvar(expt=acom2_model_run,
                                             variable='mld',
                                             session=cc_sess,
                                             frequency=freq,
                                             start_time=dt0_str,
                                             end_time=dtN_str).chunk(cnk)
                print("computing temperature derivative over time")
                dTdt    = temp.differentiate("time")
                print("generating regridding weights on the fly using xesmf via method: ", regrid_method)
                print("\nsource grid, looks like this: ",G_ERA5)
                lat_mom = temp.yt_ocean
                lon_mom = temp.xt_ocean
                G_tmp   = xr.Dataset({'lon':lon_mom,'lat':lat_mom})
                F_net   = self.compute_era5_net_atmospheric_surface_heat_flux(G_dst              = G_tmp,
                                                                              yr_str             = str(dt0.year),
                                                                              return_daily_mean  = True,
                                                                              regrid_method      = regrid_method,
                                                                              regrid_periodicity = regrid_periodicity,
                                                                              cnk_dict           = self.ERA5_chunking,
                                                                              generate_weights   = generate_weights,
                                                                              weights_file_name  = F_weights)
                F_net = F_net.rename({'nx':'xt_ocean','ny':'yt_ocean'})
                for j in range(len(F_net.time)):
                    F_qdp  = os.path.join(P_qdp,'qdp_yd{:03d}.nc'.format(j))
                    if os.path.exists(F_qdp):
                        print("File exists: {:s}\n skip writing out".format(F_qdp))
                        continue
                    print("\n\ncurrent time: ",datetime.now().time())
                    mld_j  = mld.isel(time=j).astype(np.single)
                    Fnet_j = F_net.isel(time=j).astype(np.single)
                    print("slicing temporal temperature derivative at the mixed layer depth")
                    dTdt_j = dTdt.isel(time=j).sel(st_ocean=mld_j,method='nearest').drop('st_ocean').astype(np.single)
                    print("slicing temporal temperature at the mixed layer depth")
                    temp_j = temp.isel(time=j).sel(st_ocean=mld_j,method='nearest').drop('st_ocean').astype(np.single)
                    print("slicing salinity temperature at the mixed layer depth")
                    salt_j = salt.isel(time=j).sel(st_ocean=mld_j,method='nearest').drop('st_ocean').astype(np.single)
                    print("computing ocean heat capacity at the mixed layer depth")
                    cp_j   = compute_ocn_heat_capacity_at_depth(salt_j, temp_j, mld_j, mld.yt_ocean).astype(np.single)
                    print("computing ocean density at the mixed layer depth")
                    rho_j  = compute_ocn_density_at_depth(salt_j, temp_j, mld_j, mld.yt_ocean).astype(np.single)
                    print("computing ocean heat flux at the mixed layer depth")
                    qdp_j = compute_ocn_heat_flux_at_depth(rho_j, cp_j, mld_j, Fnet_j, dTdt_j)
                    qdp_j  = qdp_j.assign_coords(time=temp.time.isel(time=j).values).expand_dims('time').astype(np.single).to_dataset(name='qdp')
                    print("writing out ocean heat flux netcdf to: ",F_qdp)
                    print("\n this what qdp_j looks like before writing:\n",qdp_j)
                    enc_dict  = {'shuffle':True,'zlib':True,'complevel':5}
                    write_job = qdp_j.to_netcdf(F_qdp, unlimited_dims=['time'], compute=False, encoding={'qdp':enc_dict})
                    with ProgressBar():
                        print(f"Writing to {F_qdp}")
                        write_job.compute()
                exit()
                sys_call = 'ncrcat {:s} {:s}/qdp_yd* {:s}'.format(ncrcat_opts,P_qdp,P_out)
                os.system(sys_call)
                if os.path.exists(P_out):
                    print("successfully created {:s}\n deleting temporary files in {:s}".format(P_out,P_qdp))
                    os.system("rm {:s}/*.nc".format(P_qdp))

    ############################################################################################
    def compute_era5_net_atmospheric_surface_heat_flux(self, G_dst='', yr_str='', return_daily_mean=True,
                                                       regrid_method='', regrid_periodicity='',
                                                       cnk_dict='', generate_weights='', weights_file_name=''):
        '''
        '''
        if not G_dst:              G_dst              = self.cice_grid_prep()
        if not regrid_method:      regrid_method      = self.ERA5_regrid_method
        if not regrid_periodicity: regrid_periodicity = self.ERA5_periodicity
        if not cnk_dict:           cnk_dict           = self.ERA5_chunking
        if not generate_weights:   generate_weights   = self.ERA5_regen_weights
        if not weights_file_name: 
            F_weights = self.F_ERA5_weights
        else:
            F_weights = weights_file_name
        if not yr_str:
            dt_obj = define_datetime_object(start_date=self.ERA5_start_date, year_offset=0)
            yr_str = year_string_from_datetime_object(dt_obj)
        G_ERA5 = self.era5_grid_prep()
        print(F_weights)
        print("generating regridding weights on the fly using xesmf via method: ", regrid_method)
        print("\nsource grid, looks like this: ",G_ERA5)
        print("\ndestination grid, looks like this: ",G_dst)
        lw    = self.era5_load_and_regrid(era5_var           = 'msdwlwrf', 
                                          yr_str             = yr_str,
                                          regrid_method      = regrid_method,
                                          regrid_periodicity = regrid_periodicity,
                                          cnk_dict           = cnk_dict,
                                          generate_weights   = generate_weights,
                                          weights_file_name  = F_weights)
        sw    = self.era5_load_and_regrid(era5_var           = 'msdwswrf', 
                                          yr_str             = yr_str,
                                          regrid_method      = regrid_method,
                                          regrid_periodicity = regrid_periodicity,
                                          cnk_dict           = cnk_dict,
                                          generate_weights   = generate_weights,
                                          weights_file_name  = F_weights)
        sh    = self.era5_load_and_regrid(era5_var           = 'msshf', 
                                          yr_str             = yr_str,
                                          regrid_method      = regrid_method,
                                          regrid_periodicity = regrid_periodicity,
                                          cnk_dict           = cnk_dict,
                                          generate_weights   = generate_weights,
                                          weights_file_name  = F_weights)
        lh    = self.era5_load_and_regrid(era5_var           = 'mslhf', 
                                          yr_str             = yr_str,
                                          regrid_method      = regrid_method,
                                          regrid_periodicity = regrid_periodicity,
                                          cnk_dict           = cnk_dict,
                                          generate_weights   = generate_weights,
                                          weights_file_name  = F_weights)
        F_net = sw.msdwswrf - lw.msdwlwrf - lh.mslhf - sh.msshf
        if return_daily_mean:
            print("computing daily mean of atmospheric heat flux")
            F_net = F_net.resample(time='1D').mean()
        return F_net