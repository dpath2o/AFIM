import numpy      as np
import xarray     as xr
import pandas     as pd
import metpy.calc as mpc
from datetime     import datetime,timedelta
from pyproj       import Geod

def days_since_1601_to_date(days_since_1601):
    """
    Converts a given number of days since January 1, 1601, to a human-readable date format.
    
    Parameters:
        days_since_1601 (int): The number of days since January 1, 1601.
    
    Output:
        Prints the calculated date in the YYYY-MM-DD format.
    """
    start_date      = date(1601, 1, 1)
    offset          = timedelta(days=days_since_1601)
    result_date     = start_date + offset
    print(result_date.strftime("%Y-%m-%d"))

############################################################################
def compute_nsdic_grid_cell_areas(ds, proj_str):
    """
    Computes the area of each grid cell in a given dataset using a specific projection string.
    
    Parameters:
        ds (xarray.Dataset): The dataset containing grid information.
        proj_str (str): The PROJ string describing the coordinate system.
        
    Returns:
        numpy.ndarray: 2D array containing the area of each grid cell.
    """
    a       = 6378273
    b       = 6356889.449
    geod    = Geod(a=a, b=b)
    x_vals  = ds['xgrid'].values
    y_vals  = ds['ygrid'].values
    dx      = ds['xgrid'].diff('x').fillna(0)
    dx_vals = np.append(dx, dx[-1])  # Repeat the last value
    dy      = ds['ygrid'].diff('y').fillna(0)
    dy_vals = np.append(dy, dy[-1])  # Repeat the last value
    proj    = pyproj.Proj(proj_str)
    A       = np.zeros((len(y_vals), len(x_vals)))
    for i in range(len(y_vals)):
        for j in range(len(x_vals)):
            x  = x_vals[j]
            y  = y_vals[i]
            dx = dx_vals[j]
            dy = dy_vals[i]
            # Define the four corners of the cell
            lon1, lat1 = proj(x, y, inverse=True)
            lon2, lat2 = proj(x + dx, y, inverse=True)
            lon3, lat3 = proj(x + dx, y + dy, inverse=True)
            lon4, lat4 = proj(x, y + dy, inverse=True)
            # Calculate the perimeter of the cell
            perimeter = geod.line_length([lon1, lon2, lon3, lon4, lon1], [lat1, lat2, lat3, lat4, lat1])
            # Calculate the area of the cell using Brahmagupta's formula
            s       = perimeter / 2
            a, _, _ = geod.inv(lon1, lat1, lon2, lat2)
            b, _, _ = geod.inv(lon2, lat2, lon3, lat3)
            c, _, _ = geod.inv(lon3, lat3, lon4, lat4)
            d, _, _ = geod.inv(lon4, lat4, lon1, lat1)
            A[i,j]  = np.sqrt((s - a) * (s - b) * (s - c) * (s - d))
    return A

############################################################################
def find_indices_within_region(LAT, LON, regn):
    """
    Find the indices of grid cells within a specified geographical region.
    
    This function takes in 2D arrays or DataArrays of latitude and longitude 
    coordinates, and a list representing the bounding box of a geographical 
    region. It returns the indices of the grid cells that are located within 
    the specified region.
    
    Parameters:
    LAT (np.array or xr.DataArray): 2D array representing the latitudes of the grid cells.
    LON (np.array or xr.DataArray): 2D array representing the longitudes of the grid cells.
    regn (list or tuple): A list or tuple containing the bounding box of the region 
          in the format [min_longitude, max_longitude, min_latitude, max_latitude].
          
    Returns:
    tuple: A tuple containing two arrays. The first array represents the row indices, 
           and the second array represents the column indices of the grid cells 
           within the specified region.
    
    Example:
    --------
    LAT  = np.array([[1, 2], [3, 4]])
    LON  = np.array([[5, 6], [7, 8]])
    regn = [5, 7, 1, 3]
    
    row_indices, col_indices = find_indices_within_region(LAT, LON, regn)
    # row_indices = [0, 1]
    # col_indices = [0, 1]
    """
    ln_mn,ln_mx,lt_mn,lt_mx = regn
    mask                    = (LAT >= lt_mn) & (LAT <= lt_mx) & (LON >= ln_mn) & (LON <= ln_mx)
    return np.where(mask)

############################################################################
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

############################################################################
def compute_ocn_sfc_slope(eta, dxdy_array, direction='x', grid_scale_factor=100):
    """
    Computes sea surface slope in a specified direction.
    
    Parameters:
        eta (float or xarray.DataArray): Sea surface height.
        dxdy_array (float or xarray.DataArray): Grid spacing array. Assumed to be in cm.
        direction (str): The direction for which to compute the sea surface slope ('x' or 'y'). Default is 'x'.
        grid_scale_factor (int): Scale factor to convert dxdy_array units to match eta. Default is 100.
        
    Returns:
        xarray.DataArray: Sea surface slope in the specified direction with appropriate units and attributes.
    """
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
def compute_ocn_heat_flux_at_depth(rho_D, cp_D, D, F_net, dTdt_D, time_unit_to_seconds=3600):
    '''
    Computes the ocean heat flux at a specific depth in W/m^2.
    
    Parameters:
        rho_D (float): Density at depth, in kg/m^3.
        cp_D (float): Heat capacity at depth, in J/(kg*C).
        D (float): Depth in meters.
        F_net (float): Atmospheric heat flux at surface, in W/m^2.
        dTdt_D (float): Time derivative of temperature at depth, in C/hr by default.
        time_unit_to_seconds (int): Conversion factor for time derivative unit to seconds. Default is 3600 seconds.
        
    Returns:
        float: Ocean heat flux at specified depth in W/m^2.
        
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
    Calculates the pressure at a specific depth in the ocean in dbars.
    
    Parameters:
        D (float): Depth in meters.
        latitude (float): Latitude in degrees.
        
    Returns:
        float: Pressure in dbars.
    
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
    Computes the Secant Bulk Modulus (K) of sea water using Equation of State 1980 (UNESCO).
    
    Parameters:
        S (float): Salinity in PSU.
        T (float): Temperature in Celsius.
        P (float): Pressure in dbar.
        
    Returns:
        float: Secant Bulk Modulus in appropriate units.
     
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
    Calculate the density of seawater at a given depth, temperature, salinity, and latitude.
    
    Parameters:
    S (float)        : Salinity (in PSU)
    T (float)        : Temperature (in degrees Celsius)
    D (float)        : Depth (in meters)
    latitude (float) : Latitude (in degrees)
    
    Returns:
    float: Density of seawater (in kg/m^3)
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
    Compute the diffuse electromagnetic radiation at the surface based on total radiation at surface and top of the atmosphere.
    
    Parameters:
    em_sfc (float) : Electromagnetic radiation at surface (W/m^2)
    em_toa (float) : Electromagnetic radiation at top of atmosphere (W/m^2)
    
    Returns:
    float: Diffuse electromagnetic radiation at surface (W/m^2)
    '''
    k_t = em_sfc / em_toa
    return 0.952 - 1.041 * np.exp(-np.exp((2.3 - 4.702*k_t)))

############################################################################
def compute_u2_from_u10(u10):
    '''
    Compute wind speed at 2 meters from reported wind speed at 10 meters.
    
    Parameters:
    u10 (float) : Wind speed at 10 meters (in m/s)
    
    Returns:
    float : Wind speed at 2 meters (in m/s)
    '''
    return (u10 * 4.87) / np.log((67.8 * 10) - 5.42)

############################################################################
def compute_sfc_airrho(t2m, d2m, sp):
    '''
    Compute air density at the surface based on temperature, dewpoint, and surface pressure.
    
    Parameters:
    t2m (float) : Temperature at 2 meters (in Kelvin)
    d2m (float) : Dewpoint at 2 meters (in Kelvin)
    sp  (float) : Surface pressure (in Pascal)
    
    Returns:
    float : Air density at the surface (in kg/m^3)
    '''
    RH  = mpc.relative_humidity_from_dewpoint(t2m,d2m)
    r   = mpc.mixing_ratio_from_relative_humidity(sp, t2m, RH)
    return mpc.density(sp, t2m, r)

############################################################################
def compute_sfc_qsat(d2m, sp):
    '''
    Compute specific humidity at 2 meters based on dewpoint and surface pressure.
    
    Parameters:
    d2m (float) : Dewpoint at 2 meters (in Kelvin)
    sp  (float) : Surface pressure (in Pascal)
    
    Returns:
    float : Specific humidity at 2 meters (dimensionless)
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
    Compute ocean temperature from potential temperature.
    
    Parameters:
    p0 (float) : Reference pressure in decibars (db)
    p  (float) : Pressure at the desired level in decibars (db)
    
    Returns:
    [Currently, this function doesn't return anything; you might want to implement this.]
    
    Note:
    This function currently calculates the pressure difference but doesn't compute the ocean temperature from potential temperature.
    '''
    dp = p0 - p
