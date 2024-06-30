import os, json, re, glob, fnmatch, dask
import xarray            as xr
import numpy             as np
import pandas            as pd
import matplotlib.pyplot as plt
import matplotlib.path   as mpath
import cartopy.crs       as ccrs
import cartopy.feature   as cft
import cmocean           as cm
from datetime            import datetime

# Function to compare a CICE run with NSDIC
def compare_cice_nsdic(P_cice, NSDIC, grid_areas, proj_str, dt0='', dtN='', threshold=0.15, cice_reG_var='aice_m'):
    """
    Compares CICE model output with NSIDC observations.
    
    Parameters:
        P_cice (str): Path to CICE model output files.
        NSDIC (xarray.Dataset): NSIDC observation dataset.
        grid_areas (numpy.ndarray): 2D array containing the area of each grid cell.
        proj_str (str): The PROJ string describing the coordinate system.
        dt0 (str): The start date for slicing time, format YYYY-MM-DD. Default is an empty string.
        dtN (str): The end date for slicing time, format YYYY-MM-DD. Default is an empty string.
        threshold (float): Threshold for masking ice concentration. Default is 0.15.
        cice_reG_var (str): CICE variable name for regridding. Default is 'aice_m'.
        
    Returns:
        xarray.DataArray: Regridded CICE sea ice area extent data.
    """
    CICE          = xr.open_mfdataset(P_cice, decode_coords=False)
    time_index    = pd.DatetimeIndex(CICE['time'].values)
    adjusted_time = time_index - pd.DateOffset(months=1)
    CICE['time']  = adjusted_time.values
    mask          = (CICE['TLAT'] < 0).compute()
    CICE          = CICE.sel(time=slice(dt0,dtN))
    CICE_SH       = CICE.where(mask, drop=True)
    #Convert xgrid and ygrid to latitude and longitude, and add the new coordinates to NSDIC Dataset
    in_proj            = pyproj.Proj(proj_str)
    out_proj           = pyproj.Proj(proj="latlong", datum="WGS84")
    transformer        = pyproj.Transformer.from_proj(in_proj, out_proj, always_xy=True)
    xx, yy             = np.meshgrid(NSDIC['xgrid'].values, NSDIC['ygrid'].values)  # Make sure x_vals and y_vals are defined
    lon_vals, lat_vals = transformer.transform(xx, yy)
    NSDIC['latitude']  = (('y', 'x'), lat_vals)
    NSDIC['longitude'] = (('y', 'x'), lon_vals)
    # Define source grid from the subsetted CICE6 data
    src_grid = xr.Dataset({'lat': (['nj', 'ni'], CICE_SH['TLAT'][0].values), 'lon': (['nj', 'ni'], CICE_SH['TLON'][0].values)})
    # Define destination grid from NSIDC
    dst_grid = xr.Dataset({'lat': (['y', 'x'], lat_vals), 'lon': (['y', 'x'], lon_vals)})
    # Create regridder object
    regridder = xe.Regridder(src_grid, dst_grid, method='bilinear', periodic=False)
    # reG
    reG_aice = regridder(CICE_SH[cice_reG_var])
    # Maskthe concentration data based on threshold for CICE6. Then computing sea ice area extent of CICE6, which requires the above cell to be run first
    mask_CICE                       = reG_aice > threshold
    reG_aice['sea_ice_area_extent'] = (mask_CICE * reG_aice * grid_areas).sum(dim=['y', 'x'])/1e7
    return reG_aice

############################################################################
def compute_cice_area(CICE='',
                      threshold=0.15,
                      hemisphere='sp',
                      monthly=True):
    """
    Computes the total ice area based on a given threshold concentration and specified hemisphere.
    
    Parameters:
        CICE (xarray.Dataset): The dataset containing CICE model variables.
        threshold (float): The threshold for ice concentration. Default is 0.15.
        hemisphere (str): Specifies the hemisphere ('sp' for South Pole, 'np' for North Pole). Default is 'sp'.
        monthly (bool): Indicates whether to use monthly ('aice_m') or general ('aice') variable. Default is True.
        
    Returns:
        float: The total area of ice in the specified hemisphere that exceeds the threshold concentration.
    """
    if hemisphere=='sp':
        lat_slice = (-90,-45)
    elif hemisphere=='np':
        nj_slice = (45,90)
    if monthly:
        var_name = 'aice_m'
    else:
        var_name = 'aice'
    lons      = CICE.TLON.isel(TLAT=slice(lat_slice))
    lats      = CICE.TLAT.isel(TLAT=slice(lat_slice))
    aice_bool = CICE[var_name].isel(time=0,TLAT=slice(lat_slice)).where(ICE1.aice_m.isel(time=0,TLAT=slice(lat_slice)) >= threshold).values.flatten()
    coords    = list(zip(lons.values.flatten(), lats.values.flatten()))
    n         = len(coords)
    coords    = [(radians(lon), radians(lat)) for lon, lat in coords]
    area      = 0
    for i in range(n):
        if aice_bool[i]:
            j          = (i + 1) % n
            x1, y1, z1 = np.cos(coords[i][1]) * np.cos(coords[i][0]), np.cos(coords[i][1]) * np.sin(coords[i][0]), np.sin(coords[i][1])
            x2, y2, z2 = np.cos(coords[j][1]) * np.cos(coords[j][0]), np.cos(coords[j][1]) * np.sin(coords[j][0]), np.sin(coords[j][1])
            area      += (x1 * y2 - x2 * y1) * (y1 + y2 + z1 + z2)
    return abs(area) / 2

############################################################################
def read_json(F_json):
    '''Read a JSON file.

    Parameters:
    -----------
    filename : str
        Path to the JSON file.

    Returns:
    --------
    dict
        Dictionary containing the parsed JSON data.
    '''
    with open(F_json, 'r') as file:
        return json.load(file)

######################################################################################################################
#########################################     ANALYSIS CLASS     #####################################################
######################################################################################################################

class analysis:
    '''
    CICE (Community Ice CodE) Analysis Class.

    Description:
    ------------
    This class is designed to load and handle parameters and configurations for analyzing CICE datasets.
    The class reads a JSON file containing necessary attributes and setups for the analysis, such as time frames,
    thresholds, directories, titles, locations, and other specific parameters.

    Attributes:
    -----------
    dt0 : str
        Start date for the analysis period.
    dtN : str
        End date for the analysis period.
    aice_thresh : float
        Threshold value for the sea ice concentration.
    FI_thresh : float
        Threshold value for some feature of interest (e.g., flux or intensity).
    aice_name : str
        Variable name for sea ice concentration within the dataset.
    proc_freq : str
        Processing frequency for the data (e.g., 'monthly').
    P_nsdic : str
        Path to some data directory or file, e.g., National Snow and Ice Data Center.
    model_run_date : str
        Date the model was run.
    P_cices : str
        Path to CICE data or files.
    titles : list of str
        List of titles for plots or analysis tasks.
    olav_locs : list of lists
        List of [longitude, latitude] locations related to the Olav region.
    maws_locs : list of lists
        List of [longitude, latitude] locations related to the Maws region.
    spacing : str
        Grid spacing for the GMT regular grid.
    search_radius : str
        Search radius for GMT's nearneighbor algorithm.
    cmap_plot_cice : str
        Color map for plotting CICE data.
    cice_labels : list of str
        Labels for CICE data plots.
    region_names : list of str
        Names of regions of interest.
    regions_info : dict
        Dictionary containing additional info about the regions of interest.
    GI_locations : list of lists
        List of geographical interest locations for analysis.

    Methods:
    --------
    Currently, this class doesn't define extra methods, but future implementations might include 
    methods for data processing, plotting, and other analysis tasks.

    Example:
    --------
    analysis = Analysis(F_json="config.json", dt0="1990-01-01", dtN="1991-01-01", aice_thresh=0.15)
    print(analysis.dt0)   # Prints the start date from the JSON configuration or the user input.

    '''
    #####################################################################################
    def __init__(self, F_json="/home/581/da1339/AFIM/src/python/AFIM/JSONs/afim_cice_analysis.json",
                 dt0=None, dtN=None, aice_thresh=None, stic_thresh=None, proc_freq=None, GI_dict=None, hemisphere=None,
                 nco_cat=None):
        '''Initialize the analysis object by reading parameters from a JSON file.'''
        data = read_json(F_json)
        # Load parameters from JSON or use user-provided values
        self.nco_cat      = nco_cat     if nco_cat     is not None else data['nco_cat']
        self.hemisphere   = hemisphere  if hemisphere  is not None else data['hemisphere']
        self.dt0          = dt0         if dt0         is not None else data['dt0']
        self.dtN          = dtN         if dtN         is not None else data['dtN']
        self.aice_thresh  = aice_thresh if aice_thresh is not None else data['aice_thresh']
        self.stic_thresh  = stic_thresh if stic_thresh is not None else data['stic_thresh']
        self.proc_freq    = proc_freq   if proc_freq   is not None else data['proc_freq']
        self.GI_dict      = GI_dict     if GI_dict     is not None else data['GI_dict']
        # Read other static attributes from JSON
        self.G_res              = data["G_res"]
        self.F_G_AFIM           = data['F_G_AFIM']
        self.G_AFIM             = xr.open_dataset(self.F_G_AFIM)
        self.D_AOM2             = data['D_AOM2']
        self.D_AOM2_local       = data['D_AOM2_local']
        self.D_AFIM_jk72        = data['D_AFIM_jk72']
        self.D_AFIM_ol01        = data['D_AFIM_ol01']
        self.D_tmp              = data['D_tmp']
        self.D_graphical        = data["D_graphical"]
        self.D_NSIDC            = data['D_NSIDC']
        self.P_NSIDC_cell_area  = data['P_NSIDC_cell_area']
        self.NSIDC_flags        = data['NSIDC_flags']
        self.F_cice_hist_dy_fmt = data['F_cice_hist_dy_fmt']
        self.F_cice_hist_mo_fmt = data['F_cice_hist_mo_fmt']
        self.cice_stic_vars     = data['cice_stic_vars']
        self.AFIM_runs          = self.construct_afim_full_path(data['AFIM_runs'])
        self.plot_var_dict      = data['plot_var_dict']
        self.grid_shorthand     = {'t': 'TLAT', 'u': 'ULAT', 'ot': 'yt_ocean', 'ou': 'yu_ocean', 'y': 'ygrid'}

    #####################################################################################
    # Function to construct the full path
    def construct_afim_full_path(self,afim_run_dict):
        for afim_run_name, value in afim_run_dict.items():
            D_base      = eval(f"self.D_AFIM_{value['project']}")
            #D_base      = getattr(self, project_key)
            value['D_'] = f"{D_base}/{value['D_']}"
        return afim_run_dict

    #####################################################################################
    def find_files(self, D_base, dt0=None, dtN=None, proc_freq=None, ds_type='afim', frcg='era5'):
        dt0       = dt0       if dt0       is not None else self.dt0
        dtN       = dtN       if dtN       is not None else self.dtN
        proc_freq = proc_freq if proc_freq is not None else self.proc_freq
        dt0_str = datetime.strptime(dt0, '%Y-%m-%d %H:%M')
        dtN_str = datetime.strptime(dtN, '%Y-%m-%d %H:%M')
        if ds_type=='afim':
            if proc_freq == 'monthly':
                search_pattern = f"{D_base}/history/monthly/iceh.*.nc"
            elif proc_freq == 'daily':
                search_pattern = f"{D_base}/history/daily/iceh.*.nc"
            else:
                raise ValueError(f"Invalid proc_freq: {proc_freq}")
            date_format = "%Y-%m" if proc_freq == 'monthly' else "%Y-%m-%d"
        elif ds_type=='aom2':
            if proc_freq == 'monthly':
                search_pattern = f"{D_base}/{frcg.upper()}/ice/{self.G_res}/monthly/iceh.*.nc"
            elif proc_freq == 'daily':
                search_pattern = f"{D_base}/{frcg.upper()}/ice/{self.G_res}/daily/iceh.*-daily.nc"
            else:
                raise ValueError(f"Invalid proc_freq: {proc_freq}")
            date_format = "%Y-%m"
        elif ds_type=='nsidc':
            if proc_freq == 'monthly':
                search_pattern = f"{D_base}/{self.hemisphere}/monthly/seaice_conc_monthly_*.nc"
            elif proc_freq == 'daily':
                search_pattern = f"{D_base}/{self.hemisphere}/daily/seaice_conc_daily_*.nc"
            else:
                raise ValueError(f"Invalid proc_freq: {proc_freq}")
            date_format = "%Y%m" if proc_freq == 'monthly' else "%Y%m%d"
        else:
            raise ValueError(f"Invalid ds_type: {ds_type}")
        F_found = sorted(glob.glob(search_pattern, recursive=True))
        print(f"Search pattern: {search_pattern}")
        print(f"Files found: {len(F_found)}")
        F_filtered = []
        for f in F_found:
            try:
                if ds_type=='nsidc':
                    file_date_str = os.path.basename(f).split('_')[4]
                elif ds_type=='aom2':
                    if proc_freq=='daily':
                        tmp_str       = os.path.basename(f).split('.')[1]
                        file_date_str = tmp_str.split('-daily')[0]
                    else:
                        file_date_str = os.path.basename(f).split('.')[1]
                else:
                    file_date_str = os.path.basename(f).split('.')[1]  # Adjust the split index as per filename format
                file_date = datetime.strptime(file_date_str, date_format)
                if dt0_str <= file_date <= dtN_str:
                    F_filtered.append(f)
            except Exception as e:
                print(f"Error processing file {f}: {e}")
        print(f'Found {len(F_filtered)} files within the date range.')
        return F_filtered
        
    #####################################################################################
    def define_hemisphere(self,hemisphere):
        if hemisphere in ['north', 'northern', 'NH', 'nh']:
            return 'nh'
        elif hemisphere in ['south', 'southern', 'SH', 'sh']:
            return 'sh'
        else:
            raise ValueError(f"Invalid hemisphere '{self.hemisphere}'. Valid options are: ['north', 'south','northern','southern','sh','nh','SH','NH']")

    #####################################################################################
    def define_date_range(self, start_date=None, stop_date=None, periods=None, frequency=None):
        if start_date is None: start_date = self.dt0
        if frequency  is None: frequency  = self.proc_freq
        if frequency=='monthly': frequency='M'
        if frequency=='daily'  : frequency='D'
        if stop_date is None:
            stop_date = self.dtN
            return pd.date_range(start=start_date, end=stop_date, freq=frequency)
        else:
            if periods is None: periods = 365
            return pd.date_range(start=start_date, end=stop_date, freq=frequency)

    #####################################################################################
    def load_dataset(self, DS_name=None, run_name=None, D_base=None,
                     dt0=None, dtN=None, proc_freq='monthly', hemisphere='south', data_type='ocean', frcg_type='era5',
                     F_search_pattern=None, variable_names=None, rtn_ds=True):
        """
        General method to load various datasets (AFIM, AOM2, NSIDC).

        Parameters:
        -----------
        DS_name : str
            Name of the dataset to load; there are three options: 'jk72', 'ol01', 'AOM2', 'NSIDC'
        run_name : str
            The name of the AFIM run
        D_base : str
            Base directory where the dataset files are located.
        dt0 : str, optional
            Start date in the format 'YYYY-MM-DD'. Default is class attribute self.dt0.
        dtN : str, optional
            End date in the format 'YYYY-MM-DD'. Default is class attribute self.dtN.
        proc_freq : str, optional
            Processing frequency ('daily', 'monthly'). Default is class attribute self.proc_freq.
        hemisphere : str, optional
            Hemisphere for NSIDC data ('north', 'south'). Default is class attribute self.hemisphere.
        data_type : str, optional, default='ocean'
            Either 'ocean' or 'ice' and this option is specific to AOM2 datasets 
        F_search_pattern : str, optional
            File pattern to match dataset files. Default patterns will be used if not provided.
        variable_names : dict, optional
            Dictionary of variable names specific to the dataset type.
        rtn_ds : bool, optional
            If True, returns the loaded dataset. Default is False.

        Returns:
        --------
        xarray.Dataset or None
            Loaded dataset if rtn_ds is True, otherwise sets the dataset as an attribute.
        """
        if dt0 is None:
            dt0 = self.dt0
        if dtN is None:
            dtN = self.dtN
        if proc_freq is None:
            proc_freq = self.proc_freq
        if hemisphere is None:
            hemisphere = self.hemisphere
        if hemisphere in ['southern', 'sh', 'south', 's', 'sthn', 'sth', 'SH']:
            hemisphere = 'south'
            hem_short_name = 'sh'
        elif hemisphere in ['northern', 'nh', 'north', 'n', 'nthn', 'nth', 'NH']:
            hemisphere = 'north'
            hem_short_name = 'nh'
        else:
            raise ValueError("Hemisphere must be 'north' or 'south'.")
        dt0 = pd.to_datetime(dt0)
        dtN = pd.to_datetime(dtN)
        if DS_name in ['jk72', 'ol01']:
            if D_base is None:
                D_base = self.AFIM_runs[run_name]['D_']
            F_search_pattern = F_search_pattern or (r"iceh\.\d{4}-\d{2}-\d{2}\.nc" if proc_freq=='daily' else r"iceh\.\d{4}-\d{2}\.nc")
            D_base = os.path.join(D_base, "history", proc_freq)
        elif DS_name == 'AOM2':
            if D_base is None:
                D_base = self.D_AOM2_local
            if frcg_type == 'era5':
                D_base = os.path.join(D_base, "ERA5")
            else:
                D_base = os.path.join(D_base, "JRA55")
            if data_type == 'ocean':
                F_search_pattern = F_search_pattern or (r"\d{4}_ocean_daily\.nc" if proc_freq=="daily" else r"\d{4}_ocean_month\.nc")
                D_base = os.path.join(D_base, "ocean", self.G_res, proc_freq)
            else:
                F_search_pattern = F_search_pattern or (r"iceh\.\d{4}-\d{2}-daily\.nc" if proc_freq=='daily' else r"iceh\.\d{4}-\d{2}\.nc")
                D_base = os.path.join(D_base, "ice", self.G_res, proc_freq)
        elif DS_name == 'NSIDC':
            if D_base is None:
                D_base = self.D_NSIDC
            F_search_pattern = F_search_pattern or (r".*\d{8}.*\.nc" if proc_freq=='daily' else r".*\d{6}.*\.nc")
            D_base = os.path.join(D_base, hemisphere, proc_freq)
        print(f"search directory {D_base} with search pattern {F_search_pattern}")
        F_full_list = [f for f in glob.glob(os.path.join(D_base, "*")) if re.match(F_search_pattern, os.path.basename(f))]
        F_filt_list = []
        for F_ in F_full_list:
            if DS_name in ['jk72','ol01']:
                if proc_freq == 'monthly':
                    F_dt_str = re.search(r"\d{4}-\d{2}", os.path.basename(F_)).group(0)
                    F_dt = pd.to_datetime(F_dt_str, format='%Y-%m')
                else:
                    F_dt_str = re.search(r"\d{4}-\d{2}-\d{2}", os.path.basename(F_)).group(0)
                    F_dt = pd.to_datetime(F_dt_str, format='%Y-%m-%d')
            elif DS_name == 'AOM2' and data_type == 'ice':
                F_dt_str = re.search(r"\d{4}-\d{2}", os.path.basename(F_)).group(0)
                F_dt = pd.to_datetime(F_dt_str, format='%Y-%m')
            elif DS_name == 'AOM2' and data_type == 'ocean':
                F_dt_str = re.search(r"\d{4}", os.path.basename(F_)).group(0)
                F_dt = pd.to_datetime(F_dt_str, format='%Y')
            elif DS_name == 'NSIDC':
                if proc_freq == 'monthly':
                    F_dt_str = re.search(r"\d{6}", os.path.basename(F_)).group(0)
                    F_dt = pd.to_datetime(F_dt_str, format='%Y%m')
                else:
                    F_dt_str = re.search(r"\d{8}", os.path.basename(F_)).group(0)
                    F_dt = pd.to_datetime(F_dt_str, format='%Y%m%d')
            if (dt0 and F_dt < dt0) or (dtN and F_dt > dtN):
                continue
            F_filt_list.append(F_)
        if not F_filt_list:
            print("No files matching the pattern and date range found.")
            return None
        print("Loading found files")
        if DS_name in ['jk72','ol01','AOM2']:
            concat_dim = 'time'
        else:
            concat_dim = 'tdim'
        with dask.config.set(scheduler='threads'):
            DS = xr.open_mfdataset(F_filt_list, concat_dim=concat_dim, combine='nested', parallel=True)
        if DS_name in ['jk72', 'ol01', 'AOM2'] and proc_freq == "monthly":
            DS['time'] = pd.to_datetime(DS['time'].values) - pd.DateOffset(months=1)
        elif DS_name in ['jk72', 'ol01', 'AOM2'] and proc_freq == 'daily':
            DS['time'] = pd.to_datetime(DS['time'].values) - pd.DateOffset(days=1)
        if DS_name in ['jk72', 'ol01']:
            DS = DS.assign_coords(tlon_i=('ni', DS['TLON'][0, :].values))
            DS = DS.assign_coords(tlat_i=('nj', DS['TLAT'][:, 0].values))
            DS = DS.assign_coords(ulon_i=('ni', DS['ULON'][0, :].values))
            DS = DS.assign_coords(ulat_i=('nj', DS['ULAT'][:, 0].values))
        if rtn_ds:
            return DS
        else:
            setattr(self, DS_name, DS.sel(time=slice(dt0, dtN)))

    #####################################################################################
    def compute_sia_or_sie(self, DA, DS_name='AFIM', sia_or_sie='sia',
                           var_name='cdr_seaice_conc_monthly',
                           area=None, aice_thresh=None, flags=None):
        """
        Compute Sea Ice Area (SIA) for the given dataset -- be sure to mask the hemisphere first!

        Parameters:
        -----------
        dataset : xarray.Dataset
        The input dataset containing sea ice concentration and cell area.
        DS_name : str, optional, DEFAULT='AFIM'
        Name of the dataset
        var_name : str, optional
        Variable name for sea ice concentration.
        area_name : str, optional
        Variable name for cell area.
        aice_thresh : float, optional
        Threshold value for the sea ice concentration.
        flags : list of float, optional
        List of flag values to be treated as NaN. Default is None.

        Notes:
        -------
            - NSIDC sea ice concentration variable name : 'cdr_seaice_conc_monthly' and 'cdr_seaice_conc'

        Returns:
        --------
        xarray.DataArray
        Computed Sea Ice Area (SIA) or Sea Ice Extent (SIE)
        """
        if aice_thresh is None:
            aice_thresh = self.aice_thresh
        if DS_name=='NSIDC':
            aice = DA[var_name]
            for flag in self.NSIDC_flags:
                aice = xr.where(aice == flag / 100, np.nan, aice)
            cell_area = xr.open_dataset(self.P_NSIDC_cell_area).cell_area
            mask = aice > aice_thresh
            if sia_or_sie=='sia':
                sia = (aice * cell_area).where(mask).sum(dim=['y', 'x'], skipna=True) / 1e12
                return sia
            else:
                sie = (mask * cell_area).sum(dim=['y', 'x'], skipna=True) / 1e12
                return sie
        else:
            aice = DA
            mask = aice > aice_thresh
            cell_area = area / 1e12  # Convert cell_area from m^2 to million km^2
            if sia_or_sie=='sia':
                sia = (aice * cell_area).where(aice > aice_thresh).sum(dim=['nj', 'ni'], skipna=True)
                return sia
            else:
                sie = (mask * cell_area).sum(dim=['nj', 'ni'], skipna=True)
                return sie

    #####################################################################################
    def mask_hemisphere(self, DS, var_name=None, grid_type='t', hemisphere='south'):
        """
        Apply a hemisphere mask to a dataset variable; *not* required for NSIDC data

        Parameters:
        -----------
        DS : xarray.Dataset
            The dataset containing the variable to be masked.
        var_name : str
            The name of the variable within the dataset to be masked.
        grid_type : str, optional
            The type of grid to use ('t', 'u', 'ot', 'ou'). Default is 't'.
        hemisphere : str, optional
            The hemisphere to mask ('north' or 'south'). Default is 'south'.

        Returns:
        --------
        xarray.DataArray
            The masked DataArray for the specified variable.
        
        Raises:
        -------
        ValueError
            If the hemisphere is not 'north' or 'south'.
        
        Example:
        --------
        masked_da = analysis.mask_hemisphere(DS, var_name='aice', grid_type='t', hemisphere='north')
        """
        DA = DS[var_name]
        if hemisphere in ['southern', 'sh', 'south', 's', 'sthn', 'sth', 'SH']:
            hemisphere = 'south'
        elif hemisphere in ['northern', 'nh', 'north', 'n', 'nthn', 'nth', 'NH']:
            hemisphere = 'north'
        else:
            raise ValueError("Hemisphere must be 'north' or 'south'.")
        op = (lambda x: x < 0) if hemisphere == 'south' else (lambda x: x > 0)
        mask = op(DS[self.grid_shorthand[grid_type]])
        return DA.where(mask)

    #####################################################################################
    def find_nearest_coordinates(self, GLAT, GLON, geo_locs):
        '''
        Identify the nearest geographic coordinates in a matrix to a set of reference coordinates.

        Parameters:
        -----------
        GLAT : numpy.ndarray
            Two-dimensional array (matrix) representing latitudes.

        GLON : numpy.ndarray
            Two-dimensional array (matrix) representing longitudes.

        geo_locs : numpy.ndarray

        Returns:
        --------
        (numpy.ndarray, numpy.ndarray)
            Two arrays representing the longitudes and latitudes in GLON and GLAT, respectively, 
            that are closest to the reference locations provided in `geo_locs`.

        Description:
        ------------
        For each reference location in `geo_locs`, the function calculates the squared 
        differences between the reference location and each coordinate pair in GLAT and GLON.
        It then identifies the matrix indices where the total squared difference (latitude and 
        longitude combined) is minimized. The function returns the GLON and GLAT values at these 
        indices as the nearest coordinates to the reference location.
        '''
        result = []
        for loc in geo_locs:
            lat_diff = (GLAT - loc[1]) ** 2
            lon_diff = (GLON - loc[0]) ** 2
            total_diff = lat_diff + lon_diff
            yi, xi = np.unravel_index(total_diff.argmin(), total_diff.shape)
            result.append([GLON[yi, xi], GLAT[yi, xi]])
        self.GI_dict = np.array(result)

    #####################################################################################
    @staticmethod
    def get_iceberg_grid_indices(iceberg_lons, iceberg_lats, grid_lons, grid_lats):
        """
        Calculate the grid indices corresponding to specified iceberg positions.

        This method determines the indices of grid cells that are closest to
        given iceberg positions based on the provided grid longitude and latitude arrays.

        Parameters:
        - iceberg_lons (array-like): Longitudes of the icebergs.
        - iceberg_lats (array-like): Latitudes of the icebergs.
        - grid_lons (2D array-like): 2D array representing the grid cell center longitudes.
        - grid_lats (2D array-like): 2D array representing the grid cell center latitudes.

        Returns:
        - list of tuples: Each tuple contains the (j, i) index corresponding to an iceberg position.
        """

        indices = []
        for lon, lat in zip(iceberg_lons, iceberg_lats):
            lon_diff = np.abs(grid_lons - lon)
            lat_diff = np.abs(grid_lats - lat)
            j, i = np.unravel_index(np.argmin(lon_diff + lat_diff), lon_diff.shape)
            print(f"Iceberg at ({lon}, {lat}) is mapped to grid cell ({j}, {i})")
            indices.append((j, i))
        return indices

    #####################################################################################
    @staticmethod
    def get_cell_boundaries(lon, lat, GLON, GLAT):
        """
        Get the boundaries of a grid cell for a given geographic point.

        This method returns the bounding longitudes and latitudes of a grid cell that 
        is closest to the provided longitude and latitude point.

        Parameters:
        - lon (float): Longitude of the given point.
        - lat (float): Latitude of the given point.
        - GLON (2D array-like): 2D array representing the grid cell center longitudes.
        - GLAT (2D array-like): 2D array representing the grid cell center latitudes.

        Returns:
        - list of floats: The boundaries of the grid cell in the format [min_lon, max_lon, min_lat, max_lat].
        """

        lon_diff         = np.abs(GLON - lon)
        lat_diff         = np.abs(GLAT - lat)
        j, i             = np.unravel_index(np.argmin(lon_diff + lat_diff), lon_diff.shape)
        min_lon, max_lon = GLON[j, i], GLON[j, i+1]
        min_lat, max_lat = GLAT[j, i], GLAT[j+1, i]
        return [min_lon, max_lon, min_lat, max_lat]

    #####################################################################################
    def construct_graphical_output_directory(self,
                                            D_base           = '',
                                            model_name       = '',
                                            coast_name       = '',
                                            mean_length_str  = '',
                                            var_name         = ''):
        '''
        Construct a directory path for graphical outputs based on model parameters and characteristics.

        Parameters:
        - D_base (str): The base directory path. Defaults to self.D_graphical if not provided.
        - model_name (str): Name of the model.
        - var_name (str): Variable name.
        - mean_length_str (str, optional): String representing the length of averaging or time-mean (e.g., "5day", "monthly").

        Sets:
        - self.D_figure (str): Constructed directory path based on the provided parameters.

        Note:
        - The function checks if the directory exists and if not, it creates one.
        '''
        if not D_base    : D_base     = self.D_graphical
        if not coast_name: coast_name = 'circumpolar'
        if mean_length_str:
            self.D_fig = os.path.join(D_base, model_name, coast_name, var_name, mean_length_str)
        else:
            self.D_fig = os.path.join(D_base, model_name, coast_name, var_name)
        if not os.path.exists(self.D_fig): os.makedirs(self.D_fig)

    #####################################################################################
    def construct_date_string_from_model(self, t_inc=None, t_search=None, var_name='', user_data=None):
        """
        Construct a date string based on the provided time index or search time and variable name.
        The date string will be stored in the instance variable `self.fig_dt_str`.

        Args:
            t_inc (int or None): Time index for which to construct the date string.
            t_search (str or None): A string representing the time to search for.
            var_name (str): Name of the variable.
            user_data (xarray.Dataset or None): Optional user-provided data. If provided,
                the time will be extracted from this data instead of the model data.

        Raises:
            ValueError: If both `t_inc` and `t_search` are None or if an error occurs while
                constructing the date string.
        """
        if t_inc is None and t_search is None:
            raise ValueError("Either t_inc or t_search must be provided.")

        try:
            if t_search:
                time_val = pd.Timestamp(t_search)
                data = user_data if user_data is not None else self.AFIM
                if 'time' in data:
                    t_index = abs(data['time'] - np.datetime64(time_val)).argmin().item()
                    time_val = data.isel(time=t_index).time.values
                else:
                    raise ValueError(f"No 'time' dimension found in the data for var_name={var_name}.")
            else:
                data = user_data if user_data is not None else self.AFIM
                time_val = data.isel(time=t_inc).time.values

            self.fig_dt_str = pd.Timestamp(time_val).to_pydatetime().strftime('%Y_%m_%d')
        except Exception as e:
            raise ValueError(f"Unable to construct date string for t_inc={t_inc}, t_search={t_search}, var_name={var_name}. Error: {e}")

    #####################################################################################
    def construct_graphical_output_filename(self,
                                            dt_str        = '',
                                            var_name      = '',
                                            model_name    = '',
                                            img_type      = '',
                                            hemisphere    = 'sh',
                                            atm_frcg_name = '',
                                            ocn_frcg_name = '',
                                            G_res         = ''):
        '''
        '''
        if not img_type      : img_type = self.img_type
        if not G_res         : G_res    = self.G_res
        if not atm_frcg_name : atm_frcg_name = self.atm_frcg_name
        if not ocn_frcg_name : ocn_frcg_name = self.ocn_frcg_name

        if hemisphere=='both':
            hemi_str = 'sh_nh'.upper()
        else:
            hemi_str = hemisphere.upper()
        self.F_fig = f'{dt_str}_{model_name}_{var_name}_{hemi_str}_{atm_frcg_name}_{ocn_frcg_name}_{G_res}.{img_type}'

    #####################################################################################
    def update_attrs(self, source_dict, var_inc, model_run_name, dt_str, 
                    tit_str=None, cbar_lab=None, cbar_units=None, cmap=None, grid_type=None):
        """
        Update attributes for a given variable.

        Args:
            source_dict (dict): Dictionary containing metadata information for various variables.
            var_inc (str): Name of the variable of interest.
            tit_str (str, optional): User-defined title string. Defaults to None.
            cbar_lab (str, optional): User-defined colorbar label. Defaults to None.
            cbar_units (str, optional): User-defined colorbar units. Defaults to None.
            cmap (str, optional): User-defined colormap. Defaults to None.
            grid_type (str, optional): User-defined grid type. Defaults to None.

        Returns:
            tuple: Descriptions, labels, units, colormaps, and grid types associated with the variable.
        """
        tit_str    = tit_str    or f"{source_dict['var_descs'][var_inc]} {model_run_name} {dt_str}"
        cbar_lab   = cbar_lab   or source_dict['var_labs'][var_inc]
        cbar_units = cbar_units or source_dict['var_units'][var_inc]
        cmap       = cmap       or source_dict['var_cmaps'][var_inc]
        grid_type  = grid_type  or source_dict['var_grids'][var_inc].upper()
        return tit_str, cbar_lab, cbar_units, cmap, grid_type

    #####################################################################################
    def spatially_trim_data_to_region(self, dat, glon_str, glat_str, region=None):
        """
        Trim spatial data to a specified region.

        Args:
            dat (xarray.DataArray or similar): The data array containing spatial information.
            glon_str (str): The longitude variable key in 'dat'.
            glat_str (str): The latitude variable key in 'dat'.
            region (tuple, optional): A tuple (lon_min, lon_max, lat_min, lat_max) defining the region to trim to.

        Returns:
            tuple: Trimmed data, latitudes, and longitudes.
        """
        if region:
            lon_key, lat_key = f"{glon_str}", f"{glat_str}"
            dat = dat.where((dat[lon_key] >= region[0]) & (dat[lon_key] <= region[1]) & (dat[lat_key] >= region[2]) & (dat[lat_key] <= region[3]))
        lon = np.ravel(dat[glon_str].values)
        lat = np.ravel(dat[glat_str].values)
        dat = np.ravel(dat.values)
        mask = ~np.isnan(dat)
        return dat[mask], lat[mask], lon[mask]

    #####################################################################################
    def set_vmin_vmax(self, dat, source_dict, var_inc):
        """
        Determine minimum and maximum values for a data array based on a source dictionary.

        Args:
            dat (array): Data array to find vmin and vmax for.
            source_dict (dict): Dictionary containing vmin and vmax information for various variables.
            var_inc (str): Name of the variable of interest.

        Returns:
            tuple: Computed or specified vmin and vmax values.
        """
        if dat.size == 0:
            print("Warning: Skipping this t_inc, as the data array is empty.")
            return None, None  # This will indicate to the caller that it should skip this iteration
        if not source_dict['var_vmins'][var_inc]:
            vmin = np.nanmin(dat)
        else:
            vmin = source_dict['var_vmins'][var_inc]
        if not source_dict['var_vmaxs'][var_inc]:
            vmax = np.nanmax(dat)
        else:
            vmax = source_dict['var_vmaxs'][var_inc]
        return vmin, vmax

    #####################################################################################
    def plot_basic_map(self, DS=None, var_name=None, model_name=None, model_forcing=None, ice_or_ocean='ice',
                       dt0=None, dtN=None, G_=None, grid_type='t', region=None, overwrite=False,
                       lon_name=None, lat_name=None, figsize=(15, 9), hemisphere=None,
                       projection=ccrs.SouthPolarStereo(), cmap=None, vmin=None, vmax=None,
                       cbar_label=None, tit_str=None, D_save=None, F_save=None, save=True, show=True):
        """
        Plot a basic map with a specified projection, color map, and region.

        Parameters:
        ----------
        DS : xarray.Dataset
            Data to plot.
        var_name : str
            Name of the variable.
        model_name : str
            Name of the model.
        model_forcing : str
            Forcing of the model.
        ice_or_ocean : str, default 'ice'
            Type of data ('ice' or 'ocean').
        dt0 : str, optional
            Start date/time string in the format YYYY-MM-DD HH:MM for the first time stamp to plot.
        dtN : str, optional
            End date/time string in the format YYYY-MM-DD HH:MM for the last time step to plot.
        G_ : xarray.Dataset, optional
            Grid dataset containing latitude and longitude coordinates.
        grid_type : str, default 't'
            Type of grid ('t' or 'u').
        region : list, optional
            Region to plot [lon_min, lon_max, lat_min, lat_max]; default is None.
        overwrite : bool, default False
            Whether to overwrite existing files.
        lon_name : str, optional
            Name of the longitude variable in the grid file.
        lat_name : str, optional
            Name of the latitude variable in the grid file.
        figsize : tuple, default (15, 9)
            Size of the figure.
        hemisphere : str, optional
            Hemisphere ('nh' or 'sh').
        projection : cartopy.crs.Projection, default ccrs.SouthPolarStereo()
            Cartopy projection to use.
        cmap : matplotlib.colors.Colormap, optional
            Colormap to use.
        vmin : float, optional
            Minimum value for colormap.
        vmax : float, optional
            Maximum value for colormap.
        cbar_label : str, optional
            Label for the colorbar.
        tit_str : str, optional
            Title of the plot.
        D_save : str, optional
            Directory to save the plot.
        F_save : str, optional
            Filename to save the plot.
        save : bool, default True
            Whether to save the plot.
        show : bool, default True
            Whether to display the plot.
        """
        if DS is None:
            raise ValueError("Dataset 'DS' must be provided")
        if var_name is None:
            raise ValueError("Variable name 'var_name' must be provided")
        if model_name is None:
            raise ValueError("Model name 'model_name' must be provided")
        if model_forcing is None:
            raise ValueError("Model forcing 'model_forcing' must be provided")
        if cbar_label is None:
            cbar_label = DS[var_name].attrs.get('units', '')
        if dt0 is not None and dtN is not None:
            DS = DS.sel(time=slice(dt0, dtN))
        elif dt0 is not None:
            DS = DS.sel(time=dt0, method='nearest')
        elif dtN is not None:
            DS = DS.sel(time=dtN, method='nearest')
        if hemisphere:
            hemisphere = self.define_hemisphere(hemisphere)
        else:
            hemisphere = 'sh'
        if hemisphere == 'sh' and region is None:
            region = [0, 360, -90, -50]
        elif hemisphere == 'nh' and region is None:
            region = [0, 360, 50, 90]
        if var_name in self.plot_var_dict.keys():
            var_dict = self.plot_var_dict[var_name]
            if G_ is None:
                if var_dict['grid'].lower() == 't':
                    ln = self.G_AFIM['tlon'].values[0, :] * 180 / np.pi
                    lt = self.G_AFIM['tlat'].values[:, 0] * 180 / np.pi
                elif var_dict['grid'].lower() == 'u':
                    ln = self.G_AFIM['ulon'].values[0, :] * 180 / np.pi
                    lt = self.G_AFIM['ulat'].values[:, 0] * 180 / np.pi
                else:
                    raise ValueError("grid_type must be either 't' or 'u'")
            else:
                if lon_name and lat_name:
                    ln = G_[lon_name]
                    lt = G_[lat_name]
                else:
                    raise ValueError("User must provide the names of geospatial coordinates in the grid file")
            if cmap is None:
                cmap = eval(f"cm.cm.{var_dict['cmap']}")
            if vmin is None:
                vmin = var_dict['vmin']
            if vmax is None:
                vmax = var_dict['vmax']
            if cbar_label is None:
                cbar_label = var_dict['units']
            for idx in range(len(DS.time)):
                if var_name=='uvel':
                    da = np.sqrt( DS['uvel'].isel(time=idx)**2 , DS['vvel'].isel(time=idx)**2 )
                elif var_name=='usurf':
                    da = np.sqrt( DS['usurf'].isel(time=idx)**2 , DS['vsurf'].isel(time=idx)**2 )
                elif var_name=='u':
                    da = np.sqrt( DS['u'].isel(time=idx)**2 , DS['v'].isel(time=idx)**2 )
                elif var_name in ['sst','surface_pot_temp']:
                    da = DS[var_name].isel(time=idx) - 273.15
                else:
                    da = DS[var_name].isel(time=idx)
                t_str   = da.time.values.astype('datetime64[D]').astype(str)
                tit_str = f"{model_name} {model_forcing} {var_dict['name']} {t_str}"
                D_save  = os.path.join(self.D_graphical, model_name, model_forcing, ice_or_ocean, var_name)
                F_save  = f"{t_str}_{model_name}_{model_forcing}_{var_name}_{hemisphere.upper()}.png"
                fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=150, facecolor="w", subplot_kw=dict(projection=projection))
                ax.set_extent(region, ccrs.PlateCarree())
                theta          = np.linspace(0, 2 * np.pi, 100)
                center, radius = [0.5, 0.5], 0.5
                verts          = np.vstack([np.sin(theta), np.cos(theta)]).T
                circle         = mpath.Path(verts * radius + center)
                ax.set_boundary(circle, transform=ax.transAxes)
                ax.add_feature(cft.LAND, color='darkgrey')
                ax.add_feature(cft.COASTLINE, linewidth=0.5)
                pcm = ax.pcolormesh(ln, lt, da, cmap=cmap, vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree(), shading='auto')
                plt.colorbar(pcm, label=cbar_label)
                plt.title(tit_str)
                ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=0.25, linestyle='--')
                if save:
                    if not os.path.exists(D_save):
                        os.makedirs(D_save)
                    P_save = os.path.join(D_save, F_save)
                    if not os.path.exists(P_save) or overwrite:
                        plt.savefig(P_save)
                if show:
                    plt.show()
                plt.close()
        else:
            print(f"{var_name} is not in {self.plot_var_dict.keys()}, can't handle variable")

    #####################################################################################
    def plot_dif_map(self, AOM2=None, AFIM=None, var_name=None, aom2_name=None, afim_name=None,
                     ice_or_ocean='ice', dt0=None, dtN=None, G_=None,
                     grid_type='t', region=None, overwrite=False, lon_name=None, lat_name=None,
                     figsize=(15, 9), hemisphere=None, projection=ccrs.SouthPolarStereo(),
                     cmap=None, vmin=None, vmax=None, cbar_label=None, tit_str=None, D_save=None,
                     F_save=None, save=True, show=True):
        """
        Plot a basic map with a specified projection, color map, and region.

        Parameters:
        ----------
        DS : xarray.Dataset
            Data to plot.
        var_name : str
            Name of the variable.
        ice_or_ocean : str, default 'ice'
            Type of data ('ice' or 'ocean').
        dt0 : str, optional
            Start date/time string in the format YYYY-MM-DD HH:MM for the first time stamp to plot.
        dtN : str, optional
            End date/time string in the format YYYY-MM-DD HH:MM for the last time step to plot.
        G_ : xarray.Dataset, optional
            Grid dataset containing latitude and longitude coordinates.
        grid_type : str, default 't'
            Type of grid ('t' or 'u').
        region : list, optional
            Region to plot [lon_min, lon_max, lat_min, lat_max]; default is None.
        overwrite : bool, default False
            Whether to overwrite existing files.
        lon_name : str, optional
            Name of the longitude variable in the grid file.
        lat_name : str, optional
            Name of the latitude variable in the grid file.
        figsize : tuple, default (15, 9)
            Size of the figure.
        hemisphere : str, optional
            Hemisphere ('nh' or 'sh').
        projection : cartopy.crs.Projection, default ccrs.SouthPolarStereo()
            Cartopy projection to use.
        cmap : matplotlib.colors.Colormap, optional
            Colormap to use.
        vmin : float, optional
            Minimum value for colormap.
        vmax : float, optional
            Maximum value for colormap.
        cbar_label : str, optional
            Label for the colorbar.
        tit_str : str, optional
            Title of the plot.
        D_save : str, optional
            Directory to save the plot.
        F_save : str, optional
            Filename to save the plot.
        save : bool, default True
            Whether to save the plot.
        show : bool, default True
            Whether to display the plot.
        """
        AFIM = AFIM.sortby('time')
        AOM2 = AOM2.sortby('time')
        if AOM2 is None:
            raise ValueError("Dataset 'DS' must be provided")
        if var_name is None:
            raise ValueError("Variable name 'var_name' must be provided")
        if hemisphere:
            hemisphere = self.define_hemisphere(hemisphere)
        else:
            hemisphere = 'sh'
        if hemisphere == 'sh' and region is None:
            region = [0, 360, -90, -50]
        elif hemisphere == 'nh' and region is None:
            region = [0, 360, 50, 90]
        if var_name in self.plot_var_dict.keys():
            var_dict = self.plot_var_dict[var_name]
            if G_ is None:
                if var_dict['grid'].lower() == 't':
                    ln = self.G_AFIM['tlon'].values[0, :] * 180 / np.pi
                    lt = self.G_AFIM['tlat'].values[:, 0] * 180 / np.pi
                elif var_dict['grid'].lower() == 'u':
                    ln = self.G_AFIM['ulon'].values[0, :] * 180 / np.pi
                    lt = self.G_AFIM['ulat'].values[:, 0] * 180 / np.pi
                else:
                    raise ValueError("grid_type must be either 't' or 'u'")
            else:
                if lon_name and lat_name:
                    ln = G_[lon_name]
                    lt = G_[lat_name]
                else:
                    raise ValueError("User must provide the names of geospatial coordinates in the grid file")
            if cmap is None:
                cmap = eval(f"cm.cm.{var_dict['cmap']}")
            if var_name in ['aice','daidtt','daidtd','iage']:
                vmin = -1
                vmax = 1
            elif var_name in ['hi','hs','uvel','strairx','strintx']:
                vmin = -2
                vmax = 2
            elif var_name in ['congel','dvidtt','dvidtd','frazil','meltt','melts','meltb','snoice','sst','sss']:
                vmin = -10
                vmax = 10
            if cbar_label is None:
                cbar_label = var_dict['units']
            for idx in range(len(AFIM.time)):
                if var_name=='uvel':
                    afim = np.sqrt( AFIM['uvel'].isel(time=idx)**2 , AFIM['vvel'].isel(time=idx)**2 )
                    aom2 = np.sqrt( AOM2['uvel'].isel(time=idx)**2 , AOM2['vvel'].isel(time=idx)**2 )
                    da = aom2 - afim
                elif var_name=='usurf':
                    afim = np.sqrt( AFIM['usurf'].isel(time=idx)**2 , AFIM['vsurf'].isel(time=idx)**2 )
                    aom2 = np.sqrt( AOM2['usurf'].isel(time=idx)**2 , AOM2['vsurf'].isel(time=idx)**2 )
                    aom2_renamed = aom2.rename({'yt_ocean': 'nj', 'xt_ocean': 'ni'})
                    aom2_aligned, afim_aligned = xr.align(aom2_renamed, afim, join='inner')
                    da = aom2_aligned - afim_aligned
                elif var_name in ['sst','surface_pot_temp']:
                    if aom2_name=='AOM2_era5':
                        aom2 = AOM2["surface_pot_temp"].isel(time=idx) - 273.15
                    else:
                        aom2 = AOM2[var_name].isel(time=idx) - 273.15
                    aom2_renamed = aom2.rename({'yt_ocean': 'nj', 'xt_ocean': 'ni'})
                    afim = AFIM[var_name].isel(time=idx)
                    aom2_aligned, afim_aligned = xr.align(aom2_renamed, afim, join='inner')
                    da = aom2_aligned - afim_aligned
                elif var_name in ['sss','surface_salt']:
                    if aom2_name=='AOM2_era5':
                        aom2 = AOM2['surface_salt'].isel(time=idx)
                    else:
                        aom2 = AOM2[var_name].isel(time=idx)
                    aom2_renamed = aom2.rename({'yt_ocean': 'nj', 'xt_ocean': 'ni'})
                    afim = AFIM[var_name].isel(time=idx)
                    aom2_aligned, afim_aligned = xr.align(aom2_renamed, afim, join='inner')
                    da = aom2_aligned - afim_aligned
                else:
                    afim = AFIM[var_name].isel(time=idx)
                    aom2 = AOM2[var_name].isel(time=idx) 
                    da = aom2 - afim
                t_str   = AFIM.isel(time=idx).time.values.astype('datetime64[D]').astype(str)
                tit_str = f"{aom2_name} minus {afim_name} {var_dict['name']} {t_str}"
                fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=150, facecolor="w", subplot_kw=dict(projection=projection))
                ax.set_extent(region, ccrs.PlateCarree())
                theta          = np.linspace(0, 2 * np.pi, 100)
                center, radius = [0.5, 0.5], 0.5
                verts          = np.vstack([np.sin(theta), np.cos(theta)]).T
                circle         = mpath.Path(verts * radius + center)
                ax.set_boundary(circle, transform=ax.transAxes)
                ax.add_feature(cft.LAND, color='darkgrey')
                ax.add_feature(cft.COASTLINE, linewidth=0.5)
                pcm = ax.pcolormesh(ln, lt, da, cmap=cm.cm.balance, vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree(), shading='auto')
                plt.colorbar(pcm, label=cbar_label)
                plt.title(tit_str)
                ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=0.25, linestyle='--')
                if save:
                    D_save = os.path.join(self.D_graphical, "dif_maps", f"{aom2_name}_minus_{afim_name}_{var_name}")
                    F_save = f"{t_str}_{aom2_name}_minus_{afim_name}_{var_name}_{hemisphere.upper()}.png"
                    if not os.path.exists(D_save):
                        os.makedirs(D_save)
                    P_save = os.path.join(D_save, F_save)
                    if not os.path.exists(P_save) or overwrite:
                        print(f"{P_save} does not exist *or* we are overwriting existing")
                        plt.savefig(P_save)
                    else:
                        print(f"{P_save} exists and not overwriting")
                if show:
                    plt.show()
                plt.close()
        else:
            print(f"{var_name} is not in {self.plot_var_dict.keys()}, can't handle variable")

    #####################################################################################
    def plot_polar_contour_panel_map(self, DS, var_name, var_dict, d_str, hemisphere=None,
                                     save=True, D_graph_base=None, overwrite=False):
        """
        Plot a panel of polar contour maps for a specified variable from a given dataset.

        Parameters:
            - DS: xarray Dataset containing the data to be plotted.
            - var_name: Name of the variable in the dataset to be plotted.
            - var_dict: Dictionary containing properties of the variable such as colormap, vmin, vmax, etc.
            - d_str: Date string used in the filename.
            - hemisphere: Specifies the hemisphere for plotting. Default is None.
            - save: Boolean flag indicating whether to save the plot. Default is True.
            - D_graph_base: Base directory where the plot will be saved. Default is None.
            - overwrite: Boolean flag indicating whether to overwrite existing plot files. Default is False.

        Returns:
            None

        Example:
            plotter = Plotter()  # Assuming Plotter is a class containing this method
            plotter.plot_polar_contour_panel_map(DS, "var_name", var_dict, "1993-01-01", hemisphere="SH")
        """
        if D_graph_base is None: D_graph_base = self.D_graphical
        lon_name       = f"{var_dict['grid']}LON"
        lat_name       = f"{var_dict['grid']}LAT"
        G_LON          = DS[lon_name].values
        G_LAT          = DS[lat_name].values
        vmin           = var_dict['vmin']
        vmax           = var_dict['vmax']
        cmap           = var_dict['cmap']
        cmap           = eval(f"cm.cm.{cmap}")
        D_save         = os.path.join(D_graph_base,var_name,hemisphere)
        F_save         = f"{d_str}_{var_name}_cice_{hemisphere}.png"
        P_save         = os.path.join(D_save,F_save)
        theta          = np.linspace(0, 2*np.pi, 100)
        center, radius = [0.5, 0.5], 0.5
        verts          = np.vstack([np.sin(theta), np.cos(theta)]).T
        circle         = mpath.Path(verts * radius + center)
        subplot_kw     = {'projection': None}
        if not os.path.exists(D_save):
            os.makedirs(D_save)
        if not overwrite and os.path.exists(P_save):
            print(f"file exists and not overwriting: {P_save}")
            return
        if hemisphere == 'SH':
            subplot_kw['projection'] = ccrs.SouthPolarStereo()
        else:
            subplot_kw['projection'] = ccrs.NorthPolarStereo()
        fig, axs = plt.subplots(nrows=3, ncols=1, subplot_kw=subplot_kw, figsize=(5, 15))
        axs = axs.flatten()
        for i, model in enumerate(DS.ocn_frcg_time.values):
            DA = DS.isel(ocn_frcg_time=i, time=0)[var_name].values
            if var_name=='uvel':
                DA2 = DS.isel(ocn_frcg_time=i, time=0)["vvel"].values
                DA  = np.sqrt( DA**2 + DA2**2 )
            elif var_name=='strairx':
                DA2 = DS.isel(ocn_frcg_time=i, time=0)["strairy"].values
                DA  = np.sqrt( DA**2 + DA2**2 )
            elif var_name=='strintx':
                DA2 = DS.isel(ocn_frcg_time=i, time=0)["strinty"].values
                DA  = np.sqrt( DA**2 + DA2**2 )
            elif var_name=='taubx':
                DA2 = DS.isel(ocn_frcg_time=i, time=0)["tauby"].values
                DA  = np.sqrt( DA**2 + DA2**2 )
            if hemisphere == 'SH':
                axs[i].set_extent([0,360,-90,-50], ccrs.PlateCarree())
            else:
                axs[i].set_extent([0,360, 50, 90], ccrs.PlateCarree())
            axs[i].set_boundary(circle, transform=axs[i].transAxes)
            axs[i].add_feature(cft.LAND, color='darkgrey')
            axs[i].add_feature(cft.COASTLINE, linewidth=.5)
            axs[i].set_title(f'ocn_frcg_time : {model}')
            cs = axs[i].pcolormesh(G_LON, G_LAT, DA, transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax, cmap=cmap )
        fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.9, wspace=0.05, hspace=0.05)
        cbar_ax = fig.add_axes([0.9, 0.2, 0.05, 0.65])
        cbar = fig.colorbar(cs, cax=cbar_ax, orientation='vertical')
        plt.suptitle(f"{var_dict['name']}, {d_str}")
        if save:
            print(f"saving: {P_save}")
            plt.savefig(P_save)
            plt.close()
        else:
            plt.show()
            plt.close()

    #####################################################################################
    def plot_polar_differences_contour_panel_map(self, DS, var_name, var_dict, d_str, hemisphere=None,
                                                 save=True, D_graph_base=None, overwrite=False):
        """
        Plot a panel of polar contour maps for a specified variable from a given dataset.

        Parameters:
            - DS: xarray Dataset containing the data to be plotted.
            - var_name: Name of the variable in the dataset to be plotted.
            - var_dict: Dictionary containing properties of the variable such as colormap, vmin, vmax, etc.
            - d_str: Date string used in the filename.
            - hemisphere: Specifies the hemisphere for plotting. Default is None.
            - save: Boolean flag indicating whether to save the plot. Default is True.
            - D_graph_base: Base directory where the plot will be saved. Default is None.
            - overwrite: Boolean flag indicating whether to overwrite existing plot files. Default is False.

        Returns:
            None

        Example:
            plotter = Plotter()  # Assuming Plotter is a class containing this method
            plotter.plot_polar_contour_panel_map(DS, "var_name", var_dict, "1993-01-01", hemisphere="SH")
        """
        if D_graph_base is None: D_graph_base = self.D_graphical
        lon_name       = f"{var_dict['grid']}LON"
        lat_name       = f"{var_dict['grid']}LAT"
        G_LON          = DS[lon_name].values
        G_LAT          = DS[lat_name].values
        vmin           = var_dict['vmin']
        vmax           = var_dict['vmax']
        cmap           = var_dict['cmap']
        D_save         = os.path.join(D_graph_base,'dif_maps',var_name,hemisphere)
        F_save         = f"{d_str}_{var_name}_cice_{hemisphere}.png"
        P_save         = os.path.join(D_save,F_save)
        theta          = np.linspace(0, 2*np.pi, 100)
        center, radius = [0.5, 0.5], 0.5
        verts          = np.vstack([np.sin(theta), np.cos(theta)]).T
        circle         = mpath.Path(verts * radius + center)
        subplot_kw     = {'projection': None}
        if not os.path.exists(D_save):
            os.makedirs(D_save)
        if not overwrite and os.path.exists(P_save):
            print(f"file exists and not overwriting: {P_save}")
            return
        if hemisphere == 'SH':
            subplot_kw['projection'] = ccrs.SouthPolarStereo()
        else:
            subplot_kw['projection'] = ccrs.NorthPolarStereo()
        fig, axs = plt.subplots(nrows=3, ncols=1, subplot_kw=subplot_kw, figsize=(5, 15))
        axs = axs.flatten()
        for i, model in enumerate(DS.ocn_frcg_time.values):
            DA = DS.isel(ocn_frcg_time=i, time=0).values
            if hemisphere == 'SH':
                axs[i].set_extent([0,360,-90,-50], ccrs.PlateCarree())
            else:
                axs[i].set_extent([0,360, 50, 90], ccrs.PlateCarree())
            axs[i].set_boundary(circle, transform=axs[i].transAxes)
            axs[i].add_feature(cft.LAND, color='darkgrey')
            axs[i].add_feature(cft.COASTLINE, linewidth=.5)
            axs[i].set_title(f'ocn_frcg_time : {model}')
            cs = axs[i].pcolormesh(G_LON, G_LAT, DA, transform=ccrs.PlateCarree(), cmap=cm.cm.balance )
        fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.9, wspace=0.05, hspace=0.05)
        cbar_ax = fig.add_axes([0.9, 0.2, 0.05, 0.65])
        cbar = fig.colorbar(cs, cax=cbar_ax, orientation='vertical')
        plt.suptitle(f"{var_dict['name']}, {d_str}")
        if save:
            print(f"saving: {P_save}")
            plt.savefig(P_save)
            plt.close()
        else:
            plt.show()
            plt.close()

    #####################################################################################
    @staticmethod
    def animate_sequenced_figures(D_seq_fig       = '',
                                  D_ani_base      = '',
                                  var_name        = '',
                                  interval        = '-delay 100',
                                  img_in_type     = 'png',
                                  img_out_type    = 'gif',
                                  mean_length_str = '',
                                  model_name      = '',
                                  model_run_date  = '',
                                  atm_frcg_name   = '',
                                  ocn_frcg_name   = '',
                                  grid_res_str    = '',
                                  overwrite       = False):
        '''
        Convert a sequence of still images into an animation format (GIF, MP4, or both).

        Parameters:
        -----------
        D_seq_fig : str
            Directory containing the sequence of still images to be animated.

        D_ani_base : str
            Base directory where the animation file will be saved.

        var_name : str
            Variable name.

        interval : str, default='-delay 100'
            Time interval between frames in the animation (specific to ImageMagick's `convert` command).

        img_in_type : str, default='png'
            Image file extension of the input still images.

        img_out_type : str, default='gif'
            Desired animation format (can be 'gif', 'mp4', or 'both').

        mean_length_str : str
            String representing the averaging length/period.

        model_name : str
            Name of the model being used.

        model_run_date : str
            Date when the model was run.

        self.atm_frcg_name : str
            Atmospheric forcing name.

        ocn_frcg_name : str
            Ocean forcing name.

        grid_res_str : str
            String indicating the grid resolution.

        overwrite : bool, default=False
            If True, will overwrite existing animation files. If False, will skip if file exists.

        Returns:
        --------
        None

        Description:
        ------------
        This method converts a sequence of still images located in `D_seq_fig` into an animation format
        based on the `img_out_type` parameter. The resulting animation will be saved in a directory constructed
        from various model parameters. If the output animation file already exists and `overwrite` is set to False,
        the animation process will be skipped.
        '''
        c0 = time.process_time()
        D_ = construct_graphical_output_directory(directory_base  = D_ani_base,
                                                    model_name      = model_name,
                                                    model_run_date  = model_run_date,
                                                    mean_length_str = mean_length_str,
                                                    var_name        = var_name)
        if img_out_type=='gif':
            F_ = construct_graphical_output_filename(atm_frcg_name   = self.atm_frcg_name,
                                            ocn_frcg_name   = ocn_frcg_name,
                                            grid_res_str    = grid_res_str,
                                            img_type        = img_out_type,
                                            dt_str          = model_run_date)
            if os.path.exists(os.path.join(D_,F_)) and not overwrite:
                print(f'{os.path.join(D_,F_)} exists and not overwriting')
                return
            print(f'converting stills in {D_seq_fig}/*.{img_in_type}\n to animation file {os.path.join(D_,F_)}\n{(time.process_time()-c0):.3f}')
            os.system(f"convert {interval} -loop 0 {D_seq_fig}/*.{img_in_type} -duplicate 1,-2-1 {os.path.join(D_,F_)}")
        elif img_out_type=='mp4':
            F_ = construct_graphical_output_filename(atm_frcg_name   = self.atm_frcg_name,
                                            ocn_frcg_name   = ocn_frcg_name,
                                            grid_res_str    = grid_res_str,
                                            img_type        = img_out_type,
                                            dt_str          = model_run_date)
            if os.path.exists(os.path.join(D_,F_)) and not overwrite:
                print(f'{os.path.join(D_,F_)} exists and not overwriting')
                return
            print(f'converting stills in {D_seq_fig}/*.{img_in_type}\n to animation file {os.path.join(D_,F_)}\n{(time.process_time()-c0):.3f}')
            os.system(f"ffmpeg -framerate 1 -pattern_type glob -i {D_seq_fig}/*.{img_in_type} -c:v libx264 -r 30 -pix_fmt yuv420p {os.path.join(D_,F_)}")
        elif img_out_type=='both':
            F_ = construct_graphical_output_filename(atm_frcg_name   = self.atm_frcg_name,
                                            ocn_frcg_name   = ocn_frcg_name,
                                            grid_res_str    = grid_res_str,
                                            img_type        = 'gif',
                                            dt_str          = model_run_date)
            if os.path.exists(os.path.join(D_,F_)) and not overwrite:
                print(f'{os.path.join(D_,F_)} exists and not overwriting')
            else:
                print(f'converting stills in {D_seq_fig}/*.{img_in_type}\n to animation file {os.path.join(D_,F_)}\n{(time.process_time()-c0):.3f}')
                os.system(f"convert {interval} -loop 0 {D_seq_fig}/*.{img_in_type} -duplicate 1,-2-1 {os.path.join(D_,F_)}")
            F_ = construct_graphical_output_filename(atm_frcg_name   = self.atm_frcg_name,
                                            ocn_frcg_name   = ocn_frcg_name,
                                            grid_res_str    = grid_res_str,
                                            img_type        = 'mp4',
                                            dt_str          = model_run_date)
            if os.path.exists(os.path.join(D_,F_)) and not overwrite:
                print(f'{os.path.join(D_,F_)} exists and not overwriting')
            else:
                print(f'converting stills in {D_seq_fig}/*.{img_in_type}\n to animation file {os.path.join(D_,F_)}\n{(time.process_time()-c0):.3f}')
                os.system(f"ffmpeg -framerate 1 -pattern_type glob -i '{D_seq_fig}/*.{img_in_type}' {os.path.join(D_,F_)}")
        print(f'finished animating {(time.process_time()-c0):.3f}')
        os.system(f"ls -lh {os.path.join(D_,F_)}")


    # def compute_NSIDC_sia(self):
    #     """
    #     """
    #     if self.proc_freq=="monthly":
    #         NSIDC_conc = self.NSIDC['cdr_seaice_conc_monthly']
    #     else:
    #         NSIDC_conc = self.NSIDC['cdr_seaice_conc']
    #     for flag in self.NSIDC_flags:
    #         NSIDC_conc      = xr.where(NSIDC_conc == flag/100, np.nan, NSIDC_conc)  # Using flag/100 because you've already converted percentages to fractions
    #     NSIDC_cell_area = xr.open_dataset(self.P_NSIDC_cell_area)
    #     NSIDC_area      = NSIDC_cell_area['cell_area'].values / 1e12 # Convert cell_area from m^2 to million km^2
    #     NSIDC_ApC       = NSIDC_conc * NSIDC_area
    #     NSIDC_mask      = NSIDC_ApC.where(NSIDC_conc > self.aice_thresh)
    #     self.NSIDC_SIA  = NSIDC_mask.sum(dim=['y','x'])

    # #####################################################################################
    # def compute_NSIDC_sie(self):
    #     """
    #     """
    #     if self.proc_freq=="monthly":
    #         NSIDC_conc = self.NSIDC['cdr_seaice_conc_monthly']
    #     else:
    #         NSIDC_conc = self.NSIDC['cdr_seaice_conc']
    #     for flag in self.NSIDC_flags:
    #         NSIDC_conc      = xr.where(NSIDC_conc == flag/100, np.nan, NSIDC_conc)  # Using flag/100 because you've already converted percentages to fractions
    #     NSIDC_cell_area = xr.open_dataset(self.P_NSIDC_cell_area)
    #     NSIDC_area      = NSIDC_cell_area['cell_area'].values / 1e12 # Convert cell_area from m^2 to million km^2
    #     NSIDC_mask      = NSIDC_conc > self.aice_thresh
    #     SIE             = (NSIDC_mask * NSIDC_area).sum(dim=['y', 'x'])
    #     self.NSIDC_SIE  = SIE



