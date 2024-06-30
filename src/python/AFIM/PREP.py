import os, sys, json
import xarray as xr
import numpy  as np
import pandas as pd

class prep:
    '''
    prep: A utility for preparing datasets for CICE6 input.

    The `prep` class facilitates the generation of input datasets for CICE6 
    from various sources, such as BRAN climatology or ERA5 climatology. The class 
    relies on a JSON configuration file to determine its parameters.

    Internal Dependencies:
    ----------------------
    It's recommended to use the Anaconda package manager to handle these dependencies. 

    Note: Exact version compatibility is not rigorously tracked, so some issues might arise from version mismatches.

    - xarray
    - metpy
    - pygmt
    - numpy
    - pandas
    - json

    Initialization:
    ---------------
    The class initialization requires a JSON configuration file with the following key fields:

    Basic Configurations:
    - CICE_ver: Version of CICE. (Default: 6)
    - start_date: Date to start data consideration. (Default: "2010-01-01")
    - n_months: Number of months for data consideration from the start date. (Default: 12)
    - n_years: Duration in years. (Default: 2)
    - G_res: Grid resolution of CICE, used in filenames and directories. (Default: "0p1")
    - regrid_interp_type: Interpolation method for regridding. (Default: "bilinear")

    Directory Paths:
    - D_data: Base directory for data output.
    - D_BRAN: Location of BRAN climatology data.
    - D_ERA5: Location of ERA5 climatology data.
    - D_modin: Used in earlier module versions, may not be necessary.
    - D_graph: Directory where plots/figures/animations are stored.
    - D_access-om2_out: Specifies which ACCESS-OM2 cycle data is considered.

    File Paths:
    - F_amo2bath: File path for bathymetry (may not be currently used, but relevant for potential plotting).
    - F_G_CICE: Contains the spherical grid for CICE.
    - F_G_BRAN: Contains the spherical grid for BRAN climatology.
    - F_G_ERA5: Contains the spherical grid for ERA5 climatology.
    - F_BRAN_weights: After generating using the `prep.esmf_generate_weights` method, this file is used for BRAN climatology regridding.
    - F_ERA5_weights: Similarly, used for ERA5 climatology regridding after generation.

    And other parameters specific to the datasets and grids involved.

    Parameters:
    -----------
    - FJSON (str): Path to the JSON configuration file.
    - **kwargs: Additional keyword arguments.
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
                            :: is used in the method prep.esmf_generate_weights(), however, that
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
                            :: prep.esmf_generate_weights this is the name of the file used in the 
                            :: method prep.regrid_climatology() for BRAN climatology
           ERA5 weights     :: "F_ERA5_weights", once created via the AFIM python module
                            :: prep.esmf_generate_weights this is the name of the file used in the 
                            :: method prep.regrid_climatology() for ERA5 climatology
        '''
        # read in the JSON file
        self.FJSON = FJSON
        with open(FJSON) as f:
            PARAMS = json.load(f)
        # logging
        self.F_log              = PARAMS['F_log']
        logging.basicConfig(filename=self.F_log, encoding='utf-8', level=logging.DEBUG)
        # miscellaneous
        self.CICE_ver           = PARAMS['CICE_ver']
        self.start_date         = PARAMS['start_date']
        self.n_months           = PARAMS['n_months']
        self.n_years            = PARAMS['n_years']
        self.ncrcat_opts_t_dim  = PARAMS["ncrcat_opts_t_dim"]
        self.ncrcat_opts_n_dim  = PARAMS["ncrcat_opts_n_dim"]
        # dataset specific
        self.acom2_model_run    = PARAMS["acom2_model_run"]
        self.acom2_start_date   = PARAMS['acom2_start_date']
        self.acom2_var_freq     = PARAMS['acom2_var_freq']
        self.acom2_chunking     = PARAMS['acom2_chunking']
        self.ERA5_chunking      = PARAMS['ERA5_chunking']
        self.ERA5_start_date    = PARAMS['ERA5_start_date']
        self.ERA5_regen_weights = PARAMS['ERA5_regen_weights']
        self.ERA5_regrid_method = PARAMS['ERA5_regrid_method']
        self.ERA5_periodicity   = PARAMS['ERA5_periodicity']
        self.BRAN_t_chunking    = PARAMS['BRAN_t_chunking']
        self.BRAN_u_chunking    = PARAMS['BRAN_u_chunking']
        self.BRAN_start_date    = PARAMS['BRAN_start_date']
        self.BRAN_regen_weights = PARAMS['BRAN_regen_weights']
        self.BRAN_regrid_method = PARAMS['BRAN_regrid_method']
        self.BRAN_periodicity   = PARAMS['BRAN_periodicity']
        self.BRAN_var_freq      = PARAMS['BRAN_var_freq']
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
        self.F_G_ACOM2       = PARAMS['F_G_ACOM2']
        self.F_G_CICE        = PARAMS['F_G_CICE']
        self.F_G_CICE_vars   = PARAMS['F_G_CICE_vars']
        self.F_G_CICE_original = PARAMS['F_G_CICE_original']
        self.F_G_BRAN        = PARAMS['F_G_BRAN']
        self.F_G_ERA5        = PARAMS['F_G_ERA5']
        self.F_ERA5_weights  = PARAMS['F_ERA5_weights']
        self.F_BRAN_t_weights= PARAMS['F_BRAN_t_weights']
        self.F_BRAN_u_weights= PARAMS['F_BRAN_u_weights']
        self.F_ERA5_reG_form = PARAMS['F_ERA5_reG_form']
        self.F_ERA5_reG_form_tmp = PARAMS["F_ERA5_reG_form_tmp"]
        self.F_BRAN_reG_form = PARAMS['F_BRAN_reG_form']
        
    ############################################################################################
    def time_series_day_month_start_or_end(self,start_date='',n_months='',n_years='',month_start=True):
        '''
        Generates a date range based on given parameters.

        Parameters:
        - start_date (str, default=''): The starting date. If not provided, will use the class's start_date.
        - n_months (str, default=''): Number of months for the period. If not provided, will use the class's n_months.
        - n_years (str, default=''): Number of years for the period. If not provided, will use the class's n_years.
        - month_start (bool, default=True): If True, will start the date range at the beginning of the month; otherwise, will start at the end of the month.

        Returns:
        - DateRange: A Pandas date range object based on the provided or default parameters.
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
    def define_datetime_object(self,start_date='',year_offset=0,full_datetime=False):
        '''
        Returns a datetime object based on the given parameters.

        Parameters:
        - start_date (str, default=''): The date to begin with. If not provided, will use the class's start_date.
        - year_offset (int, default=0): The number of years to offset from the given start_date.
        - full_datetime (bool, default=False): If True, expects start_date in the format '%Y-%m-%d %H:%M'. If False, expects '%Y-%m-%d'.

        Returns:
        - datetime_obj: A datetime object based on the provided or default parameters with the year_offset applied.
        '''
        if not start_date: start_date = self.start_date
        if full_datetime:
            return datetime.strptime(start_date,'%Y-%m-%d %H:%M').date() + pd.DateOffset(years=year_offset)
        else:   
            return datetime.strptime(start_date,'%Y-%m-%d').date() + pd.DateOffset(years=year_offset)
    
    ############################################################################################
    def year_string_from_datetime_object(self,dt_obj):
        '''
        Returns a year string from a datetime object.

        Parameters:
        - dt_obj (datetime): A datetime object.

        Returns:
        - str: A string representing the year from the provided datetime object.
        '''
        return dt_obj.strftime('%Y')
    
    ############################################################################################
    def define_era5_variable_path(self,var_name='',stt='',stp=''):
        '''
        Constructs and returns the file path for a specified ERA5 variable.

        Parameters:
        - var_name (str, default=''): The name of the ERA5 variable.
        - stt (str, default=''): The start date/time. If not provided, will use the class's time_series_day_month_start_or_end() function.
        - stp (str, default=''): The stop date/time. If not provided, will use the class's time_series_day_month_end() function.

        Returns:
        - str: A constructed file path string for the specified ERA5 variable.
        '''
        if not stt:
            stt = self.time_series_day_month_start_or_end()
        if not stp:
            ptp = self.time_series_day_month_end()
        F = '{var:s}_era5_oper_sfc_{stt:s}-{stp:s}.nc'.format(var=var_name,stt=stt,stp=stp)
        return os.path.join(self.D_ERA5,var_name,self.yr_str,F)
    
    ############################################################################################
    def era5_load_and_regrid(self,era5_var='2t',yr_str='', mo_str='', D_base='', G_dst='',
                             regrid_method='', regrid_periodicity='',cnk_dict='',
                             generate_weights=False, weights_file_name=''):
        '''
        '''
        # optional inputs
        if not G_dst:              G_dst              = self.cice_grid_prep()
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
        logging.info("Source file looks like this: %s",G_ERA5)
        logging.info("Destination file looks like this: %s",G_dst)
        if generate_weights or not os.path.exists(F_weights):
            logging.info("User has requested to generate weights file *or* %s does not exist",F_weights)
            logging.info("Creating weights can take some time ... sometimes it can take a long time!!")
            logging.info("Intialising xesmf regridder")
            logging.info("Using xesmf method %s",regrid_method)
            logging.info("Full path to weight file is: %s",F_weights)
            rg = xe.Regridder(G_ERA5, G_dst, method=regrid_method,  periodic=regrid_periodicity, filename=F_weights, reuse_weights=False)
        else:
            logging.info("REUSING EXISTING WEIGHT FILE: %s",F_weights)
            logging.info("Using xesmf method %s",regrid_method)
            rg = xe.Regridder(G_ERA5, G_dst, method=regrid_method,  periodic=regrid_periodicity, filename=F_weights, reuse_weights=True)
        dat_n = self.era5_load(era5_var=era5_var, D_base=D_ERA5, yr_str=yr_str, mo_str='', cnk_dict=cnk_dict )
        print("regridding ERA5")
        return rg(dat_n)
    
    ############################################################################################
    def era5_load(self,era5_var='2t',D_base='',yr_str='',mo_str='',cnk_dict=''):
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
        if mo_str:
            P_dat = os.path.join(D_ERA5,era5_var,yr_str,'{:s}_era5_oper_sfc_{:s}{:s}*.nc'.format(era5_var,yr_str,mo_str))
        else:
            P_dat = os.path.join(D_ERA5,era5_var,yr_str,'*.nc')
        print("loading ERA5: {:s}".format(P_dat))
        prt_str = "".join(str(key) + str(value) for key, value in cnk_dict.items())
        logging.info("chunking dictionary: %s",prt_str)
        return xr.open_mfdataset( P_dat , chunks=cnk_dict )
     
    ############################################################################################
    def bran_load_and_regrid(self, bran_var='temp', yr_str='', mo_str='', D_base='', G_dst='', var_freq='', grid_type='',
                             regrid_method='', regrid_periodicity='', cnk_dict='',
                             generate_weights=False, weights_file_name=''):
        '''
        '''
        # optional inputs
        if not grid_type:          grid_type          = 't'
        if not G_dst:              G_dst              = self.cice_grid_prep()
        if not D_base:             D_BRAN             = self.D_BRAN
        if not regrid_method:      regrid_method      = self.BRAN_regrid_method
        if not regrid_periodicity: regrid_periodicity = self.BRAN_periodicity
        if not var_freq:           var_freq           = self.BRAN_var_freq
        if not weights_file_name:
            if grid_type=='t':
                F_weights = self.F_BRAN_t_weights
            elif grid_type=='u':
                F_weights = self.F_BRAN_u_weights
        else:
            F_weights = weights_file_name
        if not yr_str:
            dt_obj = define_datetime_object(start_date=self.BRAN_start_date, year_offset=0)
            yr_str = year_string_from_datetime_object(dt_obj)
        if grid_type=='t' and not cnk_dict:
            cnk_dict = self.BRAN_t_chunking
        elif grid_type=='u' and not cnk_dict:
            cnk_dict = self.BRAN_u_chunking
        # pull-in the grids
        G_BRAN = self.bran_grid_prep()
        logging.info("Source file looks like this: %s",G_BRAN)
        logging.info("Destination file looks like this: %s",G_dst)
        if generate_weights or not os.path.exists(F_weights):
            logging.info("User has requested to generate weights file *or* %s does not exist",F_weights)
            logging.info("Creating weights can take some time ... sometimes it can take a long time!!")
            logging.info("Intialising xesmf regridder")
            logging.info("Using xesmf method %s",regrid_method)
            logging.info("Full path to weight file is: %s",F_weights)
            rg = xe.Regridder(G_BRAN, G_dst, method=regrid_method,  periodic=regrid_periodicity, filename=F_weights, reuse_weights=False)
        else:
            logging.info("REUSING EXISTING WEIGHT FILE: %s",F_weights)
            logging.info("Using xesmf method %s",regrid_method)
            rg = xe.Regridder(G_BRAN, G_dst, method=regrid_method,  periodic=regrid_periodicity, filename=F_weights, reuse_weights=True)
        dat_n = self.bran_load(bran_var=bran_var, D_base=D_BRAN, yr_str=yr_str, mo_str=mo_str, var_freq=var_freq, cnk_dict=cnk_dict , grid_type=grid_type)
        print("regridding BRAN")
        return rg(dat_n[bran_var])

    ############################################################################################
    def bran_load(self,bran_var='temp',D_base='',yr_str='',mo_str='',var_freq='',cnk_dict='',grid_type=''):
        '''
        '''
        #optional inputs
        if not D_base:
            D_BRAN = self.D_BRAN
        else:
            D_BRAN = D_base
        if not var_freq: var_freq = self.BRAN_var_freq
        if not yr_str:
            dt_obj = define_datetime_object(start_date=self.BRAN_start_date, year_offset=0)
            yr_str = year_string_from_datetime_object(dt_obj)
        if not grid_type: grid_type = 't'
        if grid_type=='t' and not cnk_dict:
            cnk_dict = self.BRAN_t_chunking
        elif grid_type=='u' and not cnk_dict:
            cnk_dict = self.BRAN_u_chunking
        if mo_str:
            P_dat = os.path.join(D_BRAN,var_freq,'ocean_{var:s}_{yr:s}_{mo:s}.nc'.format(var=bran_var,yr=yr_str,mo=mo_str))
            print("loading BRAN: {:s}".format(P_dat))
            return xr.open_dataset( P_dat )
        else:
            P_dat = os.path.join(D_BRAN,var_freq,'ocean_{var:s}_{yr:s}*.nc'.format(var=bran_var,yr=yr_str))
            print("loading BRAN: {:s}".format(P_dat))
            #prt_str = "".join(str(key) + str(value) for key, value in cnk_dict.items())
            #logging.info("chunking dictionary: %s",prt_str)
            # BRAN chunking works but sticking it into the re-gridder fails!!
            return xr.open_mfdataset( P_dat )#, chunks=cnk_dict )
        
    ############################################################################################
    def era5_grid_prep(self):
        '''
        '''
        return xr.open_dataset(self.F_G_ERA5)
    
    ############################################################################################
    def bran_grid_prep(self,grid_type=''):
        '''
        '''
        if not grid_type: grid_type = 't'
        if grid_type=='t':
            lat_name = 'yt_ocean'
            lon_name = 'xt_ocean'
        elif grid_type=='u':
            lat_name = 'yu_ocean'
            lon_name = 'xu_ocean'
        G_BRAN        = xr.open_dataset(self.F_G_BRAN)
        LN,LT         = np.meshgrid(G_BRAN[lon_name].values,G_BRAN[lat_name].values)
        G_BRAN['lon'] = (['ny','nx'],LN,{'units':'degrees_east'})
        G_BRAN['lat'] = (['ny','nx'],LT,{'units':'degrees_north'})
        G_BRAN        = G_BRAN.drop(('xt_ocean','yt_ocean','xu_ocean','yu_ocean','hu','ht','kmt',
                                     'kmu','umask','tmask','st_edges_ocean','Time','st_ocean'))
        return G_BRAN

    ############################################################################################
    def cice_grid_prep(self,F_grid='',grid_type=''):
        '''
        '''
        if not grid_type: grid_type = 't'
        if grid_type=='t':
            lat = 'tlat'
            lon = 'tlon'
        elif grid_type=='u':
            lat = 'ulat'
            lon = 'ulon'
        if F_grid:
            G_CICE  = xr.open_dataset(F_grid)
        else:
            G_CICE  = xr.open_dataset(self.F_G_CICE_original)
        G_CICE['lat'] = (['ny','nx'],G_CICE[lat].values*(180/np.pi),{'units':'degrees_north'})
        G_CICE['lon'] = (['ny','nx'],G_CICE[lon].values*(180/np.pi),{'units':'degrees_east'})
        if self.G_res=='0p1':
            G_CICE = G_CICE.drop(('ulat','ulon','tlat','tlon','clon_t','clat_t','clat_u','clon_u','angle','uarea'))
        else:
            G_CICE = G_CICE.drop(('ulat','ulon','tlat','tlon','angle','angleT','tarea','uarea'))
        G_CICE = G_CICE.assign_coords(xt_ocean=G_CICE['lon'][0,:],yt_ocean=G_CICE['lat'][:,0])
        G_CICE = G_CICE.set_index(nx='xt_ocean',ny='yt_ocean')
        return G_CICE

    ############################################################################################
    def regrid_era5_for_cice6(self,start_date=''):
        '''
        '''
        if not start_date: start_date = self.ERA5_start_date
        G_CICE    = self.cice_grid_prep()
        D_ATM     = os.path.join(self.D_modin,'ERA5')
        D_ATM_tmp = os.path.join(D_ATM,'monthly')
        for i in range(self.n_years):
            dt_obj = self.define_datetime_object(start_date=start_date, year_offset=i)
            yr_str = self.year_string_from_datetime_object(dt_obj)
            P_ATM  = os.path.join(D_ATM,'ERA5_hourly_forcing_{res:s}_{yr:s}.nc'.format(res=self.G_res,yr=yr_str))
            if os.path.exists(P_ATM):
                print('\tFile exists: {:s}\n\t ... advancing to next iteration'.format(P_ATM))
                continue
            print("\nCreating file: {:s}".format(P_ATM))
            print("\tcurrent time: ",datetime.now().time())
            t2m    = self.era5_load_and_regrid(era5_var = '2t'      , yr_str = yr_str, generate_weights=True)
            lw     = self.era5_load_and_regrid(era5_var = 'msdwlwrf', yr_str = yr_str, generate_weights=True)
            sw     = self.era5_load_and_regrid(era5_var = 'msdwswrf', yr_str = yr_str, generate_weights=True)
            mtpr   = self.era5_load_and_regrid(era5_var = 'mtpr'    , yr_str = yr_str, generate_weights=True)
            u10    = self.era5_load_and_regrid(era5_var = '10u'     , yr_str = yr_str, generate_weights=True)
            v10    = self.era5_load_and_regrid(era5_var = '10v'     , yr_str = yr_str, generate_weights=True)
            d2m    = self.era5_load_and_regrid(era5_var = '2d'      , yr_str = yr_str, generate_weights=True)
            sp     = self.era5_load_and_regrid(era5_var = 'sp'      , yr_str = yr_str, generate_weights=True)
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
                                   '_FillValue':-2e16}),
                      "dlwsfc" : (['time','nj','ni'],lw.msdwlwrf.astype(np.single).data,
                                  {'long_name':"Mean surface downward long-wave radiation flux",
                                   'units'    :"W m**-2",
                                   '_FillValue':-2e16}),
                      "glbrad" : (['time','nj','ni'],sw.msdwswrf.astype(np.single).data,
                                  {'long_name':"Mean surface downward short-wave radiation flux",
                                   'units'    :"W m**-2",
                                   '_FillValue':-2e16}),
                      "spchmd" : (['time','nj','ni'],qsat.astype(np.single).data,
                                  {'long_name':"specific humidity",
                                   'units'    :"kg/kg",
                                   '_FillValue':-2e16}),
                      "ttlpcp" : (['time','nj','ni'],mtpr.mtpr.astype(np.single).data,
                                  {'long_name':"Mean total precipitation rate",
                                   'units'    :"kg m**-2 s**-1",
                                   '_FillValue':-2e16}),
                      "wndewd" : (['time','nj','ni'],u10.u10.astype(np.single).data,
                                  {'long_name':"10 metre meridional wind component",
                                   'units'    :"m s**-1",
                                   '_FillValue':-2e16}),
                      "wndnwd" : (['time','nj','ni'],v10.v10.astype(np.single).data,
                                  {'long_name':"10 metre zonal wind component",
                                   'units'    :"m s**-1",
                                   '_FillValue':-2e16}) }
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
            mo0_dates = pd.date_range(start=dt_obj,freq='MS',periods=self.n_months)
            moN_dates = pd.date_range(start=dt_obj,freq='M',periods=self.n_months)
            #enc_dict  = {'shuffle':True,'zlib':True,'complevel':5}
            for j in mo0_dates:
                dt0_str   = j.strftime('%Y-%m-%d %H:%M')
                dtN_str   = moN_dates[cnt].strftime('%Y-%m-%d %H:%M')
                yrmo0_str = j.strftime('%Y_%m')
                P_ATM_tmp = os.path.join(D_ATM_tmp,self.F_ERA5_reG_form_tmp.format(res=self.G_res,dt_str=yrmo0_str))
                if not os.path.exists(P_ATM_tmp):
                    ATM_tmp   = ATM.sel(time=slice(dt0_str,dtN_str))
                    write_job = ATM_tmp.to_netcdf(P_ATM_tmp,unlimited_dims=['time'],compute=False)#,encoding={'airtmp':enc_dict,
                                                                                                  #          'dlwsfc':enc_dict,
                                                                                                  #          'glbrad':enc_dict,
                                                                                                  #          'spchmd':enc_dict,
                                                                                                  #          'ttlpcp':enc_dict,
                                                                                                  #          'wndewd':enc_dict,
                                                                                                  #          'wndnwd':enc_dict})
                    with ProgressBar():
                        print(f"Writing to {P_ATM_tmp}")
                        write_job.compute()
                cnt+=1
            print("concatenating monthly atmosphere forcing files in {:s}".format(D_ATM_tmp))
            print("into a year-long file {:s}".format(P_ATM))
            print("\tcurrent time: ",datetime.now().time())
            sys_call = 'ncrcat {:s} {:s}/ERA5_hourly_forcing_{:s}_{:s}* {:s}'.format(self.ncrcat_opts_n_dim,
                                                                                     D_ATM_tmp,
                                                                                     self.G_res,
                                                                                     yr_str,
                                                                                     P_ATM)
            os.system(sys_call)
                
    ############################################################################################
    def regrid_bran_for_cice6(self,start_date=''):
        '''
        '''
        if not start_date: start_date = self.BRAN_start_date
        G_CICE = self.cice_grid_prep()
        for i in range(self.n_years):
            dt_obj = self.define_datetime_object(start_date=start_date, year_offset=i)
            yr_str = self.year_string_from_datetime_object(dt_obj)
            D_ocn  = os.path.join(self.D_modin,'BRAN')
            P_ocn  = os.path.join(D_ocn,'bran_ocn_frcg_cice6_{res:s}_{yr:s}.nc'.format(res=self.G_res,yr=yr_str))
            print("\nCreating file: {:s}".format(P_ocn))
            print("\tcurrent time: ",datetime.now().time())
            if os.path.exists(P_ocn):
                print('\tFile exists: {:s}\n\t ... advancing to next iteration'.format(P_ocn))
                continue
            for j in range(12):
                mo_str = '{:02d}'.format(j+1)
                D_ocn_j= os.path.join(D_ocn,'yr_mo')
                P_ocn_j= os.path.join(D_ocn_j,'bran_ocn_frcg_cice6_{res:s}_{yr:s}_{mo:s}.nc'.format(res=self.G_res,yr=yr_str,mo=mo_str))
                print("\tCreating file: {:s}".format(P_ocn_j))
                if os.path.exists(P_ocn_j):
                    print('\tFile exists: {:s}\n\t ... advancing to next iteration'.format(P_ocn_j))
                else:
                    temp   = self.bran_load_and_regrid( bran_var = 'temp' , yr_str = yr_str, mo_str = mo_str, grid_type = 't')
                    salt   = self.bran_load_and_regrid( bran_var = 'salt' , yr_str = yr_str, mo_str = mo_str, grid_type = 't')
                    mld    = self.bran_load_and_regrid( bran_var = 'mld'  , yr_str = yr_str, mo_str = mo_str, grid_type = 't')
                    print("\tcomputing temperature derivative over time")
                    dTdt   = temp.differentiate('Time')
                    print("\tcomputing atmospheric surface heat flux")
                    F_net  = self.compute_era5_net_atmospheric_surface_heat_flux(yr_str = yr_str, mo_str=mo_str)
                    D_qdp  = os.path.join(D_ocn,'qdp')
                    D_qdp_k = os.path.join(D_qdp,'md')
                    P_qdp  = os.path.join(D_qdp,'bran_qdp_{:s}_{:s}.nc'.format(yr_str,mo_str))
                    print("\tCreating file: {:s}".format(P_qdp))
                    if os.path.exists(P_qdp):
                        print("\tFile exists: {:s}\n ... advancing to next iteration".format(P_qdp))
                    else:
                        for k in range(len(temp.Time)):
                            P_qdp_k = os.path.join(D_qdp_k,'qdp_md{:02d}.nc'.format(k))
                            if os.path.exists(P_qdp_k):
                                print("\t\tFile exists: {:s}\n  ... advancing to next iteration".format(P_qdp_k))
                            else:
                                print("\t\tcurrent time: ",datetime.now().time())
                                mld_k  = mld.isel(Time=k).astype(np.single).load()
                                Fnet_k = F_net.isel(time=k).astype(np.single).load()
                                print("\t\tslicing temporal temperature derivative at the mixed layer depth")
                                dTdt_k = dTdt.isel(Time=j).sel(st_ocean=mld_k,method='nearest').drop('st_ocean').astype(np.single).load()
                                print("\t\tslicing temporal temperature at the mixed layer depth")
                                temp_k = temp.isel(Time=j).sel(st_ocean=mld_k,method='nearest').drop('st_ocean').astype(np.single).load()
                                print("\t\tslicing salinity temperature at the mixed layer depth")
                                salt_k = salt.isel(Time=j).sel(st_ocean=mld_k,method='nearest').drop('st_ocean').astype(np.single).load()
                                print("\t\tcomputing ocean heat capacity at the mixed layer depth")
                                cp_k = compute_ocn_heat_capacity_at_depth(salt_k, temp_k, mld_k, mld.ny).astype(np.single)
                                cp_k = cp_k.load()
                                print("\t\tcomputing ocean density at the mixed layer depth")
                                rho_k = compute_ocn_density_at_depth(salt_k, temp_k, mld_k, mld.ny).astype(np.single)
                                rho_k = rho_k.load()
                                print("\t\tcomputing ocean heat flux at the mixed layer depth")
                                qdp_k = compute_ocn_heat_flux_at_depth(rho_k, cp_k, mld_k, Fnet_k, dTdt_k)
                                qdp_k = qdp_k.load()
                                qdp_k = qdp_k.assign_coords(time=temp.Time.isel(Time=k).values).expand_dims('time').astype(np.single).to_dataset(name='qdp')
                                print("\t\twriting out ocean heat flux netcdf to: ",P_qdp_k)
                                logging.info("his what qdp_k looks like before writing: %s",qdp_k)
                                enc_dict  = {'shuffle':True,'zlib':True,'complevel':5}
                                write_qdp_k = qdp_k.to_netcdf(P_qdp_k, unlimited_dims=['time'], compute=False, encoding={'qdp':enc_dict})
                                with ProgressBar():
                                    print(f"\t\tWriting to {P_qdp_k}")
                                    write_qdp_k.compute()
                                print("\t\tcompleted writing qdp, current time: ",datetime.now().time())
                        #time.sleep(60)
                        print("\tconcatenating qdp for one month, files in {:s}".format(D_qdp_k))
                        print("\tcurrent time: ",datetime.now().time())
                        sys_call = 'ncrcat {:s} {:s}/qdp_md* {:s}'.format(self.ncrcat_opts_t_dim,D_qdp_k,P_qdp)
                        os.system(sys_call)
                        print("\tcompleted concatenating qdp for one month, current time: ",datetime.now().time())
                        if os.path.exists(P_qdp):
                            print("\tsuccessfully created {:s}\n\tdeleting temporary files in {:s}".format(P_qdp,D_qdp_k))
                            os.system("rm {:s}/*.nc".format(D_qdp_k))
                    print("\tloading in qdp one-month file {:s} ... and persisting in memory".format(P_qdp))
                    qdp  = xr.open_dataset(P_qdp).load()
                    print("\tslicing SST out of full ocean temperatures and keeping data in memory")
                    sst  = temp.sel(st_ocean=0,method='nearest').load()
                    print("\tslicing SSS out of full ocean salinity and keeping data in memory")
                    sss  = salt.sel(st_ocean=0,method='nearest').load()
                    uocn = self.bran_load_and_regrid( bran_var = 'u', yr_str = yr_str, mo_str = mo_str, grid_type = 'u')
                    vocn = self.bran_load_and_regrid( bran_var = 'v', yr_str = yr_str, mo_str = mo_str, grid_type = 'u')
                    print("\tslicing surface zonal velocities out of full ocean zonal velocities and keeping data in memory")
                    uocn = uocn.sel(st_ocean=0,method="nearest").load()
                    print("\tslicing surface meridional velocities out of full ocean meridional velocities and keeping data in memory")
                    vocn = vocn.sel(st_ocean=0,method="nearest").load()
                    eta  = self.bran_load_and_regrid( bran_var = 'eta_t', yr_str = yr_str, mo_str = mo_str, grid_type = 'u')
                    eta  = eta.load()
                    print("\tcomputing ocean surface slopes in the zonal direction and keeping data in memory")
                    dhdx = compute_ocn_sfc_slope(eta, G_CICE.hte.data, direction='x')
                    dhdx = dhdx.load()
                    print("\tcomputing ocean surface slopes in the meridional direction and keeping data in memory")
                    dhdy = compute_ocn_sfc_slope(eta, G_CICE.htn.data, direction='y')
                    dhdy = dhdy.load()
                    print("\tloading and keeping in memory, geographic ocean grid for writing into ocean forcing file {:s}".format(self.F_G_ACOM2))
                    #LN   = xr.open_dataset(self.F_G_ACOM2).geolon_t.load()
                    #LT   = xr.open_dataset(self.F_G_ACOM2).geolat_t.load()
                    print("\tcreating ocean forcing file data structure")
                    data_vars = {'qdp' : (['time','nj','ni'],qdp.qdp.data,{'units':'W/m^2','long_name':'deep ocean heat flux'}),
                                 'T'   : (['time','nj','ni'],sst.data,{'units':'C','long_name':'sea surface temperature'}),
                                 'S'   : (['time','nj','ni'],sss.data,{'units':'psu','long_name':'sea surface salinity'}),
                                 'hblt': (['time','nj','ni'],mld.data,{'units':'m','long_name':'mixed layer depth'}),
                                 'U'   : (['time','nj','ni'],uocn.data,{'units':'m/s','long_name':'meridional surface current'}),
                                 'V'   : (['time','nj','ni'],vocn.data,{'units':'m/s','long_name':'zonal surface current'}),
                                 'dhdx': (['time','nj','ni'],dhdx.data,{'units':'m','long_name':'zonal sea surface slope'}),
                                 'dhdy': (['time','nj','ni'],dhdy.data,{'units':'m','long_name':'meridional sea surface slope'})}
                    coords = {#'xc'   : (['nj','ni'],LN.data,{'units':'degrees_east'}),
                              #'yc'   : (['nj','ni'],LT.data,{'units':'degrees_north'}),
                              'time' : (['time'],sst.Time.data,{'long_name':'time','cartesian_axis':'T','calendar_type':'GREGORIAN'})}
                    attrs = {'creation_date': datetime.now().strftime('%Y-%m-%d %H'),
                             'conventions'  : 'CCSM data model domain description -- for CICE6 standalone "NCAR" ocean option',
                             'title'        : 'daily averaged ocean forcing from BRAN2020 output for CICE6 standalone ocean forcing',
                             'source'       : 'BRAN2020',
                             'comment'      : 'source files found on gadi, /g/data/gb6/BRAN/BRAN2020',
                             'calendar'     : 'standard',
                             'note1'        : 'grid is ACCESS-OM2 t-grid',
                             'note2'        : 'u,v, and eta (i.e. dhdx and dhdy) are interpolated to t-grid',
                             'author'       : 'Daniel P Atwater',
                             'email'        : 'daniel.atwater@utas.edu.au'}
                    print("\twriting one month ocean forcing file {:s}".format(P_ocn_j))
                    print("\tcurrent time: ",datetime.now().time())
                    OCN       = xr.Dataset(data_vars=data_vars,coords=coords,attrs=attrs)
                    enc_dict  = {'shuffle':True,'zlib':True,'complevel':5}
                    write_job = OCN.to_netcdf(P_ocn_j, unlimited_dims=['time'], compute=False,
                                            encoding={'qdp':enc_dict, 'T':enc_dict, 'S':enc_dict, 'hblt':enc_dict,
                                            'U':enc_dict, 'V':enc_dict, 'dhdx':enc_dict, 'dhdy':enc_dict})
                    with ProgressBar():
                        print(f'Writing to {P_ocn_j}')
                        write_job.compute()
                    print("\tcompleted writing one month ocean forcing file, current time: ",datetime.now().time())
            print("concatenating monthly ocean forcing files in {:s}".format(D_ocn_j))
            print("into a year-long file {:s}".format(P_ocn))
            print("\tcurrent time: ",datetime.now().time())
            sys_call = 'ncrcat {:s} {:s}/bran_ocn_frcg_cice6_{:s}_{:s}* {:s}'.format(self.ncrcat_opts_n_dim,D_ocn_j,self.G_res,yr_str,P_ocn)
            os.system(sys_call)
            print("\tcompleted writing year-long ocean forcing file, current time: ",datetime.now().time())

    ###########################################################################################
    def acom2_coallate_for_cice_forcing(self, start_date='',
                                        acom2_model_run='', cnk_dict='', var_freq='',
                                        regrid_method='', regrid_periodicity='', generate_weights='',
                                        weights_file_name=''):
        '''
        '''
        self.F_G_ACOM2 = '/g/data/ik11/grids/ocean_grid_01.nc'
        self.acom2_coallate_for_cice_forcing_has_been_called = True
        if not var_freq:           var_freq           = self.acom2_var_freq 
        if not acom2_model_run:    acom2_model_run    = self.acom2_model_run
        if not start_date:         start_date         = self.acom2_start_date
        if not cnk_dict:           cnk_dict           = self.acom2_chunking
        if not regrid_method:      regrid_method      = self.ERA5_regrid_method
        if not regrid_periodicity: regrid_periodicity = self.ERA5_periodicity
        if not generate_weights:   generate_weights   = self.ERA5_regen_weights
        if not weights_file_name: 
            F_weights = self.F_ERA5_weights
        else:
            F_weights = weight_file_name
        G_CICE  = self.cice_grid_prep()
        G_ERA5  = self.era5_grid_prep()
        for i in range(self.n_years):
            dt_obj = self.define_datetime_object(start_date=start_date, year_offset=i, full_datetime=True)
            yr_str = self.year_string_from_datetime_object(dt_obj)
            D_ocn  = os.path.join(self.D_modin,'ACOM2')
            P_ocn  = os.path.join(self.D_modin,'acom2_ocn_frcg_cice6_{res:s}_{yr:s}.nc'.format(res=self.G_res,yr=yr_str))
            D_ocn_j= os.path.join(D_ocn,'yr_mo')
            if os.path.exists(P_ocn):
                print('File exists: {:s}\n SKIPPING concatentation')
            else:
                for j in range(12):
                    if i>0:
                        dt0 = datetime.strptime(start_date,'%Y-%m-%d %H:%M') + timedelta(days=(i*365)) + relativedelta(months=j)
                        dtN = datetime.strptime(start_date,'%Y-%m-%d %H:%M') + timedelta(days=(i*364*2)) + relativedelta(months=j+1) - timedelta(days=1) 
                    else:
                        dt0 = datetime.strptime(start_date,'%Y-%m-%d %H:%M') + relativedelta(months=j)
                        dtN = datetime.strptime(start_date,'%Y-%m-%d %H:%M') + relativedelta(months=j+1) - timedelta(days=1)
                    dt0_str = dt0.strftime('%Y-%m-%d %H:%M')
                    dtN_str = dtN.strftime('%Y-%m-%d %H:%M')
                    mo_str = '{:02d}'.format(j+1)
                    P_ocn_j= os.path.join(D_ocn_j,'acom2_ocn_frcg_cice6_{res:s}_{yr:s}_{mo:s}.nc'.format(res=self.G_res,yr=yr_str,mo=mo_str))
                    print("\tCreating file: {:s}".format(P_ocn_j))
                    if os.path.exists(P_ocn_j):
                        print('\tFile exists: {:s}\n\t ... advancing to next iteration'.format(P_ocn_j))
                    else:
                        cc_sess = cc.database.create_session()
                        print("extracting ac-om2 data for QDP computation, from model run: ",acom2_model_run," from start dt: ",dt0_str," to stop dt: ",dtN_str)
                        temp = cc.querying.getvar(expt=acom2_model_run,variable='temp',session=cc_sess,frequency= var_freq,start_time=dt0_str,end_time=dtN_str).load()
                        salt = cc.querying.getvar(expt=acom2_model_run,variable='salt',session=cc_sess,frequency= var_freq,start_time=dt0_str,end_time=dtN_str).load()
                        mld  = cc.querying.getvar(expt=acom2_model_run,variable='mld',session=cc_sess,frequency= var_freq,start_time=dt0_str,end_time=dtN_str).load()
                        print("computing temperature derivative over time")
                        dTdt = temp.differentiate("time")
                        print("\tcomputing atmospheric surface heat flux")
                        lat_mom = temp.yt_ocean
                        lon_mom = temp.xt_ocean
                        G_tmp   = xr.Dataset({'lon':lon_mom,'lat':lat_mom})
                        F_net   = self.compute_era5_net_atmospheric_surface_heat_flux(G_dst              = G_tmp,
                                                                                      yr_str             = str(dt0.year),
                                                                                      mean_type          = 'daily',
                                                                                      regrid_method      = regrid_method,
                                                                                      regrid_periodicity = regrid_periodicity,
                                                                                      cnk_dict           = self.ERA5_chunking,
                                                                                      generate_weights   = generate_weights,
                                                                                      weights_file_name  = F_weights)
                        D_qdp   = os.path.join(D_ocn,'qdp')
                        D_qdp_k = os.path.join(D_qdp,'md')
                        P_qdp   = os.path.join(D_qdp,'acom2_qdp_{:s}_{:s}.nc'.format(yr_str,mo_str))
                        print("\tCreating file: {:s}".format(P_qdp))
                        if os.path.exists(P_qdp):
                            print("\tFile exists: {:s}\n ... advancing to next iteration".format(P_qdp))
                        else:
                            for k in range(len(temp.time)):
                                P_qdp_k = os.path.join(D_qdp_k,'qdp_md{:02d}.nc'.format(k))
                                if os.path.exists(P_qdp_k):
                                    print("\t\tFile exists: {:s}\n  ... advancing to next iteration".format(P_qdp_k))
                                else:
                                    print("\t\tcurrent time: ",datetime.now().time())
                                    mld_k  = mld.isel(time=k).astype(np.single).load()
                                    Fnet_k = F_net.isel(time=k).astype(np.single).load()
                                    print("\t\tslicing temporal temperature derivative at the mixed layer depth")
                                    dTdt_k = dTdt.isel(time=k).sel(st_ocean=mld_k,method='nearest').drop('st_ocean').astype(np.single).load()
                                    print("\t\tslicing temporal temperature at the mixed layer depth")
                                    temp_k = temp.isel(time=k).sel(st_ocean=mld_k,method='nearest').drop('st_ocean').astype(np.single).load()
                                    print("\t\tslicing salinity temperature at the mixed layer depth")
                                    salt_k = salt.isel(time=k).sel(st_ocean=mld_k,method='nearest').drop('st_ocean').astype(np.single).load()
                                    print("\t\tcomputing ocean heat capacity at the mixed layer depth")
                                    cp_k = compute_ocn_heat_capacity_at_depth(salt_k, temp_k, mld_k, mld.yt_ocean).astype(np.single)
                                    cp_k = cp_k.load()
                                    print("\t\tcomputing ocean density at the mixed layer depth")
                                    rho_k = compute_ocn_density_at_depth(salt_k, temp_k, mld_k, mld.yt_ocean).astype(np.single)
                                    rho_k = rho_k.load()
                                    print("\t\tcomputing ocean heat flux at the mixed layer depth")
                                    qdp_k = compute_ocn_heat_flux_at_depth(rho_k, cp_k, mld_k, Fnet_k, dTdt_k)
                                    qdp_k = qdp_k.load()
                                    qdp_k = qdp_k.assign_coords(time=temp.time.isel(time=k).values).expand_dims('time').astype(np.single).to_dataset(name='qdp')
                                    print("\t\twriting out ocean heat flux netcdf to: ",P_qdp_k)
                                    print("his what qdp_k looks like before writing: %s",qdp_k)
                                    enc_dict    = {'shuffle':True,'zlib':True,'complevel':5}
                                    write_qdp_k = qdp_k.to_netcdf(P_qdp_k, unlimited_dims=['time'], compute=False, encoding={'qdp':enc_dict})
                                    with ProgressBar():
                                        print(f"\t\tWriting to {P_qdp_k}")
                                        write_qdp_k.compute()
                                        print("\t\tcompleted writing qdp, current time: ",datetime.now().time())
                            print("\tconcatenating qdp for one month, files in {:s}".format(D_qdp_k))
                            print("\tcurrent time: ",datetime.now().time())
                            sys_call = 'ncrcat {:s} {:s}/qdp_md* {:s}'.format(self.ncrcat_opts_t_dim,D_qdp_k,P_qdp)
                            os.system(sys_call)
                            print("\tcompleted concatenating qdp for one month, current time: ",datetime.now().time())
                            if os.path.exists(P_qdp):
                                print("\tsuccessfully created {:s}\n\tdeleting temporary files in {:s}".format(P_qdp,D_qdp_k))
                                os.system("rm {:s}/*.nc".format(D_qdp_k))
                        print("\tloading in qdp one-month file {:s} ... and persisting in memory".format(P_qdp))
                        qdp  = xr.open_dataset(P_qdp).load()
                        print("extracting ac-om2 data from model run: ",acom2_model_run," from start dt: ",dt0_str," to stop dt: ",dtN_str)
                        sst  = temp.isel(st_ocean=0).drop('st_ocean').load()
                        sss  = salt.isel(st_ocean=0).drop('st_ocean').load()                    
                        uocn = cc.querying.getvar(expt=acom2_model_run,variable='u',session=cc_sess,frequency=var_freq,start_time=dt0_str,end_time=dtN_str)
                        uocn = uocn.isel(st_ocean=0).drop('st_ocean').load()
                        vocn = cc.querying.getvar(expt=acom2_model_run,variable='v',session=cc_sess,frequency=var_freq,start_time=dt0_str,end_time=dtN_str)
                        vocn = vocn.isel(st_ocean=0).drop('st_ocean').load()
                        mld  = mld.sel(time=slice(dt0_str,dtN_str)).load()
                        eta  = cc.querying.getvar(expt=acom2_model_run,variable='sea_level',session=cc_sess,frequency=var_freq,start_time=dt0_str,end_time=dtN_str)
                        eta  = eta.sel(time=slice(dt0_str,dtN_str)).load()
                        print("\tcomputing ocean surface slopes in the zonal direction and keeping data in memory")
                        dhdx = compute_ocn_sfc_slope(eta, G_CICE.hte.data, direction='x').astype(np.single)
                        dhdx = dhdx.load()
                        print("\tcomputing ocean surface slopes in the meridional direction and keeping data in memory")
                        dhdy = compute_ocn_sfc_slope(eta, G_CICE.htn.data, direction='y').astype(np.single)
                        dhdy = dhdy.load()
                        #print("\tloading and keeping in memory, geographic ocean grid for writing into ocean forcing file {:s}".format(self.F_G_ACOM2))
                        #LN   = xr.open_dataset(self.F_G_ACOM2).geolon_t.load()
                        #LT   = xr.open_dataset(self.F_G_ACOM2).geolat_t.load()
                        print("\tcreating ocean forcing file data structure")
                        data_vars = {'qdp' : (['time','nj','ni'],qdp.qdp.data,{'units':'W/m^2','long_name':'deep ocean heat flux','_FillValue':9.96921e+36}),
                                     'T'   : (['time','nj','ni'],sst.data-273.15,{'units':'degC','long_name':'sea surface temperature','_FillValue':9.96921e+36}),
                                     'S'   : (['time','nj','ni'],sss.data,{'units':'ppt','long_name':'sea surface salinity','_FillValue':9.96921e+36}),
                                     'hblt': (['time','nj','ni'],mld.data,{'units':'m','long_name':'mixed layer depth','_FillValue':9.96921e+36}),
                                     'U'   : (['time','nj','ni'],uocn.data,{'units':'m/s','long_name':'meridional surface current','_FillValue':9.96921e+36}),
                                     'V'   : (['time','nj','ni'],vocn.data,{'units':'m/s','long_name':'zonal surface current','_FillValue':9.96921e+36}),
                                     'dhdx': (['time','nj','ni'],dhdx.data,{'units':'m','long_name':'zonal sea surface slope','_FillValue':9.96921e+36}),
                                     'dhdy': (['time','nj','ni'],dhdy.data,{'units':'m','long_name':'meridional sea surface slope','_FillValue':9.96921e+36})}
                        coords = {#'xc'   : (['ni'],LN.data[0,:],{'units':'degrees_east'}),
                                  #'yc'   : (['nj'],LT.data[:,0],{'units':'degrees_north'}),
                                  'time' : (['time'],sst.time.data,{'long_name':'observation time','cartesian_axis':'T','calendar_type':'GREGORIAN'})}
                        attrs = {'creation_date': datetime.now().strftime('%Y-%m-%d %H'),
                                 'conventions'  : 'CCSM data model domain description -- for CICE6 standalone "NCAR" ocean option',
                                 'title'        : 'daily averaged ocean forcing from ACCESS-OM2-01 outputn for CICE6 standalone ocean forcing',
                                 'source'       : 'ACCESS-OM2-01',
                                 'comment'      : 'source files found using COSIMA Cookbook',
                                 'calendar'     : 'standard',
                                 'note1'        : 'grid is ACCESS-OM2 t-grid',
                                 'note2'        : 'u,v, and eta (i.e. dhdx and dhdy) are interpolated to t-grid',
                                 'note3'        : 'dhdx,dhdy are described in ACCESS-OM2 output as effective sea level (eta_t + patm/(rho0*g)) on T cells',
                                 'author'       : 'Daniel P Atwater',
                                 'email'        : 'daniel.atwater@utas.edu.au'}
                        print("\twriting one month ocean forcing file {:s}".format(P_ocn_j))
                        print("\tcurrent time: ",datetime.now().time())
                        OCN       = xr.Dataset(data_vars=data_vars,coords=coords,attrs=attrs)
                        enc_dict  = {'shuffle':True,'zlib':True,'complevel':5}
                        write_job = OCN.to_netcdf(P_ocn_j, unlimited_dims=['time'], compute=False,
                                                encoding={'qdp':enc_dict, 'T':enc_dict, 'S':enc_dict, 'hblt':enc_dict,
                                                'U':enc_dict, 'V':enc_dict, 'dhdx':enc_dict, 'dhdy':enc_dict})
                        with ProgressBar():
                            print(f'Writing to {P_ocn_j}')
                            write_job.compute()
                        print("\tcompleted writing one month ocean forcing file, current time: ",datetime.now().time())
            print("concatenating monthly ocean forcing files in {:s}".format(D_ocn_j))
            print("into a year-long file {:s}".format(P_ocn))
            print("\tcurrent time: ",datetime.now().time())
            sys_call = 'ncrcat {:s} {:s}/acom2_ocn_frcg_cice6_{:s}_{:s}* {:s}'.format(self.ncrcat_opts_n_dim,D_ocn_j,self.G_res,yr_str,P_ocn)
            os.system(sys_call)
            print("\tcompleted writing year-long ocean forcing file, current time: ",datetime.now().time())

    ############################################################################################
    def compute_era5_net_atmospheric_surface_heat_flux(self, G_dst='', yr_str='', mo_str='', mean_type='daily',
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
        logging.info("generating regridding weights on the fly using xesmf via method: %s",regrid_method)
        lw    = self.era5_load_and_regrid(era5_var           = 'msdwlwrf',
                                          G_dst              = G_dst,
                                          yr_str             = yr_str,
                                          mo_str             = mo_str,
                                          regrid_method      = regrid_method,
                                          regrid_periodicity = regrid_periodicity,
                                          cnk_dict           = cnk_dict,
                                          generate_weights   = generate_weights,
                                          weights_file_name  = F_weights)
        sw    = self.era5_load_and_regrid(era5_var           = 'msdwswrf',
                                          G_dst              = G_dst,
                                          yr_str             = yr_str,
                                          mo_str             = mo_str,
                                          regrid_method      = regrid_method,
                                          regrid_periodicity = regrid_periodicity,
                                          cnk_dict           = cnk_dict,
                                          generate_weights   = generate_weights,
                                          weights_file_name  = F_weights)
        sh    = self.era5_load_and_regrid(era5_var           = 'msshf',
                                          G_dst              = G_dst,
                                          yr_str             = yr_str,
                                          mo_str             = mo_str,
                                          regrid_method      = regrid_method,
                                          regrid_periodicity = regrid_periodicity,
                                          cnk_dict           = cnk_dict,
                                          generate_weights   = generate_weights,
                                          weights_file_name  = F_weights)
        lh    = self.era5_load_and_regrid(era5_var           = 'mslhf',
                                          G_dst              = G_dst,
                                          yr_str             = yr_str,
                                          mo_str             = mo_str,
                                          regrid_method      = regrid_method,
                                          regrid_periodicity = regrid_periodicity,
                                          cnk_dict           = cnk_dict,
                                          generate_weights   = generate_weights,
                                          weights_file_name  = F_weights)
        F_net = sw.msdwswrf - lw.msdwlwrf - lh.mslhf - sh.msshf
        if mean_type=='daily':
            print("computing daily mean of atmospheric heat flux")
            F_net = F_net.resample(time='1D').mean()
        elif mean_type=='monthly':
            print("computing monthly mean of atmospheric heat flux")
            F_net = F_net.resample(time='1M').mean()
        return F_net

