import os, json, shutil
import xarray      as xr
import numpy       as np
import pandas      as pd
import geopandas   as gpd
from scipy.spatial import cKDTree
from pathlib       import Path
from datetime      import datetime

class GroundedIcebergProcessor:
    """
    Class to preprocess grounded iceberg data for use in sea ice model analysis.

    This class loads grounded iceberg shapefiles, converts them to the model grid,
    applies probabilistic thinning (if requested), and generates a modified landmask
    (KMT array) for grounded iceberg-affected grid cells. These modifications are
    used when analyzing fast ice and pack ice fields in CICE model output.

    Parameters
    ----------
    config     : dict
                 Configuration dictionary (typically loaded from afim_cice_analysis.json).
    sim_name   : str
                 Simulation name (used for locating paths and logs).
    P_gi_shp   : str or Path, optional
                 Path to grounded iceberg shapefile (overrides config).
    P_KMT      : str or Path, optional
                 Path to KMT landmask NetCDF file (overrides config).
    thin_prob  : float, optional
                 Probability threshold for thinning grounded iceberg cells (0–1 range).
                 If None, no thinning is applied.

    Attributes
    ----------
    G_t           : xarray.Dataset
                    Model T-grid loaded from NetCDF (with lon/lat/area).
    G_u           : xarray.Dataset
                    Model U-grid loaded from NetCDF (with lon/lat).
    GI_mask       : xarray.DataArray
                    Binary mask (1 where grounded iceberg exists).
    GI_lons       : np.ndarray
                    Longitude of grounded iceberg centroids (for plotting).
    GI_lats       : np.ndarray
                    Latitude of grounded iceberg centroids (for plotting).
    GI_lon_cells  : np.ndarray
                    Longitude of affected model grid cells (for KMT plotting).
    GI_lat_cells  : np.ndarray
                    Latitude of affected model grid cells (for KMT plotting).
    KMT_path      : str
                    Path to landmask (KMT) file.
    total_area    : float
                    Total area (in m²) affected by grounded icebergs.
    use_gi        : bool
                    True if grounded iceberg masking is active.
    logger        : logging.Logger
                    Logger instance for status and debugging messages.

    Examples
    --------
    Load grounded iceberg and landmask information:

    >>> from grounded_iceberg_processor import GroundedIcebergProcessor
    >>> proc = GroundedIcebergProcessor(config, sim_name="Rothrock")
    >>> proc.load_GI()
    >>> gi_mask = proc.GI_mask

    Plot modified landmask with PyGMT:

    >>> fig.plot(x=proc.GI_lon_cells, y=proc.GI_lat_cells, style="c0.1c", fill="purple")
    >>> fig.plot(x=proc.GI_lons, y=proc.GI_lats, style="s0.1c", fill="black")

    References
    ----------
    Full AFIM repository: https://github.com/dpath2o/AFIM
    JSON config:          AFIM/JSONs/afim_cice_analysis.json
    Shapefile source:     SCAR Antarctic Digital Database (ADD) or similar
    """

    def __init__(self, P_json=None, sim_name=None, overwrite=False, GI_raw_version=None):
        """
        Initialize the GroundedIcebergProcessor with configuration data and simulation name.

        Parameters:
        -----------
        config : dict, optional
            Dictionary loaded from a JSON config file.
        sim_name : str, optional
            The name of the simulation to process.
        overwrite : bool, optional
            If True, allows overwriting of existing NetCDF outputs. Defaults to False.
        """
        if P_json is None:
            P_json = '/home/581/da1339/AFIM/src/AFIM/src/JSONs/afim_cice_analysis.json'
        with open(P_json, 'r') as f:
            self.config = json.load(f)
        self.overwrite          = overwrite
        self.ice_shelves_loaded = False
        self.CICE_dict          = self.config['CICE_dict']
        self.CICE_P_G           = self.CICE_dict['P_G']
        self.spatial_dims       = self.CICE_dict.get('spatial_dims', ('nj', 'ni'))
        self.GI_dict            = self.config['GI_dict']
        self.P_KMT_org          = Path(self.GI_dict["D_GI_thin"],self.GI_dict['KMT_org_fmt'])
        self.FIC_scale          = self.config.get("FIC_scale", 1e9)
        if sim_name is not None:
            self.sim_name       = sim_name
            self.sim_config     = self.config['sim_dict'][sim_name]
            self.GI_thin_fact   = self.sim_config.get('GI_thin_fact')
            self.GI_version     = self.sim_config.get('GI_version', 1.5)
            self.GI_version_str = f"{self.GI_version:.2f}".replace('.', 'p')
        else:
            print("Simulation name not provided. Class assumes you are creating new landmask for sea ice modelling.")
            self.GI_version     = GI_raw_version if GI_raw_version is not None else 1.5
            self.GI_version_str = f"{self.GI_version:.2f}".replace('.', 'p')
            self.GI_thin_fact   = None
        if self.GI_thin_fact is not None:
            GI_thin_str      = f"{self.GI_thin_fact:0.2f}".replace('.', 'p')
            self.GI_P_counts = os.path.join(self.GI_dict['D_GI_thin'],
                                            self.GI_dict['GI_thin_fmt'].format(GI_thin_fact=GI_thin_str, version=self.GI_version_str))
            self.P_KMT_mod   = os.path.join(self.GI_dict['D_GI_thin'],
                                            self.GI_dict['KMT_mod_fmt'].format(GI_thin_fact=GI_thin_str, version=self.GI_version_str))
            self.use_gi      = True
        else:
            self.GI_P_counts = ''
            self.P_KMT_mod   = self.P_KMT_org 
            self.use_gi      = False
            self.total_area  = 0
        self.bgrid_loaded              = False
        self.modified_landmask_aligned = False

    def radians_to_degrees(self, da):
        return (da * 180) / np.pi

    def normalise_longitudes(self,lon):
        return (lon + 360) % 360

    def build_grid_dict(self, lat, lon):
        ny, nx               = lat.shape
        lat_b                = np.zeros((ny + 1, nx + 1))
        lon_b                = np.zeros((ny + 1, nx + 1))
        lat_b[1:-1, 1:-1]    = 0.25 * (lat[:-1, :-1] + lat[:-1, 1:] + lat[1:, :-1] + lat[1:, 1:])
        lon_b[1:-1, 1:-1]    = 0.25 * (lon[:-1, :-1] + lon[:-1, 1:] + lon[1:, :-1] + lon[1:, 1:])
        lat_b[[0, -1], 1:-1] = lat[[0, -1], :-1]
        lon_b[[0, -1], 1:-1] = lon[[0, -1], :-1]
        lat_b[1:-1, [0, -1]] = lat[:-1, [0, -1]]
        lon_b[1:-1, [0, -1]] = lon[:-1, [0, -1]]
        lat_b[0, 0]          = lat[0, 0]
        lat_b[0, -1]         = lat[0, -1]
        lat_b[-1, 0]         = lat[-1, 0]
        lat_b[-1, -1]        = lat[-1, -1]
        lon_b[0, 0]          = lon[0, 0]
        lon_b[0, -1]         = lon[0, -1]
        lon_b[-1, 0]         = lon[-1, 0]
        lon_b[-1, -1]        = lon[-1, -1]        
        return {"lat"  : lat,
                "lon"  : self.normalise_longitudes(lon),
                "lat_b": lat_b,
                "lon_b": self.normalise_longitudes(lon_b)}            

    def load_bgrid(self):
        """
        Load the B-grid (t-grid and u-grid) and create cooordinates which are converted to degrees and standardized to [-180, 180].
        Sets attributes `self.G_t`, `self.G_u`.
        """
        G        = xr.open_dataset(self.CICE_P_G)
        KMT_org  = xr.open_dataset(self.P_KMT_org).kmt.data
        if self.use_gi:
            KMT_mod  = xr.open_dataset(self.P_KMT_mod).kmt.data
        else:
            KMT_mod = KMT_org
        TLAT     = self.radians_to_degrees(G['tlat'].data)
        TLON     = self.radians_to_degrees(G['tlon'].data)
        ULAT     = self.radians_to_degrees(G['ulat'].data)
        ULON     = self.radians_to_degrees(G['ulon'].data)
        TCOORDS  = self.build_grid_dict( TLAT , TLON )
        UCOORDS  = self.build_grid_dict( ULAT , ULON )
        T_ANGLE  = self.radians_to_degrees(G['angleT'].data)
        U_ANGLE  = self.radians_to_degrees(G['angle'].data)
        TAREA    = G['tarea'].data
        UAREA    = G['uarea'].data
        # Dimensions
        nj, ni           = TLAT.shape
        nj_b, ni_b       = nj + 1, ni + 1
        native_dim_names = ('nj','ni')
        native_dims      = (nj, ni)
        extend_dim_names = ('nj_b','ni_b')
        extend_dims      = (nj_b, ni_b)
        self.G_t = xr.Dataset(data_vars = { 'lat'     : (native_dim_names, TCOORDS['lat']   , {'units': 'degrees'}),
                                            'lat_b'   : (extend_dim_names, TCOORDS['lat_b'] , {'units': 'degrees'}),
                                            'lon'     : (native_dim_names, TCOORDS['lon']   , {'units': 'degrees'}),
                                            'lon_b'   : (extend_dim_names, TCOORDS['lon_b'] , {'units': 'degrees'}),
                                            'angle'   : (native_dim_names, T_ANGLE          , {'units': 'degrees'}),
                                            'area'    : (native_dim_names, TAREA            , {'units': 'm^2'}),
                                            'kmt_org' : (native_dim_names, KMT_org          , {'units'      : 'binary',
                                                                                                'description': '1=land, 0=ocean',
                                                                                                'long_name'  : 'original landmask on t-grid'}),
                                            'kmt_mod' : (native_dim_names, KMT_mod          , {'units'      : 'binary',
                                                                                                'description': '1=land, 0=ocean',
                                                                                                'long_name'  : 'modified t-grid-landmask to simulate grounded icebergs'})},
                               coords   = { 'nj'   : np.arange(nj),
                                            'ni'   : np.arange(ni),
                                            'nj_b' : np.arange(nj_b),
                                            'ni_b' : np.arange(ni_b)})
        self.G_u = xr.Dataset(data_vars = { 'lat'     : (native_dim_names, UCOORDS['lat']   , {'units': 'degrees'}),
                                            'lat_b'   : (extend_dim_names, UCOORDS['lat_b'] , {'units': 'degrees'}),
                                            'lon'     : (native_dim_names, UCOORDS['lon']   , {'units': 'degrees'}),
                                            'lon_b'   : (extend_dim_names, UCOORDS['lon_b'] , {'units': 'degrees'}),
                                            'angle'   : (native_dim_names, U_ANGLE          , {'units': 'degrees'}),
                                            'area'    : (native_dim_names, UAREA            , {'units': 'm^2'}),
                                            'kmt_org' : (native_dim_names, KMT_org          , {'units'      : 'binary',
                                                                                                'description': '1=land, 0=ocean',
                                                                                                'long_name'  : 'original landmask on t-grid'}),
                                            'kmt_mod' : (native_dim_names, KMT_mod          , {'units'      : 'binary',
                                                                                                'description': '1=land, 0=ocean',
                                                                                                'long_name'  : 'modified t-grid-landmask to simulate grounded icebergs'})},
                              coords    = { 'nj'   : np.arange(nj),
                                            'ni'   : np.arange(ni),
                                            'nj_b' : np.arange(nj_b),
                                            'ni_b' : np.arange(ni_b)})
        self.bgrid_loaded = True

    def align_modified_landmask(self):
        """
        Compute the grounded iceberg mask and coordinates from the difference
        between original and modified KMT arrays. Store in self.G_t.
        """
        if not self.bgrid_loaded:
            self.load_bgrid()
        kmt_mod = self.G_t['kmt_mod'].data
        kmt_org = self.G_t['kmt_org'].data
        lat     = self.G_t['lat'].data
        lon     = self.G_t['lon'].data
        area    = self.G_t['area'].data
        # Difference: grounded icebergs are cells that changed from ocean (1) to land (0)
        diff_mask           = (kmt_org == 1) & (kmt_mod == 0)
        self.G_t['GI_mask'] = (('nj', 'ni'), diff_mask)
        # Get coordinates of affected cells (shifted west by one ni index to match B-grid layout)
        nj_idx, ni_idx = np.where(diff_mask)
        ni_idx_shifted = ni_idx - 1
        valid          = ni_idx_shifted >= 0
        nj_idx         = nj_idx[valid]
        ni_idx         = ni_idx_shifted[valid]
        GI_lat         = lat[nj_idx, ni_idx]
        GI_lon          = lon[nj_idx, ni_idx]
        # Save 1D arrays with iceberg IDs
        iceberg_id         = np.arange(len(GI_lat))
        self.G_t['GI_lat'] = (('iceberg_id',), GI_lat)
        self.G_t['GI_lon'] = (('iceberg_id',), GI_lon)
        self.G_t['GI_lat'].attrs.update({'long_name': 'Grounded iceberg latitude'})
        self.G_t['GI_lon'].attrs.update({'long_name': 'Grounded iceberg longitude'})
        print(f"{len(iceberg_id)} circumpolar grounded icebergs associated with {self.sim_name}")
        self.modified_landmask_aligned = True

    def compute_grounded_iceberg_area(self, region=None):
        '''
        Compute grounded iceberg area from modified KMT relative to original KMT.

        Parameters
        ----------
        region : dict, optional
            Dictionary with keys as region names and values as sub-dictionaries:
            {'REG_NAME': {'ext': [lon_min, lon_max, lat_min, lat_max]}}

        Returns
        -------
        float or dict
            Total grounded iceberg area in m^2 if region is None.
            Dictionary of region-specific grounded iceberg areas otherwise.
        '''
        if not self.modified_landmask_aligned:
            self.align_modified_landmask()
        area = self.G_t['area'].values
        mask = self.G_t['GI_mask']
        lon  = self.G_t['lon'].values
        lat  = self.G_t['lat'].values
        total_area = 0
        if region is not None:
            area_dict = {}
            for reg_name, reg_cfg in region.items():
                lon_min, lon_max, lat_min, lat_max = reg_cfg['ext']
                region_mask                        = ( (lon >= lon_min) &
                                                       (lon <= lon_max) &
                                                       (lat >= lat_min) &
                                                       (lat <= lat_max) )
                total                              = np.sum(area[mask & region_mask])
                area_dict[reg_name]                = total
            self.G_t['GI_regional_area']           = area_dict
            return area_dict
        else:
            total_area                = np.sum(area[mask])
            self.G_t['GI_total_area'] = total_area
            print(f"{total_area:0.2f} m^2 total circumpolar grounded iceberg area for {self.sim_name}")
            return total_area

    def load_existing_thinned(self):
        """
        Load previously saved thinned grounded iceberg data and modified KMT file.

        Sets the attributes `self.GI_thin_da`, `self.GI_thin_mask`, and `self.GI_KMT`
        if their respective NetCDF files exist.
        """
        if os.path.exists(self.GI_P_counts):
            self.GI_thin_da = xr.open_dataarray(self.GI_P_counts)
            self.GI_thin_mask = self.GI_thin_da.values > 0
        if os.path.exists(self.P_KMT_mod):
            self.GI_KMT = xr.open_dataarray(self.P_KMT_mod)

    def process_raw_grounded_iceberg_database(self, GI_P_raw=None):
        """
        Load and spatially map raw grounded iceberg locations from CSV to the model grid.

        Uses a KD-tree nearest-neighbor search to assign iceberg coordinates to grid cells.
        Sets `self.GI_clean` with attached (nj, ni) indices.
        """
        self.load_bgrid()
        GI_P_raw                                        = GI_P_raw if GI_P_raw is not None else self.GI_dict['P_raw']
        GI_raw                                          = pd.read_csv(GI_P_raw)
        GI_clean                                        = GI_raw.dropna(subset=['X', 'Y']).copy()
        G_pts                                           = np.column_stack((self.G_t['lon'].values.ravel(), self.G_t['lat'].values.ravel()))
        tree                                            = cKDTree(G_pts)
        GI_locs                                         = np.column_stack((GI_clean['X'].values, GI_clean['Y'].values))
        _, nrst_i                                       = tree.query(GI_locs)
        nj, ni                                          = self.G_t['lon'].shape
        nj_i, ni_i                                      = np.unravel_index(nrst_i, (nj, ni))
        GI_clean['nj_i']                                = nj_i
        GI_clean['ni_i']                                = ni_i
        self.GI_clean                                   = GI_clean
        GI_cnts                                         = np.zeros((nj, ni), dtype=int)
        GI_cnts_df                                      = self.GI_clean.groupby(['nj_i', 'ni_i']).size().reset_index(name='count')
        GI_cnts[GI_cnts_df['nj_i'], GI_cnts_df['ni_i']] = GI_cnts_df['count']
        GI_cnts_da                                      = xr.DataArray(GI_cnts, dims=self.spatial_dims, name='GI_counts')
        self.GI_cnts                                    = GI_cnts_da

    def modify_landmask_with_grounded_icebergs(self, p_min=0.1, p_max=0.9, scaling_exponent=1.5):
        """
        Apply probabilistic thinning to the grounded iceberg mask based on iceberg counts.

        Parameters:
        -----------
        p_min : float
            Minimum probability of retaining a cell with low counts.
        p_max : float
            Maximum probability of retaining a cell with high counts.
        scaling_exponent : float
            Controls steepness of probability scaling with count.

        Sets:
        ------
        self.GI_thin_da : xarray.DataArray
        self.GI_thin_mask : ndarray (boolean)
        self.GI_thin_fact : float (fraction thinned out)
        """
        kmt                   = self.KMT_mod.values.copy()
        mask                  = self.GI_cnts.values >= 1
        counts                = self.GI_cnts.values
        thin_mask, thin_fact  = self._thin_mask(mask, counts, p_min, p_max, scaling_exponent)
        GI_thinned            = counts * thin_mask
        self.GI_thin_da       = xr.DataArray(data   = GI_thinned,
                                            dims   = self.spatial_dims,
                                            coords = {'lat': (self.spatial_dims, self.G_t['lat'].values),
                                                      'lon': (self.spatial_dims, self.G_t['lon'].values)},
                                            name   = 'GI_counts')
        self.GI_thin_mask     = thin_mask
        self.GI_thin_fact     = thin_fact
        thin_nj, thin_ni      = np.where(thin_mask)
        kmt[thin_nj, thin_ni] = 0
        self.GI_KMT           = xr.DataArray(data   = kmt,
                                             dims   = self.spatial_dims,
                                             coords = {'lat': (self.spatial_dims, self.G_t['lat'].values),
                                                       'lon': (self.spatial_dims, self.G_t['lon'].values)},
                                             name   = 'kmt' )

    def write_modified_landmask_and_counts_to_disk(self, write=True, overwrite=None):
        """
        Save the thinned grounded iceberg counts and modified KMT landmask to NetCDF files.

        File paths are generated from the config template using GI thinning factor and version.
        """
        overwrite        = overwrite if overwrite is not None else self.overwrite
        GI_thin_fact_str = f"{self.GI_thin_fact:.2f}".replace('.', 'p')
        F_thin           = self.GI_dict['GI_thin_fmt'].format(GI_thin_fact=GI_thin_fact_str, version=self.GI_version_str)
        F_KMT_mod        = self.GI_dict['KMT_mod_fmt'].format(GI_thin_fact=GI_thin_fact_str, version=self.GI_version_str)
        self.GI_P_counts = os.path.join(self.GI_dict['D_GI_thin'], F_thin)
        self.P_KMT_mod      = os.path.join(self.GI_dict['D_GI_thin'], F_KMT_mod)
        print(f"modified KMT     : {self.GI_P_counts}")
        print(f"thinned GI count : {self.P_KMT_mod}")
        if write:
            if os.path.exists(self.GI_P_counts) or os.path.exists(self.P_KMT_mod):
                if not overwrite:
                    print(f"GI counts file and/or KMT file already exists:\n  - {self.GI_P_counts}\n  - {self.P_KMT_mod}\n"
                          "These files may have been used in model simulations. Set `overwrite=True` in the constructor to allow overwriting." )
                else:
                    date_tag   = datetime.now().strftime("%Y-%m-%d")
                    base_gi    = os.path.splitext(self.GI_P_counts)[0]
                    base_kmt   = os.path.splitext(self.P_KMT_mod)[0]
                    backup_gi  = f"{base_gi}_{date_tag}.nc"
                    backup_kmt = f"{base_kmt}_{date_tag}.nc"
                    shutil.copy2(self.GI_P_counts, backup_gi)
                    shutil.copy2(self.P_KMT_mod, backup_kmt)
                    print("*** WARNING OVER-WRITING EXISTING FILES ***")
                    print(f"Backups created:\n - {backup_gi}\n - {backup_kmt}")
                    print("Saving thinned GI count and modified KMT files")
                    self.GI_thin_da.to_netcdf(self.GI_P_counts, mode='w')
                    self.GI_KMT.to_netcdf(self.P_KMT_mod, mode='w')
            else:
                print("Modified KMT and GI-count-file do not exist. Writing out the above files")
                self.GI_thin_da.to_netcdf(self.GI_P_counts, mode='w')
                self.GI_KMT.to_netcdf(self.P_KMT_mod, mode='w')

    def _is_isolated(self, y, x, mask):
        """
        Determine if a given grid cell is isolated (no grounded icebergs in 8-connected neighbors).

        Parameters:
        -----------
        y, x : int
            Grid coordinates to test.
        mask : ndarray
            Boolean mask where True = grounded iceberg.

        Returns:
        --------
        bool
            True if isolated, False otherwise.
        """
        ny, nx = mask.shape
        adjacent_coords = [(y-1, x-1), (y-1, x)  , (y-1, x+1),
                           (y,   x-1), (y,   x+1),
                           (y+1, x-1), (y+1, x)  , (y+1, x+1) ]
        for adj_y, adj_x in adjacent_coords:
            if 0 <= adj_y < ny and 0 <= adj_x < nx:
                if mask[adj_y, adj_x]:
                    return False
        return True

    def _thin_mask(self, mask, counts, p_min=0.2, p_max=0.8, scaling_exponent=2.0):
        """
        Perform probabilistic thinning of grounded iceberg cells based on their counts.

        Parameters:
        -----------
        mask : ndarray
            Boolean array of original grounded iceberg presence.
        counts : ndarray
            Corresponding iceberg counts.
        p_min : float
            Minimum probability of retaining a cell.
        p_max : float
            Maximum probability of retaining a cell.
        scaling_exponent : float
            Controls how sharply retention probability scales with count.

        Returns:
        --------
        thinned_mask : ndarray (bool)
            Mask after thinning.
        GI_thinning_factor : float
            Fraction of cells that were removed.
        """
        ny, nx               = mask.shape
        thinned_mask         = np.zeros_like(mask, dtype=bool)
        count_min, count_max = np.nanmin(counts), np.nanmax(counts)
        if count_max > count_min:
            norm_counts = (counts - count_min) / (count_max - count_min)
            norm_counts = norm_counts ** scaling_exponent
        else:
            norm_counts = np.ones_like(counts)
        retain_probs   = p_min + norm_counts * (p_max - p_min)
        retained_count = 0
        total_count    = np.count_nonzero(mask)
        for y in range(ny):
            for x in range(nx):
                if mask[y, x]:
                    retain_prob = retain_probs[y, x]
                    if self._is_isolated(y, x, mask) or np.random.rand() < retain_prob:
                        thinned_mask[y, x] =  True
                        retained_count     += 1
        GI_thinning_factor = 1 - (retained_count / total_count) if total_count > 0 else 0
        return thinned_mask, GI_thinning_factor
