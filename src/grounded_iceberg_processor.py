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

    def load_AFIM_GI(self):
        """
        Load pre-generated grounded iceberg mask and modified KMT information.

        This method does not create new files or apply thinning. It reads existing
        outputs from NetCDF files based on the configured simulation dictionary.
        """
        if not self.use_gi:
            return None
        self.load_grid_and_landmask()
        GI_locs             = xr.open_dataset(self.GI_P_counts)
        G_pts               = np.column_stack((self.G_t["lon"].values.ravel(), self.G_t["lat"].values.ravel()))
        GI_pts              = np.column_stack((GI_locs["lon"].values.ravel(), GI_locs["lat"].values.ravel()))
        tree                = cKDTree(G_pts)
        _, nrst_i           = tree.query(GI_pts)
        nj, ni              = self.G_t['lon'].shape
        nj_i, ni_i          = np.unravel_index(nrst_i, (nj, ni))
        GI_mask             = np.zeros((nj, ni), dtype=bool)
        GI_mask[nj_i, ni_i] = True
        GI_total_area       = np.sum(self.G_t["area"].values[GI_mask]) / self.FIC_scale
        KMT_org             = xr.open_dataset(self.P_KMT_org)
        KMT_mod             = xr.open_dataset(self.P_KMT_mod)
        kmt_orig_arr        = KMT_org['kmt'].values
        kmt_mod_arr         = KMT_mod['kmt'].values
        lat_arr             = KMT_mod['lat'].values
        lon_arr             = KMT_mod['lon'].values
        lon_arr_wrapped     = (lon_arr + 360) % 360
        diff_mask           = (kmt_orig_arr == 1) & (kmt_mod_arr == 0)
        nj, ni              = np.where(diff_mask)
        lon_grid_cells      = lon_arr_wrapped[nj, ni]
        lat_grid_cells      = lat_arr[nj, ni]
        self.GI_mask        = GI_mask
        self.total_area     = GI_total_area
        self.GI_lons        = GI_locs["lon"].values.ravel()
        self.GI_lats        = GI_locs["lat"].values.ravel()
        self.GI_lon_cells   = lon_grid_cells
        self.GI_lat_cells   = lat_grid_cells
        self.nj_i           = nj_i
        self.ni_i           = ni_i

    def load_grid_and_landmask(self):
        """
        Load the CICE grid (t-grid and u-grid) and the landmask (KMT).

        Coordinates are converted to degrees and standardized to [-180, 180].
        Sets attributes `self.G_t`, `self.G_u`, and `self.KMT`.
        """
        grid_ds = xr.open_dataset(self.CICE_P_G)
        G_t     = xr.Dataset({'lat'  : (self.spatial_dims, grid_ds.tlat.data * 180 / np.pi),
                              'lon'  : (self.spatial_dims, ((grid_ds.tlon.data * 180 / np.pi + 180) % 360) - 180),
                              'area' : (self.spatial_dims, grid_ds.tarea.data)})
        G_u     = xr.Dataset({'lat'  : (self.spatial_dims, grid_ds.ulat.data * 180 / np.pi),
                              'lon'  : (self.spatial_dims, ((grid_ds.ulon.data * 180 / np.pi + 180) % 360) - 180),
                              'area' : (self.spatial_dims, grid_ds.uarea.data)})
        # original KMT
        KMT_org = xr.DataArray(data   = xr.open_dataset(self.P_KMT_org).kmt.values,
                               dims   = self.spatial_dims,
                               coords = {'lat': (self.spatial_dims, G_t['lat'].values),
                                         'lon': (self.spatial_dims, G_t['lon'].values)},
                               name   = 'kmt' )
        # modified KMT
        KMT_mod = xr.open_dataset(self.P_KMT_mod).kmt
        # just in case someone passes the original KMT as the modified KMT
        if not self.use_gi and 'ny' in KMT_mod.dims and 'nx' in KMT_mod.dims:
            KMT_mod = KMT_mod.rename({'ny': 'nj', 'nx': 'ni'})
        self.G_t     = G_t
        self.G_u     = G_u
        self.KMT_org = KMT_org
        self.KMT_mod = KMT_mod

    def process_raw_grounded_iceberg_database(self, GI_P_raw=None):
        """
        Load and spatially map raw grounded iceberg locations from CSV to the model grid.

        Uses a KD-tree nearest-neighbor search to assign iceberg coordinates to grid cells.
        Sets `self.GI_clean` with attached (nj, ni) indices.
        """
        self.load_grid_and_landmask()
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
