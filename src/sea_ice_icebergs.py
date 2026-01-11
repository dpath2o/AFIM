import os, shutil
import xarray      as xr
import numpy       as np
import pandas      as pd
import geopandas   as gpd
from scipy.spatial import cKDTree
from pathlib       import Path
from datetime      import datetime

class SeaIceIcebergs:
    
    def __init__(self, **kwargs):
        return
    
    def load_GI_lon_lats(self):
        self.load_cice_grid(slice_hem=True)

        lon  = self.G_t["lon"].reset_coords(drop=True)
        lat  = self.G_t["lat"].reset_coords(drop=True)
        mask = self.G_GI["mask"].reset_coords(drop=True).astype(bool)

        return {
            "lon": lon.data[mask.data].ravel(),
            "lat": lat.data[mask.data].ravel(),
        }
    
    def _region_lon_mask(self, lon2d, lon_min, lon_max):
        """
        Build a longitude mask that works for either:
        - non-wrapping bounds (lon_min <= lon_max)
        - seam-crossing bounds (lon_min > lon_max) e.g. 330..10 on a 0..360 grid
        """
        if lon_min <= lon_max:
            return (lon2d >= lon_min) & (lon2d <= lon_max)
        else:
            # crosses seam
            return (lon2d >= lon_min) | (lon2d <= lon_max)    
    
    def compute_grounded_iceberg_area(self, region=None, scale=1e6):
        """
        Returns:
        - dict of km^2 by region if region is provided
        - total GI area (m^2) if region is None (or km^2 if you divide by scale)
        """
        self.load_cice_grid(slice_hem=True)

        # IMPORTANT: use FULL-grid lon/lat for region tests (no NaNs)
        lon2d = self.G_t["lon"].values
        lat2d = self.G_t["lat"].values

        # GI area field already has NaNs outside GI in your loader
        area2d = self.G_GI["area"].values

        # IMPORTANT: convert mask back to boolean for indexing
        gi_mask = self.G_GI["mask"].values.astype(bool)

        # Logging (nan-safe)
        self.logger.info(f"GI-area_calc: grid lon min/max: {np.nanmin(lon2d):.2f} / {np.nanmax(lon2d):.2f}")
        self.logger.info(f"GI-area_calc: grid lat min/max: {np.nanmin(lat2d):.2f} / {np.nanmax(lat2d):.2f}")

        if region is None:
            # total GI area (m^2)
            return np.nansum(area2d[gi_mask])

        # Determine grid lon convention
        grid_is_0360 = (np.nanmin(lon2d) >= 0.0) and (np.nanmax(lon2d) > 180.0)

        area_dict = {}
        for reg_name, reg_cfg in region.items():
            lon_min, lon_max, lat_min, lat_max = reg_cfg["plot_region"]

            # Convert region bounds if grid is 0..360 but bounds are negative / -180..180 style
            if grid_is_0360:
                lon_min_c = self.normalise_longitudes(lon_min)
                lon_max_c = self.normalise_longitudes(lon_max)
            else:
                lon_min_c, lon_max_c = lon_min, lon_max

            lon_mask = self._region_lon_mask(lon2d, lon_min_c, lon_max_c)
            lat_mask = (lat2d >= lat_min) & (lat2d <= lat_max)
            region_mask = lon_mask & lat_mask

            # Sum GI areas within region
            area_km2 = np.nansum(area2d[gi_mask & region_mask]) / scale
            area_dict[reg_name] = area_km2

            self.logger.info(f"GI-area_calc: region {reg_name} GI area = {area_km2:0.2f} km^2")

        return area_dict

    def check_GI_coverage(self, da, varname="GI_counts", lon_name="lon", lat_name="lat"):
        """
        Assess grounded iceberg (GI) or landmask coverage across Antarctica.

        Parameters
        ----------
        da : xr.DataArray
            The GI or KMT data.
        varname : str
            Determines masking logic: 'GI_counts' uses >0, anything else uses >1.
        lon_name : str, default 'lon'
            Name of the longitude coordinate.
        lat_name : str, default 'lat'
            Name of the latitude coordinate.

        Returns
        -------
        report : dict
            Dictionary of summary statistics and Antarctic sector coverage.
        """
        report = {}
        lon = da[lon_name]
        lat = da[lat_name]
        lon_vals = lon.values
        lat_vals = lat.values
        # Handle longitude wrapping
        if lon_vals.max() <= 180:
            lon_wrapped = np.mod(lon_vals + 360, 360)
            if lon_vals.ndim == 1:
                da = da.assign_coords({lon_name: (da[lon_name].dims, lon_wrapped)})
                da = da.sortby(lon_name)
            elif lon_vals.ndim == 2:
                da = da.assign_coords({lon_name: (da[lon_name].dims, lon_wrapped)})
                self.logger.warning("Skipping sortby(): longitude is 2D.")
            else:
                raise ValueError("Longitude data must be 1D or 2D.")
        # Construct binary coverage mask
        if varname == "GI_counts":
            mask = da.values > 0
        else:
            mask = da.values > 1
        # Basic spatial stats
        report['lon_range'] = (float(np.nanmin(lon_vals)), float(np.nanmax(lon_vals)))
        report['lat_range'] = (float(np.nanmin(lat_vals)), float(np.nanmax(lat_vals)))
        report['n_valid_cells'] = int(np.count_nonzero(mask))
        report['n_total_cells'] = int(np.prod(mask.shape))
        report['coverage_fraction'] = report['n_valid_cells'] / report['n_total_cells']
        report['nonfinite_count'] = int(np.sum(~np.isfinite(da.values)))
        # Value stats
        report['mean_value'] = float(np.nanmean(da.values)) if np.any(mask) else float('nan')
        report['max_value'] = float(np.nanmax(da.values)) if np.any(mask) else float('nan')
        report['min_nonzero'] = float(np.nanmin(da.values[da.values > 0])) if np.any(da.values > 0) else float('nan')
        # Sector-based coverage check
        if lon_vals.ndim == 2 and lat_vals.ndim == 2:
            try:
                sector_edges = np.linspace(0, 360, 9)  # 8 sectors
                presence_by_sector = []
                for i in range(len(sector_edges) - 1):
                    lon_min, lon_max = sector_edges[i], sector_edges[i + 1]
                    sector_mask = (
                        (lon_vals >= lon_min) & (lon_vals < lon_max) &
                        (lat_vals < -60)  # Limit to Antarctic region
                    )
                    if np.any(sector_mask):
                        sector_data = da.values[sector_mask]
                        presence_by_sector.append(np.any(sector_data > 0))
                    else:
                        presence_by_sector.append(False)
                report['sector_gaps'] = [i for i, present in enumerate(presence_by_sector) if not present]
                report['sector_gap_count'] = len(report['sector_gaps'])
            except Exception as e:
                report['sector_check_failed'] = str(e)
        return report

    def load_existing_thinned(self):
        """
        Load previously saved thinned grounded iceberg data and modified KMT file.

        Sets the attributes `self.GI_thin_da`, `self.GI_thin_mask`, and `self.GI_KMT`
        if their respective NetCDF files exist.
        """
        if os.path.exists(self.P_GI_thin):
            self.GI_thin_da = xr.open_dataarray(self.P_GI_thin)
            self.logger.info(f"GROUNDED ICEBERG COUNTS PER GRID CELL REPORT {self.P_GI_thin}")
            self.logger.info(self.check_GI_coverage(self.GI_thin_da))
            self.GI_thin_mask = self.GI_thin_da.values > 0
        if os.path.exists(self.P_KMT_mod):
            self.GI_KMT = xr.open_dataarray(self.P_KMT_mod)
            self.logger.info(f"MODIFIED LANDMASK REPORT {self.P_KMT_mod}")
            self.logger.info(self.check_GI_coverage(self.GI_KMT, varname="kmt"))

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
        kmt                   = xr.open_dataset(self.P_KMT_org).kmt.values
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

    def regrid_gi_to_fip_bool(self,
                              gi_mask_da: xr.DataArray,
                                gi_lon2d: xr.DataArray,
                                gi_lat2d: xr.DataArray,
                                fip_lon2d: xr.DataArray,
                                fip_lat2d: xr.DataArray,
                                dilate_pixels: int = 1):   # set 0 to disable adjacency expansion
        """
        Regrids a boolean GI mask from the model grid onto the FIP grid via nearest-neighbour
        on great-circle distance (xesmf 'nearest_s2d'). Returns a boolean mask on the FIP grid.
        Optionally dilates by N pixels on the target grid to include coast-adjacent cells.
        """
        from scipy.ndimage import binary_dilation
        import xesmf as xe
        # xESMF expects float fields; cast True->1, False->0
        gi_float = gi_mask_da.astype(np.float32)

        # Build source/target 'grid' dictionaries from 2D lon/lat
        src = {'lon': gi_lon2d.values, 'lat': gi_lat2d.values}
        dst = {'lon': fip_lon2d.values, 'lat': fip_lat2d.values}

        # Regridder (nearest source-to-destination). reuse_weights=True if you call repeatedly.
        regridder = xe.Regridder(src, dst, method='nearest_s2d', locstream_in=False, locstream_out=False)

        gi_on_fip_float = regridder(gi_float.values)  # ndarray on target shape
        gi_on_fip = xr.DataArray(gi_on_fip_float > 0.5, coords=fip_lon2d.coords, dims=fip_lon2d.dims)

        # Optional: 1-pixel dilation to catch coast-adjacent cells next to GI
        if dilate_pixels and dilate_pixels > 0:
            # ndimage expects numpy array; preserve NaNs as False for mask ops
            mask_np = np.asarray(gi_on_fip.fillna(False).values, dtype=bool)
            structure = np.ones((1 + 2*dilate_pixels, 1 + 2*dilate_pixels), dtype=bool)
            mask_np = binary_dilation(mask_np, structure=structure)
            gi_on_fip = xr.DataArray(mask_np, coords=fip_lon2d.coords, dims=fip_lon2d.dims)

        gi_on_fip.name = "GI_mask_on_FIP"
        gi_on_fip.attrs.update(dict(long_name="Grounded iceberg mask (regridded to FIP)", comment="nearest_s2d + optional dilation"))
        return gi_on_fip