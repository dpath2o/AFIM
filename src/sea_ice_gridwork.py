import os
import xarray as xr
import numpy  as np
import pandas as pd
from pathlib  import Path
__all__ = ["SeaIceGridWork"]
class SeaIceGridWork:
    """
    Grid utilities and CICE grid loaders for AFIM sea-ice workflows.

    `SeaIceGridWork` centralizes common spatial operations used across the AFIM
    toolbox, including longitude normalization, landmask application, construction
    of cell-corner and face-coordinate geometry, and loading CICE grid metadata
    (T-grid and U-grid) from the model grid file.

    The class is intended to be mixed into higher-level analysis classes (e.g.,
    toolbox/plotter/regridder) and therefore typically relies on attributes defined
    elsewhere, such as:
    - `self.CICE_dict` : configuration dictionary (paths, dimension names, etc.)
    - `self.logger`   : logger instance
    - `self.use_gi`   : bool; whether grounded-iceberg modified landmask is active
    - `self.P_KMT_org`, `self.P_KMT_mod` : paths to original/modified landmask files
    - `self.slice_hemisphere(...)` : optional helper to subset to SH/NH
    - `self.radians_to_degrees(...)` : helper for coordinate conversion
    - state flags: `self.grid_loaded`, `self.bgrid_loaded`, `self.cgrid_loaded`

    Notes
    -----
    - CICE grids are curvilinear; longitudes require seam-aware handling. The class
    provides robust longitude wrapping and a circular mean for longitudes.
    - `load_cice_grid` sets multiple datasets as **attributes**; it is a loader with
    side effects and does not return the grids explicitly.
    """

    def __init__(self, **kwargs):
        """
        Initialize the gridwork mixin.

        This initializer is intentionally minimal. In typical AFIM usage, configuration
        and state (e.g., `self.CICE_dict`, logging, hemisphere settings, and file paths)
        are injected via `**kwargs` and/or by parent classes.

        Parameters
        ----------
        **kwargs
            Arbitrary keyword arguments forwarded from a parent class initializer.

        Notes
        -----
        This method does not validate configuration. Callers should ensure required
        attributes (e.g., `self.CICE_dict`, `self.logger`) are set prior to invoking
        grid-dependent methods.
        """
        return

    def normalise_longitudes(self, lon, to="0-360", eps=1e-12):
        """
        Normalize longitudes to a consistent wrap convention.

        This routine wraps longitude values either to:
        - ``"0-360"``     : [0, 360)
        - ``"-180-180"``  : [-180, 180)

        It supports scalars, NumPy arrays, and xarray objects (DataArray / Dataset
        variables). NaNs are preserved.

        Parameters
        ----------
        lon : scalar, numpy.ndarray, or xarray.DataArray
            Longitude values in degrees east. May be any numeric dtype; NaNs pass through.
        to : {"0-360", "-180-180"}, default "0-360"
            Target longitude convention.
        eps : float, default 1e-12
            Numerical tolerance used to collapse boundary values (e.g., 360→0, 180→-180)
            when floating-point rounding produces near-boundary values.

        Returns
        -------
        same type as `lon`
            Longitude values wrapped to the requested convention.

        Raises
        ------
        ValueError
            If `to` is not one of {"0-360", "-180-180"}.

        Notes
        -----
        - For the "-180-180" convention, values exactly equal to +180 (within `eps`) are
        mapped to -180 for consistency with a half-open interval [-180, 180).
        """
        # First get [0, 360)
        lon_wrapped = ((lon % 360) + 360) % 360  # safe for negatives, NaNs pass through
        if to == "0-360":
            # Collapse values extremely close to 360 back to 0
            if isinstance(lon_wrapped, xr.DataArray):
                lon_wrapped = xr.where(np.isclose(lon_wrapped, 360.0, atol=eps), 0.0, lon_wrapped)
            else:
                lon_wrapped = np.where(np.isclose(lon_wrapped, 360.0, atol=eps), 0.0, lon_wrapped)
            return lon_wrapped
        elif to == "-180-180":
            lon_180 = ((lon_wrapped + 180.0) % 360.0) - 180.0  # -> (-180, 180]
            # Prefer [-180, 180) by mapping exactly 180 to -180
            if isinstance(lon_180, xr.DataArray):
                lon_180 = xr.where(np.isclose(lon_180, 180.0, atol=eps), -180.0, lon_180)
            else:
                lon_180 = np.where(np.isclose(lon_180, 180.0, atol=eps), -180.0, lon_180)
            return lon_180
        else:
            raise ValueError("to must be '0-360' or '-180-180'")
        
    def _xy_to_lonlat(self, x, y):
        """
        Convert EPSG:3031 projected coordinates to geographic lon/lat (EPSG:4326).

        Parameters
        ----------
        x, y : array-like
            Projected coordinates in EPSG:3031 meters.

        Returns
        -------
        (lon, lat) : tuple of array-like
            Longitude and latitude in degrees.
        """
        from pyproj import Transformer
        T = Transformer.from_crs(3031, 4326, always_xy=True)
        lon, lat = T.transform(x, y)
        return lon, lat
    
    def _match_lon_range(self, lon_1d: np.ndarray, grid_lon_2d: np.ndarray) -> np.ndarray:
        """
        Map swath longitudes onto the same wrap convention as a target grid.

        Given a 1-D longitude array (often from swath samples) and a target grid longitude
        field, this chooses between:
        - [-180, 180) wrapping, or
        - [0, 360) wrapping,
        based on the numeric range of the target grid longitudes.

        Parameters
        ----------
        lon_1d : numpy.ndarray
            1-D longitudes of swath samples (degrees east).
        grid_lon_2d : numpy.ndarray
            2-D longitudes of the target grid (degrees east).

        Returns
        -------
        numpy.ndarray
            Longitudes mapped to the target grid convention.
        """
        gmin, gmax = np.nanmin(grid_lon_2d), np.nanmax(grid_lon_2d)
        grid_is_180 = (gmin >= -180.0) and (gmax <= 180.0)
        lon_a = ((lon_1d + 180.0) % 360.0) - 180.0   # [-180, 180)
        lon_b = lon_1d % 360.0                       # [0, 360)
        lon_b[lon_b < 0] += 360.0
        # pick the representation whose range best matches the grid
        if grid_is_180:
            return lon_a
        else:
            return lon_b

    def reapply_landmask(self, DS, apply_unmodified=False):
        """
        Apply the CICE land mask to all spatial variables in a dataset.

        This method masks land points by setting values to NaN for all variables that
        contain both spatial dimensions ``("nj", "ni")``. The mask is constructed from
        either the original landmask (kmt_org) or the grounded-iceberg modified landmask
        (kmt_mod), depending on configuration.

        Parameters
        ----------
        DS : xarray.Dataset
            Dataset containing sea-ice fields. Any variable whose dimensions include
            both "nj" and "ni" will be masked.
        apply_unmodified : bool, default False
            If True, always apply the **original** landmask (`P_KMT_org`) even when
            `self.use_gi` is True. If False (default) and `self.use_gi` is True, apply
            the modified landmask (`P_KMT_mod`).

        Returns
        -------
        xarray.Dataset
            The same dataset instance with land cells set to NaN for spatial variables.

        Side Effects
        ------------
        - Opens the landmask dataset from disk via `xarray.open_dataset`.
        - Modifies `DS` in place by assigning masked variables.

        Notes
        -----
        - The landmask is assumed to be encoded such that ocean == 1 and land == 0
        (or equivalent). The mask applied is ``kmt == 1``.
        - Non-spatial variables (without both "nj" and "ni") are left unchanged.
        """
        if self.use_gi and not apply_unmodified:
            kmt_mask = xr.open_dataset(self.P_KMT_mod)['kmt'] == 1
        else:
            kmt_mask = xr.open_dataset(self.P_KMT_org)['kmt'] == 1
        for var in DS.data_vars:
            da = DS[var]
            if {"nj", "ni"}.issubset(da.dims):  # Only apply to spatial fields
                self.logger.debug(f"Masking land for variable: {var}")
                DS[var] = da.where(kmt_mask)
        self.logger.info("Applied landmask to rolled dataset")
        return DS
    
    def _circular_mean_lon_deg(self, *lon_deg_arrays):
        """
        Compute a seam-safe circular mean of longitudes in degrees.

        This helper computes the mean longitude using a circular (vector) average:
        1) wrap all longitudes into [0, 360)
        2) convert to radians
        3) average cos/sin components
        4) convert back to degrees and re-wrap to [0, 360)

        Parameters
        ----------
        *lon_deg_arrays : array-like
            One or more longitude arrays in degrees. Inputs must be broadcastable to a
            common shape (e.g., multiple (nj, ni) arrays).

        Returns
        -------
        numpy.ndarray
            Circular mean longitude in degrees wrapped to [0, 360), with the broadcasted
            shape of the input arrays.

        Notes
        -----
        This is preferred over arithmetic averaging near the dateline (e.g., averaging
        359° and 1° should yield 0°, not 180°).
        """
        lon = np.stack(lon_deg_arrays, axis=0)  # (n, ...)
        lon = self.normalise_longitudes(lon, to="0-360")
        lonr = np.deg2rad(lon)
        x = np.cos(lonr).mean(axis=0)
        y = np.sin(lonr).mean(axis=0)
        out = np.rad2deg(np.arctan2(y, x))
        return self.normalise_longitudes(out, to="0-360")

    def build_grid_faces(self, lon, lat, source_in_radians=False):
        """
        Construct corner and face-coordinate geometry from cell-centre lon/lat fields.

        Given 2-D cell-centre longitudes and latitudes (typically the CICE T-grid
        centres), this routine constructs:
        - cell corner coordinates (``lon_b``, ``lat_b``) with shape (nj+1, ni+1)
        - vertical face-centre coordinates (``lon_e``, ``lat_e``) with shape (nj, ni+1)
        - horizontal face-centre coordinates (``lon_n``, ``lat_n``) with shape (nj+1, ni)

        Longitudes are handled using a circular mean to avoid dateline artifacts.

        Parameters
        ----------
        lon, lat : array-like
            2-D arrays of cell-centre coordinates with shape (nj, ni). Units are degrees
            unless `source_in_radians=True`.
        source_in_radians : bool, default False
            If True, inputs are in radians and will be converted to degrees before
            geometry construction.

        Returns
        -------
        lon_b, lat_b : numpy.ndarray
            Corner coordinates with shape (nj+1, ni+1) in degrees.
        lon_e, lat_e : numpy.ndarray
            Vertical face-centre coordinates with shape (nj, ni+1) in degrees.
        lon_n, lat_n : numpy.ndarray
            Horizontal face-centre coordinates with shape (nj+1, ni) in degrees.

        Notes
        -----
        - Interior corners are computed as the mean of the four surrounding cell centres.
        Longitude uses a circular mean; latitude uses an arithmetic mean.
        - Edge corners are computed from adjacent two-cell averages (more consistent than
        copying boundary cells).
        - The four outermost corners are set to the nearest cell-centre values.
        """
        if source_in_radians:
            lon = self.radians_to_degrees(lon)
            lat = self.radians_to_degrees(lat)
        lon = self.normalise_longitudes(lon, to="0-360")
        lat = np.asarray(lat)

        nj, ni = lat.shape

        # --- corners (nj+1, ni+1)
        lon_b = np.full((nj + 1, ni + 1), np.nan, dtype=float)
        lat_b = np.full((nj + 1, ni + 1), np.nan, dtype=float)

        # Interior corners: mean of surrounding 4 cell centres
        lat_b[1:-1, 1:-1] = 0.25 * (lat[:-1, :-1] + lat[:-1, 1:] + lat[1:, :-1] + lat[1:, 1:])
        lon_b[1:-1, 1:-1] = self._circular_mean_lon_deg(lon[:-1, :-1], lon[:-1, 1:], lon[1:, :-1], lon[1:, 1:])

        # Edges: mean of adjacent 2 cell centres (more consistent than copying)
        # top edge (j=0): between cells (0, i-1) and (0, i)
        lat_b[0, 1:-1] = 0.5 * (lat[0, :-1] + lat[0, 1:])
        lon_b[0, 1:-1] = self._circular_mean_lon_deg(lon[0, :-1], lon[0, 1:])

        # bottom edge (j=nj): between cells (nj-1, i-1) and (nj-1, i)
        lat_b[-1, 1:-1] = 0.5 * (lat[-1, :-1] + lat[-1, 1:])
        lon_b[-1, 1:-1] = self._circular_mean_lon_deg(lon[-1, :-1], lon[-1, 1:])

        # left edge (i=0): between cells (j-1,0) and (j,0)
        lat_b[1:-1, 0] = 0.5 * (lat[:-1, 0] + lat[1:, 0])
        lon_b[1:-1, 0] = self._circular_mean_lon_deg(lon[:-1, 0], lon[1:, 0])

        # right edge (i=ni): between cells (j-1,ni-1) and (j,ni-1)
        lat_b[1:-1, -1] = 0.5 * (lat[:-1, -1] + lat[1:, -1])
        lon_b[1:-1, -1] = self._circular_mean_lon_deg(lon[:-1, -1], lon[1:, -1])

        # Four outer corners: nearest cell centre
        lat_b[0, 0] = lat[0, 0]
        lat_b[0, -1] = lat[0, -1]
        lat_b[-1, 0] = lat[-1, 0]
        lat_b[-1, -1] = lat[-1, -1]

        lon_b[0, 0] = lon[0, 0]
        lon_b[0, -1] = lon[0, -1]
        lon_b[-1, 0] = lon[-1, 0]
        lon_b[-1, -1] = lon[-1, -1]

        lon_b = self.normalise_longitudes(lon_b, to="0-360")

        # --- face centres from corners
        # Vertical faces (E/W): average corners above/below along a vertical edge
        # shape (nj, ni+1) = (nj, ni_b)
        lat_e = 0.5 * (lat_b[:-1, :] + lat_b[1:, :])
        lon_e = self._circular_mean_lon_deg(lon_b[:-1, :], lon_b[1:, :])

        # Horizontal faces (N/S): average corners left/right along a horizontal edge
        # shape (nj+1, ni) = (nj_b, ni)
        lat_n = 0.5 * (lat_b[:, :-1] + lat_b[:, 1:])
        lon_n = self._circular_mean_lon_deg(lon_b[:, :-1], lon_b[:, 1:])

        return lon_b, lat_b, lon_e, lat_e, lon_n, lat_n


    def build_grid_corners(self, lat_rads, lon_rads, grid_res=0.25, source_in_radians=True):
        """
        Construct corner coordinates from cell-centre coordinates using simple averaging.

        This routine builds (ny+1, nx+1) corner arrays from (ny, nx) cell-centre arrays
        via arithmetic averaging. It is a legacy/simple alternative to `build_grid_faces`.

        Parameters
        ----------
        lat_rads, lon_rads : array-like
            2-D arrays of cell-centre lat/lon. Interpreted as radians when
            `source_in_radians=True`, else as degrees.
        grid_res : float, default 0.25
            Scaling factor applied to the interior arithmetic mean. Historically, this
            function used `grid_res` as a multiplier; for conventional corner averaging
            this should be 0.25.
        source_in_radians : bool, default True
            If True, interpret inputs as radians and convert to degrees.

        Returns
        -------
        lon_b, lat_b : numpy.ndarray
            Corner longitudes/latitudes with shape (ny+1, nx+1) in degrees.

        Notes
        -----
        - Longitudes are normalized using `normalise_longitudes` after construction.
        - This method uses arithmetic averaging for longitude and is not seam-safe near
        the dateline. Prefer `build_grid_faces` when longitude seam handling matters.
        """
        if source_in_radians:
            lon = self.radians_to_degrees(lon_rads)
            lat = self.radians_to_degrees(lat_rads)
        else:
            lon = lon_rads
            lat = lat_rads
        ny, nx               = lat.shape
        lat_b                = np.zeros((ny + 1, nx + 1))
        lon_b                = np.zeros((ny + 1, nx + 1))
        lat_b[1:-1, 1:-1]    = grid_res * (lat[:-1, :-1] + lat[:-1, 1:] + lat[1:, :-1] + lat[1:, 1:])
        lon_b[1:-1, 1:-1]    = grid_res * (lon[:-1, :-1] + lon[:-1, 1:] + lon[1:, :-1] + lon[1:, 1:])
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
        lon_b                = self.normalise_longitudes(lon_b)
        return lon_b, lat_b

    def load_cice_grid(self, slice_hem=False, build_faces=True):
        """
        Load CICE grid metadata and construct T-grid/U-grid datasets and (optionally) face grids.

        This unified loader reads the CICE grid file and associated landmask files, then
        constructs and stores as attributes:

        - `self.G_t` : T-grid (cell centres) dataset with dims (nj, ni)
        - `self.G_u` : U-grid (B-grid velocity points) dataset with dims:
                    centres (nj, ni) and corners (nj_b, ni_b) = (nj+1, ni+1)
        - `self.G_e` : C-grid vertical face coordinates derived from T-grid corners,
                    with dims (nj, ni_b) = (nj, ni+1) when `build_faces=True`
        - `self.G_n` : C-grid horizontal face coordinates derived from T-grid corners,
                    with dims (nj_b, ni) = (nj+1, ni) when `build_faces=True`

        Landmask handling:
        - `self.kmt_org` is always loaded (original landmask).
        - If `self.use_gi` is True, `self.kmt_mod` is loaded (modified landmask including
        grounded icebergs), and an additional dataset `self.G_GI` is constructed to
        describe grounded-iceberg cells and their properties.

        Parameters
        ----------
        slice_hem : bool, default False
            If True, apply `self.slice_hemisphere(...)` to all constructed datasets
            (grids, masks, and grounded-iceberg dataset) to subset to the configured
            hemisphere.
        build_faces : bool, default True
            If True, compute face-coordinate datasets `self.G_e` and `self.G_n` using
            `build_grid_faces(...)`. If False, these are set to None and `self.cgrid_loaded`
            is False.

        Returns
        -------
        None
            Grids and masks are stored on the instance as attributes.

        Side Effects
        ------------
        - Opens multiple NetCDF files from disk using xarray:
        - the main CICE grid file (`self.CICE_dict["P_G"]`)
        - the original landmask (`self.P_KMT_org`)
        - the modified landmask (`self.P_KMT_mod`) if `self.use_gi` is True
        - Sets internal state flags:
        - `self.bgrid_loaded = True`
        - `self.cgrid_loaded = bool(build_faces)`
        - `self.grid_loaded  = True`

        Attributes Set
        --------------
        G_t : xarray.Dataset
            Variables: lat, lon, angle, area on (nj, ni).
        G_u : xarray.Dataset
            Variables: lat, lon, lat_b, lon_b, angle, area with centre dims (nj, ni) and
            corner dims (nj_b, ni_b) = (nj+1, ni+1).
        G_e, G_n : xarray.Dataset or None
            Face-coordinate datasets created when `build_faces=True`.
        kmt_org : xarray.Dataset
            Binary mask dataset on (nj, ni) describing the original ocean/land mask.
        kmt_mod : xarray.Dataset, optional
            Binary mask dataset on (nj, ni) describing the modified mask including
            grounded icebergs (only when `self.use_gi` is True).
        G_GI : xarray.Dataset, optional
            Grounded-iceberg diagnostics dataset, including:
            - `mask(nj,ni)` boolean mask of grounded-iceberg cells
            - `area/lat/lon(nj,ni)` with NaN outside GI cells
            - 1-D lists `area_nGI`, `lat_nGI`, `lon_nGI` for GI cell locations
            Only present when `self.use_gi` is True.

        Notes
        -----
        - Longitudes are normalized to 0–360 degrees to avoid seam ambiguity.
        - Corner/face coordinates are derived from T-grid centres to enable consistent
        C-grid geometry for plotting and regridding.
        - The grounded-iceberg 1-D coordinate lists include an index shift in `ni` (westward)
        to better align with B-grid layout conventions used elsewhere in the toolbox.
        """
        if self.grid_loaded:
            self.logger.info("CICE grid previously loaded ... returning")
            return
        nat_dim = self.CICE_dict["spatial_dims"]      # e.g. ("nj","ni")
        ext_dim = tuple(f"{d}_b" for d in nat_dim)    # e.g. ("nj_b","ni_b")
        with xr.open_dataset(self.CICE_dict["P_G"]) as G:
            # land masks
            kmt_org = xr.open_dataset(self.P_KMT_org).kmt.data
            kmt_mod = xr.open_dataset(self.P_KMT_mod).kmt.data if self.use_gi else kmt_org
            # --- T grid (centres)
            TLAT = self.radians_to_degrees(G["tlat"].data)
            TLON = self.normalise_longitudes(self.radians_to_degrees(G["tlon"].data), to="0-360")
            TAREA = G["tarea"].data
            T_ANGLE = self.radians_to_degrees(G["angleT"].data)
            # corners + faces from t-grid centres
            ULON_b, ULAT_b, lon_e, lat_e, lon_n, lat_n = self.build_grid_faces(TLON, TLAT, source_in_radians=False)
            # --- U grid (legacy / B-grid velocity points)
            ULAT = self.radians_to_degrees(G["ulat"].data)
            ULON = self.normalise_longitudes(self.radians_to_degrees(G["ulon"].data), to="0-360")
            UAREA = G["uarea"].data
            U_ANGLE = self.radians_to_degrees(G["angle"].data)
        # dims
        nj, ni = TLAT.shape
        coords_nat = {nat_dim[0]: np.arange(nj),
                      nat_dim[1]: np.arange(ni)}
        coords_ext = {nat_dim[0]: np.arange(nj),
                      nat_dim[1]: np.arange(ni),
                      ext_dim[0]: np.arange(nj + 1),
                      ext_dim[1]: np.arange(ni + 1)}
        # landmasks
        self.kmt_org = xr.Dataset(data_vars = {"kmt_org": (nat_dim, kmt_org, {"units"    : "binary", 
                                                                              "long_name": "original landmask on t-grid"})},
                                  coords    = coords_nat)
        # Grounded icebergs
        if self.use_gi:
            self.kmt_mod = xr.Dataset(data_vars = {"kmt_mod": (nat_dim, kmt_mod, {"units"    : "binary", 
                                                                                "long_name": "modified landmask on t-grid (grounded icebergs)"})},
                                    coords    = coords_nat)
            # Difference: grounded icebergs are cells that changed from ocean (1) to land (0)
            GI_mask = (kmt_org == 1) & (kmt_mod == 0)
            # 2-D fields with same shape as GI_mask (NaN outside GI)
            GI_area = np.where(GI_mask, TAREA, np.nan).astype("float32")
            GI_lat  = np.where(GI_mask, TLAT,  np.nan).astype("float32")
            GI_lon  = np.where(GI_mask, TLON,  np.nan).astype("float32")
            self.G_GI = xr.Dataset(data_vars={"mask": (nat_dim, GI_mask.astype("uint8"), {"units": "binary",
                                                                                          "long_name": "grounded iceberg mask on t-grid"}),
                                              "area": (nat_dim, GI_area, {"units": "m2",
                                                                          "long_name": "t-cell area at grounded iceberg locations (NaN elsewhere)"}),
                                              "lat":  (nat_dim, GI_lat,  {"units": "degrees_north",
                                                                          "long_name": "latitude at grounded iceberg locations (NaN elsewhere)"}),
                                              "lon":  (nat_dim, GI_lon,  {"units": "degrees_east",
                                                                          "long_name": "longitude at grounded iceberg locations (NaN elsewhere)"})},
                                   coords = coords_nat)
            # Get coordinates of affected cells (shifted west by one ni index to match B-grid layout)
            nj_idx, ni_idx = np.where(GI_mask)
            ni_idx_shifted = ni_idx - 1
            valid          = ni_idx_shifted >= 0
            nj_idx         = nj_idx[valid]
            ni_idx         = ni_idx_shifted[valid]
            # Save 1D arrays with iceberg IDs
            self.G_GI             = self.G_GI.assign_coords(nGI=np.arange(nj_idx.size))
            self.G_GI["area_nGI"] = (("nGI",), TAREA[nj_idx, ni_idx].astype("float32"))
            self.G_GI["lat_nGI"]  = (("nGI",), TLAT[nj_idx, ni_idx].astype("float32"))
            self.G_GI["lon_nGI"]  = (("nGI",), TLON[nj_idx, ni_idx].astype("float32"))
            self.G_GI["area_nGI"].attrs.update({"units": "m2", "long_name": "t-cell area (1D list) at GI cells"})
            self.G_GI["lat_nGI"].attrs.update({"units": "degrees_north", "long_name": "GI latitude (1D list)"})
            self.G_GI["lon_nGI"].attrs.update({"units": "degrees_east",  "long_name": "GI longitude (1D list)"})
        # --- build datasets
        self.G_t = xr.Dataset(data_vars = {"lat"    : (nat_dim, TLAT, {"units": "degrees"}),
                                           "lon"    : (nat_dim, TLON, {"units": "degrees"}),
                                           "angle"  : (nat_dim, T_ANGLE, {"units": "degrees"}),
                                           "area"   : (nat_dim, TAREA, {"units": "m^2"})},
                            coords    = coords_nat)
        self.G_u = xr.Dataset(data_vars = {"lat"    : (nat_dim, ULAT, {"units": "degrees"}),
                                           "lon"    : (nat_dim, ULON, {"units": "degrees"}),
                                           "lat_b"  : (ext_dim, ULAT_b, {"units": "degrees"}),
                                           "lon_b"  : (ext_dim, ULON_b, {"units": "degrees"}),
                                           "angle"  : (nat_dim, U_ANGLE, {"units": "degrees"}),
                                           "area"   : (nat_dim, UAREA, {"units": "m^2"})},
                            coords    = coords_ext)
        if build_faces:
            # Use ext_dim/nat_dim naming so shapes are explicit:
            # lon_e/lat_e: (nj, ni_b) == (nat_dim[0], ext_dim[1])
            # lon_n/lat_n: (nj_b, ni) == (ext_dim[0], nat_dim[1])
            self.G_e = xr.Dataset(data_vars = {"lat": ((nat_dim[0], ext_dim[1]), lat_e, {"units": "degrees"}),
                                               "lon": ((nat_dim[0], ext_dim[1]), lon_e, {"units": "degrees"})},
                                  coords    = coords_ext)
            self.G_n = xr.Dataset(data_vars = {"lat": ((ext_dim[0], nat_dim[1]), lat_n, {"units": "degrees"}),
                                               "lon": ((ext_dim[0], nat_dim[1]), lon_n, {"units": "degrees"})},
                                  coords    = coords_ext)
        else:
            self.G_e = None
            self.G_n = None
        if slice_hem:
            self.G_t = self.slice_hemisphere(self.G_t)
            self.G_u = self.slice_hemisphere(self.G_u)
            if self.G_e is not None:
                self.G_e = self.slice_hemisphere(self.G_e)
            if self.G_n is not None:
                self.G_n = self.slice_hemisphere(self.G_n)
            if self.use_gi and hasattr(self, "G_GI") and (self.G_GI is not None):
                self.G_GI = self.slice_hemisphere(self.G_GI)
            if hasattr(self, "kmt_org") and (self.kmt_org is not None):
                self.kmt_org = self.slice_hemisphere(self.kmt_org)
            if self.use_gi and hasattr(self, "kmt_mod") and (self.kmt_mod is not None):
                self.kmt_mod = self.slice_hemisphere(self.kmt_mod)
        # flags for downstream
        self.bgrid_loaded = True
        self.cgrid_loaded = bool(build_faces)
        self.grid_loaded  = True
        
    def compute_regular_grid_area(self, da):
        """
        Compute cell areas (km^2) for a regular lat/lon grid (1D coords).

        Assumes 1D, monotonically increasing latitude and longitude coordinates with
        uniform longitudinal spacing. Uses spherical geometry to integrate the area
        between latitude edges for each latitude band, broadcasting across longitudes.

        Parameters
        ----------
        da : xr.DataArray
            Any DataArray that carries 1D coordinates named ``'lat'`` and ``'lon'``.
            Values are not used; only the coordinates are read.

        Returns
        -------
        xr.DataArray
            2D area array with dims (``lat``, ``lon``) in **km²**, aligned to the
            provided coordinate vectors.

        Notes
        -----
        - Earth radius is fixed at 6,371,000 m.
        - For non-uniform lon spacing or curvilinear grids, use model grid areas
        (e.g., CICE `tarea`) rather than this helper.
        """
        R   = 6371000.0  # Earth radius in meters
        lat = np.deg2rad(da['lat'].values)
        lon = np.deg2rad(da['lon'].values)
        # Latitude edges
        lat_edges       = np.zeros(len(lat)+1)
        lat_edges[1:-1] = (lat[:-1] + lat[1:]) / 2
        lat_edges[0]    = lat[0] - (lat[1]-lat[0])/2
        lat_edges[-1]   = lat[-1] + (lat[-1]-lat[-2])/2
        # Longitude spacing
        dlon = lon[1] - lon[0]
        # 1D cell area per latitude
        dA_lat = (R**2) * dlon * (np.sin(lat_edges[1:]) - np.sin(lat_edges[:-1]))
        # Broadcast to 2D
        area_2d = np.tile(dA_lat[:, np.newaxis], (1, len(lon)))
        # Convert to km^2
        area_2d /= 1e6
        area_da = xr.DataArray(area_2d, dims=("lat","lon"), coords={"lat":da['lat'], "lon":da['lon']})
        return area_da
    
    ######################################################################
    #                           FORM FACTORS                             #
    ######################################################################
    def _infer_deg_from_grid_units(self, arr, name="var"):
        """
        Infer whether an angular array is in radians and convert to degrees if so.

        This helper applies a simple magnitude heuristic to distinguish radians from
        degrees for grid longitude/latitude/angle fields. If the maximum absolute
        finite value is approximately within 2π, the input is treated as radians and
        converted to degrees; otherwise it is returned unchanged.

        Parameters
        ----------
        arr : array-like
            Input values representing an angular quantity (e.g., lon, lat).
        name : str, default="var"
            Variable name used for logging/debug context (no logging performed here,
            but retained for consistency and future use).

        Returns
        -------
        out : numpy.ndarray
            Array of dtype float64. If inferred radians, returned in degrees; else
            returned in original units.
        """
        a = np.asarray(arr, dtype="float64")
        finite = np.isfinite(a)
        if not finite.any():
            return a
        amax = np.nanmax(np.abs(a[finite]))
        # Heuristic: radians for lon/lat/angle typically within ~2*pi
        if amax <= (2.0 * np.pi + 1e-6):
            return np.rad2deg(a)
        return a

    def _to_meters(self, arr, units_attr: str | None, name=""):
        """
        Convert a length-like array to meters using units metadata when available.

        Supported explicit units are centimeters ("cm") and meters ("m"). If units are
        missing/unknown, a magnitude heuristic is used:
        - median > 1e3  : assume meters
        - median > 10   : ambiguous; assume meters and emit a warning
        - otherwise     : suspicious; assume meters and emit a warning

        Parameters
        ----------
        arr : array-like
            Input length values.
        units_attr : str or None
            Units attribute string (e.g., from `DataArray.attrs.get("units")`).
        name : str, default=""
            Variable name for warning messages.

        Returns
        -------
        out : numpy.ndarray
            Values converted to meters (float64).

        Notes
        -----
        The heuristic is intentionally conservative (defaults to meters) because grid
        metrics for Antarctic climate configurations are typically O(1e4) meters.
        If your grid uses kilometers or other conventions, supply a correct units
        attribute to avoid ambiguity.
        """
        a = np.asarray(arr, dtype="float64")
        u = (units_attr or "").strip().lower()
        if u in ("cm", "centimeter", "centimeters"):
            return a * 0.01
        if u in ("m", "meter", "meters"):
            return a
        # fallback heuristic by magnitude (typical dx ~ O(1e4) m for 0.25° near Antarctica)
        med = np.nanmedian(a[np.isfinite(a)])
        if med > 1e3:   # likely meters already
            return a
        if med > 10:    # ambiguous; could be km or something odd
            self.logger.warning(f"{name} units ambiguous (median={med}); assuming meters.")
            return a
        # if it's ~250 (and you expect 25 km), could be km*10? etc.
        self.logger.warning(f"{name} units suspicious (median={med}); please verify units.")
        return a

    def _open_cice_cgrid_for_F2(self, P_grid=None):
        """
        Open a CICE C-grid file and extract geometry required for F2 computations.

        This routine loads a CICE grid NetCDF and returns:
        - T-cell center longitude/latitude in degrees
        - T-cell local grid rotation angle (anglet) in radians
        - T-cell metric lengths dx, dy in meters

        Variable names are detected robustly across common naming conventions.

        Parameters
        ----------
        P_grid : str or pathlib.Path, optional
            Path to the CICE grid file. If None, uses `self.CICE_dict["P_G"]`.

        Returns
        -------
        tlon_deg : numpy.ndarray
            T-cell longitude in degrees east, shape (nj, ni).
        tlat_deg : numpy.ndarray
            T-cell latitude in degrees north, shape (nj, ni).
        anglet_rad : numpy.ndarray
            Local rotation angle of the T-grid in radians, shape (nj, ni). If the
            source variable appears to be in degrees (max > 2π and <= 360), it is
            converted to radians.
        dx_m : numpy.ndarray
            Grid metric dx in meters, shape (nj, ni).
        dy_m : numpy.ndarray
            Grid metric dy in meters, shape (nj, ni).
        ds : xarray.Dataset
            Opened dataset handle (caller may close if desired).

        Raises
        ------
        ValueError
            If `P_grid` is not provided and `self.CICE_dict["P_G"]` is unset.
        KeyError
            If required variables cannot be found in the dataset.

        Notes
        -----
        Expected for ACCESS-OM3 C-grid (common case):
        - tlon, tlat in radians
        - anglet in radians
        - hte, htn in centimeters

        Longitudes are normalised to [-180, 180] to support Antarctic projection
        transforms and to reduce antimeridian issues.
        """
        if P_grid is None:
            P_grid = self.CICE_dict.get("P_G", None)
        if P_grid is None:
            raise ValueError("P_grid is None and self.CICE_dict['P_G'] is not set.")
        self.logger.info(f"Opening CICE grid: {P_grid}")
        ds = xr.open_dataset(P_grid, decode_times=False)
        # ---- locate vars robustly
        def pick(*names):
            for n in names:
                if n in ds.variables:
                    return n
            return None
        v_tlon = pick("tlon", "TLON", "t_lon", "lon_t")
        v_tlat = pick("tlat", "TLAT", "t_lat", "lat_t")
        v_angt = pick("anglet", "angleT", "ANGLET", "angle_t")
        v_hte  = pick("hte", "HTE", "dxT", "dxt")
        v_htn  = pick("htn", "HTN", "dyT", "dyt")
        for v, nm in [(v_tlon, "tlon"), (v_tlat, "tlat"), (v_angt, "anglet"), (v_hte, "hte"), (v_htn, "htn")]:
            if v is None:
                raise KeyError(f"Could not find required grid variable '{nm}' in {P_grid}")
        tlon = ds[v_tlon].values
        tlat = ds[v_tlat].values
        angt = ds[v_angt].values
        hte  = ds[v_hte].values
        htn  = ds[v_htn].values
        # lon/lat -> degrees, anglet stays radians
        tlon_deg = self._infer_deg_from_grid_units(tlon, "tlon")
        tlat_deg = self._infer_deg_from_grid_units(tlat, "tlat")
        # ensure lon is suitable for Antarctic projection transforms
        tlon_deg = self.normalise_longitudes(tlon_deg, to="-180-180")
        # anglet: if it looks like degrees, convert to rad
        angt_arr = np.asarray(angt, dtype="float64")
        finite = np.isfinite(angt_arr)
        if finite.any():
            amax = np.nanmax(np.abs(angt_arr[finite]))
            # if looks like degrees (e.g. up to ~180), assume degrees
            if amax > (2.0 * np.pi + 1e-6) and amax <= 360.0:
                anglet_rad = np.deg2rad(angt_arr)
            else:
                anglet_rad = angt_arr
        else:
            anglet_rad = angt_arr
        # hte/htn expected in cm -> m
        dx_m = self._to_meters(hte, ds[v_hte].attrs.get("units"), name=v_hte)
        dy_m = self._to_meters(htn, ds[v_htn].attrs.get("units"), name=v_htn)
        return tlon_deg, tlat_deg, anglet_rad, dx_m, dy_m, ds

    def _iter_exterior_rings_lonlat(self, P_shp,
                                    target_crs    : str   = "EPSG:4326",
                                    keep_surfaces : tuple = ("land", "ice shelf", "ice tongue", "rumple"),
                                    dissolve      : bool  = True,
                                    # robustness knobs
                                    fix_invalid   : bool  = True,
                                    union_batch   : float = 2000,
                                    precision_m   : float = None):   # e.g. 1.0 or 5.0 meters; None disables snapping
        """
        Stream exterior polygon rings from a coastline/landmask file as lon/lat arrays.

        The shapefile (or any vector dataset readable by GeoPandas) is read, filtered
        to polygonal geometries, optionally filtered by a `surface` attribute, and
        optionally dissolved (unioned) in the native CRS to remove internal boundaries
        (e.g., grounding-line edges between land and ice shelf polygons). Each returned
        ring is the exterior boundary of a Polygon in the requested CRS.

        Parameters
        ----------
        P_shp : str or pathlib.Path
            Input polygon dataset (e.g., Antarctic coastline/landmask) typically in
            EPSG:3031.
        target_crs : str, default="EPSG:4326"
            CRS for returned coordinates. EPSG:4326 returns lon/lat degrees.
        keep_surfaces : tuple of str, default=("land","ice shelf","ice tongue","rumple")
            If a `surface` column exists, only these classes are retained. If None,
            no filtering is applied.
        dissolve : bool, default=True
            If True, union all retained polygons in the native CRS before reprojecting,
            which removes interior boundaries and yields a cleaner coastline.
        fix_invalid : bool, default=True
            Attempt to repair invalid geometries before unioning (via `make_valid` if
            available, else `buffer(0)`).
        union_batch : float, default=2000
            Batch size for union operations. Large datasets can be unioned in chunks
            to limit memory and improve robustness.
        precision_m : float, optional
            If provided and supported by Shapely, snap coordinates to this precision
            (in meters) before unioning to reduce slivers/spikes. `None` disables
            snapping.

        Yields
        ------
        lon_deg : numpy.ndarray
            Exterior ring longitudes in degrees (float64), normalised to [-180, 180].
        lat_deg : numpy.ndarray
            Exterior ring latitudes in degrees (float64).

        Notes
        -----
        - Non-polygonal geometries (lines/points/empties) are ignored.
        - Dissolving is performed before reprojection to avoid antimeridian and wrap
        artefacts.
        """
        import geopandas as gpd
        from shapely.geometry import Polygon, MultiPolygon
        # shapely 2.x: make_valid + set_precision are available; fall back where needed
        try:
            from shapely.validation import make_valid
        except Exception:
            make_valid = None
        try:
            from shapely import set_precision
        except Exception:
            set_precision = None
        from shapely.ops import unary_union
        # internal functions
        def _polygonal_only(g):
            """Keep polygonal parts only; drop empties/lines/points."""
            if g is None or g.is_empty:
                return None
            # Repair if invalid
            if fix_invalid and hasattr(g, "is_valid") and (not g.is_valid):
                if make_valid is not None:
                    g = make_valid(g)
                else:
                    # classic fallback; fixes many self-intersections
                    g = g.buffer(0)
            if g is None or g.is_empty:
                return None
            # make_valid can return GeometryCollection; keep only polygonal components
            gt = getattr(g, "geom_type", "")
            if gt == "Polygon" or gt == "MultiPolygon":
                return g
            if gt == "GeometryCollection" and hasattr(g, "geoms"):
                polys = []
                for gg in g.geoms:
                    if gg.geom_type in ("Polygon", "MultiPolygon"):
                        polys.append(gg)
                if not polys:
                    return None
                return unary_union(polys)
            return None
        def _apply_precision(g):
            if precision_m is None or set_precision is None:
                return g
            # keep_collapsed avoids creating invalid spikes; drop tiny remnants
            return set_precision(g, float(precision_m), mode="keep_collapsed")
        # ---- read (pyogrio avoids Fiona; you said geopandas works)
        gdf = gpd.read_file(P_shp, engine="pyogrio")
        # ensure CRS (your file is EPSG:3031)
        if gdf.crs is None:
            gdf = gdf.set_crs("EPSG:3031")
        # filter by surface classes if present
        if ("surface" in gdf.columns) and (keep_surfaces is not None):
            gdf = gdf[gdf["surface"].isin(list(keep_surfaces))].copy()
        geoms = [_polygonal_only(g) for g in gdf.geometry.values]
        geoms = [g for g in geoms if (g is not None and not g.is_empty)]
        if len(geoms) == 0:
            return
        # optional snapping (do before union)
        if precision_m is not None and set_precision is not None:
            geoms = [_apply_precision(g) for g in geoms]
            geoms = [g for g in geoms if (g is not None and not g.is_empty)]
        # ---- dissolve in native CRS (EPSG:3031) to remove grounding-line internal edges
        if dissolve:
            parts = []
            for k in range(0, len(geoms), int(union_batch)):
                chunk = geoms[k : k + int(union_batch)]
                if not chunk:
                    continue
                try:
                    parts.append(unary_union(chunk))
                except Exception:
                    # last-chance: buffer(0) on the chunk then union
                    chunk2 = [gg.buffer(0) for gg in chunk]
                    parts.append(unary_union(chunk2))
            # union the already-unioned parts
            try:
                geom_u = unary_union(parts)
            except Exception:
                # last-chance: buffer(0) on parts
                geom_u = unary_union([p.buffer(0) for p in parts])

            geom_list = [geom_u]
        else:
            geom_list = geoms
        # reproject AFTER dissolve (safer + avoids antimeridian / lon wrap weirdness)
        gs = gpd.GeoSeries(geom_list, crs=gdf.crs).to_crs(target_crs)
        for geom in gs.values:
            if geom is None or geom.is_empty:
                continue
            if isinstance(geom, Polygon):
                polys = [geom]
            elif isinstance(geom, MultiPolygon):
                polys = list(geom.geoms)
            else:
                # ignore other types
                continue
            for poly in polys:
                x, y = poly.exterior.coords.xy
                lon = np.asarray(x, dtype="float64")
                lat = np.asarray(y, dtype="float64")
                lon = self.normalise_longitudes(lon, to="-180-180")
                yield lon, lat
     
    def _infer_ocean_mask_from_grid_ds(self, ds):
        """
        Infer a T-grid wet/ocean mask from a grid dataset using common conventions.

        The method searches for typical mask variables (e.g., tmask, wet, kmt, tarea)
        and applies a simple interpretation:
        - kmt-like variables: ocean where kmt > 0
        - mask/wet/tarea-like variables: ocean where finite and > 0

        Parameters
        ----------
        ds : xarray.Dataset
            Grid dataset containing potential mask variables.

        Returns
        -------
        ocean_mask : numpy.ndarray or None
            Boolean array of shape (nj, ni) where True indicates ocean/wet cells.
            Returns None if no suitable variable is found.

        Notes
        -----
        This is heuristic and intended for restricting search domains (e.g., building
        a KDTree over coastal-ocean cells). If your grid uses a different mask
        convention, provide an explicit mask or extend the candidate list.
        """
        # common names across CICE/MOM grids
        candidates = ["tmask", "TMASK", "maskT", "MASKT", "wet", "WET",
                      "kmt", "KMT", "kmt_t", "KMT_T",
                      "tarea", "TAREA"]
        for v in candidates:
            if v in ds.variables:
                a = ds[v].values
                a = np.asarray(a)
                # heuristic:
                # - kmt: 0 land, >0 ocean
                # - masks: 0 land, 1 ocean
                # - tarea: 0 or NaN on land (not always)
                if np.issubdtype(a.dtype, np.integer) or np.issubdtype(a.dtype, np.floating):
                    if v.lower().startswith("kmt"):
                        m = a > 0
                    else:
                        # treat finite positive as ocean
                        m = np.isfinite(a) & (a > 0)
                    if m.ndim == 2 and m.any():
                        return m.astype(bool)
        return None

    def build_F2_form_factors_from_high_res_coast(self,
                                                  P_shp                     : str = None,
                                                  P_grid                    : str = None,
                                                  P_out                     : str = None,
                                                  proj_crs                  : str = "EPSG:3031",
                                                  chunk_segments            : int = 2_000_000,
                                                  coast_write_stride        : int = 25,
                                                  lat_subset_max            : float = -30.0,
                                                  netcdf_compression        : int  = 4,
                                                  use_coastal_ocean_kdtree  : bool = True,
                                                  coast_buffer_cells        : int  = 1,
                                                  max_assign_km             : float = 50.0):
        """
        Compute Liu et al. (2022) coastal form factors (F2x, F2y) on the CICE T-grid
        from a high-resolution polygon coastline.

        The cell-based form factors (Liu et al. 2022, Eqs. 9–10) are computed as:
            F2x(i,j) = Σ_n | l_n cos(θ_n) | / dx(i,j)
            F2y(i,j) = Σ_n | l_n sin(θ_n) | / dy(i,j)

        where each coastline segment n has geodesic length l_n (WGS84) and bearing
        projected into the local model axes using:
            θ_n = α_n - anglet(i,j),
        with α_n the segment angle in Earth-referenced east-north coordinates and
        anglet the local model grid rotation (radians).

        Parameters
        ----------
        P_shp : str, optional
            Path to the high-resolution coastline/land polygons. If None, uses
            `self.CICE_dict["P_high_res_coast"]`.
        P_grid : str, optional
            Path to the CICE grid file. If None, uses `self.CICE_dict["P_G"]`.
        P_out : str, optional
            Output NetCDF path. If None, uses `self.CICE_dict["P_F2_coast"]`.
        proj_crs : str, default="EPSG:3031"
            Projected CRS used to assign each coastline segment midpoint to the
            nearest T-cell via KDTree (meters).
        chunk_segments : int, default=2_000_000
            Process coastline segments in chunks to limit memory use.
        coast_write_stride : int, default=25
            Subsampling stride for storing coastline vertices in the output dataset
            for provenance/plotting. Set <=0 or None to disable.
        lat_subset_max : float, default=-30.0
            Only T-cells with latitude <= this threshold are included in the KDTree
            (a performance and relevance filter for Antarctic use).
        netcdf_compression : int, default=4
            NetCDF zlib compression level for output variables.
        use_coastal_ocean_kdtree : bool, default=True
            If True and an ocean mask can be inferred from the grid file, restrict the
            KDTree to ocean cells within `coast_buffer_cells` of land (to avoid
            mapping segments to interior ice-shelf/grounding-line regions).
        coast_buffer_cells : int, default=1
            Buffer (in grid cells) used to define the coastal-ocean band via binary
            dilation of the land mask.
        max_assign_km : float, default=50.0
            Reject segment-to-cell assignments whose KDTree distance exceeds this
            threshold (in km). Use None to disable.

        Returns
        -------
        ds_out : xarray.Dataset
            Dataset containing:
            - F2x(nj,ni), F2y(nj,ni) : float32, unitless
            - lon(ncoast), lat(ncoast) : thinned coastline vertices (float32)
            with provenance attributes describing inputs and settings.

        Raises
        ------
        ValueError
            If required paths are not provided via arguments or configuration.
        RuntimeError
            If no grid cells satisfy the latitude subset filter (likely units issue).

        Notes
        -----
        - This computes the *cell-based* integrals (Liu f2^u and f2^v). Conversion to
        C-grid u/v points (Liu Eqs. 11–12) is intentionally deferred (e.g., Fortran).
        - Segment lengths and azimuths are computed geodesically on WGS84 (pyproj.Geod).
        - Nearest-cell assignment uses segment midpoints in the projected CRS.
        """
        if P_shp is None:
            P_shp = self.CICE_dict.get("P_high_res_coast", None)
        if P_out is None:
            P_out = self.CICE_dict.get("P_F2_coast", None)
        if P_shp is None:
            raise ValueError("P_shp is None and self.CICE_dict['P_high_res_coast'] is not set.")
        if P_out is None:
            raise ValueError("P_out is None and self.CICE_dict['P_F2_coast'] is not set.")
        from pyproj import CRS, Transformer, Geod
        from scipy.spatial import cKDTree
        from scipy.ndimage import binary_dilation
        tlon_deg, tlat_deg, anglet_rad, dx_m, dy_m, ds_grid = self._open_cice_cgrid_for_F2(P_grid=P_grid)
        nj, ni = tlon_deg.shape
        ncell = nj * ni
        # Base lat subset
        mask = np.isfinite(tlat_deg) & (tlat_deg <= lat_subset_max)
        if not mask.any():
            raise RuntimeError(f"No grid cells found with tlat <= {lat_subset_max}. Check grid/units.")
        # Optional: restrict KDTree to *coastal-ocean* cells to exclude grounding line
        if use_coastal_ocean_kdtree:
            ocean_mask = self._infer_ocean_mask_from_grid_ds(ds_grid)
            if ocean_mask is None:
                self.logger.warning("Could not infer ocean mask from grid dataset; falling back to lat-only KDTree.")
            else:
                land_mask = ~ocean_mask
                coast_ocean = ocean_mask & binary_dilation(land_mask, iterations=int(coast_buffer_cells))
                if coast_ocean.any():
                    mask = mask & coast_ocean
                    self.logger.info(
                        f"KDTree restricted to coastal-ocean cells: buffer_cells={coast_buffer_cells}, "
                        f"kept={mask.sum():,} cells"
                    )
                else:
                    self.logger.warning("Coastal-ocean mask is empty; falling back to lat-only KDTree.")
        flat_idx = np.arange(ncell, dtype="int64")
        sub_flat = flat_idx[mask.ravel()]
        # Project T-cell centers for KDTree
        tfm_to_proj = Transformer.from_crs(CRS.from_epsg(4326), CRS.from_user_input(proj_crs), always_xy=True)
        x_sub, y_sub = tfm_to_proj.transform(tlon_deg.ravel()[sub_flat], tlat_deg.ravel()[sub_flat])
        pts = np.column_stack([np.asarray(x_sub, dtype="float64"), np.asarray(y_sub, dtype="float64")])
        self.logger.info(f"Building KDTree on {pts.shape[0]:,} T-cells (masked).")
        tree = cKDTree(pts)
        # Accumulators
        Sx = np.zeros(ncell, dtype="float64")
        Sy = np.zeros(ncell, dtype="float64")
        # Geodesic for segment length/azimuth
        geod = Geod(ellps="WGS84")
        # For output coastline provenance (thin vertices)
        coast_lon_out = []
        coast_lat_out = []
        self.logger.info(f"Streaming coastline rings from: {P_shp}")
        ring_count = 0
        seg_total = 0
        for lon_ring, lat_ring in self._iter_exterior_rings_lonlat(P_shp):
            ring_count += 1
            # store thinned vertices for output
            if coast_write_stride is not None and coast_write_stride > 0:
                coast_lon_out.append(lon_ring[::coast_write_stride])
                coast_lat_out.append(lat_ring[::coast_write_stride])
                # separator
                coast_lon_out.append(np.asarray([np.nan], dtype="float64"))
                coast_lat_out.append(np.asarray([np.nan], dtype="float64"))
            # segments from consecutive vertices
            if len(lon_ring) < 2:
                continue
            lon0_all = lon_ring[:-1]
            lat0_all = lat_ring[:-1]
            lon1_all = lon_ring[1:]
            lat1_all = lat_ring[1:]
            nseg     = lon0_all.size
            if nseg == 0:
                continue
            # chunk over segments for memory safety
            for s0 in range(0, nseg, int(chunk_segments)):
                s1   = min(nseg, s0 + int(chunk_segments))
                lon0 = lon0_all[s0:s1]
                lat0 = lat0_all[s0:s1]
                lon1 = lon1_all[s0:s1]
                lat1 = lat1_all[s0:s1]
                # geodesic azimuth (deg from north) and length (m)
                az12, _, dist_m = geod.inv(lon0, lat0, lon1, lat1)
                dist_m          = np.asarray(dist_m, dtype="float64")
                # drop zero/NaN segments
                good = np.isfinite(dist_m) & (dist_m > 0.0) & np.isfinite(az12)
                if not good.any():
                    continue
                lon0   = lon0[good]; lat0 = lat0[good]
                lon1   = lon1[good]; lat1 = lat1[good]
                az12   = np.asarray(az12, dtype="float64")[good]
                dist_m = dist_m[good]
                # midpoint in projected CRS for nearest-cell assignment
                x0, y0    = tfm_to_proj.transform(lon0, lat0)
                x1, y1    = tfm_to_proj.transform(lon1, lat1)
                xm        = 0.5 * (np.asarray(x0, dtype="float64") + np.asarray(x1, dtype="float64"))
                ym        = 0.5 * (np.asarray(y0, dtype="float64") + np.asarray(y1, dtype="float64"))
                dist_nn, nn = tree.query(np.column_stack([xm, ym]), workers=-1)
                if ring_count == 1 and s0 == 0:
                    self.logger.info(f"KDTree mapping distances (m): min={np.nanmin(dist_nn):.1f}, "
                                    f"p50={np.nanmedian(dist_nn):.1f}, max={np.nanmax(dist_nn):.1f}")
                good_nn = np.isfinite(dist_nn)
                if (max_assign_km is not None) and np.isfinite(max_assign_km):
                    good_nn &= (dist_nn <= (1000.0 * float(max_assign_km)))
                if not good_nn.any():
                    continue
                # apply good_nn to everything
                dist_m    = dist_m[good_nn]
                az12      = az12[good_nn]
                cell_flat = sub_flat[nn[good_nn]]
                seg_ang_e = np.deg2rad(90.0 - az12)
                ang_local = anglet_rad.ravel()[cell_flat]
                ang_local = np.where(np.isfinite(ang_local), ang_local, 0.0)
                theta = seg_ang_e - ang_local
                sx = np.abs(dist_m * np.cos(theta))
                sy = np.abs(dist_m * np.sin(theta))
                np.add.at(Sx, cell_flat, sx)
                np.add.at(Sy, cell_flat, sy)
                seg_total += int(dist_m.size)
            if ring_count % 50 == 0:
                self.logger.info(f"Processed {ring_count} rings; ~{seg_total:,} segments accumulated so far.")
        self.logger.info(f"Finished coastline pass: rings={ring_count:,}, segments={seg_total:,}")
        # create F2x/y arrays 
        dx           = dx_m.ravel()
        dy           = dy_m.ravel()
        good_dx      = np.isfinite(dx) & (dx > 0)
        good_dy      = np.isfinite(dy) & (dy > 0)
        F2x          = np.zeros(ncell, dtype="float64")
        F2y          = np.zeros(ncell, dtype="float64")
        F2x[good_dx] = Sx[good_dx] / dx[good_dx]
        F2y[good_dy] = Sy[good_dy] / dy[good_dy]
        F2x          = F2x.reshape((nj, ni)).astype("float32")
        F2y          = F2y.reshape((nj, ni)).astype("float32")
        # diagnostic
        # j,i  = np.unravel_index(np.nanargmax(F2x.values), F2x.shape)
        # Lx_m = float(F2x.values[j,i] * dx_m[j,i])
        # Ly_m = float(F2y.values[j,i] * dy_m[j,i])
        # print(j,i, F2x.values[j,i], dx_m[j,i], Lx_m/1000)
        # coastline output arrays
        if coast_lon_out:
            coast_lon = np.concatenate(coast_lon_out).astype("float32")
            coast_lat = np.concatenate(coast_lat_out).astype("float32")
        else:
            coast_lon = np.asarray([], dtype="float32")
            coast_lat = np.asarray([], dtype="float32")
        ds_out = xr.Dataset(data_vars=dict(F2x = (("nj", "ni"), F2x, dict(long_name   = "Liu et al. (2022) F2 form factor, x-projection (cell-based)",
                                                                          units       = "1",
                                                                          description = "sum_n |l_n cos(theta_n)| / dx; theta measured relative to local model x-axis")),
                                           F2y = (("nj", "ni"), F2y, dict(long_name   = "Liu et al. (2022) F2 form factor, y-projection (cell-based)",
                                                                          units       = "1",
                                                                          description = "sum_n |l_n sin(theta_n)| / dy; theta measured relative to local model x-axis")),
                                           lon = (("ncoast",), coast_lon, dict(long_name = "coast longitude",
                                                                               units     = "degrees_east")),
                                           lat = (("ncoast",), coast_lat, dict(long_name = "coast latitude",
                                                                               units     = "degrees_north"))),
                            coords = dict(ni     = np.arange(ni, dtype="int32"),
                                          nj     = np.arange(nj, dtype="int32"),
                                          ncoast = np.arange(coast_lon.size, dtype="int64")),
                            attrs  = dict(title              = "High-resolution coastline-derived F2 form factors for CICE",
                                          references         = "Liu et al. (2022) JGR Oceans, doi:10.1029/2022JC018413 (Eqs 9-10)",
                                          coastline_source   = str(P_shp),
                                          grid_source        = str(P_grid if P_grid is not None else self.CICE_dict.get("P_G")),
                                          proj_crs           = str(proj_crs),
                                          coast_write_stride = int(coast_write_stride),
                                          lat_subset_max     = float(lat_subset_max),
                                          created_by         = "SeaIceToolbox.SeaIceGridWork.build_F2_form_factors_from_high_res_coast"))
        # write netcdf
        self.logger.info(f"Writing F2 NetCDF: {P_out}")
        enc = {"F2x": {"zlib": True, "complevel": int(netcdf_compression), "dtype": "float32"},
               "F2y": {"zlib": True, "complevel": int(netcdf_compression), "dtype": "float32"},
               "lon": {"zlib": True, "complevel": int(netcdf_compression), "dtype": "float32"},
               "lat": {"zlib": True, "complevel": int(netcdf_compression), "dtype": "float32"}}
        ds_out.to_netcdf(P_out, encoding=enc, mode='w')
        return ds_out

    def read_grounded_iceberg_csv(self, P_csv,
                                  lon_col          : str = "lon",
                                  lat_col          : str = "lat",
                                  normalise_lon_to : str = "0-360"):
        """
        Read a grounded-iceberg (GI) point CSV and standardise columns.

        The CSV must contain longitude and latitude columns (degrees). Optional
        columns are preserved when present (e.g., area_m2, C_gi, orient_deg, t_start,
        t_end). Longitudes are normalised to the chosen convention.

        Parameters
        ----------
        P_csv : str or pathlib.Path
            Path to the CSV file.
        lon_col : str, default="lon"
            Name of the longitude column in the CSV.
        lat_col : str, default="lat"
            Name of the latitude column in the CSV.
        normalise_lon_to : {"0-360","-180-180"}, default="0-360"
            Longitude convention applied after reading.

        Returns
        -------
        df : pandas.DataFrame
            DataFrame with at least columns:
            - lon : float (degrees)
            - lat : float (degrees)
            and any recognised optional columns preserved.

        Raises
        ------
        ValueError
            If required columns are missing.

        Notes
        -----
        Rows with missing or non-numeric lon/lat are dropped. Latitudes are bounded
        to [-90, 90] as a basic sanity check.
        """
        df = pd.read_csv(P_csv)
        if lon_col not in df.columns or lat_col not in df.columns:
            raise ValueError(f"CSV must contain columns '{lon_col}' and '{lat_col}'. "
                             f"Found: {list(df.columns)}")
        df = df.rename(columns={lon_col: "lon", lat_col: "lat"}).copy()
        df["lon"] = pd.to_numeric(df["lon"], errors="coerce")
        df["lat"] = pd.to_numeric(df["lat"], errors="coerce")
        df = df.dropna(subset=["lon", "lat"]).reset_index(drop=True)
        # normalise lon to match CICE grid convention
        if normalise_lon_to == "0-360":
            df["lon"] = (df["lon"] % 360.0 + 360.0) % 360.0
        elif normalise_lon_to == "-180-180":
            df["lon"] = ((df["lon"] + 180.0) % 360.0) - 180.0
        # optional columns (not *yet* required)
        for opt in ["area_m2", "C_gi", "orient_deg", "t_start", "t_end"]:
            if opt in df.columns:
                # keep as-is; coerce numeric where appropriate
                if opt in ["area_m2", "C_gi", "orient_deg"]:
                    df[opt] = pd.to_numeric(df[opt], errors="coerce")
        # sanity check ...
        df = df[(df["lat"] >= -90.0) & (df["lat"] <= 90.0)].reset_index(drop=True)
        if getattr(self, "logger", None) is not None:
            self.logger.info(f"Loaded GI CSV: {P_csv} (n={len(df)})")
        return df
    
    def read_grounded_iceberg_gpkg(self,
                                    P_gpkg           : str,
                                    layer            : str = "grounded_icebergs",
                                    dedup_uid        : bool = True,
                                    uid_col          : str = "Global_UID",
                                    x_col            : str = "Easting_3031",
                                    y_col            : str = "Northing_3031",
                                    lon_col          : str = "Longitude",
                                    lat_col          : str = "Latitude",
                                    area_km2_col     : str = "Area_Mean_km2",
                                    use_attr_area    : bool = True,
                                    normalise_lon_to : str = "0-360"):
        """
        Read a grounded-iceberg polygon GeoPackage and derive mapping-ready attributes.

        The GeoPackage is expected to contain grounded iceberg polygons (typically in
        EPSG:3031). This method computes perimeter and area from geometry, optionally
        prefers an attribute-provided mean area (e.g., Area_Mean_km2), and returns a
        tabular representation suitable for mapping to the model grid.

        Parameters
        ----------
        P_gpkg : str or pathlib.Path
            Path to the GeoPackage.
        layer : str, default="grounded_icebergs"
            Layer name to read.
        dedup_uid : bool, default=True
            If True, aggregate multiple detections sharing `uid_col` into one record
            (mean position; median size metrics).
        uid_col : str, default="Global_UID"
            Unique identifier column for deduplication.
        x_col, y_col : str
            Column names containing projected coordinates in EPSG:3031 meters.
        lon_col, lat_col : str
            Optional geographic lon/lat columns (degrees) for metadata/debug.
        area_km2_col : str, default="Area_Mean_km2"
            Attribute column containing mean iceberg area in km^2.
        use_attr_area : bool, default=True
            If True and `area_km2_col` exists, use it (falling back to geometry area
            when missing/invalid). If False, use geometry area only.
        normalise_lon_to : {"0-360","-180-180"}, default="0-360"
            Longitude convention applied to the lon column (if present).

        Returns
        -------
        df : pandas.DataFrame
            Table with columns:
            - uid, x_m, y_m, lon, lat
            - area_m2, perim_m
            - bed_depth (if present), plus additional metadata columns when available

        Notes
        -----
        Rows without finite x/y coordinates are dropped. Perimeter and geometry area
        are computed in the native CRS units (meters for EPSG:3031).
        """
        import geopandas as gpd
        gdf = gpd.read_file(P_gpkg, layer=layer, engine="pyogrio")
        if gdf.crs is None:
            gdf = gdf.set_crs("EPSG:3031")
        # metric perimeter/area from geometry in EPSG:3031
        perim_m      = gdf.geometry.length.astype("float64").to_numpy()
        geom_area_m2 = gdf.geometry.area.astype("float64").to_numpy()
        # preferred area from attribute (mean across detections), fallback to geometry
        if use_attr_area and (area_km2_col in gdf.columns):
            area_m2      = pd.to_numeric(gdf[area_km2_col], errors="coerce").to_numpy(dtype="float64") * 1e6
            bad          = ~np.isfinite(area_m2) | (area_m2 <= 0)
            area_m2[bad] = geom_area_m2[bad]
        else:
            area_m2 = geom_area_m2
        # coordinates for mapping (x/y in metres)
        x = pd.to_numeric(gdf.get(x_col, np.nan), errors="coerce").to_numpy(dtype="float64")
        y = pd.to_numeric(gdf.get(y_col, np.nan), errors="coerce").to_numpy(dtype="float64")
        # lon/lat for metadata (optional)
        lon = pd.to_numeric(gdf.get(lon_col, np.nan), errors="coerce").to_numpy(dtype="float64")
        lat = pd.to_numeric(gdf.get(lat_col, np.nan), errors="coerce").to_numpy(dtype="float64")
        lon = self.normalise_longitudes(lon, to = normalise_lon_to)
        df = pd.DataFrame({"uid"         : gdf[uid_col].values if uid_col in gdf.columns else np.arange(len(gdf)),
                            "x_m"        : x,
                            "y_m"        : y,
                            "lon"        : lon,
                            "lat"        : lat,
                            "area_m2"    : area_m2,
                            "perim_m"    : perim_m,
                            "bed_depth"  : pd.to_numeric(gdf.get("Bed_Depth", np.nan), errors="coerce"),
                            "timestamp"  : gdf.get("Timestamp", None),
                            "date_range" : gdf.get("Date_Range", None),
                            "orbit"      : gdf.get("Orbit", None),
                            "swath_mode" : gdf.get("Swath_Mode", None)})
        # drop unusable rows
        df = df[np.isfinite(df["x_m"]) & np.isfinite(df["y_m"])].reset_index(drop=True)
        # deduplicate by uid (recommended)
        if dedup_uid and "uid" in df.columns:
            # keep one row per uid; use median for sizes, mean for position
            agg = {"x_m"      : "mean",   "y_m"    : "mean",
                   "lon"      : "mean",   "lat"    : "mean",
                   "area_m2"  : "median", "perim_m": "median",
                   "bed_depth": "median"}
            df = df.groupby("uid", as_index=False).agg(agg)
        self.logger.info(f"Loaded GI GPKG: {P_gpkg} layer={layer} n={len(df)} (dedup_uid={dedup_uid})")
        return df
    
    def map_xy_to_tgrid(self, x_m, y_m,
                        proj_crs    : str   = "EPSG:3031", 
                        max_dist_km : float = None):
        """
        Map projected x/y points to the nearest CICE T-cell using a cached KDTree.

        Parameters
        ----------
        x_m, y_m : array-like
            Point coordinates in meters in the CRS specified by `proj_crs`.
        proj_crs : str, default="EPSG:3031"
            CRS of input x/y coordinates. Must match the CRS used to build the KDTree.
        max_dist_km : float, optional
            If provided, points farther than this distance (km) from the nearest
            T-cell center are flagged invalid.

        Returns
        -------
        ii : numpy.ndarray
            Nearest-cell i indices (int64).
        jj : numpy.ndarray
            Nearest-cell j indices (int64).
        dist_m : numpy.ndarray
            Distance (meters) from each point to the nearest T-cell center.
        valid : numpy.ndarray
            Boolean mask indicating finite distances and (if provided) max_dist_km
            compliance.

        Notes
        -----
        Requires `self.G_t` to be loaded; if not present, `load_cice_grid` is called
        with `build_faces=True` to ensure consistent grid metadata.
        """
        if (not hasattr(self, "G_t")) or (self.G_t is None):
            self.load_cice_grid(slice_hem=False, build_faces=True)
        tree, ij = self._get_tgrid_kdtree(proj_crs=proj_crs, force=False)
        q        = np.column_stack([np.asarray(x_m, dtype="float64"),
                                    np.asarray(y_m, dtype="float64")])
        dist_m, idx = tree.query(q, k=1)
        ii          = ij[idx, 0].astype("int64")
        jj          = ij[idx, 1].astype("int64")
        valid       = np.isfinite(dist_m)
        if max_dist_km is not None:
            valid &= (dist_m <= (1000.0 * float(max_dist_km)))
        return ii, jj, dist_m, valid

    def compute_tcell_dxdy_m(self, 
                             force = False,
                             cache = True):
        """
        Compute approximate T-cell metric lengths dx and dy (meters).

        The method estimates T-cell sizes using distances between adjacent face-center
        coordinates derived from grid corner/face metadata (self.G_e, self.G_n). Great-
        circle (ellipsoidal) distances are computed on the WGS84 ellipsoid using
        `pyproj.Geod.inv`.

        Parameters
        ----------
        force : bool, default=False
            If True, recompute dx/dy even if cached values exist.
        cache : bool, default=True
            If True, store results in `self._dxdy_cache` and reuse on subsequent calls.

        Returns
        -------
        dx_m : numpy.ndarray
            Zonal metric length between adjacent vertical-face centers, shape (nj, ni),
            in meters (float64).
        dy_m : numpy.ndarray
            Meridional metric length between adjacent horizontal-face centers, shape
            (nj, ni), in meters (float64).

        Raises
        ------
        ImportError
            If `pyproj` is not available.

        Notes
        -----
        - Distances are returned as absolute values; non-finite or non-positive values
        are set to NaN.
        - This is an approximate metric consistent with face-center spacing; it is
        appropriate for normalising form-factor sums and related diagnostics.
        """
        try:
            from pyproj import Geod
        except Exception as e:
            raise ImportError("pyproj.Geod is required to compute dx/dy. Install pyproj.") from e
        if (not force) and hasattr(self, "_dxdy_cache") and cache:
            return self._dxdy_cache["dx_m"], self._dxdy_cache["dy_m"]
        # ensure grid is loaded with faces
        if (not hasattr(self, "G_t")) or (self.G_t is None) or (not hasattr(self, "G_e")) or (self.G_e is None):
            self.load_cice_grid(slice_hem=False, build_faces=True)
        lon_e = self.G_e["lon"].values  # shape (nj, ni+1)
        lat_e = self.G_e["lat"].values
        lon_n = self.G_n["lon"].values  # shape (nj+1, ni)
        lat_n = self.G_n["lat"].values
        geod = Geod(ellps="WGS84")
        # dx: distance between adjacent vertical-face centers: (j,i) edge -> (j,i+1) edge
        lon1       = self.normalise_longitudes(lon_e[:, :-1], to="-180-180")
        lat1       = lat_e[:, :-1]
        lon2       = self.normalise_longitudes(lon_e[:, 1:], to="-180-180")
        lat2       = lat_e[:, 1:]
        _, _, dx_m = geod.inv(lon1, lat1, lon2, lat2)  # returns meters, shape (nj, ni)
        # dy: distance between adjacent horizontal-face centers: (j,i) edge -> (j+1,i) edge
        lon1       = self.normalise_longitudes(lon_n[:-1, :], to="-180-180")
        lat1       = lat_n[:-1, :]
        lon2       = self.normalise_longitudes(lon_n[1:, :], to="-180-180")
        lat2       = lat_n[1:, :]
        _, _, dy_m = geod.inv(lon1, lat1, lon2, lat2)  # shape (nj, ni)
        # |dx| & |dy|, and convert to float64 
        dx_m = np.abs(dx_m).astype("float64")
        dy_m = np.abs(dy_m).astype("float64")
        # avoid zeros ... statistically unsafe
        dx_m[~np.isfinite(dx_m)] = np.nan
        dy_m[~np.isfinite(dy_m)] = np.nan
        dx_m[dx_m <= 0]          = np.nan
        dy_m[dy_m <= 0]          = np.nan
        if cache:
            self._dxdy_cache = {"dx_m": dx_m, "dy_m": dy_m}
        return dx_m, dy_m

    def _get_tgrid_kdtree(self, 
                          proj_crs = "EPSG:3031",
                          force     = False):
        """
        Build (or retrieve) a cached KDTree of projected T-cell centers.

        The KDTree enables fast nearest-neighbor mapping from projected coordinates
        (x,y) to model T-grid indices. Lon/lat are projected from EPSG:4326 into
        `proj_crs` using pyproj.Transformer.

        Parameters
        ----------
        proj_crs : str, default="EPSG:3031"
            Projected CRS used to represent T-cell centers in meters.
        force : bool, default=False
            If True, rebuild the KDTree even if a cached tree exists.

        Returns
        -------
        tree : scipy.spatial.cKDTree
            KDTree built over projected (x,y) coordinates for all T-cells.
        ij : numpy.ndarray
            Integer index array of shape (N, 2) mapping KDTree point order to
            (i, j) indices: ij[:,0] = i, ij[:,1] = j.

        Notes
        -----
        The cache key includes the CRS and the grid shape. If the grid changes or a
        different CRS is requested, the tree is rebuilt.
        """
        from scipy.spatial import cKDTree
        from pyproj import Transformer
        if (not force) and hasattr(self, "_tgrid_tree_cache"):
            key = self._tgrid_tree_cache.get("key", None)
            if key == (proj_crs, tuple(self.G_t["lat"].shape)):
                return self._tgrid_tree_cache["tree"], self._tgrid_tree_cache["ij"]
        lon = self.G_t["lon"].values
        lat = self.G_t["lat"].values
        # Project lon/lat -> x/y (meters)
        tr     = Transformer.from_crs("EPSG:4326", proj_crs, always_xy=True)
        x, y   = tr.transform(lon.astype("float64"), lat.astype("float64"))
        nj, ni = lon.shape
        ij     = np.stack(np.meshgrid(np.arange(ni), np.arange(nj)), axis=-1).reshape(-1, 2)  # (N,2): [i,j]
        pts    = np.column_stack([x.reshape(-1), y.reshape(-1)])  # (N,2)
        tree   = cKDTree(pts)
        # strore in an underscore ...
        self._tgrid_tree_cache = {"key": (proj_crs, (nj, ni)), "tree": tree, "ij": ij}
        return tree, ij

    def map_points_to_tgrid(self, lon_deg, lat_deg,
                            proj_crs    : str ="EPSG:3031",
                            max_dist_km : float = None):
        """
        Map geographic lon/lat points to the nearest CICE T-cell.

        Points are projected into `proj_crs` and mapped to the nearest projected
        T-cell center using a cached KDTree.

        Parameters
        ----------
        lon_deg, lat_deg : array-like
            Longitudes and latitudes in degrees (EPSG:4326).
        proj_crs : str, default="EPSG:3031"
            Projected CRS used for KDTree search (meters).
        max_dist_km : float, optional
            If provided, points farther than this distance (km) from the nearest
            T-cell center are flagged invalid.

        Returns
        -------
        ii : numpy.ndarray
            Nearest-cell i indices (int64).
        jj : numpy.ndarray
            Nearest-cell j indices (int64).
        dist_m : numpy.ndarray
            Distance (meters) from each point to the nearest T-cell center.
        valid : numpy.ndarray
            Boolean mask indicating finite distances and (if provided) max_dist_km
            compliance.

        Notes
        -----
        Requires `self.G_t` to be loaded; if absent, `load_cice_grid` is called.
        """
        from pyproj import Transformer
        if (not hasattr(self, "G_t")) or (self.G_t is None):
            self.load_cice_grid(slice_hem=False, build_faces=True)
        tree, ij = self._get_tgrid_kdtree(proj_crs=proj_crs, force=False)
        tr       = Transformer.from_crs("EPSG:4326", proj_crs, always_xy=True)
        x, y     = tr.transform(np.asarray(lon_deg, dtype="float64"), np.asarray(lat_deg, dtype="float64"))
        q        = np.column_stack([x, y])
        # compute distances to each t-cell
        dist_m, idx = tree.query(q, k=1)
        ii          = ij[idx, 0].astype("int64")
        jj          = ij[idx, 1].astype("int64")
        # only use valid points
        valid = np.isfinite(dist_m)
        if max_dist_km is not None:
            valid = valid & (dist_m <= (1000.0 * float(max_dist_km)))
        return ii, jj, dist_m, valid

    def _dbscan_like_clusters(self, x, y,
                              eps_m       : float =15000.0,
                              min_samples : int =2 ):
        """
        Cluster points using a minimal DBSCAN-like algorithm (no sklearn dependency).

        Connectivity is defined by an epsilon-neighborhood in Euclidean projected
        space. Points with at least `min_samples` neighbors (including themselves)
        are treated as core points; clusters are grown by breadth-first search from
        core points, and border points are assigned to a cluster if reachable.
        Points not reachable from any core point are labeled as noise (-1).

        Parameters
        ----------
        x, y : array-like
            Point coordinates in meters (projected planar coordinates).
        eps_m : float, default=15000.0
            Neighborhood radius in meters.
        min_samples : int, default=2
            Minimum neighbor count required for a point to be considered a core point.

        Returns
        -------
        labels : numpy.ndarray
            Integer labels of shape (N,). Cluster IDs are 0..K-1; noise points are -1.

        Notes
        -----
        This implementation is intended for modest point counts and avoids external
        dependencies. For very large N, consider sklearn.cluster.DBSCAN for improved
        performance and additional options.
        """
        from scipy.spatial import cKDTree
        from collections import deque
        x    = np.asarray(x, dtype="float64")
        y    = np.asarray(y, dtype="float64")
        n    = x.size
        pts  = np.column_stack([x, y])
        tree = cKDTree(pts)
        # get neighbours from tree and count the number 
        neigh        = tree.query_ball_point(pts, r=float(eps_m))
        neigh_counts = np.array([len(v) for v in neigh], dtype="int64")
        is_core      = neigh_counts >= int(min_samples)
        labels       = -np.ones(n, dtype="int64")
        cid          = 0
        # loop over lons to create a label of good ('1') and bad ('-1') points
        for i in range(n):
            if labels[i] != -1 or (not is_core[i]):
                continue
            # Start a new cluster
            labels[i] = cid
            q = deque([i])
            while q:
                p = q.popleft()
                for j in neigh[p]:
                    if labels[j] == -1:
                        labels[j] = cid
                        if is_core[j]:
                            q.append(j)
            cid += 1
        return labels
    
    def build_F2_GI_from_df(self, df,
                            method                : str = "simple-geometry",
                            length_scale          : str = "perimeter",        # {"perimeter","area"}
                            base_area_m2          : float = None,
                            C_gi                  : float = 1.0,
                            proj_crs              : str = "EPSG:3031",
                            max_map_dist_km       : float = 50.0,
                            eps_cluster_km        : float = 15.0,
                            min_cluster_size      : int   = 3,
                            cluster_buffer_m      : float = None,
                            cluster_amplification : float = 0.0,
                            weight_by             : str   = "perim"):          # for cluster-axis share: {"equal","perim","area"}
        """
        Build grounded-iceberg (GI) contributions to F2 form factors on the CICE T-grid.

        This method maps grounded iceberg locations (and optional size metrics) onto
        the model T-grid and accumulates additional form-factor terms (F2x_gi, F2y_gi)
        intended to represent sub-grid-scale coastal/obstacle form drag associated
        with grounded icebergs.

        Two approaches are supported:

        1) method="simple-geometry"
        Each GI contributes an isotropic projected length scale Lproj normalised by
        local grid metrics:
            F2x_gi += C_gi * (Lproj / dx)
            F2y_gi += C_gi * (Lproj / dy)
        where Lproj is derived from either perimeter (preferred) or area.

        2) method="cluster-axis"
        Nearby GI points are clustered (DBSCAN-like). Each cluster is represented
        by an oriented ellipse/rectangle described by principal axes in projected
        space. Cluster-aligned projected lengths are rotated into the local model
        coordinate frame using the grid angle, and then distributed back to member
        points with optional weighting.

        Parameters
        ----------
        df : pandas.DataFrame
            Input table of GI features. Expected columns:
            - lon, lat (degrees) and/or x_m, y_m (meters in `proj_crs`)
            Optional columns:
            - area_m2 : feature area in m^2
            - perim_m : feature perimeter in m
            - C_gi    : per-feature scaling coefficient
        method : {"simple-geometry","cluster-axis"}, default="simple-geometry"
            Parameterisation for converting GI features into form factors.
        length_scale : {"perimeter","area"}, default="perimeter"
            Length-scale basis for "simple-geometry". If perimeter is unavailable,
            falls back to area.
        base_area_m2 : float, optional
            Default area used when df lacks area_m2. If None, uses median T-cell area.
        C_gi : float, default=1.0
            Default scaling coefficient applied when df lacks per-feature C_gi.
        proj_crs : str, default="EPSG:3031"
            CRS for x/y mapping and clustering operations (meters).
        max_map_dist_km : float, default=50.0
            Maximum allowed distance between a GI point and its nearest T-cell center.
            Features beyond this are dropped.
        eps_cluster_km : float, default=15.0
            Clustering neighborhood radius (km) for "cluster-axis".
        min_cluster_size : int, default=3
            Minimum number of points required to form a cluster in "cluster-axis".
        cluster_buffer_m : float, optional
            Buffer added to cluster axis lengths (meters). If None, defaults to an
            equivalent-radius buffer derived from `base_area_m2`.
        cluster_amplification : float, default=0.0
            Optional amplification factor for cluster lengths:
                amp = 1 + cluster_amplification * log1p(n_cluster)
            (applied to both principal axes).
        weight_by : {"equal","perim","area"}, default="perim"
            Weighting used to distribute cluster-scale form factor back to member
            points/cells. Falls back to equal weights when required metrics are missing.

        Returns
        -------
        ds : xarray.Dataset
            Dataset containing:
            - F2x_gi(nj,ni), F2y_gi(nj,ni) : float32, unitless
            - gi_count(nj,ni) : int32 count of mapped GI features per cell
            plus diagnostic vectors indexed by "ngi" (mapped features):
            - gi_lon, gi_lat, gi_i, gi_j, gi_cluster_id, gi_map_dist_m,
                gi_area_m2, gi_perim_m

        Raises
        ------
        ValueError
            If `method` is unsupported.

        Notes
        -----
        - Grid metrics dx/dy are computed (and cached) via `compute_tcell_dxdy_m`.
        - Mapping is performed by nearest-neighbor search in projected space using
        a cached KDTree of T-cell centers.
        """
        from pyproj import Transformer
        if (not hasattr(self, "G_t")) or (self.G_t is None):
            self.load_cice_grid(slice_hem=False, build_faces=True)
        nat_dim = self.CICE_dict.get("spatial_dims", ("nj", "ni"))
        nj, ni  = self.G_t["lat"].shape
        if base_area_m2 is None:
            base_area_m2 = float(np.nanmedian(self.G_t["area"].values))
        # get coords for mapping -- allows for either lon/lat or x/y points
        has_xy = ("x_m" in df.columns) and ("y_m" in df.columns) and np.isfinite(df["x_m"]).any()
        if has_xy:
            ii, jj, dist_m, valid = self.map_xy_to_tgrid(df["x_m"].values, df["y_m"].values,
                                                        proj_crs=proj_crs, max_dist_km=max_map_dist_km)
        else:
            ii, jj, dist_m, valid = self.map_points_to_tgrid(df["lon"].values, df["lat"].values,
                                                            proj_crs=proj_crs, max_dist_km=max_map_dist_km)
        # update dataframe with valid-only points
        df     = df.loc[valid].reset_index(drop=True)
        ii     = ii[valid]
        jj     = jj[valid]
        dist_m = dist_m[valid]
        # dx/dy on T-grid
        dx_m, dy_m = self.compute_tcell_dxdy_m(force=False, cache=True)
        # flatten indexing for vectorised accumulation
        flat = (jj.astype("int64") * int(ni) + ii.astype("int64"))
        dx   = dx_m.ravel()[flat]
        dy   = dy_m.ravel()[flat]
        # ensure dx/dy are finite and >0
        ok = np.isfinite(dx) & (dx > 0) & np.isfinite(dy) & (dy > 0)
        # areas/perimeters
        area_m2 = df["area_m2"].values.astype("float64") if "area_m2" in df.columns else np.full(len(df), base_area_m2)
        perim_m = df["perim_m"].values.astype("float64") if "perim_m" in df.columns else np.full(len(df), np.nan)
        # C_gi per point if present
        if "C_gi" in df.columns:
            Cvals = np.asarray(df["C_gi"].fillna(C_gi).values, dtype="float64")
        else:
            Cvals = np.full(len(df), float(C_gi), dtype="float64")
        # 'isotropic' (my local definition) projected length scale (for simple-geometry)
        if length_scale == "perimeter" and np.isfinite(perim_m).any():
            Lproj = (2.0 / np.pi) * perim_m
            # fallback to area for any missing perimeters
            bad = ~np.isfinite(Lproj) | (Lproj <= 0)
            if bad.any():
                Lproj[bad] = 4.0 * np.sqrt(np.maximum(area_m2[bad], 0.0) / np.pi)
        else:
            Lproj = 4.0 * np.sqrt(np.maximum(area_m2, 0.0) / np.pi)
        # initialise x/y form factors
        F2x_gi   = np.zeros((nj, ni), dtype="float64")
        F2y_gi   = np.zeros((nj, ni), dtype="float64")
        gi_count = np.zeros((nj, ni), dtype="int32")
        if method == "simple-geometry":
            # vectorised add-at
            fx = np.zeros_like(dx); fy = np.zeros_like(dy)
            fx[ok] = Cvals[ok] * (Lproj[ok] / dx[ok])
            fy[ok] = Cvals[ok] * (Lproj[ok] / dy[ok])
            np.add.at(F2x_gi.ravel(), flat[ok], fx[ok])
            np.add.at(F2y_gi.ravel(), flat[ok], fy[ok])
            np.add.at(gi_count.ravel(), flat[ok], 1)
            labels = np.full(len(df), -1, dtype="int64")
        elif method == "cluster-axis":
            # need projected x/y for clustering; if not present, derive from lon/lat
            if has_xy:
                x_gi = df["x_m"].values.astype("float64")
                y_gi = df["y_m"].values.astype("float64")
            else:
                tr         = Transformer.from_crs("EPSG:4326", proj_crs, always_xy=True)
                x_gi, y_gi = tr.transform(df["lon"].values.astype("float64"), df["lat"].values.astype("float64"))
                x_gi       = np.asarray(x_gi, dtype="float64")
                y_gi       = np.asarray(y_gi, dtype="float64")
            # make sure inputs are inputs into this function
            labels = self._dbscan_like_clusters(x_gi, y_gi,
                                                eps_m       = float(eps_cluster_km) * 1000.0,
                                                min_samples = int(min_cluster_size))
            # local grid angle (rad)
            ang_rad = np.deg2rad(self.G_t["angle"].values.astype("float64"))
            # default buffer
            if cluster_buffer_m is None:
                cluster_buffer_m = float(np.sqrt(base_area_m2 / np.pi))
            uniq = sorted([c for c in np.unique(labels) if c >= 0])
            for cid in uniq:
                idx = np.where(labels == cid)[0]
                npt = idx.size
                if npt < int(min_cluster_size):
                    continue
                # get eigenvecs and normalise for PCA part of computation
                X                = np.column_stack([x_gi[idx], y_gi[idx]])
                Xc               = X - X.mean(axis=0)
                C                = np.cov(Xc.T)
                eigvals, eigvecs = np.linalg.eigh(C)
                v                = eigvecs[:, np.argmax(eigvals)]
                v                = v / np.linalg.norm(v)
                u                = np.array([-v[1], v[0]])
                s1               = Xc @ v
                s2               = Xc @ u
                L                = (np.nanmax(s1) - np.nanmin(s1)) + 2.0 * cluster_buffer_m
                W                = (np.nanmax(s2) - np.nanmin(s2)) + 2.0 * cluster_buffer_m
                # cluster_amplication should never be 0, but ...
                if cluster_amplification and cluster_amplification != 0.0:
                    amp = 1.0 + float(cluster_amplification) * np.log1p(float(npt))
                    L *= amp; W *= amp
                phi = np.arctan2(v[1], v[0])  # rad from +x (east) in projected plane
                # cluster share weights
                if weight_by == "perim" and "perim_m" in df.columns and np.isfinite(perim_m[idx]).all():
                    w = perim_m[idx].astype("float64")
                elif weight_by == "area" and "area_m2" in df.columns:
                    w = area_m2[idx].astype("float64")
                else:
                    w = np.ones(npt, dtype="float64")
                w    = np.maximum(w, 0.0)
                wsum = w.sum()
                if wsum <= 0:
                    w = np.ones(npt, dtype="float64"); wsum = float(npt)
                share = w / wsum
                # create form factors for every grid cell based on underlying grounded berg cluster-isation (!) and axis
                for kk, k in enumerate(idx):
                    i = int(ii[k]); j = int(jj[k])
                    if not (0 <= i < ni and 0 <= j < nj):
                        continue
                    dxk = dx_m[j, i]; dyk = dy_m[j, i]
                    if (not np.isfinite(dxk)) or (not np.isfinite(dyk)) or dxk <= 0 or dyk <= 0:
                        continue
                    theta           = phi - ang_rad[j, i]
                    c               = np.abs(np.cos(theta)); s = np.abs(np.sin(theta))
                    Lx_tot          = 2.0 * L * c + 2.0 * W * s
                    Ly_tot          = 2.0 * L * s + 2.0 * W * c
                    F2x_gi[j, i]   += Cvals[k] * share[kk] * (Lx_tot / dxk)
                    F2y_gi[j, i]   += Cvals[k] * share[kk] * (Ly_tot / dyk)
                    gi_count[j, i] += 1
        else:
            raise ValueError("method must be 'simple-geometry' or 'cluster-axis'")
        ds = xr.Dataset(data_vars = {"F2x_gi"   : (nat_dim, F2x_gi.astype("float32"), {"long_name"   : "Grounded-iceberg additional form factor, x-projection (T-cell)",
                                                                                       "units"       : "1", "method": method, "C_gi_default": float(C_gi),
                                                                                       "length_scale": str(length_scale), "weight_by": str(weight_by)}),
                                     "F2y_gi"   : (nat_dim, F2y_gi.astype("float32"), {"long_name"   : "Grounded-iceberg additional form factor, y-projection (T-cell)",
                                                                                       "units"       : "1", "method": method, "C_gi_default": float(C_gi),
                                                                                       "length_scale": str(length_scale), "weight_by": str(weight_by)}),
                                     "gi_count" : (nat_dim, gi_count.astype("int32"), {"long_name"   : "Number of GI polygons mapped into each T-cell",
                                                                                       "units"       : "count"})},
                        coords    = {nat_dim[0]: np.arange(nj),
                                     nat_dim[1]: np.arange(ni)},
                        attrs     = {"proj_crs": str(proj_crs)})
        # diagnostics
        ds["gi_lon"]        = ("ngi", df["lon"].values.astype("float32") if "lon" in df.columns else np.full(len(df), np.nan, "float32"))
        ds["gi_lat"]        = ("ngi", df["lat"].values.astype("float32") if "lat" in df.columns else np.full(len(df), np.nan, "float32"))
        ds["gi_i"]          = ("ngi", ii.astype("int32"))
        ds["gi_j"]          = ("ngi", jj.astype("int32"))
        ds["gi_cluster_id"] = ("ngi", labels.astype("int32"))
        ds["gi_map_dist_m"] = ("ngi", dist_m.astype("float32"))
        ds["gi_area_m2"]    = ("ngi", df["area_m2"].values.astype("float32") if "area_m2" in df.columns else np.full(len(df), np.nan, "float32"))
        ds["gi_perim_m"]    = ("ngi", df["perim_m"].values.astype("float32") if "perim_m" in df.columns else np.full(len(df), np.nan, "float32"))
        self.logger.info(f"Built GI form factors ({method}): n_gi={len(df)}, nonzero_cells={(gi_count>0).sum()}")
        return ds
    
    # wrapper function for build_F2_GI_from_df
    def build_F2_GI_from_gpkg(self,
                                P_gpkg          : str,
                                layer           : str = "grounded_icebergs",
                                method          : str = "simple-geometry",
                                length_scale    : str = "perimeter",
                                C_gi            : float = 1.0,
                                proj_crs        : str = "EPSG:3031",
                                max_map_dist_km : float = 50.0,
                                dedup_uid       : bool = True, **kwargs):
        """
        Convenience wrapper to build GI form factors directly from a GeoPackage.

        This reads grounded iceberg polygons from a GeoPackage layer using
        `read_grounded_iceberg_gpkg`, then delegates to `build_F2_GI_from_df` for
        mapping and accumulation on the CICE T-grid.

        Parameters
        ----------
        P_gpkg : str or pathlib.Path
            GeoPackage path.
        layer : str, default="grounded_icebergs"
            Layer name containing GI features.
        method : {"simple-geometry","cluster-axis"}, default="simple-geometry"
            GI form-factor parameterisation.
        length_scale : {"perimeter","area"}, default="perimeter"
            Length-scale option forwarded to `build_F2_GI_from_df`.
        C_gi : float, default=1.0
            Default scaling coefficient forwarded to `build_F2_GI_from_df`.
        proj_crs : str, default="EPSG:3031"
            CRS used for mapping/clustering.
        max_map_dist_km : float, default=50.0
            Maximum mapping distance threshold (km).
        dedup_uid : bool, default=True
            Whether to deduplicate GI features by UID prior to mapping.
        **kwargs
            Additional keyword arguments forwarded to `build_F2_GI_from_df` (e.g.,
            clustering controls).

        Returns
        -------
        ds : xarray.Dataset
            Output dataset from `build_F2_GI_from_df`.
        """
        df = self.read_grounded_iceberg_gpkg(P_gpkg,
                                             layer            = layer,
                                             dedup_uid        = dedup_uid,
                                             normalise_lon_to = "0-360")
        return self.build_F2_GI_from_df(df,
                                        method          = method,
                                        length_scale    = length_scale,
                                        C_gi            = C_gi,
                                        proj_crs        = proj_crs,
                                        max_map_dist_km = max_map_dist_km,
                                        **kwargs)
        
    def _nc_attr(self, v):
        # netCDF4 can't store bool attrs
        if isinstance(v, (bool, np.bool_)):
            return int(v)
        # numpy scalars -> python scalars
        if isinstance(v, np.generic):
            return v.item()
        # paths -> str
        if isinstance(v, Path):
            return str(v)
        # None -> omit (or empty string if you prefer)
        if v is None:
            return None
        # dict/list/tuple -> JSON string (netCDF attrs must be scalar-ish)
        if isinstance(v, (dict, list, tuple)):
            return json.dumps(v)
        return v
    
    # wrapper and writer 
    def write_F2_with_GI(self,
                        P_F2_coast            : str   = None,
                        P_GI                  : str   = None,   # NEW: can be .csv or .gpkg
                        P_GI_CSV              : str   = None,   # OPTIONAL: backward compatibility
                        P_out                 : str   = None,
                        method                : str   = "simple-geometry",
                        base_area_m2          : float = None,
                        C_gi                  : float = 1.0,
                        # pass-through clustering controls
                        eps_cluster_km        : float = 15.0,
                        min_cluster_size      : int   = 3,
                        cluster_buffer_m      : float = None,
                        cluster_amplification : float = 0.0,
                        weight_by             : str   = 'perim',
                        # new optional knobs for GPKG workflow
                        length_scale          : str   = "perimeter",   # {"perimeter","area"} used by build_F2_GI_from_df/gpkg
                        dedup_uid             : bool  = True,
                        overwrite             : bool  = False):
        """
        Combine coast-only F2 form factors with grounded-iceberg (GI) contributions
        and write a new NetCDF.

        This method loads an existing coastline-derived F2 file (containing `F2x` and
        `F2y`), computes GI contributions from either:
        - a CSV of GI points (.csv), or
        - a GeoPackage of GI polygons (.gpkg),
        then writes a combined dataset containing:
        - total:     F2x, F2y
        - components F2x_coast, F2y_coast, F2x_gi, F2y_gi
        plus any available coastline vertex vectors and GI diagnostic vectors.

        Parameters
        ----------
        P_F2_coast : str, optional
            Input coast-only F2 NetCDF path. If None, uses `self.CICE_dict["P_F2_coast"]`.
        P_GI : str, optional
            Grounded iceberg input path. Supported extensions are ".csv" and ".gpkg".
            If None, falls back to `P_GI_CSV` (legacy) or `self.GI_dict.get("P_KJ")`.
        P_GI_CSV : str, optional
            Legacy CSV input path (kept for backward compatibility). Ignored if `P_GI`
            is provided.
        P_out : str, optional
            Output NetCDF path. If None, uses `self.CICE_dict["P_F2_GI"]`.
        method : {"simple-geometry","cluster-axis"}, default="simple-geometry"
            GI parameterisation used when constructing GI form factors.
        base_area_m2 : float, optional
            Default area used when GI input lacks feature area. Passed through to GI
            builders where supported.
        C_gi : float, default=1.0
            Default GI scaling coefficient.
        eps_cluster_km : float, default=15.0
            Clustering radius (km) used for "cluster-axis" GI method.
        min_cluster_size : int, default=3
            Minimum cluster size used for "cluster-axis".
        cluster_buffer_m : float, optional
            Buffer added to derived cluster axis lengths (meters).
        cluster_amplification : float, default=0.0
            Optional amplification of cluster axis lengths as a function of cluster size.
        length_scale : {"perimeter","area"}, default="perimeter"
            Length-scale selection for the GI "simple-geometry" method.
        dedup_uid : bool, default=True
            For GeoPackage inputs, whether to deduplicate features by UID prior to mapping.
        overwrite : bool, default=False
            If False and `P_out` exists, raise FileExistsError. If True, overwrite.

        Returns
        -------
        ds_out : xarray.Dataset
            Combined dataset written to disk.

        Raises
        ------
        ValueError
            If required inputs are missing or GI file type is unsupported.
        FileExistsError
            If `overwrite=False` and output path already exists.

        Notes
        -----
        Dim ordering is aligned to the coast file. GI arrays are transposed if needed
        to match the coast file (typically (nj, ni)).
        """
        P_F2_coast = P_F2_coast or self.CICE_dict.get("P_F2_coast", None)
        P_out      = P_out      or self.CICE_dict.get("P_F2_GI", None)
        # Backward compatible resolution of GI path
        if P_GI is None:
            P_GI = P_GI_CSV or self.GI_dict.get("P_KJ", None)
        if P_GI is None:
            raise ValueError("No grounded iceberg input provided. Set P_GI (or legacy P_GI_CSV) or self.GI_dict['P_raw'].")
        if getattr(self, "logger", None) is not None:
            self.logger.info(f"GI input: {P_GI}")
        if (not overwrite) and os.path.exists(P_out):
            raise FileExistsError(f"Output already exists: {P_out} (set overwrite=True to replace)")
        ds_coast = xr.open_dataset(P_F2_coast)
        if "F2x" not in ds_coast.variables or "F2y" not in ds_coast.variables:
            raise ValueError(f"{P_F2_coast} must contain variables 'F2x' and 'F2y'.")
        # -----------------------------
        # Dispatch: CSV vs GPKG
        # -----------------------------
        ext = os.path.splitext(str(P_GI))[1].lower()
        if ext == ".csv":
            ds_gi = self.build_F2_GI_from_csv(P_GI,
                                                method                = method,
                                                base_area_m2          = base_area_m2,
                                                C_gi                  = C_gi,
                                                eps_cluster_km        = eps_cluster_km,
                                                min_cluster_size      = min_cluster_size,
                                                cluster_buffer_m      = cluster_buffer_m,
                                                cluster_amplification = cluster_amplification)

        elif ext == ".gpkg":
            # Prefer the new GPKG path; base_area_m2 can still be passed through via kwargs if you use it inside.
            ds_gi = self.build_F2_GI_from_gpkg(P_GI,
                                                method                = method,
                                                length_scale          = length_scale,
                                                C_gi                  = C_gi,
                                                proj_crs              = "EPSG:3031",
                                                max_map_dist_km       = 50.0,
                                                dedup_uid             = dedup_uid,
                                                # clustering controls forwarded
                                                eps_cluster_km        = eps_cluster_km,
                                                min_cluster_size      = min_cluster_size,
                                                cluster_buffer_m      = cluster_buffer_m,
                                                cluster_amplification = cluster_amplification,
                                                weight_by             = weight_by,
                                                # if build_F2_GI_from_df supports base_area_m2, include it:
                                                base_area_m2          = base_area_m2)
        else:
            raise ValueError(f"Unsupported GI input type '{ext}'. Expected .csv or .gpkg. Path: {P_GI}")
        # Align dims (coast file uses (nj,ni) order)
        F2x_coast = ds_coast["F2x"].astype("float32")
        F2y_coast = ds_coast["F2y"].astype("float32")
        gi_x      = ds_gi["F2x_gi"]
        gi_y      = ds_gi["F2y_gi"]
        if tuple(F2x_coast.dims) != tuple(gi_x.dims):
            gi_x = gi_x.transpose(*F2x_coast.dims)
            gi_y = gi_y.transpose(*F2y_coast.dims)
        F2x_total = (F2x_coast + gi_x).astype("float32")
        F2y_total = (F2y_coast + gi_y).astype("float32")
        # Build output dataset
        ds_out = xr.Dataset()
        ds_out["F2x"] = F2x_total
        ds_out["F2y"] = F2y_total
        ds_out["F2x"].attrs.update({"long_name": "Combined coastal + grounded-iceberg form factor, x-projection (cell-based)",
                                    "units": "1"})
        ds_out["F2y"].attrs.update({"long_name": "Combined coastal + grounded-iceberg form factor, y-projection (cell-based)",
                                    "units": "1",})
        # Components
        ds_out["F2x_coast"] = F2x_coast
        ds_out["F2y_coast"] = F2y_coast
        ds_out["F2x_gi"]    = gi_x.astype("float32")
        ds_out["F2y_gi"]    = gi_y.astype("float32")
        # Carry through coastline lon/lat vectors if present
        for v in ["lon", "lat", "ncoast"]:
            if v in ds_coast.variables:
                ds_out[v] = ds_coast[v]
        # Add GI metadata vectors (debug/traceability)
        for v in ["gi_lon", "gi_lat", "gi_i", "gi_j", "gi_cluster_id", "gi_map_dist_m", "gi_area_m2", "gi_perim_m"]:
            if v in ds_gi.variables:
                ds_out[v] = ds_gi[v]
        ds_out.attrs.update({"P_F2_coast_in"        : str(P_F2_coast),
                            "P_GI_in"               : str(P_GI),
                            "method_gi"             : str(method),
                            "C_gi_default"          : float(C_gi),
                            "base_area_m2"          : float(base_area_m2) if base_area_m2 is not None else np.nan,
                            "length_scale"          : str(length_scale),
                            "dedup_uid"             : bool(dedup_uid),
                            "eps_cluster_km"        : float(eps_cluster_km),
                            "min_cluster_size"      : int(min_cluster_size),
                            "cluster_buffer_m"      : float(cluster_buffer_m) if cluster_buffer_m is not None else np.nan,
                            "cluster_amplification" : float(cluster_amplification)})
        if getattr(self, "logger", None) is not None:
            self.logger.info(f"Writing combined F2 (coast+GI): {P_out}")
        clean = {}
        for k, v in ds_out.attrs.items():
            vv = self._nc_attr(v)
            if vv is not None:
                clean[str(k)] = vv
        ds_out.attrs = clean
        ds_out.to_netcdf(P_out)
        ds_coast.close()
        return ds_out
    
    def define_regular_G(self, grid_res,
                         region            = [0,360,-90,0],
                         spatial_dim_names = ("nj","ni")):
        """
        Create a regular (lat, lon) destination grid as an xarray.Dataset.

        Parameters
        ----------
        grid_res : float
            Grid resolution in degrees.
        region : list[float], default [0, 360, -90, 0]
            [lon_min, lon_max, lat_min, lat_max] in degrees.
        spatial_dim_names : tuple[str, str], default ("nj","ni")
            Output dimension names.

        Returns
        -------
        xr.Dataset
            Dataset containing 2D lon/lat coordinates on a regular grid.

        Notes
        -----
        - Longitudes are generated in the same numeric range as `region`.
        - Intended as a destination grid for xESMF regridding of persistence-style maps.
        """
        lon_min, lon_max, lat_min, lat_max = region
        lon_regular = np.arange(lon_min, lon_max + grid_res, grid_res)
        lat_regular = np.arange(lat_min, lat_max + grid_res, grid_res)
        LON, LAT    = np.meshgrid(lon_regular,lat_regular)
        return xr.Dataset({self.CICE_dict['lon_coord_name'] : (spatial_dim_names, LON),
                           self.CICE_dict['lat_coord_name'] : (spatial_dim_names, LAT)})
