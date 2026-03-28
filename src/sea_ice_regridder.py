from __future__ import annotations
import os
import xarray as xr
import numpy  as np
import pandas as pd
from pathlib       import Path
from scipy.spatial import cKDTree
from scipy         import sparse
from scipy.sparse  import save_npz, load_npz

__all__ = ["SeaIceRegridder"]

class SeaIceRegridder:
    """
    Regridding and geometric utilities for Antarctic sea-ice analysis workflows.

    This class provides two broad groups of functionality:

    1) CICE grid-to-grid remapping (B-grid → T-grid)
       - xESMF-based remapping using persistent weight files (recommended when you need
         a consistent, reproducible mapping operator).
       - a lightweight, Dask-safe 2×2 corner averaging operator as a fast alternative
         when an exact xESMF regrid is unnecessary.

    2) Swath-to-grid remapping and EPSG:3031 helpers
       - projection of lat/lon swaths to Antarctic polar stereographic (EPSG:3031),
       - extent unioning/snap-to-grid utilities for building a common analysis grid,
       - nearest-neighbour resampling of swath data onto an AreaDefinition grid,
       - convenience functions for adding lon/lat coordinates to EPSG:3031 grids,
       - geographic subsetting of curvilinear model grids using lon/lat masks that
         handle dateline/seam crossing.

    Expected external configuration
    ------------------------------
    This class is designed to sit inside the broader AFIM stack. The following
    attributes/methods are expected to exist on `self` (typically injected by kwargs
    or provided by a base class/toolbox manager):

    Attributes
    - logger : logging.Logger
        For progress/debug logging.
    - CICE_dict : dict
        Must contain keys used by the B→T routines, typically:
          * "x_dim", "y_dim"
          * "x_dim_length", "y_dim_length"
          * "bcoord_names" (e.g., ["ULON","ULAT"])
          * "P_reG_u2t_weights" (path to xESMF weight file)
    - G_u, G_t : dict-like or xr.Dataset
        CICE source and target grid definitions used by xESMF.
    - kmt_org : xr.Dataset or dict-like
        Must include "kmt_org" used as the target grid mask for xESMF weight generation.

    Methods
    - load_cice_grid(...)
        Must populate `self.G_u` and `self.G_t` (and typically `self.kmt_org`).
    - define_reG_weights()
        Creates `self.reG` (xESMF regridder) and sets `self.reG_weights_defined`.
    - normalise_longitudes(lon, wrap=...)
        Used by EPSG:3031 and geographic subsetting utilities.

    Notes
    -----
    - xESMF weight generation is usually expensive; this class is structured to
      reuse persistent weight files where possible.
    - Longitudes are explicitly wrapped before building swath definitions and before
      geographic masking to avoid seam artefacts.
    """
    
    def __init__(self, **kwargs):
        """
        Initialise the regridding helper.

        Parameters
        ----------
        **kwargs
            Optional configuration injected into `self` by the surrounding workflow.
            Common keys include `logger`, `CICE_dict`, and file paths.

        Notes
        -----
        - The shown implementation is a no-op; in production you typically assign kwargs
          onto `self` (as you do in your other classes) or inherit from a shared base class.
        """
        return

    def _ensure_reG_defined(self):
        """
        Ensure that an xESMF regridder is available on `self`.

        This method checks for an existing `self.reG` operator and, if missing,
        calls `define_reG_weights()` to create (or reuse) xESMF weights and instantiate
        the regridder.

        Returns
        -------
        None

        Notes
        -----
        - Intended as a lightweight guard at the top of functions that depend on `self.reG`.
        - Weight reuse behaviour is controlled inside `define_reG_weights()`.
        """
        if getattr(self, "reG", None) is not None and getattr(self, "reG", None) is not None:
            return
        self.define_reG_weights()

    def define_reG_weights(self):
        """
        Define and store an xESMF regridder to remap CICE B-grid (U-point) data to the T-grid.

        This method constructs/loads the source (U-point) and target (T-point) CICE grids,
        attaches a land/sea mask to the target grid, and then instantiates an xESMF
        `Regridder` object. If a weight file already exists, it is reused; otherwise, new
        weights are generated and written to disk.

        Regridding configuration
        ------------------------
        - Source grid : `self.G_u` (B-grid / U-point style grid definition)
        - Target grid : `self.G_t` (T-grid / cell-center grid definition)
        - Target mask : `self.kmt_org["kmt_org"]` stored as `G_t["mask"]`
        - Method      : "bilinear"
        - periodic    : True (assumes global periodicity in longitude)
        - ignore_degenerate : True (skip degenerate cells)
        - extrap_method : "nearest_s2d" (nearest-neighbour extrapolation source→dest)
        - Weight file  : `self.CICE_dict["P_reG_u2t_weights"]`
        - reuse_weights: True if the file exists, else False (weights created)

        Side effects
        ------------
        - Sets `self.reG` to the instantiated xESMF regridder.
        - Sets `self.reG_weights_defined = True` upon success.

        Returns
        -------
        None

        Raises
        ------
        KeyError
            If required keys are missing from `self.CICE_dict` (e.g., weight path).
        AttributeError
            If `load_cice_grid()` does not populate `self.G_u`, `self.G_t`, or `self.kmt_org`.
        ImportError
            If `xesmf` is not available in the runtime environment.
        Exception
            Propagates xESMF errors arising from grid definitions or weight I/O.

        Notes
        -----
        - The grid definitions must include valid "lon"/"lat" (and ideally corner) information
          for the selected method.
        - For conservative methods you would typically require corner arrays; for bilinear,
          centers are sufficient, but corner metadata may still improve robustness.
        """
        import xesmf as xe
        self.load_cice_grid()
        G_u           = self.G_u
        G_t           = self.G_t
        G_t['mask']   = self.kmt_org['kmt_org']
        F_weights     = self.CICE_dict["P_reG_u2t_weights"]
        weights_exist = os.path.exists(F_weights)
        self.logger.info(f"{'Reusing' if weights_exist else 'Creating'} regrid weights: {F_weights}")
        self.reG = xe.Regridder(G_u, G_t,
                                method            = "bilinear",
                                periodic          = True,
                                ignore_degenerate = True,
                                extrap_method     = "nearest_s2d",
                                reuse_weights     = weights_exist,
                                filename          = F_weights)
        self.reG_weights_defined = True

    def reG_bgrid_to_tgrid_xesmf(self, da, coord_names=None):
        """
        Regrid a B-grid DataArray to the T-grid using the pre-defined xESMF regridder.

        The input DataArray must carry longitude/latitude coordinate variables. If
        `coord_names` is not provided, this method uses the configured coordinate names
        in `self.CICE_dict["bcoord_names"]` (commonly ["ULON","ULAT"]).

        Parameters
        ----------
        da : xr.DataArray
            B-grid variable to regrid. Must include coordinate variables corresponding
            to longitude and latitude (either 1D or 2D, depending on your grid definition).
        coord_names : list[str], optional
            Names of the coordinates on `da` that represent lon/lat. If omitted, defaults
            to `self.CICE_dict["bcoord_names"]`.

        Returns
        -------
        xr.DataArray or None
            Regridded DataArray on the T-grid. Returns None if required coordinates are
            missing or if regridding fails.

        Raises
        ------
        RuntimeError
            If `self.reG` is not defined. Call `_ensure_reG_defined()` or `define_reG_weights()`
            first (or ensure your workflow does so).

        Notes
        -----
        - The method renames the provided coordinate variables to the xESMF-conventional
          names "lon" and "lat" before applying `self.reG(...)`.
        - Error handling is intentionally conservative: failures are logged and `None`
          is returned to allow calling workflows to skip problematic variables gracefully.
        """
        coord_names = coord_names if coord_names is not None else self.CICE_dict["bcoord_names"]
        if not set(coord_names).issubset(set(da.coords)):
            self.logger.error(f"Cannot regrid: as {coord_names} not found in coordinates.")
            return None
        coord_map = {}
        for name in coord_names:
            if "LAT" in name.upper():
                coord_map[name] = "lat"
            elif "LON" in name.upper():
                coord_map[name] = "lon"
        if set(coord_map.values()) != {"lat", "lon"}:
            self.logger.error(f"Could not identify lat/lon from coord_names: {coord_names}")
            return None
        da_tmp = da.rename(coord_map)
        try:
            da_reG = self.reG(da_tmp)
        except Exception as e:
            self.logger.error(f"Regridding failed: {e}")
            return None
        return da_reG
    
    def define_reG_regular_weights(self, da,
                                    G_res             = 0.15,
                                    region            = [0,360,-90,0],
                                    AF2020            = False,
                                    variable_name     = None,
                                    lon_coord_name    = None,
                                    lat_coord_name    = None,
                                    time_coord_name   = None,
                                    spatial_dim_names = None,
                                    reG_method        = "bilinear",
                                    periodic          = False,
                                    reuse_weights     = True,
                                    P_weights         = None):
        """
        Define an xESMF regridder from a curvilinear source grid to a regular lat/lon grid.

        Parameters
        ----------
        da : xr.Dataset or xr.DataArray
            Source object from which lat/lon coordinates are extracted.
        G_res : float, default 0.15
            Regular grid resolution (degrees).
        region : list[float], default [0,360,-90,0]
            Regular grid bounding box [lon_min, lon_max, lat_min, lat_max].
        AF2020 : bool, default False
            If True, use AF2020 naming and default weights path; else CICE/model.
        reG_method : str, default "bilinear"
            xESMF method (e.g., "bilinear", "nearest_s2d", "conservative").
        periodic : bool, default False
            Set True for global periodic longitude grids.
        reuse_weights : bool, default True
            If True, reuse an existing weights file if present.
        P_weights : str or Path, optional
            Override weights filename. Defaults to AF2020/CICE configured weights path.

        Returns
        -------
        xesmf.Regridder
            Configured regridder from source to regular grid.

        Notes
        -----
        - Source/target grids are constructed from 2D lat/lon arrays.
        - Conservative regridding requires grid corner information; this helper builds
          only center coordinates unless you extend it.
        """
        import xesmf as xe
        coords  = self.define_fast_ice_coordinates(da,
                                                    AF2020            = AF2020,
                                                    variable_name     = variable_name,
                                                    lon_coord_name    = lon_coord_name,
                                                    lat_coord_name    = lat_coord_name,
                                                    time_coord_name   = time_coord_name,
                                                    spatial_dim_names = spatial_dim_names)
        if AF2020: 
            P_weights = P_weights if P_weights is not None else self.AF_FI_dict["P_reG_reg_weights"]
        else:
            P_weights = P_weights if P_weights is not None else self.CICE_dict["P_reG_reg_weights"]
        G_src = xr.Dataset({self.CICE_dict['lat_coord_name'] : (coords['names']['spatial_dims'], coords['latitudes']),
                            self.CICE_dict['lon_coord_name'] : (coords['names']['spatial_dims'], coords['longitudes'])})
        G_dst = self.define_regular_G(G_res, region=region, spatial_dim_names=coords['names']['spatial_dims'])
        return xe.Regridder(G_src, G_dst,
                            method            = reG_method,
                            periodic          = periodic,
                            ignore_degenerate = True,
                            reuse_weights     = reuse_weights,
                            filename          = P_weights)

    def simple_spatial_averaging_bgrid_to_tgrid(self, var):
        """
        Compute a Dask-safe 2×2 unweighted average from B-grid corner points to T-grid centers.

        This provides a lightweight alternative to xESMF remapping when a simple local
        averaging is sufficient. The operation is performed by slicing four corner-shifted
        views (v00, v01, v10, v11), averaging them, padding back to the original grid size,
        and applying a cyclic wrap to the last x column (last equals first).

        Parameters
        ----------
        var : xr.DataArray
            Input array on the B-grid. May be 2D (nj, ni) or 3D (time, nj, ni) depending on
            your configuration. Dimension names are taken from `self.CICE_dict["y_dim"]` and
            `self.CICE_dict["x_dim"]`.

        Returns
        -------
        xr.DataArray
            Averaged field on the T-grid with the same nominal shape as the configured
            target sizes (`y_dim_length`, `x_dim_length`). NaNs are used where padding is
            required.

        Raises
        ------
        KeyError
            If required dimension names/lengths are missing from `self.CICE_dict`.
        AssertionError
            If the final array shape does not match the configured target sizes.

        Notes
        -----
        - This is not a conservative remap; it is a local arithmetic average.
        - The last x column is set equal to the first to maintain periodicity.
        - Indexes for spatial dims are dropped to avoid downstream surprises when mixing
          integer and coordinate-based selection.
        """
        x_dim = self.CICE_dict["x_dim"]
        y_dim = self.CICE_dict["y_dim"]
        x_len = self.CICE_dict["x_dim_length"]
        y_len = self.CICE_dict["y_dim_length"]
        self.logger.info(f"input shape to spatial averaging: {var.shape}")
        self.logger.info("  → Slicing corner points for averaging...")
        v00 = var.isel({y_dim: slice(None, -1), x_dim: slice(None, -1)})
        v01 = var.isel({y_dim: slice(None, -1), x_dim: slice(1, None)})
        v10 = var.isel({y_dim: slice(1, None), x_dim: slice(None, -1)})
        v11 = var.isel({y_dim: slice(1, None), x_dim: slice(1, None)})
        self.logger.info("  → Computing mean of four corners...")
        avg = (v00 + v01 + v10 + v11) / 4.0
        self.logger.info("  → Padding with NaNs to restore original grid size...")
        pad_y = max(y_len - avg.sizes.get(y_dim, 0), 0)
        pad_x = max(x_len - avg.sizes.get(x_dim, 0), 0)
        avg = avg.pad({y_dim: (0, pad_y), x_dim: (0, pad_x)}, constant_values=np.nan)
        self.logger.info("  → Applying cyclic wrap for last column...")
        if avg.sizes.get(x_dim, 0) > 1:
            avg[{x_dim: -1}] = avg.isel({x_dim: 0})
        # Force re-slicing to expected grid size to ensure consistency
        avg = avg.isel({y_dim: slice(0, y_len), x_dim: slice(0, x_len)})
        if "time" in var.coords:
            avg = avg.assign_coords(time=var["time"])
            self.logger.info("  → Time coordinate restored.")
        assert avg.sizes[y_dim] == y_len, f"{y_dim} mismatch: got {avg.sizes[y_dim]}, expected {y_len}"
        assert avg.sizes[x_dim] == x_len, f"{x_dim} mismatch: got {avg.sizes[x_dim]}, expected {x_len}"
        for dim in [y_dim, x_dim]:
            if dim in avg.indexes:
                avg = avg.drop_indexes(dim)
        return avg

    def pygmt_regrid(self, da, lon, lat, grid_res=None, region=None, search_radius="200k"):
        """
        Regrid a 2D data array using PyGMT's nearneighbor interpolation.

        This method applies PyGMT's `nearneighbor` algorithm to interpolate scattered
        data values (`da`) onto a regular grid based on specified longitude and latitude 
        arrays. The input is masked to ignore NaNs or non-finite values.

        Parameters
        ----------
        da : xarray.DataArray
            2D array of data values to interpolate (e.g., sea ice thickness).
        lon : xarray.DataArray or np.ndarray
            Longitude values corresponding to `da`, same shape.
        lat : xarray.DataArray or np.ndarray
            Latitude values corresponding to `da`, same shape.
        grid_res : str or float, optional
            Grid spacing for the output grid (e.g., "0.5", "10m"). Required by PyGMT.
        region : list or tuple, optional
            Bounding box for the output grid in the form [west, east, south, north].
        search_radius : str or float, default "200k"
            Search radius for PyGMT's nearneighbor (e.g., "100k" for 100 km).

        Returns
        -------
        gridded : xarray.DataArray
            Gridded output with interpolated values over the defined region.

        Notes
        -----
        - All non-finite values in `da` are excluded prior to interpolation.
        - PyGMT must be properly installed and configured with GMT for this to work.
        """
        import pygmt
        mask = np.isfinite(da)
        df   = pd.DataFrame({"longitude": lon.values[mask].ravel(), 
                             "latitude" : lat.values[mask].ravel(),
                             "z"        : da.values[mask].ravel()})
        return pygmt.nearneighbor(data          = df,
                                  spacing       = grid_res,
                                  region        = region,
                                  search_radius = search_radius)

    ##############################################################
    #                         PYRESAMPLE                         #
    ##############################################################

    def to_3031_extent(self, lat2d, lon2d, buffer_m=20_000):
        """
        Project a swath's latitude/longitude coordinates to EPSG:3031 and return an extent.

        This helper converts 2D lat/lon arrays to Antarctic polar stereographic coordinates
        (EPSG:3031), then returns an axis-aligned bounding box:
            [xmin, ymin, xmax, ymax]
        expanded by an optional buffer (meters).

        Parameters
        ----------
        lat2d, lon2d : array-like
            2D arrays of latitude and longitude (degrees). Must be broadcast-compatible.
        buffer_m : float, default 20000
            Buffer added to each side of the extent (meters).

        Returns
        -------
        list[float]
            Extent [xmin, ymin, xmax, ymax] in EPSG:3031 meters (including buffer).

        Notes
        -----
        - Longitudes are wrapped to [-180, 180) prior to projection to avoid dateline issues.
        - Non-finite projected points are excluded when computing min/max.
        """
        from pyproj import Transformer
        lon2d = self.normalise_longitudes(lon2d, to="-180-180")
        transformer = Transformer.from_crs("EPSG:4326", "EPSG:3031", always_xy=True)
        x, y = transformer.transform(lon2d.ravel(), lat2d.ravel())
        x = np.asarray(x); y = np.asarray(y)
        finite = np.isfinite(x) & np.isfinite(y)
        xmin, xmax = x[finite].min(), x[finite].max()
        ymin, ymax = y[finite].min(), y[finite].max()
        return [xmin - buffer_m, ymin - buffer_m, xmax + buffer_m, ymax + buffer_m]

    def union_extents(self, extents):
        """
        Compute the union (bounding) extent of multiple EPSG:3031 extents.

        Parameters
        ----------
        extents : list of extents
            Each extent must be [xmin, ymin, xmax, ymax] in meters.

        Returns
        -------
        list[float]
            Union extent [xmin, ymin, xmax, ymax] spanning all input extents.
        """
        xs = [e[0] for e in extents] + [e[2] for e in extents]
        ys = [e[1] for e in extents] + [e[3] for e in extents]
        return [min(xs), min(ys), max(xs), max(ys)]

    def snap_extent_to_grid(self, extent, pixel_size):
        """
        Snap an extent to a regular grid defined by `pixel_size`.

        The snapped extent is expanded (never shrunk) such that:
          - xmin, ymin are floored to the nearest multiple of pixel_size
          - xmax, ymax are ceiled  to the nearest multiple of pixel_size

        Parameters
        ----------
        extent : list[float]
            [xmin, ymin, xmax, ymax] in meters.
        pixel_size : float
            Grid resolution (meters).

        Returns
        -------
        list[float]
            Snapped extent [xmin, ymin, xmax, ymax] in meters.
        """
        xmin, ymin, xmax, ymax = extent
        xmin = np.floor(xmin / pixel_size) * pixel_size
        ymin = np.floor(ymin / pixel_size) * pixel_size
        xmax = np.ceil (xmax / pixel_size) * pixel_size
        ymax = np.ceil (ymax / pixel_size) * pixel_size
        return [xmin, ymin, xmax, ymax]

    def make_area_definition(self, extent, pixel_size=5_000, area_id="epsg3031_5km_union"):
        """
        Create a PyResample AreaDefinition for a regular EPSG:3031 grid.

        Parameters
        ----------
        extent : list[float]
            [xmin, ymin, xmax, ymax] in EPSG:3031 meters. Typically produced by
            `to_3031_extent()`, `union_extents()`, and `snap_extent_to_grid()`.
        pixel_size : float, default 5000
            Grid resolution (meters).
        area_id : str, default "epsg3031_5km_union"
            Identifier string for the AreaDefinition.

        Returns
        -------
        pyresample.geometry.AreaDefinition
            Regular projected grid definition suitable for PyResample resampling.

        Notes
        -----
        - Width/height are computed from the extent and rounded to integer pixel counts.
        - The returned extent is adjusted so that xmax/ymax align exactly with width/height.
        """
        from pyresample.geometry import AreaDefinition
        xmin, ymin, xmax, ymax = extent
        width  = int(round((xmax - xmin) / pixel_size))
        height = int(round((ymax - ymin) / pixel_size))
        xmax = xmin + width  * pixel_size
        ymax = ymin + height * pixel_size
        area_def = AreaDefinition(
            area_id=area_id,
            description="Common 5 km EPSG:3031 grid (union of inputs)",
            proj_id="epsg3031",
            projection="EPSG:3031",
            width=width,
            height=height,
            area_extent=(xmin, ymin, xmax, ymax),
        )
        return area_def

    def grid_coords_from_area(self, area_def, pixel_size=5_000):
        """
        Construct 1D x/y coordinate arrays (cell centers) from a PyResample AreaDefinition.

        Parameters
        ----------
        area_def : pyresample.geometry.AreaDefinition
            Target grid definition in EPSG:3031.
        pixel_size : float, default 5000
            Grid resolution (meters).

        Returns
        -------
        (x, y) : tuple[np.ndarray, np.ndarray]
            1D arrays of cell-center coordinates (meters). `x` increases eastward.
            `y` decreases from top to bottom (north → south) consistent with the
            area extent definition.
        """
        xmin, ymin, xmax, ymax = area_def.area_extent
        width, height = area_def.width, area_def.height
        # Cell centers
        x = xmin + (np.arange(width) + 0.5) * pixel_size
        y = ymax - (np.arange(height) + 0.5) * pixel_size  # top->down (north->south)
        return x, y

    def resample_swath_to_area(self, src_da, lat2d, lon2d, area_def, radius=10_000, fill_value=np.nan, pixel_size=5_000):
        """
        Nearest-neighbour resample a 2D swath to an EPSG:3031 AreaDefinition grid.

        Parameters
        ----------
        src_da : xr.DataArray
            2D swath data array (e.g., satellite swath variable). Must be aligned with
            `lat2d`/`lon2d` in shape.
        lat2d, lon2d : array-like
            2D latitude/longitude arrays (degrees) describing the swath geolocation.
        area_def : pyresample.geometry.AreaDefinition
            Target grid definition (EPSG:3031).
        radius : float, default 10000
            Radius of influence in meters used for nearest-neighbour resampling.
        fill_value : scalar, default np.nan
            Fill value assigned where no source points fall within `radius`.
        pixel_size : float, default 5000
            Target grid resolution used to construct the output x/y coordinates.

        Returns
        -------
        xr.DataArray
            Resampled 2D field on the target grid with dims ("y","x") and coordinates
            "x" and "y" in meters, plus metadata indicating EPSG:3031.

        Notes
        -----
        - Longitudes are wrapped to [-180, 180) before constructing the SwathDefinition.
          This is critical to avoid dateline discontinuities in the KD-tree search.
        - `nprocs=0` uses serial execution; set >0 to parallelise if appropriate.
        """
        from pyresample.geometry import SwathDefinition
        from pyresample.kd_tree  import resample_nearest
        lon2d = self.normalise_longitudes(lon2d, to="-180-180")  # << key fix: wrap before building the swath
        swath = SwathDefinition(lons=lon2d, lats=lat2d)
        out2d = resample_nearest(source_geo_def     = swath,
                                 data               = src_da.values,
                                 target_geo_def     = area_def,
                                 radius_of_influence= radius,
                                 fill_value         = fill_value,
                                 nprocs             = 0,            # set >0 to parallelise
                                 reduce_data        = True)
        x, y = self.grid_coords_from_area(area_def, pixel_size=pixel_size)
        da_out = xr.DataArray(out2d,
                              dims   = ("y", "x"),
                              coords = {"x": ("x", x, {"units": "m", "standard_name": "projection_x_coordinate"}),
                                        "y": ("y", y, {"units": "m", "standard_name": "projection_y_coordinate"})},
                              name   = src_da.name,
                              attrs  = {"crs": "EPSG:3031", "grid_mapping": "spstereo", "res": float(pixel_size), **src_da.attrs})
        return da_out

    def add_lonlat_from_epsg3031(self, ds, 
                                x_name    = "x",
                                y_name    = "y",
                                wrap      = "0..360",       # or "-180..180"
                                out_dtype = "float32"):
        """
        Add 2D lon/lat coordinate fields to a dataset defined on an EPSG:3031 grid.

        This function broadcasts the 1D projected x/y coordinates to 2D, converts them
        to lon/lat using EPSG:3031 → EPSG:4326, wraps longitudes to the requested
        convention, and attaches the results as dataset coordinates.

        Parameters
        ----------
        ds : xr.Dataset
            Dataset with 1D projected coordinate dimensions `x_name` and `y_name`.
        x_name, y_name : str, default ("x","y")
            Names of the projected coordinate dimensions in `ds`.
        wrap : {"0..360","-180..180"}, default "0..360"
            Longitude wrapping convention for the output coordinate.
        out_dtype : str, default "float32"
            Output dtype for lon/lat coordinates. Use None to keep float64.

        Returns
        -------
        xr.Dataset
            Dataset with added coordinates:
              - lon(y, x)
              - lat(y, x)

        Raises
        ------
        ValueError
            If `x_name` or `y_name` are not present as dataset dimensions.

        Notes
        -----
        - Uses `xr.apply_ufunc(..., dask="parallelized")` so it can operate on Dask-backed
          x/y coordinates without materialising large intermediate arrays.
        """
        if x_name not in ds.dims or y_name not in ds.dims:
            raise ValueError(f"Expected dims '{y_name}', '{x_name}' in dataset.")
        # broadcast 1-D x/y -> 2-D (y,x)
        X2D, Y2D = xr.broadcast(ds[x_name], ds[y_name])  # shapes (y,x)
        lon, lat = xr.apply_ufunc(self._xy_to_lonlat, X2D, Y2D,
                                input_core_dims=[[y_name, x_name], [y_name, x_name]],
                                output_core_dims=[[y_name, x_name], [y_name, x_name]],
                                dask="parallelized",
                                vectorize=False,
                                output_dtypes=[np.float64, np.float64],)
        # wrap longitudes & cast 
        if wrap == "0..360":
            lon = lon % 360.0
        else:
            lon = ((lon + 180.0) % 360.0) - 180.0
        if out_dtype:
            lon = lon.astype(out_dtype)
            lat = lat.astype(out_dtype)
        # attach as coordinates (on same (y,x) dims)
        return ds.assign_coords(lon=lon, lat=lat)

    def subset_by_lonlat_box(self, da: xr.DataArray, lon_range, lat_range,
                            lon_name="TLON",
                            lat_name="TLAT",
                            jdim="nj",
                            idim="ni",
                            wrap="-180-180",
                            crop=True):
        """
        Subset a curvilinear grid DataArray by a geographic lon/lat bounding box.

        This method constructs a boolean mask from 2D lon/lat coordinates and applies it
        to `da` using `.where(mask)`. It supports bounding boxes that cross the longitude
        seam/dateline by interpreting a range where lon_min > lon_max as a wrapped interval.

        Parameters
        ----------
        da : xr.DataArray
            DataArray on a curvilinear grid, typically with dims (time, nj, ni) or (nj, ni).
            Must include 2D coordinate fields `lon_name` and `lat_name`.
        lon_range : tuple[float, float]
            (lon_min, lon_max) in degrees. If the range crosses the seam, use a wrapped
            specification (e.g., (350, 20) for 0..360 convention).
        lat_range : tuple[float, float]
            (lat_min, lat_max) in degrees.
        lon_name, lat_name : str, default ("TLON","TLAT")
            Names of 2D longitude/latitude coordinate variables in `da.coords`.
        jdim, idim : str, default ("nj","ni")
            Names of the two spatial dimensions.
        wrap : {"0-360","-180-180"}, default "-180-180"
            Longitude wrap convention to apply prior to masking.
        crop : bool, default True
            If True, crops the output to the minimal bounding index box that contains
            at least one valid cell (reduces storage/plotting cost). Cropping is done by
            finding any-True rows/columns in the 2D mask.

        Returns
        -------
        xr.DataArray
            Subsetted DataArray. Values outside the lon/lat box are masked (NaN).
            If `crop=True`, spatial dims are also index-cropped to a tight bounding box.

        Notes
        -----
        - The masking step (`da.where(mask)`) is Dask-friendly.
        - The cropping step currently uses `.values` on 1D any-masks; for very large grids
          this may trigger compute. In typical Antarctic subsetting use-cases, this is
          acceptable because the mask reduction is cheap relative to the full field.
        - Requires `self.normalise_longitudes(...)` to exist and to accept the `wrap` argument.
        """
        # coords (2-D)
        TLON = self.normalise_longitudes(da.coords[lon_name], wrap)
        TLAT = da.coords[lat_name]
        lon_min, lon_max = lon_range
        lat_min, lat_max = lat_range
        # Longitude mask (handle seam)
        if wrap == "0-360":
            lon_min, lon_max = lon_min % 360, lon_max % 360
        else:
            lon_min = ((lon_min + 180) % 360) - 180
            lon_max = ((lon_max + 180) % 360) - 180
        if lon_min <= lon_max:
            mask_lon = (TLON >= lon_min) & (TLON <= lon_max)
        else:
            # crosses dateline / wrap seam
            mask_lon = (TLON >= lon_min) | (TLON <= lon_max)
        mask_lat = (TLAT >= lat_min) & (TLAT <= lat_max)
        mask     = mask_lon & mask_lat                      # (nj, ni)
        out = da.where(mask)                            # broadcast over time
        if crop:
            j_any = mask.any(dim=idim)
            i_any = mask.any(dim=jdim)
            j_idx = np.where(j_any.values)[0]
            i_idx = np.where(i_any.values)[0]
            if j_idx.size and i_idx.size:
                j0, j1 = j_idx.min(), j_idx.max() + 1
                i0, i1 = i_idx.min(), i_idx.max() + 1
                out    = out.isel({jdim: slice(j0, j1), idim: slice(i0, i1)})
        return out

    #------------------------------------------------------------------------------------
    #                                        CAWCR
    #------------------------------------------------------------------------------------
    def _lonlat_to_unit_sphere_xyz(self, lon, lat):
        """
        Convert lon/lat in degrees to unit-sphere Cartesian coordinates.
        """
        lon  = np.asarray(lon, dtype=float)
        lat  = np.asarray(lat, dtype=float)
        lonr = np.deg2rad(((lon + 180.0) % 360.0) - 180.0)
        latr = np.deg2rad(lat)
        x    = np.cos(latr) * np.cos(lonr)
        y    = np.cos(latr) * np.sin(lonr)
        z    = np.sin(latr)
        return np.column_stack([x.ravel(), y.ravel(), z.ravel()])

    def _angular_radius_to_chord(self, radius_km, earth_radius_km=6371.0):
        """
        Convert great-circle radius in km to chord distance on the unit sphere.
        """
        if radius_km is None:
            return None
        ang = radius_km / earth_radius_km
        return 2.0 * np.sin(0.5 * ang)

    def build_station_to_curvilinear_sparse_weights(self, src_lon, src_lat, tgt_lon, tgt_lat, p_weights,
                                                    target_mask = None,
                                                    k           = 1,
                                                    power       = 1.0,
                                                    radius_km   = None,
                                                    overwrite   = False,
                                                    chunk_size  = 100_000):
        """
        Build or load a sparse IDW remapping operator from station points to a
        curvilinear target grid.

        Parameters
        ----------
        src_lon, src_lat : 1D arrays
            Station coordinates.
        tgt_lon, tgt_lat : 2D arrays
            Target grid coordinates on (nj, ni).
        p_weights : str or Path
            Path to a sparse .npz weight matrix.
        target_mask : 2D array-like, optional
            1 for valid ocean points, 0 for masked points.
        k : int, default 8
            Number of nearest stations used in IDW.
        power : float, default 2.0
            IDW exponent.
        radius_km : float, optional
            Optional cutoff radius; if no station lies within radius, the nearest
            station is still used as fallback.
        overwrite : bool, default False
            Rebuild even if weight file already exists.
        chunk_size : int, default 100000
            Number of target points per KDTree query chunk.

        Returns
        -------
        scipy.sparse.csr_matrix
            Sparse interpolation matrix of shape (n_target, n_station).
        """
        p_weights = Path(p_weights)
        p_weights.parent.mkdir(parents=True, exist_ok=True)
        if p_weights.exists() and not overwrite:
            self.logger.info(f"Reusing wave sparse weights: {p_weights}")
            return load_npz(p_weights)
        self.logger.info(f"Creating wave sparse weights: {p_weights}")
        src_xyz      = self._lonlat_to_unit_sphere_xyz(src_lon, src_lat)
        tgt_xyz      = self._lonlat_to_unit_sphere_xyz(tgt_lon, tgt_lat)
        n_station    = src_xyz.shape[0]
        n_target     = tgt_xyz.shape[0]
        k            = int(min(k, n_station))
        radius_chord = self._angular_radius_to_chord(radius_km)
        tree         = cKDTree(src_xyz)
        if target_mask is None:
            active = np.ones(n_target, dtype=bool)
        else:
            active = np.asarray(target_mask).ravel().astype(bool)
        row_parts  = []
        col_parts  = []
        data_parts = []
        active_idx = np.where(active)[0]
        chunk_size = int(chunk_size)
        n_chunks   = (active_idx.size + chunk_size - 1) // chunk_size
        for ichunk, i0 in enumerate(range(0, active_idx.size, chunk_size), start=1):
            ii = active_idx[i0:i0 + chunk_size]
            self.logger.info(f"Wave weight chunk {ichunk}/{n_chunks}: {ii.size} target cells")
            d, ind = tree.query(tgt_xyz[ii], k=k)
            d, ind = tree.query(tgt_xyz[ii], k=k)
            if k == 1:
                d   = d[:, None]
                ind = ind[:, None]
            # Optional cutoff radius
            if radius_chord is not None:
                within = d <= radius_chord
            else:
                within = np.ones_like(d, dtype=bool)
            # Ensure at least nearest neighbour is used
            no_valid            = ~within.any(axis=1)
            within[no_valid, 0] = True
            # Inverse-distance weights
            d_safe = np.maximum(d, 1.0e-12)
            w      = np.where(within, d_safe ** (-power), 0.0)
            # Exact hit -> give full weight to nearest station
            exact  = d[:, 0] < 1.0e-12
            if np.any(exact):
                w[exact, :] = 0.0
                w[exact, 0] = 1.0
            wsum = w.sum(axis=1, keepdims=True)
            w    = np.where(wsum > 0.0, w / wsum, 0.0)
            rr   = np.repeat(ii, k)
            cc   = ind.reshape(-1)
            vv   = w.reshape(-1)
            keep = vv > 0.0
            row_parts.append(rr[keep])
            col_parts.append(cc[keep])
            data_parts.append(vv[keep])
        rows = np.concatenate(row_parts) if row_parts else np.array([], dtype=np.int64)
        cols = np.concatenate(col_parts) if col_parts else np.array([], dtype=np.int64)
        vals = np.concatenate(data_parts) if data_parts else np.array([], dtype=float)
        vals = vals.astype(np.float32, copy=False)
        W    = sparse.csr_matrix((vals, (rows, cols)), shape=(n_target, n_station), dtype=np.float32)
        save_npz(p_weights, W)
        self.logger.info(f"weight file saved to {p_weights}")
        return W

    def build_or_load_station_to_curvilinear_sparse_weights(self, src_lon, src_lat, tgt_lon, tgt_lat, p_weights,
                                                            target_mask = None,
                                                            k           = 8,
                                                            power       = 2.0,
                                                            radius_km   = None,
                                                            overwrite   = False,
                                                            chunk_size  = 100_000):
        """
        Thin wrapper for a more readable calling path from SeaIceWaves.
        """
        return self.build_station_to_curvilinear_sparse_weights(src_lon     = src_lon,
                                                                src_lat     = src_lat,
                                                                tgt_lon     = tgt_lon,
                                                                tgt_lat     = tgt_lat,
                                                                p_weights   = p_weights,
                                                                target_mask = target_mask,
                                                                k           = k,
                                                                power       = power,
                                                                radius_km   = radius_km,
                                                                overwrite   = overwrite,
                                                                chunk_size  = chunk_size)

    def apply_sparse_station_regridder(
        self,
        values,
        weights,
        target_shape,
        fill_value=np.nan,
        time_chunk=24,
    ):
        """
        Apply a sparse station->grid operator to values with shape:

            (time, station, frequency)

        Works with NumPy, xarray, or dask-backed arrays.
        Processes in time chunks to avoid forcing a full remote/lazy load.
        """
        import numpy as np
        import xarray as xr

        # Normalize input
        if isinstance(values, xr.DataArray):
            da = values.transpose("time", "station", "frequency")
            arr = da.data
            time_coord = da["time"].values
        else:
            arr = values
            time_coord = None

        shape = getattr(arr, "shape", None)
        if shape is None or len(shape) != 3:
            raise ValueError(f"Expected values with ndim=3 (time,station,frequency), got {shape}")

        n_time, n_station, n_freq = shape
        n_target = target_shape[0] * target_shape[1]

        if weights.shape != (n_target, n_station):
            raise ValueError(
                f"Weight matrix shape {weights.shape} incompatible with "
                f"target/station sizes {(n_target, n_station)}"
            )

        self.logger.info(
            f"Applying sparse station regridder in time chunks: "
            f"n_time={n_time}, n_station={n_station}, n_freq={n_freq}, "
            f"time_chunk={time_chunk}"
        )

        out_chunks = []

        for t0 in range(0, n_time, time_chunk):
            t1 = min(n_time, t0 + time_chunk)
            self.logger.info(f"Regridding wave spectra time chunk {t0}:{t1}")

            # This is where the remote/lazy load happens, but only for one chunk
            block = np.asarray(arr[t0:t1, :, :], dtype=np.float32)   # (tb, station, freq)

            tb = block.shape[0]

            # Recast to (station, tb*freq)
            vals2 = np.moveaxis(block, 1, 0).reshape(n_station, tb * n_freq)

            # Sparse apply
            out2 = weights @ vals2                                   # (n_target, tb*freq)

            # Reshape back to (tb, nj, ni, freq)
            out = out2.reshape(n_target, tb, n_freq)
            out = out.reshape(target_shape[0], target_shape[1], tb, n_freq)
            out = np.moveaxis(out, 2, 0)

            out = np.where(np.isfinite(out), out, fill_value).astype(np.float32, copy=False)
            out_chunks.append(out)

        return np.concatenate(out_chunks, axis=0)
