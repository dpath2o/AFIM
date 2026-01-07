import gc, re, dask, os
import xarray            as xr
import numpy             as np
import pandas            as pd
import xesmf             as xe
from pathlib             import Path
from pyproj              import Transformer
from pyresample.geometry import AreaDefinition, SwathDefinition
from pyresample.kd_tree  import resample_nearest

__all__ = ["SeaIceRegridder"]

class SeaIceRegridder:
    def __init__(self, **kwargs):
        return

    def _ensure_reG_defined(self):
        """Create/reuse xESMF regridder located sea_ice_regridder.py module"""
        if getattr(self, "reG", None) is not None and getattr(self, "reG", None) is not None:
            return
        self.define_reG_weights()

    def define_reG_weights(self):
        """
        Define and store an xESMF regridder to remap CICE B-grid (U-point) data to the T-grid.

        This method constructs two CICE grids—one at the U-point (B-grid) and one at the T-point (centered)—using 
        internally defined grid definitions. It then either reuses existing xESMF regridding weights or creates 
        them if they do not already exist.

        The resulting regridder is stored as `self.reG`, and a flag `self.reG_weights_defined` is set to True 
        upon successful creation or reuse of the regridder.

        Grid details:
        - The source grid (`G_u`) is built using U-point coordinates (`grid_type='u'`).
        - The target grid (`G_t`) is built using T-point coordinates (`grid_type='t'`) and includes a land-sea mask.
        - Both grids include corner information to support conservative or bilinear regridding.

        Regridding parameters:
        - Method: "bilinear"
        - Periodic: True (assumes global grid)
        - Degenerate cells: ignored
        - Extrapolation: nearest-neighbor (source to destination)
        - Weight reuse: enabled if existing file is found at path `self.CICE_dict["P_reG_u2t_weights"]`

        Logging:
        - Logs whether existing weights are reused or new weights are created.

        Returns
        -------
        None
            The regridder is stored as `self.reG` and is not returned explicitly.
        """
        self.define_cice_grid()
        G_u           = self.G_u
        G_t           = self.G_t
        G_t['mask']   = G_t['kmt_mod']
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
        Regrid a single B-grid DataArray to the T-grid using pre-defined xESMF regridder.

        If coord_names are not provided then assumes coordinate names provided in JSON configuration file,
        which are like ["ULON","ULAT"]

        INPUTS:
        da : xr.DataArray; b-grid variable with coordinates coord_names

        OUTPUTS:
        xr.DataArray; re-gridded DataArray on the T-grid with dimensions (time, nj, ni).
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

    def simple_spatial_averaging_bgrid_to_tgrid(self, var):
        """
        Dask-safe 4-point unweighted average from B-grid to T-grid.

        Uses efficient array shifting and avoids costly concatenation over new dimensions.

        Parameters
        ----------
        var : xr.DataArray
            2D or 3D (time, nj, ni) array on B-grid.

        Returns
        -------
        xr.DataArray
            Averaged field on T-grid with shape (time, nj, ni).
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
        Project a swath's lat/lon to EPSG:3031 and return [xmin, ymin, xmax, ymax] (+buffer).
        Wrap longitudes to [-180, 180) first to avoid dateline issues.
        """
        lon2d = self.normalise_longitudes(lon2d, wrap="-180-180")
        transformer = Transformer.from_crs("EPSG:4326", "EPSG:3031", always_xy=True)
        x, y = transformer.transform(lon2d.ravel(), lat2d.ravel())
        x = np.asarray(x); y = np.asarray(y)
        finite = np.isfinite(x) & np.isfinite(y)
        xmin, xmax = x[finite].min(), x[finite].max()
        ymin, ymax = y[finite].min(), y[finite].max()
        return [xmin - buffer_m, ymin - buffer_m, xmax + buffer_m, ymax + buffer_m]

    def union_extents(self, extents):
        xs = [e[0] for e in extents] + [e[2] for e in extents]
        ys = [e[1] for e in extents] + [e[3] for e in extents]
        return [min(xs), min(ys), max(xs), max(ys)]

    def snap_extent_to_grid(self, extent, pixel_size):
        xmin, ymin, xmax, ymax = extent
        xmin = np.floor(xmin / pixel_size) * pixel_size
        ymin = np.floor(ymin / pixel_size) * pixel_size
        xmax = np.ceil (xmax / pixel_size) * pixel_size
        ymax = np.ceil (ymax / pixel_size) * pixel_size
        return [xmin, ymin, xmax, ymax]

    def make_area_definition(self, extent, pixel_size=5_000, area_id="epsg3031_5km_union"):
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
        xmin, ymin, xmax, ymax = area_def.area_extent
        width, height = area_def.width, area_def.height
        # Cell centers
        x = xmin + (np.arange(width) + 0.5) * pixel_size
        y = ymax - (np.arange(height) + 0.5) * pixel_size  # top->down (north->south)
        return x, y

    def resample_swath_to_area(self, src_da, lat2d, lon2d, area_def, radius=10_000, fill_value=np.nan, pixel_size=5_000):
        """Nearest-neighbour resample a 2D swath (lat2d, lon2d) to an AreaDefinition grid."""
        lon2d = self.normalise_longitudes(lon2d, wrap="-180-180")  # << key fix: wrap before building the swath
        swath = SwathDefinition(lons=lon2d, lats=lat2d)
        out2d = resample_nearest(
            source_geo_def=swath,
            data=src_da.values,
            target_geo_def=area_def,
            radius_of_influence=radius,
            fill_value=fill_value,
            nprocs=0,            # set >0 to parallelise
            reduce_data=True,
        )
        x, y = grid_coords_from_area(area_def, pixel_size=pixel_size)
        da_out = xr.DataArray(
            out2d,
            dims=("y", "x"),
            coords={"x": ("x", x, {"units": "m", "standard_name": "projection_x_coordinate"}),
                    "y": ("y", y, {"units": "m", "standard_name": "projection_y_coordinate"})},
            name=src_da.name,
            attrs={"crs": "EPSG:3031", "grid_mapping": "spstereo", "res": float(pixel_size), **src_da.attrs},
        )
        return da_out

    def _xy_to_lonlat(self, x, y):
        T = Transformer.from_crs(3031, 4326, always_xy=True)
        lon, lat = T.transform(x, y)
        return lon, lat

    def add_lonlat_from_epsg3031(self, ds, 
                                x_name    = "x",
                                y_name    = "y",
                                wrap      = "0..360",       # or "-180..180"
                                out_dtype = "float32"):
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
        Subset a curvilinear grid by a geographic box. Works with dask.
        If the lon_range crosses the seam (e.g. (350, 20) in 0..360), it handles it.
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