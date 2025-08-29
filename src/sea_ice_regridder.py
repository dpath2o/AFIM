import gc, re, dask
import xarray            as xr
import numpy             as np
import pandas            as pd
from pathlib             import Path
from pyproj              import Transformer
from pyresample.geometry import AreaDefinition, SwathDefinition
from pyresample.kd_tree  import resample_nearest

__all__ = ["SeaIceRegridder"]

class SeaIceRegridder:
    def __init__(self, **kwargs):
        return

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

    def resample_swath_to_area(self, src_da, lat2d, lon2d, area_def, radius=10_000, fill_value=np.nan):
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
        x, y = grid_coords_from_area(area_def, pixel_size=PIXEL)
        da_out = xr.DataArray(
            out2d,
            dims=("y", "x"),
            coords={"x": ("x", x, {"units": "m", "standard_name": "projection_x_coordinate"}),
                    "y": ("y", y, {"units": "m", "standard_name": "projection_y_coordinate"})},
            name=src_da.name,
            attrs={"crs": "EPSG:3031", "grid_mapping": "spstereo", "res": float(PIXEL), **src_da.attrs},
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
        lon, lat = xr.apply_ufunc(_xy_to_lonlat, X2D, Y2D,
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
        TLON = SI_tools.normalise_longitudes(da.coords[lon_name], wrap)
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