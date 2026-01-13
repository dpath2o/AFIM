import os, gc, re, dask
import xarray            as xr
import numpy             as np
import pandas            as pd
import xesmf             as xe
from pathlib             import Path
from pyproj              import Transformer
from pyresample.geometry import AreaDefinition, SwathDefinition
from pyresample.kd_tree  import resample_nearest

__all__ = ["SeaIceGridWork"]

class SeaIceGridWork:

    def __init__(self, **kwargs):
        return

    def normalise_longitudes(self, lon, to="0-360", eps=1e-12):
        """
        Wrap longitudes to either [0, 360) or (-180, 180].
        Works with numpy arrays, scalars, and xarray objects.
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

    def reapply_landmask(self, DS, apply_unmodified=False):
        """

        Apply landmask to all spatial variables in the dataset.

        Uses the modified bathymetry (`kmt_mod`) field from the model grid to mask out land cells.
        Applies this mask to all variables with dimensions `("nj", "ni")`, ensuring that land areas
        are excluded from any subsequent analysis or output.

        INPUTS:
           DS : xarray.Dataset; Dataset containing sea ice fields to be masked.

        OUTPUTS:
           xarray.Dataset; Same dataset with land cells set to NaN for spatial variables.

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
        Circular mean for longitudes in degrees. Inputs must be broadcastable to same shape.
        Returns degrees in [0, 360).
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
        Build C-grid geometry from *cell-centre* lon/lat arrays.

        Parameters
        ----------
        lon, lat : array-like
            2D arrays of cell-centre lon/lat, either in degrees or radians.
            Shape (nj, ni).
        source_in_radians : bool
            If True, lon/lat are radians and will be converted to degrees.

        Returns
        -------
        lon_b, lat_b : ndarray
            Corner coordinates, shape (nj+1, ni+1)
        lon_e, lat_e : ndarray
            Vertical face centres (E/W faces), shape (nj, ni+1)
            (dims can be interpreted as (nj, ni_b))
        lon_n, lat_n : ndarray
            Horizontal face centres (N/S faces), shape (nj+1, ni)
            (dims can be interpreted as (nj_b, ni))
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
        Unified grid loader.

        Sets class attributes (no return):
        - self.G_t : t-grid dataset (centres + corners, mask, angle, area)
        - self.G_u : legacy u-grid dataset (centres + corners, mask, angle, area)
        - self.G_e : C-grid vertical-face coords derived from t-grid corners (lon_e/lat_e)
        - self.G_n : C-grid horizontal-face coords derived from t-grid corners (lon_n/lat_n)

        Parameters
        ----------
        slice_hem : bool
            If True, applies self.slice_hemisphere() to all constructed grid datasets.
        build_faces : bool
            If True, compute C-grid face coordinate datasets (G_e, G_n).
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
        Convert radians->degrees if values look like radians.
        Returns numpy array (float64).
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

    def _open_cice_cgrid_for_F2(self, P_grid=None):
        """
        Open a CICE C-grid NetCDF and return (tlon_deg, tlat_deg, anglet_rad, dx_m, dy_m).
        Handles common variable name variants.

        Expected for ACCESS-OM3 C-grid:
          tlon, tlat in radians
          anglet in radians
          hte, htn in cm
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
        dx_m = np.asarray(hte, dtype="float64") * 0.01
        dy_m = np.asarray(htn, dtype="float64") * 0.01
        return tlon_deg, tlat_deg, anglet_rad, dx_m, dy_m, ds

    def _iter_exterior_rings_lonlat(self, P_shp, target_crs="EPSG:4326"):
        """
        Stream exterior rings from a Polygon/MultiPolygon shapefile.

        Yields:
            (lon_deg_1d, lat_deg_1d) for each exterior ring, in target_crs degrees.
        """
        try:
            import fiona
            from pyproj import CRS, Transformer
        except Exception as e:
            raise ImportError("Requires fiona and pyproj for streaming shapefile rings.") from e
        with fiona.open(P_shp) as src:
            # Fiona may provide crs or crs_wkt
            src_crs = None
            if src.crs_wkt:
                src_crs = CRS.from_wkt(src.crs_wkt)
            elif src.crs:
                src_crs = CRS.from_user_input(src.crs)
            else:
                # assume already lon/lat
                src_crs = CRS.from_epsg(4326)
            dst_crs = CRS.from_user_input(target_crs)
            tfm = Transformer.from_crs(src_crs, dst_crs, always_xy=True)
            for feat in src:
                geom = feat.get("geometry", None)
                if geom is None:
                    continue
                gtype = geom.get("type", None)
                coords = geom.get("coordinates", None)
                if coords is None:
                    continue
                # Polygon: coords = [ring0, ring1, ...], ring0 is exterior
                # MultiPolygon: coords = [[ring0, ...], [ring0, ...], ...]
                if gtype == "Polygon":
                    polys = [coords]
                elif gtype == "MultiPolygon":
                    polys = coords
                else:
                    # ignore non-polygons
                    continue
                for poly in polys:
                    if not poly:
                        continue
                    exterior = poly[0]
                    if exterior is None or len(exterior) < 2:
                        continue
                    x = np.asarray([p[0] for p in exterior], dtype="float64")
                    y = np.asarray([p[1] for p in exterior], dtype="float64")
                    lon, lat = tfm.transform(x, y)
                    lon = np.asarray(lon, dtype="float64")
                    lat = np.asarray(lat, dtype="float64")
                    # normalise for downstream projection
                    lon = self.normalise_longitudes(lon, to="-180-180")
                    yield lon, lat

    def build_F2_form_factors_from_high_res_coast(self,
                                                  P_shp              = None,
                                                  P_grid             = None,
                                                  P_out              = None,
                                                  proj_crs           = "EPSG:3031",
                                                  chunk_segments     = 2_000_000,
                                                  coast_write_stride = 25,
                                                  lat_subset_max     = -30.0,
                                                  netcdf_compression = 4):
        """
        Compute Liu et al. (2022) F2 form factors (Equations 9–10) on the CICE T-grid:

            F2x(i,j) = sum_n |l_n cos(theta_n)| / dx(i,j)
            F2y(i,j) = sum_n |l_n sin(theta_n)| / dy(i,j)

        where theta_n is the angle between the coastline segment and the local model x-axis.

        Output NetCDF contains:
          - F2x(nj,ni), F2y(nj,ni)
          - lon(ncoast), lat(ncoast)  [thinned coastline vertices for provenance/plotting]

        Notes:
          * This produces the *cell-based* integrals (Liu f^u_2 and f^v_2). The C-grid u/v-point
            combination (their Eqs 11–12) is intentionally left for the Fortran side.
          * Coastline segment length + azimuth are computed geodesically (WGS84).
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
        tlon_deg, tlat_deg, anglet_rad, dx_m, dy_m, ds_grid = self._open_cice_cgrid_for_F2(P_grid=P_grid)
        nj, ni = tlon_deg.shape
        ncell = nj * ni
        # Subset to SH band for KDTree efficiency (Antarctic coastline only)
        mask = np.isfinite(tlat_deg) & (tlat_deg <= lat_subset_max)
        if not mask.any():
            raise RuntimeError(f"No grid cells found with tlat <= {lat_subset_max}. Check grid/units.")
        flat_idx = np.arange(ncell, dtype="int64")
        sub_flat = flat_idx[mask.ravel()]
        # Project T-cell centers for KDTree
        tfm_to_proj = Transformer.from_crs(CRS.from_epsg(4326), CRS.from_user_input(proj_crs), always_xy=True)
        x_sub, y_sub = tfm_to_proj.transform(tlon_deg.ravel()[sub_flat], tlat_deg.ravel()[sub_flat])
        pts = np.column_stack([np.asarray(x_sub, dtype="float64"), np.asarray(y_sub, dtype="float64")])
        self.logger.info(f"Building KDTree on {pts.shape[0]:,} T-cells (subset for lat <= {lat_subset_max}).")
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
        for lon_ring, lat_ring in self._iter_exterior_rings_lonlat(P_shp, target_crs="EPSG:4326"):
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
                _, nn     = tree.query(np.column_stack([xm, ym]), workers=-1)
                cell_flat = sub_flat[nn]  # map KDTree index -> full-grid flat index
                # theta = angle between segment direction and local model x-axis
                # az12 is degrees clockwise from north; convert to radians from east:
                # angle_from_east = 90deg - az_from_north
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
        ds_out.to_netcdf(P_out, encoding=enc)
        return ds_out

    def read_grounded_iceberg_csv(self, P_csv,
                                  lon_col          : str = "lon",
                                  lat_col          : str = "lat",
                                  normalise_lon_to : str = "0-360"):
        """
        Read a simple GI CSV with columns lon,lat (degrees). Future-proofed for
        optional columns such as area_m2, C_gi, orient_deg, t_start, t_end.

        Returns
        -------
        df : pandas.DataFrame with at least columns ["lon","lat"] in degrees.
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

    def compute_tcell_dxdy_m(self, force=False, cache=True):
        """
        Compute approximate T-cell dx,dy (meters) from face-center coordinates
        derived from t-grid corners (self.G_e, self.G_n).

        Uses Geod.inv on WGS84 ellipsoid.

        Returns
        -------
        dx_m, dy_m : np.ndarray, shape (nj, ni)
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

    def _get_tgrid_kdtree(self, proj_crs="EPSG:3031", force=False):
        """
        Build/cache a KDTree of projected T-cell centres for fast point->cell mapping.
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
        Map lon/lat points to nearest T-cell.

        Returns
        -------
        ii, jj : int arrays (same length as input)
        dist_m : float array (same length)
        valid  : bool array (same length) applying max_dist_km if provided
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
        Minimal DBSCAN-like clustering without sklearn dependency.
        Returns labels where -1 indicates noise.

        Notes
        -----
        - Connectivity is defined by neighbors within eps_m.
        - Core points: neighbor count >= min_samples.
        - Clusters grown from core points; border points are included if reachable.
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

    def build_F2_GI_from_csv(self, P_GI_CSV,
                             method                : str ="simple-geometry",
                             base_area_m2          : float = None,
                             C_gi                  : float = 1.0,
                             proj_crs              : str = "EPSG:3031",
                             max_map_dist_km       : float = 50.0,
                             # cluster-axis controls
                             eps_cluster_km        : float = 15.0,
                             min_cluster_size      : int = 3,
                             cluster_buffer_m      : float = None,
                             cluster_amplification : float = 0.0):
        """
        Build a GI-only (additional) form factor field on the T-grid.

        Parameters
        ----------
        method : {"simple-geometry","cluster-axis"}
        - simple-geometry: isotropic circle per GI point
        - cluster-axis: DBSCAN-like clusters + PCA along a principal axis -> oriented rectangle drag
        base_area_m2 : float
        Default iceberg area if CSV doesn't provide 'area_m2'. If None, uses median T-cell area.
        C_gi : float
        Base intrinsic coefficient applied to GI drag contributions (dimensionless).
        CSV may override per-point if it provides a 'C_gi' column.
        max_map_dist_km : float
        Reject GI points that map too far from the nearest T-cell centre.
        eps_cluster_km, min_cluster_size : clustering parameters (cluster-axis).
        cluster_buffer_m : optional buffer added to inferred cluster L,W. If None, uses circle radius from base_area_m2.
        cluster_amplification : float
        Optional multiplicative amplification for clusters:
            L,W *= (1 + cluster_amplification * log1p(n_points_in_cluster))

        Returns
        -------
        ds_gi : xr.Dataset with variables:
        - F2x_gi(nj,ni), F2y_gi(nj,ni)
        - gi_count(nj,ni)
        - plus 1D metadata arrays: gi_lon, gi_lat, gi_cluster_id, gi_i, gi_j (useful for debugging)
        """
        from pyproj import Transformer
        if (not hasattr(self, "G_t")) or (self.G_t is None):
            self.load_cice_grid(slice_hem=False, build_faces=True)
        nat_dim = self.CICE_dict.get("spatial_dims", ("nj", "ni"))
        nj, ni  = self.G_t["lat"].shape
        df = self.read_grounded_iceberg_csv(P_GI_CSV, normalise_lon_to="0-360")
        if base_area_m2 is None:
            base_area_m2 = float(np.nanmedian(self.G_t["area"].values))
        # Map points to grid
        ii, jj, dist_m, valid = self.map_points_to_tgrid(df["lon"].values, df["lat"].values,
                                                         proj_crs    = proj_crs,
                                                         max_dist_km = max_map_dist_km)
        df = df.loc[valid].reset_index(drop=True)
        ii = ii[valid]
        jj = jj[valid]
        # dx/dy (meters) on T-grid
        dx_m, dy_m = self.compute_tcell_dxdy_m(force=False, cache=True)
        # local grid rotation (degrees) from east->x axis; convert to radians
        ang_deg = self.G_t["angle"].values.astype("float64")
        ang_rad = np.deg2rad(ang_deg)
        # GI ... outputs on T-grid
        F2x_gi   = np.zeros((nj, ni), dtype="float64")
        F2y_gi   = np.zeros((nj, ni), dtype="float64")
        gi_count = np.zeros((nj, ni), dtype="int32")
        # per-point overrides ... future framework
        if "area_m2" in df.columns:
            area_m2 = df["area_m2"].fillna(base_area_m2).values.astype("float64")
        else:
            area_m2 = np.full(len(df), float(base_area_m2), dtype="float64")
        if "C_gi" in df.columns:
            Cvals = df["C_gi"].fillna(C_gi).values.astype("float64")
        else:
            Cvals = np.full(len(df), float(C_gi), dtype="float64")
        # Project GI points for clustering/orientation
        tr = Transformer.from_crs("EPSG:4326", proj_crs, always_xy=True)
        x_gi, y_gi = tr.transform(df["lon"].values.astype("float64"), df["lat"].values.astype("float64"))
        # ------------------------------------------------------------------
        # Method A: simple-geometry (isotropic circles)
        if method == "simple-geometry":
            # circle radius per GI
            r = np.sqrt(np.maximum(area_m2, 0.0) / np.pi)  # meters
            Lproj = 4.0 * r  # integral of |cos| around circle gives 4r
            for k in range(len(df)):
                i = int(ii[k])
                j = int(jj[k])
                if not (0 <= i < ni and 0 <= j < nj):
                    continue
                dx = dx_m[j, i]
                dy = dy_m[j, i]
                if (not np.isfinite(dx)) or (not np.isfinite(dy)):
                    continue
                # identical projections before dx/dy scaling
                F2x_gi[j, i] += Cvals[k] * (Lproj[k] / dx)
                F2y_gi[j, i] += Cvals[k] * (Lproj[k] / dy)
                gi_count[j, i] += 1
            labels = np.full(len(df), -1, dtype="int64")
        # ------------------------------------------------------------------
        # Method B: cluster-axis (anisotropic oriented rectangles)
        elif method == "cluster-axis":
            labels = self._dbscan_like_clusters(x_gi, y_gi,
                                                eps_m       = float(eps_cluster_km) * 1000.0,
                                                min_samples = int(min_cluster_size))
            # For each cluster: infer PCA axis + extents
            uniq = sorted([c for c in np.unique(labels) if c >= 0])
            # Default cluster buffer: circle radius from base_area_m2 (or median)
            if cluster_buffer_m is None:
                cluster_buffer_m = float(np.sqrt(base_area_m2 / np.pi))
            for cid in uniq:
                m = labels == cid
                idx = np.where(m)[0]
                npt = idx.size
                if npt < int(min_cluster_size):
                    continue
                X = np.column_stack([x_gi[idx], y_gi[idx]])  # (npt,2)
                # PCA principal axis
                Xc = X - X.mean(axis=0)
                C = np.cov(Xc.T)
                eigvals, eigvecs = np.linalg.eigh(C)
                v = eigvecs[:, np.argmax(eigvals)]  # principal axis unit vector in projected coords
                v = v / np.linalg.norm(v)
                # perpendicular
                u = np.array([-v[1], v[0]])
                s1 = Xc @ v
                s2 = Xc @ u
                L = (np.nanmax(s1) - np.nanmin(s1)) + 2.0 * cluster_buffer_m
                W = (np.nanmax(s2) - np.nanmin(s2)) + 2.0 * cluster_buffer_m
                # Optional amplification with cluster size
                if cluster_amplification and cluster_amplification != 0.0:
                    amp = 1.0 + float(cluster_amplification) * np.log1p(float(npt))
                    L *= amp
                    W *= amp
                # Orientation angle (from east, radians) in geographic tangent plane ~ projected plane
                phi = np.arctan2(v[1], v[0])  # rad, relative to +x=east
                # Distribute the cluster rectangle perimeter "budget" evenly across its points
                # Rectangle perimeter contributions:
                # x-projection total: 2L|cos(theta)| + 2W|sin(theta)|
                # y-projection total: 2L|sin(theta)| + 2W|cos(theta)|
                for k in idx:
                    i = int(ii[k])
                    j = int(jj[k])
                    if not (0 <= i < ni and 0 <= j < nj):
                        continue
                    dx = dx_m[j, i]
                    dy = dy_m[j, i]
                    if (not np.isfinite(dx)) or (not np.isfinite(dy)):
                        continue
                    # convert cluster axis to local model x-axis coordinates:
                    theta = phi - ang_rad[j, i]  # rad relative to local model x-axis
                    c = np.abs(np.cos(theta))
                    s = np.abs(np.sin(theta))
                    # total projected "length" budget for cluster, then per-point share
                    Lx_tot = 2.0 * L * c + 2.0 * W * s
                    Ly_tot = 2.0 * L * s + 2.0 * W * c
                    share = 1.0 / float(npt)
                    F2x_gi[j, i] += Cvals[k] * share * (Lx_tot / dx)
                    F2y_gi[j, i] += Cvals[k] * share * (Ly_tot / dy)
                    gi_count[j, i] += 1
        else:
            raise ValueError("method must be 'simple-geometry' or 'cluster-axis'")
        # Package result dataset
        ds = xr.Dataset(data_vars = {"F2x_gi": (nat_dim, F2x_gi.astype("float32"), {"long_name": "Grounded-iceberg additional form factor, x-projection (T-cell)",
                                                                                    "units": "1",
                                                                                    "method": method,
                                                                                    "C_gi_default": float(C_gi),
                                                                                    "base_area_m2": float(base_area_m2)}),
                                     "F2y_gi": (nat_dim, F2y_gi.astype("float32"), {"long_name": "Grounded-iceberg additional form factor, y-projection (T-cell)",
                                                                                    "units": "1",
                                                                                    "method": method,
                                                                                    "C_gi_default": float(C_gi),
                                                                                    "base_area_m2": float(base_area_m2)}),
                                     "gi_count" : (nat_dim, gi_count.astype("int32"), {"long_name": "Number of GI points mapped into each T-cell",
                                                                                       "units": "count"})},
                        coords = {nat_dim[0]: np.arange(nj),
                                  nat_dim[1]: np.arange(ni)},
                        attrs  = {"P_GI_CSV": str(P_GI_CSV),
                                  "proj_crs": str(proj_crs)})
        # 1D metadata ... diagnostics
        ds["gi_lon"]        = ("ngi", df["lon"].values.astype("float32"))
        ds["gi_lat"]        = ("ngi", df["lat"].values.astype("float32"))
        ds["gi_i"]          = ("ngi", ii.astype("int32"))
        ds["gi_j"]          = ("ngi", jj.astype("int32"))
        ds["gi_cluster_id"] = ("ngi", labels.astype("int32"))
        ds["gi_map_dist_m"] = ("ngi", dist_m[valid].astype("float32"))
        if getattr(self, "logger", None) is not None:
            self.logger.info(f"Built GI form factors ({method}): "
                             f"n_gi={len(df)}, nonzero_cells={(ds['gi_count'].values>0).sum()}")
        return ds

    def write_F2_with_GI(self,
                         P_F2_coast            : str = None,
                         P_GI_CSV              : str = None,
                         P_out                 : str = None,
                         method                : str = "simple-geometry",
                         base_area_m2          : float = None,
                         C_gi                  : float = 1.0,
                         # pass-through clustering controls
                         eps_cluster_km        : float = 15.0,
                         min_cluster_size      : int = 3,
                         cluster_buffer_m      : float = None,
                         cluster_amplification : float = 0.0,
                         overwrite             : bool = False):
        """
        Load an existing coast-only F2 file (contains F2x/F2y and lon/lat coast points),
        add GI contributions, and write an updated F2 NetCDF.

        Output conventions
        ------------------
        - F2x/F2y in output are the *combined* (coast + GI) fields (same names as coast file),
        so CICE can keep reading varnames 'F2x','F2y' unchanged.
        - Also writes:
            F2x_coast, F2y_coast  (original)
            F2x_gi,    F2y_gi     (GI-only)
        - Copies through coastline lon/lat variables if present in input file.

        Returns
        -------
        ds_out : xr.Dataset (also written to disk)
        """
        P_F2_coast = P_F2_coast or self.CICE_dict.get("P_F2_coast", None)
        P_out      = P_out      or self.CICE_dict.get("P_F2_GI", None)
        P_GI_CSV   = P_GI_CSV   or self.GI_dict.get("P_raw", None)
        self.logger.info(P_GI_CSV)
        if (not overwrite) and os.path.exists(P_out):
            raise FileExistsError(f"Output already exists: {P_out} (set overwrite=True to replace)")
        ds_coast = xr.open_dataset(P_F2_coast)
        if "F2x" not in ds_coast.variables or "F2y" not in ds_coast.variables:
            raise ValueError(f"{P_F2_coast} must contain variables 'F2x' and 'F2y'.")
        # Build GI-only increments on T-grid
        ds_gi = self.build_F2_GI_from_csv(P_GI_CSV,
                                          method=method,
                                          base_area_m2=base_area_m2,
                                          C_gi=C_gi,
                                          eps_cluster_km=eps_cluster_km,
                                          min_cluster_size=min_cluster_size,
                                          cluster_buffer_m=cluster_buffer_m,
                                          cluster_amplification=cluster_amplification)
        # Align dims (coast file uses (nj,ni) order)
        F2x_coast = ds_coast["F2x"].astype("float32")
        F2y_coast = ds_coast["F2y"].astype("float32")
        # ds_gi uses nat_dim from config; ensure same ordering as ds_coast
        # If ds_coast dims are ('nj','ni') and ds_gi dims are same, this is trivial.
        # Otherwise, we try to transpose ds_gi to match ds_coast.
        gi_x = ds_gi["F2x_gi"]
        gi_y = ds_gi["F2y_gi"]
        if tuple(F2x_coast.dims) != tuple(gi_x.dims):
            gi_x = gi_x.transpose(*F2x_coast.dims)
            gi_y = gi_y.transpose(*F2y_coast.dims)
        F2x_total = (F2x_coast + gi_x).astype("float32")
        F2y_total = (F2y_coast + gi_y).astype("float32")
        # Build output dataset
        ds_out = xr.Dataset()
        # Combined (primary)
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
        for v in ["gi_lon", "gi_lat", "gi_i", "gi_j", "gi_cluster_id", "gi_map_dist_m"]:
            if v in ds_gi.variables:
                ds_out[v] = ds_gi[v]
        ds_out.attrs.update({"P_F2_coast_in"         : str(P_F2_coast),
                             "P_GI_CSV"              : str(P_GI_CSV),
                             "method_gi"             : str(method),
                             "C_gi_default"          : float(C_gi),
                             "base_area_m2"          : float(base_area_m2) if base_area_m2 is not None else np.nan,
                             "eps_cluster_km"        : float(eps_cluster_km),
                             "min_cluster_size"      : int(min_cluster_size),
                             "cluster_buffer_m"      : float(cluster_buffer_m) if cluster_buffer_m is not None else np.nan,
                             "cluster_amplification" : float(cluster_amplification)})
        if getattr(self, "logger", None) is not None:
            self.logger.info(f"Writing combined F2 (coast+GI): {P_out}")
        ds_out.to_netcdf(P_out)
        ds_coast.close()
        return ds_out
