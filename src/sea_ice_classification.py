"""
SeaIceClassification
====================

Class utilities for classifying Antarctic fast ice using sea-ice concentration
(`aice`) and speed thresholds on multiple grids:

- B  : native CICE B-grid (cell centers) magnitude computed from (uvel, vvel)
- Ta : 2×2 B-grid box average (area-like) magnitude
- Tx : B-grid regridded to T-grid magnitude via precomputed weights

High-level API
--------------
- classify_fast_ice(...): end-to-end daily fast-ice mask, binary-day mask, and
  optional rolling-mean variant.
- classify_binary_days_fast_ice(...): “persistence” classification over a sliding
  window (N days) requiring ≥M fast-ice days within the window.

Assumptions / required attributes
---------------------------------
Before calling the high-level methods, the following attributes/methods must exist:

Attributes (populated externally or via your toolbox manager):
- self.CICE_dict: dict with keys
    "y_dim_length", "x_dim_length", "time_dim", "FI_chunks",
    "spatial_dims", "drop_coords"
- self.icon_thresh: float   (sea-ice concentration threshold, e.g., 0.15)
- self.ispd_thresh: float   (speed threshold for fast ice, e.g., 1e-4 m/s)
- self.bin_win_days: int    (persistence window length, e.g., 31)
- self.bin_min_days: int    (min fast-ice count in window, e.g., 16)
- self.mean_period: int     (rolling mean length in days)
- self.BT_composite_grids: iterable of {"B","Ta","Tx"} to average into composite
- self.G_t, self.G_u: dicts providing T- and U-grid lon/lat arrays
- self.dt_range: pandas.DatetimeIndex used by np3d_to_xr3d

Methods (provided elsewhere in your stack):
- self.load_bgrid()
- self.define_reG_weights()
- self.load_cice_zarr(...)
- self.slice_hemisphere(xr.DataArray/xr.Dataset)
- self.define_datetime_vars(dt0_str, dtN_str)
- self.reG(xr.DataArray)   (B→T regridding using precomputed weights)
"""

import json, os, shutil, logging, re, dask, gc
import xarray    as xr
import pandas    as pd
import numpy     as np
import xesmf     as xe
from collections import defaultdict
from pathlib     import Path
from datetime    import datetime, timedelta
import concurrent.futures
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

class SeaIceClassification:
    """
    Classify fast ice from CICE outputs using concentration and speed thresholds,
    with support for daily, persistence (binary-day), and rolling-mean variants.

    Parameters
    ----------
    sim_name : str, optional
        Name/identifier of the simulation (used for logs/paths).
    logger : logging.Logger, optional
        Logger instance. If None, uses a module-level logger.
    **kwargs
        Arbitrary configuration attributes set directly on `self`. Common keys
        include thresholds, window lengths, and dictionaries noted in the module
        docstring.

    Notes
    -----
    - The “daily” fast-ice mask is defined where:
        (aice > icon_thresh) AND (0 < ispd <= ispd_thresh)
      after composing `ispd` from one or more of {B, Ta, Tx}.
    - The binary-day mask applies a centered sliding window that marks a day as
      “persistent fast ice” if at least `bin_min_days` within `bin_win_days` are
      fast-ice days.
    - Several operations use Dask via `xr.apply_ufunc(..., dask="parallelized")`,
      but key intermediates are explicitly `.compute()`’d for determinism. Ensure
      available memory/cluster config is appropriate for your domain size.
    """
    def __init__(self, sim_name=None, logger=None, **kwargs):
        self.sim_name = sim_name
        self.logger = logger or logging.getLogger(__name__)
        for k, v in kwargs.items():
            setattr(self, k, v)

    def _rolling_mean(self, 
                      np_mag      : np.ndarray,
                      mean_period : int = None) -> np.ndarray:
        """
        Apply a centered 1-D uniform (boxcar) filter along the time axis.

        Parameters
        ----------
        np_mag : np.ndarray
            Array of shape (time, y, x) containing speed magnitudes on a common grid.
        mean_period : int, optional
            Window length in days. Defaults to `self.mean_period`.

        Returns
        -------
        np.ndarray
            Array with the same shape as `np_mag`, boxcar-smoothed along axis 0.

        Notes
        -----
        Uses `scipy.ndimage.uniform_filter1d(..., mode="nearest")`, which is
        approximately centered. For odd window sizes, alignment is symmetric.
        """
        mean_period = mean_period if mean_period is not None else self.mean_period 
        from scipy.ndimage import uniform_filter1d
        return uniform_filter1d(np_mag, size=mean_period, axis=0, mode="nearest")

    def B_mag_to_Ta_mag_numpy(self, B_mag_np: np.ndarray) -> np.ndarray:
        """
        Compute a 2×2 area-average (“Ta”) from a B-grid (cell-centered) magnitude.

        Parameters
        ----------
        B_mag_np : np.ndarray
            B-grid magnitude with shape (time, y, x).

        Returns
        -------
        np.ndarray
            Area-averaged array with shape (time, y_len, x_len) where (y_len, x_len)
            are taken from `self.CICE_dict["y_dim_length"]` and `["x_dim_length"]`.

        Algorithm
        ---------
        Averages the four neighbors:
            v00 = B[:, :-1, :-1], v01 = B[:, :-1,  1:],
            v10 = B[:,  1:, :-1], v11 = B[:,  1:,  1:]
        Pads the last row/col as needed and wraps the final x-column to match the
        first (cyclic longitude), then trims to (y_len, x_len).

        Notes
        -----
        - Input must be (t, y, x) and contiguous enough for slicing.
        - Longitude wrapping assumes periodicity in the x-direction.
        """
        y_len = self.CICE_dict["y_dim_length"]
        x_len = self.CICE_dict["x_dim_length"]
        self.logger.debug(f"input shape to Ta spatial averaging: {B_mag_np.shape}")
        # Compute average of 4 neighboring cells
        v00 = B_mag_np[:, :-1, :-1]
        v01 = B_mag_np[:, :-1, 1:]
        v10 = B_mag_np[:, 1:, :-1]
        v11 = B_mag_np[:, 1:, 1:]
        avg = (v00 + v01 + v10 + v11) / 4.0
        # Pad to restore shape (time, y_len, x_len)
        pad_y = max(y_len - avg.shape[1], 0)
        pad_x = max(x_len - avg.shape[2], 0)
        avg = np.pad(avg, pad_width=((0, 0), (0, pad_y), (0, pad_x)), constant_values=np.nan)
        # Wrap x-dimension if necessary
        if avg.shape[2] > 1:
            avg[:, :, -1] = avg[:, :, 0]
        # Final trim in case padding overshoots
        return avg[:, :y_len, :x_len]

    def np3d_to_xr3d(self, np_da, tgrid=True, name="ispd_B"):
        """
        Wrap a (t, y, x) NumPy array as an `xarray.DataArray` with lon/lat coords.

        Parameters
        ----------
        np_da : np.ndarray
            Array of shape (time, nj, ni).
        tgrid : bool, default True
            If True, use T-grid lon/lat from `self.G_t`; otherwise use U-grid lon/lat
            from `self.G_u`.
        name : str, default "ispd_B"
            Name for the resulting DataArray.

        Returns
        -------
        xarray.DataArray
            DataArray with dims ("time", "nj", "ni") and coords {"lon","lat","time"}.

        Notes
        -----
        Requires `self.dt_range` to be a `DatetimeIndex` aligned with `np_da.shape[0]`.
        Ensure `self.G_t`/`self.G_u` provide 2-D lon/lat arrays matching (nj, ni).
        """
        if tgrid:
            lon = self.G_t["lon"].values
            lat = self.G_t["lon"].values
        else:
            lon = self.G_u["lon"].values
            lat = self.G_u["lon"].values
        return xr.DataArray(data   = np_da,
                            dims   = ("time", "nj", "ni"),
                            coords = {"lon"  : (("nj", "ni"), lon),
                                      "lat"  : (("nj", "ni"), lat),
                                      "time" : self.dt_range},
                            name   = "ispd_B" )

    def B_mag_to_Tx_mag(self, np_B):
        """
        Regrid a B-grid magnitude cube to the T-grid using precomputed weights.

        Parameters
        ----------
        np_B : np.ndarray
            B-grid magnitude with shape (time, nj, ni).

        Returns
        -------
        xarray.DataArray
            T-grid magnitude with shape (time, nj, ni), produced by `self.reG(...)`.

        Requirements
        ------------
        - `self.define_reG_weights()` must have been called.
        - `self.reG(xr.DataArray)` must exist and implement the B→T mapping.

        Notes
        -----
        The input is first wrapped with `np3d_to_xr3d(..., tgrid=False)` to attach
        U-grid coordinates, then passed to `self.reG`. The result is computed lazily
        unless downstream `.compute()` is invoked.
        """

        xr_B = self.np3d_to_xr3d(np_B, tgrid=False)
        self.logger.debug(f"input shape to Tx spatial averaging: {np.shape(xr_B)}")
        return self.reG(xr_B)

    def _binary_days(self, mask):
        """
        Internal NumPy implementation of the binary-day persistence filter.

        Parameters
        ----------
        mask : np.ndarray
            Boolean/0-1 mask of shape (time, nj, ni) for daily fast-ice presence.

        Returns
        -------
        np.ndarray
            UInt8 mask (0/1) of shape (time, nj, ni) after applying a centered
            sliding window of length `self.bin_win_days` and requiring
            `self.bin_min_days` fast-ice days within the window.

        Notes
        -----
        - Uses `numpy.lib.stride_tricks.sliding_window_view` for efficiency.
        - Output is padded to match the original time length by centering.
        """
        from numpy.lib.stride_tricks import sliding_window_view
        time_len = mask.shape[0]
        self.logger.debug(">> mask shape:", mask.shape)
        windowed = sliding_window_view(mask, window_shape=self.bin_win_days, axis=0)
        self.logger.debug(">> windowed shape:", windowed.shape)  # Should be (time - win + 1, nj, ni, win)
        summed = windowed.sum(axis=-1)
        self.logger.debug(">> summed shape:", summed.shape)
        meets_thresh = (summed >= self.bin_min_days).astype(np.uint8)
        self.logger.debug(">> min_thesh shape:", meets_thresh.shape)
        pad_before = (time_len - meets_thresh.shape[0]) // 2
        self.logger.debug(">> pad_before:", pad_before)
        pad_after  = time_len - meets_thresh.shape[0] - pad_before
        self.logger.debug(">> pad_after:", pad_after)
        padded = np.pad(meets_thresh,
                        pad_width=((pad_before, pad_after), (0, 0), (0, 0)),
                        constant_values=0)
        self.logger.debug(">> after padding shape:",np.shape(padded))
        return padded

    def classify_binary_days_fast_ice(self, FI_mask,
                                      bin_win_days = None,
                                      bin_min_days = None,
                                      time_dim     = None):
        """
        Create a binary-day (persistence) fast-ice mask from a daily mask.

        Parameters
        ----------
        FI_mask : xarray.DataArray or xarray.Dataset
            Daily fast-ice mask with boolean/0-1 values over dims (time, y, x).
            If a Dataset is provided, it must contain a variable named "FI_mask".
        bin_win_days : int, optional
            Sliding-window length (days). Defaults to `self.bin_win_days`.
        bin_min_days : int, optional
            Minimum number of fast-ice days in the window. Defaults to `self.bin_min_days`.
        time_dim : str, optional
            Name of the time dimension. Defaults to `self.CICE_dict["time_dim"]`.

        Returns
        -------
        xarray.Dataset
            Dataset with a single variable `"FI_mask"` (dtype=uint8) holding the
            persistence classification, chunked per `self.CICE_dict["FI_chunks"]`.

        Raises
        ------
        ValueError
            If two spatial dimensions cannot be inferred from the input.
        """
        time_dim     = time_dim     or self.CICE_dict["time_dim"]
        bin_win_days = bin_win_days or self.bin_win_days
        bin_min_days = bin_min_days or self.bin_min_days
        spatial_dims = [d for d in FI_mask.dims if d != time_dim]
        if len(spatial_dims) != 2:
            raise ValueError(f"Expected two spatial dimensions, got: {spatial_dims}")
        dim_time, dim_y, dim_x = time_dim, spatial_dims[0], spatial_dims[1]
        self.logger.debug(f"FI_mask shape before apply_ufunc: {FI_mask.shape}")
        self.logger.debug(f"FI_mask dimensions before apply_ufunc: {FI_mask.dims}")
        da_FI_bin    = xr.apply_ufunc(self._binary_days,
                                      FI_mask["FI_mask"].data if isinstance(FI_mask, xr.Dataset) else FI_mask,
                                      input_core_dims=[[time_dim, dim_y, dim_x]],
                                      output_core_dims=[[time_dim, dim_y, dim_x]],
                                      dask = "parallelized", output_dtypes=[np.uint8],)
        return da_FI_bin.rename("FI_mask").to_dataset().chunk(self.CICE_dict["FI_chunks"])

    def _compute_fast_ice_mask(self, ispd_B, aice, dt0_str, dtN_str, label=""):
        """
        Internal: build the daily fast-ice mask on the analysis hemisphere.

        Parameters
        ----------
        ispd_B : xarray.DataArray
            B-grid speed magnitude (t, y, x) derived from (uvel, vvel).
        aice : xarray.DataArray
            Sea-ice concentration (t, y, x).
        dt0_str, dtN_str : str
            Requested start/end (ISO yyyy-mm-dd). Used for logs and cropping.
        label : str
            Tag for logging (e.g., "daily", "rolling-mean").

        Returns
        -------
        xarray.Dataset
            Dataset with boolean `"FI_mask"` over the selected hemisphere and time
            span, chunked per `self.CICE_dict["FI_chunks"]`.

        Steps
        -----
        1) Compute Ta via 2×2 average, Tx via regridding, and composite them per
           `self.BT_composite_grids` by mean over a concat dimension (skipna=True).
        2) Slice both speed and concentration to the target hemisphere.
        3) Apply thresholds:
             (aice > icon_thresh) & (0 < ispd <= ispd_thresh)
        4) Package as Dataset, drop spatial coords if they duplicate dims.

        Notes
        -----
        Some intermediates are `.compute()`’d to avoid deferred graph growth.
        """
        ispd_Ta    = xr.apply_ufunc(self.B_mag_to_Ta_mag_numpy, ispd_B, dask="parallelized", output_dtypes=[np.float32])
        ispd_Tx   = self.B_mag_to_Tx_mag(ispd_B).compute()
        xr_arrays = []
        valid_ispds = {"B": lambda: self.np3d_to_xr3d(ispd_B, tgrid=False),
                       "Ta": lambda: self.np3d_to_xr3d(ispd_Ta, name="ispd_Ta"),
                       "Tx": lambda: self.np3d_to_xr3d(ispd_Tx, name="ispd_Tx")}
        self.logger.info("creating composite sea ice speed 3D array ('ispd_BT') from:")
        for key in self.BT_composite_grids:
            xr_obj = valid_ispds[key]()
            self.logger.info(f"   {key}")
            self.logger.debug(f"   [{label}] shape to {key}: {xr_obj.shape}")
            xr_arrays.append(xr_obj)
        ispd_BT  = xr.concat(xr_arrays, dim="__concat_dim__").mean(dim="__concat_dim__", skipna=True).compute()
        ispd_hem = self.slice_hemisphere(ispd_BT)
        aice_hem = self.slice_hemisphere(aice)
        self.logger.info(f"FAST ICE (3D DATA ARRAY) MASK CREATION:")
        self.logger.info(f"   1. masking sea ice concentration: {self.icon_thresh:0.2f} < 'aice_hem'")
        self.logger.info(f"   2. masking sea ice speed        : 0 < 'ispd_hem' <= {self.ispd_thresh:.1e}")
        BT_mask = (aice_hem > self.icon_thresh) & (ispd_hem > 0) & (ispd_hem <= self.ispd_thresh)
        self.logger.info(f"   3. mask array to_dataset(), chunk() and compute()")
        FI_raw   = BT_mask.rename("FI_mask").to_dataset().chunk(self.CICE_dict["FI_chunks"]).compute()
        self.logger.info(f"   4. drop any spatial coordinates that are spatial dimensions")
        for dim in self.CICE_dict['spatial_dims']:
            if dim in FI_raw.coords:
                FI_raw = FI_raw.drop_vars(dim)
        return FI_raw

    def classify_fast_ice(self,
                          dt0_str               : str  = None,
                          dtN_str               : str  = None,
                          bin_win_days          : int  = None,
                          bin_min_days          : int  = None,
                          enable_rolling_output : bool = False,
                          roll_win_days         : int  = None,
                          time_dim_name         : str  = None):
        """
        End-to-end fast-ice classification over a date range.

        Parameters
        ----------
        dt0_str, dtN_str : str
            Start/end dates (inclusive) in "YYYY-mm-dd". The method will load a
            slightly extended range to accommodate centered windows (rolling or
            binary-day), then crop back to the requested interval.
        bin_win_days : int, optional
            Persistence window length (days). Defaults to `self.bin_win_days`.
        bin_min_days : int, optional
            Minimum number of fast-ice days within the window. Defaults to `self.bin_min_days`.
        enable_rolling_output : bool, default False
            If True, additionally compute a rolling-mean classification (speed
            smoothed with `roll_win_days`) alongside daily and binary-day products.
        roll_win_days : int, optional
            Rolling window size for speed smoothing. Defaults to `self.mean_period`.
        time_dim_name : str, optional
            Name of the time dimension. Defaults to `self.CICE_dict["time_dim"]`.

        Returns
        -------
        tuple
            (FI_dly, FI_bin, FI_roll)
            where:
              - FI_dly : xarray.Dataset with boolean `"FI_mask"` (daily fast-ice)
              - FI_bin : xarray.Dataset with uint8   `"FI_mask"` (binary-day persistence)
              - FI_roll: xarray.Dataset (as FI_dly) if `enable_rolling_output`,
                         otherwise `None`.
            All outputs are cropped to [dt0_str, dtN_str] and chunked per
            `self.CICE_dict["FI_chunks"]`.

        Workflow
        --------
        1) `load_bgrid()` and `define_reG_weights()` to ensure grids/weights ready.
        2) Extend the requested time range by half the relevant window length:
           - rolling enabled  → extend by `mean_period // 2`
           - else (binary-day)→ extend by `bin_win_days // 2`
        3) `load_cice_zarr(...)` for ['aice','uvel','vvel'] within the extended range,
           then clamp to actually available data (logs will note if truncated).
        4) Compute base `ispd_B = hypot(uvel, vvel)` (float32).
        5) Build daily mask via `_compute_fast_ice_mask(...)`.
        6) Build binary-day mask via `classify_binary_days_fast_ice(...)`.
        7) If enabled, smooth `ispd_B` with `_rolling_mean(...)` and re-classify.

        Performance
        -----------
        - `ispd_B`, the composite speed, and the daily mask are `.compute()`’d, which
          materializes arrays. Use appropriate Dask cluster settings for your domain.
        - Drop of `self.CICE_dict["drop_coords"]` minimizes coordinate bloat early.

        See Also
        --------
        classify_binary_days_fast_ice : Standalone persistence classification.
        """
        time_dim      = time_dim_name or self.CICE_dict["time_dim"]
        bin_win_days  = bin_win_days  or self.bin_win_days
        bin_min_days  = bin_min_days  or self.bin_min_days
        roll_win_days = roll_win_days or self.mean_period
        self.load_bgrid()
        self.define_reG_weights()
        # Define extended period based on method
        ext_per = (self.mean_period // 2) if enable_rolling_output else (self.bin_win_days // 2)
        dt0_ext = pd.to_datetime(dt0_str) - timedelta(days=ext_per)
        dtN_ext = pd.to_datetime(dtN_str) + timedelta(days=ext_per)
        # Load model data (clamped internally)
        self.logger.info(f"loading model data between {dt0_ext.strftime('%Y-%m-%d')} and {dtN_ext.strftime('%Y-%m-%d')}")
        CICE_all = self.load_cice_zarr( slice_hem = False,
                                        variables = ['aice', 'uvel', 'vvel'],
                                        dt0_str   = dt0_ext.strftime('%Y-%m-%d'),
                                        dtN_str   = dtN_ext.strftime('%Y-%m-%d'))
        # Clamp to actual available data
        dt0_avail = pd.to_datetime(CICE_all.time.values[0]).normalize()
        dtN_avail = pd.to_datetime(CICE_all.time.values[-1]).normalize()
        dt0_str   = pd.to_datetime(dt0_str).normalize()
        dtN_str   = pd.to_datetime(dtN_str).normalize()
        dt0_ext   = pd.to_datetime(dt0_ext).normalize()
        dtN_ext   = pd.to_datetime(dtN_ext).normalize()
        self.logger.debug(f"REQUESTED RANGE : {dt0_str} to {dtN_str}")
        self.logger.debug(f"EXTENDED  RANGE : {dt0_ext} to {dtN_ext}")
        self.logger.debug(f"AVAILABLE RANGE : {dt0_avail} to {dtN_avail}")
        dt0_load = dt0_avail if dt0_avail > dt0_ext else dt0_ext
        dtN_load = dtN_avail if dtN_avail < dtN_ext else dtN_ext
        if dt0_load != dt0_ext:
            self.logger.info(f"** EARLIEST AVAILABLE DATA IS {dt0_avail.date()} — differs from extended start date {dt0_ext.date()}")
        if dtN_load != dtN_ext:
            self.logger.info(f"** MOST RECENT AVAILABLE DATA IS {dtN_avail.date()} — differs from extended end date {dtN_ext.date()}")
        # Define the time range for internal calculations
        self.define_datetime_vars(dt0_str=dt0_load, dtN_str=dtN_load)
        # Drop unnecessary variables early
        aice = CICE_all['aice'].drop_vars(self.CICE_dict["drop_coords"])
        uvel = CICE_all['uvel'].drop_vars(self.CICE_dict["drop_coords"])
        vvel = CICE_all['vvel'].drop_vars(self.CICE_dict["drop_coords"])
        # Compute base ispd_B
        self.logger.info("computing ispd_B from uvel and vvel")
        ispd_B = xr.apply_ufunc(np.hypot, uvel, vvel, dask="parallelized", output_dtypes=[np.float32]).compute()
        # Prepare base daily fast ice mask
        FI_dly = self._compute_fast_ice_mask(ispd_B, aice, dt0_str, dtN_str, label="daily")
        # Prepare binary-day mask
        self.logger.info("classifying binary-day fast ice mask")
        FI_bin = self.classify_binary_days_fast_ice(FI_dly["FI_mask"],
                                                    bin_win_days = bin_win_days,
                                                    bin_min_days = bin_min_days,
                                                    time_dim     = time_dim).sel(time=slice(dt0_str, dtN_str))
        FI_dly = FI_dly.sel(time=slice(dt0_str, dtN_str))
        # Optional: prepare rolling classification
        if enable_rolling_output:
            self.logger.info(f"applying rolling mean (period={roll_win_days})")
            ispd_B_roll = self._rolling_mean(ispd_B, mean_period=roll_win_days)
            FI_roll     = self._compute_fast_ice_mask(ispd_B_roll, aice, dt0_str, dtN_str, label="rolling-mean").sel(time=slice(dt0_str, dtN_str))
        else:
            FI_roll = None
        return FI_dly, FI_bin, FI_roll