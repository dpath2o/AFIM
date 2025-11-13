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
- self.BT_composite_grids: iterable of {"B","Ta","Tb","Tx"} to average into composite
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
from numpy.lib.stride_tricks import sliding_window_view
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

    def Bavg_methA(self, B_component: np.ndarray, y_len, x_len) -> np.ndarray:
        """
        """
        y_len = self.CICE_dict["y_dim_length"]
        x_len = self.CICE_dict["x_dim_length"]
        self.logger.debug(f"input shape to Ta spatial averaging: {B_component.shape}")
        # Compute average of 4 neighboring cells
        v00 = B_component[:, :-1, :-1]
        v01 = B_component[:, :-1, 1:]
        v10 = B_component[:, 1:, :-1]
        v11 = B_component[:, 1:, 1:]
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

    def B2Ta(self, uB, vB, y_len, x_len):
        uT = self.Bavg_methA(uB, y_len, x_len)
        vT = self.Bavg_methA(vB, y_len, x_len)
        sT = np.sqrt(uT**2 + vT**2)
        return sT

    def Bavg_methB(self, B_component, y_len, x_len, wrap_x):
        """
        Average corner-staggered (B-grid) components onto T points by simple 2x2 mean with no-slip 
        at land corners (NaNs -> 0) and denom=4. Shapes: uB,vB: (t, y, x) on corners; 
        return: (uT, y_len, x_len)
        """
        v00 = np.nan_to_num(B_component[:, :-1, :-1], nan=0.0)
        v01 = np.nan_to_num(B_component[:, :-1,  1:], nan=0.0)
        v10 = np.nan_to_num(B_component[:,  1:, :-1], nan=0.0)
        v11 = np.nan_to_num(B_component[:,  1:,  1:], nan=0.0)
        out = 0.25 * (v00 + v01 + v10 + v11)
        pad_y = max(y_len - out.shape[1], 0)
        pad_x = max(x_len - out.shape[2], 0)
        if pad_y or pad_x:
            out = np.pad(out, ((0, 0), (0, pad_y), (0, pad_x)), constant_values=np.nan)
        if wrap_x and out.shape[2] > 1:
            out[:, :, -1] = out[:, :, 0]
        return out[:, :y_len, :x_len].astype(np.float32)

    def B2Tb(self, uB, vB, y_len, x_len, wrap_x):
        uT = self.Bavg_methB(uB, y_len, x_len, wrap_x)
        vT = self.Bavg_methB(vB, y_len, x_len, wrap_x)
        sT = np.sqrt(uT**2 + vT**2)
        return sT#, uT, vT

    def B2Tx(self, uB_da, vB_da):
        """Bilinear components B->T, then hypot."""
        self._ensure_reG_defined()
        uTx = self.reG(uB_da.fillna(0.0))  # no-slip: NaNs->0
        vTx = self.reG(vB_da.fillna(0.0))
        return xr.apply_ufunc(np.hypot, uTx, vTx, dask="parallelized", output_dtypes=[np.float32])

    def _binary_days(self, mask, 
                     win     : int, 
                     min_days: int,
                     centered: bool = True):
        """
        NumPy implementation of the binary-day persistence filter.

        Parameters
        ----------
        mask : np.ndarray
            0/1 (or bool) array shaped (time, nj, ni).
        win : int
            Sliding window length in days.
        min_days : int
            Minimum number of '1' days inside the window.
        centered : bool
            If True, window is centered on t; else trailing (window ends at t).

        Returns
        -------
        np.ndarray (uint8)
            0/1 array (time, nj, ni) after the windowed threshold.
        """
        m = np.asarray(mask)
        if m.dtype != np.uint8:
            # convert {NaN,0,1} → {0,0,1}
            m = np.where(np.isfinite(m), m, 0).astype(np.uint8)
        win = int(win); min_days = int(min_days)
        if min_days > win:
            min_days = win
        T = m.shape[0]
        wv = sliding_window_view(m, window_shape=win, axis=0)  # (T-win+1, ny, nx, win)
        meets = (wv.sum(axis=-1) >= min_days).astype(np.uint8)
        if centered:
            pad_before = (win - 1) // 2
            pad_after  = (win - 1) - pad_before
        else:  # trailing window
            pad_before, pad_after = win - 1, 0
        return np.pad(meets, ((pad_before, pad_after), (0, 0), (0, 0)), constant_values=0)

    def classify_binary_days_fast_ice(self, FI_mask,
                                      bin_win_days : int  = None,
                                      bin_min_days : int  = None,
                                      time_dim     : str  = None,
                                      centered     : bool = True):
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
        if isinstance(FI_mask, xr.Dataset):
            if "FI_mask" not in FI_mask:
                raise ValueError("Dataset must contain variable 'FI_mask'.")
            da = FI_mask["FI_mask"]
        else:
            da = FI_mask
        # infer spatial dims
        spatial_dims = [d for d in da.dims if d != time_dim]
        if len(spatial_dims) != 2:
            raise ValueError(f"Expected two spatial dimensions, got: {spatial_dims}")
        dim_y, dim_x = spatial_dims
        # make sure time is a single chunk to avoid window breaks at chunk edges
        dai = da.fillna(0).astype("i1").chunk({time_dim: -1})
        self.logger.debug(f"   FI_mask dims before rolling: {dai.dims}, shape={tuple(dai.shape)}")
        self.logger.info(f"   {bin_min_days:d} out of {bin_win_days:d} days must be immobile sea ice")
        da_FI_bin = xr.apply_ufunc(self._binary_days, dai,
                                   input_core_dims  = [[time_dim, dim_y, dim_x]],
                                   output_core_dims = [[time_dim, dim_y, dim_x]],
                                   dask             = "parallelized",
                                   output_dtypes    = [np.uint8],
                                   kwargs           = {"win"      : bin_win_days,
                                                       "min_days" : bin_min_days, 
                                                       "centered" : centered})
        out = da_FI_bin.rename("FI_mask").to_dataset()
        # record parameters on the variable (not just Dataset attrs)
        out["FI_mask"].attrs.update({"long_name"    : "binary-day fast ice mask",
                                     "bin_win_days" : bin_win_days,
                                     "bin_min_days" : bin_min_days,
                                     "window_type"  : "centered" if centered else "trailing",})
        return out.chunk(self.CICE_dict["FI_chunks"])

    def _compose_T_speed(self,
                         uvel  : xr.DataArray,
                         vvel  : xr.DataArray, 
                         label : str  = "") -> xr.DataArray:
        """
        """
        y_len             = self.CICE_dict["y_dim_length"]
        x_len             = self.CICE_dict["x_dim_length"]
        wrap_x            = self.CICE_dict.get("wrap_x", True)
        members           = []
        if "Ta" in self.BT_composite_grids:
            if hasattr(self, "B2Ta"):
                ispd    = self.B2Ta(uvel, vvel, y_len, x_len).astype(np.float32)
                ispd_Ta = xr.DataArray(ispd, dims=(self.CICE_dict["three_dims"]), name="ispd_Ta")
                members.append(ispd_Ta)
            else:
                self.logger.warning("Ta requested in BT_composite_grids but B2Ta not defined. Skipping Ta.")
        if "Tb" in self.BT_composite_grids:
            if hasattr(self, "B2Tb"):
                ispd    = self.B2Tb(uvel, vvel, y_len, x_len, wrap_x=wrap_x).astype(np.float32)
                ispd_Tb = xr.DataArray(ispd, dims=(self.CICE_dict["three_dims"]), name="ispd_Tb")
                members.append(ispd_Tb)
            else:
                self.logger.warning("Tb requested in BT_composite_grids but B2Tb not defined. Skipping Tb.")
        if "Tx" in self.BT_composite_grids:
            if hasattr(self, "B2Tx"):
                ispd_Tx = self.B2Tx(uvel, vvel).astype(np.float32).rename("ispd_Tx")
                members.append(ispd_Tx)
            else:
                self.logger.warning("Tx requested in BT_composite_grids but B2Tx not defined. Skipping Tx.")
        if len(members) == 1:
            return members[0]
        return xr.concat(members, dim="__concat_dim__").mean("__concat_dim__", skipna=True).astype(np.float32)

    def _apply_thresholds(self, ispd_T: xr.DataArray, aice: xr.DataArray, label: str) -> xr.Dataset:
        """
        Apply hemisphere slice and thresholds, package as a Dataset with proper chunking.
        """
        ispd_hem = self.slice_hemisphere(ispd_T)
        aice_hem = self.slice_hemisphere(aice)
        self.logger.info(f"   thresholding ({label}): aice>{self.icon_thresh:0.2f} and 0<speed<={self.ispd_thresh:.3e}")
        BT_mask = (aice_hem > self.icon_thresh) & (ispd_hem > 0) & (ispd_hem <= self.ispd_thresh)
        FI      = BT_mask.rename("FI_mask").to_dataset().chunk(self.CICE_dict["FI_chunks"]).compute()
        # Drop any spatial coords that duplicate dims
        for dim in self.CICE_dict['spatial_dims']:
            if dim in FI.coords:
                FI = FI.drop_vars(dim)
        return FI

    def _compute_fast_ice_mask(self,
                               uvel    : xr.DataArray,
                               vvel    : xr.DataArray,
                               aice    : xr.DataArray,
                               dt0_str : str,
                               dtN_str : str,
                               label   : str = "") -> xr.Dataset:
        """
        Build the daily fast-ice mask on the analysis hemisphere from B-grid components.

        Parameters
        ----------
        uvel, vvel : xarray.DataArray
            Corner-staggered velocity components (t, y, x) on the B grid.
        aice : xarray.DataArray
            Sea-ice concentration (t, y, x) on the T grid.
        dt0_str, dtN_str : str
            Requested start/end dates (ISO yyyy-mm-dd). Used for logs and slicing.
        label : str
            Tag for logging (e.g., "daily", "rolling-mean").

        Returns
        -------
        xarray.Dataset
            Dataset with boolean "FI_mask" over the selected hemisphere and time span.
        """
        self.logger.info(f"FAST ICE MASK [{label}]:")
        # Compose T-point speed correctly from components (Tb, optionally TxC)
        ispd_T = self._compose_T_speed(uvel, vvel, label=label).compute()
        # Thresholds
        FI = self._apply_thresholds(ispd_T, aice, label=label)
        # Final time crop (safety)
        return FI.sel(time=slice(pd.to_datetime(dt0_str), pd.to_datetime(dtN_str)))

    def classify_fast_ice(self,
                          dt0_str               : str  = None,
                          dtN_str               : str  = None,
                          bin_win_days          : int  = None,
                          bin_min_days          : int  = None,
                          enable_rolling_output : bool = False,
                          roll_win_days         : int  = None,
                          time_dim_name         : str  = None):
        """
        End-to-end fast-ice classification over a date range, corrected to use
        component-safe B->T interpolation (no-slip zeros, denom=4) for speed.
        """
        time_dim      = time_dim_name or self.CICE_dict["time_dim"]
        bin_win_days  = bin_win_days  or self.bin_win_days
        bin_min_days  = bin_min_days  or self.bin_min_days
        roll_win_days = roll_win_days or self.mean_period
        self.load_bgrid()
        self.define_reG_weights()
        # Extend the requested time range for rolling/binary windows
        ext_per = (self.mean_period // 2) if enable_rolling_output else (self.bin_win_days // 2)
        dt0_ext = pd.to_datetime(dt0_str) - timedelta(days=ext_per)
        dtN_ext = pd.to_datetime(dtN_str) + timedelta(days=ext_per)
        # Load required variables
        self.logger.info(f"loading model data between {dt0_ext.strftime('%Y-%m-%d')} and {dtN_ext.strftime('%Y-%m-%d')}")
        CICE_all = self.load_cice_zarr(slice_hem=False,
                                    variables=['aice', 'uvel', 'vvel'],
                                    dt0_str=dt0_ext.strftime('%Y-%m-%d'),
                                    dtN_str=dtN_ext.strftime('%Y-%m-%d'))
        # Clamp to available range
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
        # Define internal calc window
        self.define_datetime_vars(dt0_str=dt0_load, dtN_str=dtN_load)
        # Drop coord bloat early
        aice = CICE_all['aice'].drop_vars(self.CICE_dict["drop_coords"])
        uvel = CICE_all['uvel'].drop_vars(self.CICE_dict["drop_coords"])
        vvel = CICE_all['vvel'].drop_vars(self.CICE_dict["drop_coords"])
        # Daily mask (component-safe interpolation inside)
        FI_dly = self._compute_fast_ice_mask(uvel, vvel, aice, dt0_str, dtN_str, label="daily")
        FI_dly = FI_dly.sel(time=slice(dt0_str, dtN_str))
        # Binary-day mask (on daily FI)
        self.logger.info("classifying binary-day fast ice mask")
        FI_bin = self.classify_binary_days_fast_ice(FI_dly["FI_mask"],
                                                    bin_win_days = bin_win_days,
                                                    bin_min_days = bin_min_days,
                                                    time_dim     = time_dim).sel(time=slice(dt0_str, dtN_str))
        # Optional rolling classification on T-grid speed
        FI_roll = None
        if enable_rolling_output:
            self.logger.info(f"applying rolling mean on T-grid speed (period={roll_win_days})")
            ispd_T      = self._compose_T_speed(uvel, vvel, label="rolling-pre")  # lazy
            ispd_T_roll = self._rolling_mean(ispd_T, mean_period=roll_win_days)
            ispd_T_da   = xr.DataArray(ispd_T_roll, dims=(self.CICE_dict["three_dims"]), name="rolling-mean")
            FI_roll     = self._apply_thresholds(ispd_T_da, aice, label="rolling-mean")
            FI_roll     = FI_roll.sel(time=slice(dt0_str, dtN_str))
        return FI_dly, FI_bin, FI_roll
