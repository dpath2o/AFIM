"""
SeaIceClassification
====================

Class utilities for classifying Antarctic fast ice using sea-ice concentration
(`aice`) and speed thresholds on multiple grids:

- B  : native CICE B-grid (cell centers) magnitude computed from (uvel, vvel)
- Ta : 2×2 B-grid box average (area-like) magnitude
- Tb : 2×2 B-grid box average (grid-like) magnitude
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
- self.B2T_type: some combination of {"B","C","Ta","Tb","Tc","Tx"} to average into composite
- self.G_t, self.G_u: dicts providing T- and U-grid lon/lat arrays
- self.dt_range: pandas.DatetimeIndex used by np3d_to_xr3d

Methods (provided elsewhere in your stack):
- self.load_cice_grid()
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
      after composing `ispd` from one or more of {B, Ta, Tb, Tx}.
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

    def _rolling_mean(self, mag, mean_period: int | None = None):
        """
        Centered rolling mean along time.
        - If `mag` is an xarray.DataArray, keep coords/dims and use xarray rolling (Dask-parallel).
        - If `mag` is a NumPy array, fall back to scipy.uniform_filter1d.
        """
        mean_period = int(mean_period or self.mean_period)

        if isinstance(mag, xr.DataArray):
            time_dim = self.CICE_dict["time_dim"]
            # ensure time chunks aren't larger than the window (helps parallel rolling)
            tlen = mag.sizes.get(time_dim, 0)
            if tlen > 0:
                tch = max(8, min(mean_period, tlen))
                if mag.chunks is None or mag.chunksizes[time_dim][0] > tch:
                    mag = mag.chunk({time_dim: tch})

            # centered rolling mean; require full window to avoid edge bias
            out = (
                mag.rolling({time_dim: mean_period}, center=True, min_periods=mean_period)
                .mean()
                .astype(np.float32)
            )
            return out

        else:
            # NumPy fallback
            from scipy.ndimage import uniform_filter1d
            return uniform_filter1d(np.asarray(mag), size=mean_period, axis=0, mode="nearest").astype(np.float32)

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

    def B2Ta(self, uB, vB, y_len, x_len):
        uT = self.Bavg_methA(uB, y_len, x_len)
        vT = self.Bavg_methA(vB, y_len, x_len)
        sT = np.sqrt(uT**2 + vT**2)
        return sT

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

    def classify_binary_days_mask(self, I_mask: xr.DataArray,
                                 mask_name   : str,
                                 bin_win_days: int | None = None,
                                 bin_min_days: int | None = None,
                                 time_dim    : str | None = None,
                                 centered    : bool = True) -> xr.Dataset:
        """
        Generic binary-day persistence: in each window of length W, mark ice as 1 if
        at least M days are True.

        Parameters
        ----------
        I_mask : xr.DataArray
            Boolean (or 0/1) mask with dimensions (time, y, x).
        mask_name : str
            Name for the returned mask variable (e.g., 'FI_mask', 'PI_mask').
        bin_win_days, bin_min_days : int, optional
            Window length and minimum count. Defaults to `self.bin_win_days` and
            `self.bin_min_days` when omitted.
        time_dim : str, optional
            Name of time dimension. Defaults to `self.CICE_dict["time_dim"]`.
        centered : bool, default True
            Center the rolling window.

        Returns
        -------
        xr.Dataset
            Dataset containing `{mask_name}` as uint8 (0/1).
        """
        time_dim = time_dim or self.CICE_dict["time_dim"]
        W        = int(bin_win_days or self.bin_win_days)
        M        = int(bin_min_days or self.bin_min_days)

        if I_mask.sizes.get(time_dim, 0) == 0:
            self.logger.warning(f"   {mask_name} has zero time length for this slice; returning empty dataset.")
            return I_mask.astype(np.uint8).rename(mask_name).to_dataset()

        tlen = I_mask.sizes[time_dim]
        tch  = max(8, min(W, tlen))
        if I_mask.chunks is None or I_mask.chunksizes[time_dim][0] in (0, None) or I_mask.chunksizes[time_dim][0] > tch:
            I_mask = I_mask.chunk({time_dim: tch})

        counts = I_mask.astype("uint8").rolling({time_dim: W}, center=centered, min_periods=M).sum()
        da_bin = (counts >= M).astype("uint8").rename(mask_name)
        return da_bin.to_dataset()

    def classify_binary_days_fast_ice(self, FI_mask: xr.DataArray,
                                      bin_win_days : int | None = None,
                                      bin_min_days : int | None = None,
                                      time_dim     : str | None = None,
                                      centered     : bool = True) -> xr.Dataset:
        """
        Backwards-compatible wrapper for binary-day persistence on fast ice.
        Returns a Dataset with variable name 'FI_mask'.
        """
        return self.classify_binary_days_mask(FI_mask,
                                             mask_name    = "FI_mask",
                                             bin_win_days = bin_win_days,
                                             bin_min_days = bin_min_days,
                                             time_dim     = time_dim,
                                             centered     = centered)

    def classify_binary_days_pack_ice(self, PI_mask: xr.DataArray,
                                      bin_win_days : int | None = None,
                                      bin_min_days : int | None = None,
                                      time_dim     : str | None = None,
                                      centered     : bool = True) -> xr.Dataset:
        """
        Binary-day persistence on pack ice (converse of fast ice).
        Returns a Dataset with variable name 'PI_mask'.
        """
        return self.classify_binary_days_mask(PI_mask,
                                             mask_name    = "PI_mask",
                                             bin_win_days = bin_win_days,
                                             bin_min_days = bin_min_days,
                                             time_dim     = time_dim,
                                             centered     = centered)

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
        if "Ta" in self.B2T_type:
            if hasattr(self, "B2Ta"):
                ispd    = self.B2Ta(uvel, vvel, y_len, x_len).astype(np.float32)
                ispd_Ta = xr.DataArray(ispd, dims=(self.CICE_dict["three_dims"]), name="ispd_Ta")
                members.append(ispd_Ta)
            else:
                self.logger.warning("Ta requested in B2T_type but B2Ta not defined. Skipping Ta.")
        if "Tb" in self.B2T_type:
            if hasattr(self, "B2Tb"):
                ispd    = self.B2Tb(uvel, vvel, y_len, x_len, wrap_x=wrap_x).astype(np.float32)
                ispd_Tb = xr.DataArray(ispd, dims=(self.CICE_dict["three_dims"]), name="ispd_Tb")
                members.append(ispd_Tb)
            else:
                self.logger.warning("Tb requested in B2T_type but B2Tb not defined. Skipping Tb.")
        if "Tx" in self.B2T_type:
            if hasattr(self, "B2Tx"):
                ispd_Tx = self.B2Tx(uvel, vvel).astype(np.float32).rename("ispd_Tx")
                members.append(ispd_Tx)
            else:
                self.logger.warning("Tx requested in B2T_type but B2Tx not defined. Skipping Tx.")
        if len(members) == 1:
            return members[0]
        return xr.concat(members, dim="__concat_dim__").mean("__concat_dim__", skipna=True).astype(np.float32)

    def _apply_thresholds(self,
                          ispd_T   : xr.DataArray,
                          aice     : xr.DataArray,
                          label    : str,
                          ice_type : str = "FI",
                          mask_name: str | None = None) -> xr.Dataset:
        """
        Apply hemisphere slice and thresholds, package as a Dataset with proper chunking.

        Parameters
        ----------
        ispd_T : xr.DataArray
            T-grid ice-speed magnitude (time, y, x).
        aice : xr.DataArray
            Sea-ice concentration (time, y, x) on the same grid/dims as `ispd_T`.
        label : str
            Human-readable label for logging (e.g., 'daily', 'rolling-mean').
        ice_type : {'FI','PI'}, default 'FI'
            Ice class to compute:
              - FI: fast ice (aice > icon_thresh AND 0 < speed <= ispd_thresh)
              - PI: pack ice (aice > icon_thresh AND speed > ispd_thresh)
        mask_name : str, optional
            Output variable name. Defaults to '{ice_type}_mask'.

        Returns
        -------
        xr.Dataset
            Dataset with one variable: `{mask_name}`.
        """
        ice_type_u = (ice_type or "FI").upper()
        mask_name  = mask_name or f"{ice_type_u}_mask"

        ispd_hem = self.slice_hemisphere(ispd_T)
        aice_hem = self.slice_hemisphere(aice)

        if ice_type_u == "FI":
            self.logger.info(
                f"   thresholding ({label}): aice>{self.icon_thresh:0.2f} and 0<speed<={self.ispd_thresh:.3e}"
            )
            I_mask = (aice_hem > self.icon_thresh) & (ispd_hem > 0) & (ispd_hem <= self.ispd_thresh)
        elif ice_type_u == "PI":
            self.logger.info(
                f"   thresholding ({label}): aice>{self.icon_thresh:0.2f} and speed>{self.ispd_thresh:.3e}"
            )
            I_mask = (aice_hem > self.icon_thresh) & (ispd_hem > self.ispd_thresh)
        else:
            raise ValueError(f"Unsupported ice_type='{ice_type}'. Expected 'FI' or 'PI'.")

        DS = I_mask.rename(mask_name).to_dataset().chunk(self.CICE_dict["FI_chunks"]).compute()

        # Drop any spatial coords that duplicate dims (keeps output compact/consistent)
        for dim in self.CICE_dict.get('spatial_dims', []):
            if dim in DS.coords:
                DS = DS.drop_vars(dim)
        return DS
    
    def _with_T_coords(self, da: xr.DataArray, name: str | None = None) -> xr.DataArray:
        """
        Attach T-grid lon/lat (and optionally *_b) coords from self.G_t to `da`.
        Assumes `da` has dims (time, y, x) matching self.CICE_dict["three_dims"].
        """
        ydim, xdim = self.CICE_dict["spatial_dims"]
        da         = da.rename(name or da.name)
        da         = da.assign_coords({"lon"  : ((ydim, xdim), self.G_t["lon"].data),
                                       "lat"  : ((ydim, xdim), self.G_t["lat"].data),
                                       "angle": ((ydim, xdim), self.G_t["angle"].data),
                                       "area" : ((ydim, xdim), self.G_t["area"].data)})
        return da

    def classify_fast_ice(self,
                        dt0_str               : str  = None,
                        dtN_str               : str  = None,
                        bin_win_days          : int  = None,
                        bin_min_days          : int  = None,
                        enable_rolling_output : bool = False,
                        roll_win_days         : int  = None,
                        time_dim_name         : str  = None):
        time_dim      = time_dim_name or self.CICE_dict["time_dim"]
        bin_win_days  = int(bin_win_days  or self.bin_win_days)
        bin_min_days  = int(bin_min_days  or self.bin_min_days)
        roll_win_days = int(roll_win_days or self.mean_period)
        self.load_cice_grid()
        self.define_reG_weights()
        ispd_out = {'daily' : None,
                    'rolly' : None}
        # extend enough time for BOTH methods
        ext_per = max(bin_win_days // 2, roll_win_days // 2 if enable_rolling_output else 0)
        dt0_ext = pd.to_datetime(dt0_str) - pd.Timedelta(days=ext_per)
        dtN_ext = pd.to_datetime(dtN_str) + pd.Timedelta(days=ext_per)
        # Load (extended) data
        self.logger.info(f"loading model data between {dt0_ext:%Y-%m-%d} and {dtN_ext:%Y-%m-%d}")
        CICE_all = self.load_cice_zarr(slice_hem = False,
                                       variables = ['aice','uvel','vvel'],
                                       dt0_str   = f"{dt0_ext:%Y-%m-%d}",
                                       dtN_str   = f"{dtN_ext:%Y-%m-%d}")
        # Compact coords & set friendly chunks
        aice   = CICE_all['aice'].drop_vars(self.CICE_dict["drop_coords"]).astype(np.float32)
        uvel   = CICE_all['uvel'].drop_vars(self.CICE_dict["drop_coords"]).astype(np.float32)
        vvel   = CICE_all['vvel'].drop_vars(self.CICE_dict["drop_coords"]).astype(np.float32)
        yck    = self.CICE_dict.get("y_chunk", 540)
        xck    = self.CICE_dict.get("x_chunk", 1440)
        tch    = max(8, min(bin_win_days, roll_win_days))  # small time chunks help rolling parallelism
        chunks = {self.CICE_dict["time_dim"]: tch, self.CICE_dict["y_dim"]: yck, self.CICE_dict["x_dim"]: xck}
        aice, uvel, vvel = aice.chunk(chunks), uvel.chunk(chunks), vvel.chunk(chunks)
        # Build T-grid speed ONCE over the extended span
        ispd_T = self._compose_T_speed(uvel, vvel, label="base").persist()
        ispd_T = self._with_T_coords(ispd_T, name="ispd_T")
        ispd_out['daily'] = ispd_T
        # DAILY mask on the EXTENDED range, THEN crop
        FI_dly_ext = self._apply_thresholds(ispd_T, aice, label="daily")
        FI_dly     = FI_dly_ext.sel(time=slice(dt0_str, dtN_str))
        # BINARY-DAYS on the EXTENDED DAILY mask, THEN crop  <-- this fixes your all-zero issue
        self.logger.info("classifying binary-day fast ice mask")
        FI_bin_ext = self.classify_binary_days_fast_ice(FI_dly_ext["FI_mask"],  # use extended
                                                        bin_win_days = bin_win_days,
                                                        bin_min_days = bin_min_days,
                                                        time_dim     = time_dim,
                                                        centered     = True)
        FI_bin = FI_bin_ext.sel(time=slice(dt0_str, dtN_str))
        # Optional ROLLING on the extended speed, THEN crop (unchanged in spirit)
        FI_roll = None
        if enable_rolling_output:
            self.logger.info(f"applying rolling mean on T-grid speed (period={roll_win_days})")
            ispd_T_roll = self._rolling_mean(ispd_T, mean_period=roll_win_days)
            ispd_T_roll = self._with_T_coords(ispd_T_roll, name="ispd_T_roll")
            FI_roll_ext = self._apply_thresholds(ispd_T_roll, aice, label="rolling-mean")
            FI_roll     = FI_roll_ext.sel(time=slice(dt0_str, dtN_str))
            ispd_out['rolly'] = ispd_T_roll
        return FI_dly, FI_bin, FI_roll, ispd_out

    def classify_pack_ice(self,
                          dt0_str               : str  = None,
                          dtN_str               : str  = None,
                          bin_win_days          : int  = None,
                          bin_min_days          : int  = None,
                          enable_rolling_output : bool = False,
                          roll_win_days         : int  = None,
                          time_dim_name         : str  = None):
        """
        Classify pack ice (PI) as the converse of fast ice (FI) under the same
        concentration threshold:

            PI := (aice > icon_thresh) AND (speed > ispd_thresh)

        Returns daily, binary-day (persistence), and optional rolling-mean PI masks.

        Notes
        -----
        - Pack ice is defined only where `aice` exceeds the concentration cutoff.
        - No special handling is applied for the zero-speed case; by definition
          `speed > ispd_thresh` implies `speed > 0`.
        """
        time_dim      = time_dim_name or self.CICE_dict["time_dim"]
        bin_win_days  = int(bin_win_days  or self.bin_win_days)
        bin_min_days  = int(bin_min_days  or self.bin_min_days)
        roll_win_days = int(roll_win_days or self.mean_period)

        self.load_cice_grid()
        self.define_reG_weights()

        ispd_out = {'daily' : None,
                    'rolly' : None}

        # extend enough time for BOTH methods
        ext_per = max(bin_win_days // 2, roll_win_days // 2 if enable_rolling_output else 0)
        dt0_ext = pd.to_datetime(dt0_str) - pd.Timedelta(days=ext_per)
        dtN_ext = pd.to_datetime(dtN_str) + pd.Timedelta(days=ext_per)

        # Load (extended) data
        self.logger.info(f"loading model data between {dt0_ext:%Y-%m-%d} and {dtN_ext:%Y-%m-%d}")
        CICE_all = self.load_cice_zarr(slice_hem = False,
                                       variables = ['aice','uvel','vvel'],
                                       dt0_str   = f"{dt0_ext:%Y-%m-%d}",
                                       dtN_str   = f"{dtN_ext:%Y-%m-%d}")

        # Compact coords & set friendly chunks
        aice   = CICE_all['aice'].drop_vars(self.CICE_dict["drop_coords"]).astype(np.float32)
        uvel   = CICE_all['uvel'].drop_vars(self.CICE_dict["drop_coords"]).astype(np.float32)
        vvel   = CICE_all['vvel'].drop_vars(self.CICE_dict["drop_coords"]).astype(np.float32)

        yck    = self.CICE_dict.get("y_chunk", 540)
        xck    = self.CICE_dict.get("x_chunk", 1440)
        tch    = max(8, min(bin_win_days, roll_win_days))
        chunks = {self.CICE_dict["time_dim"]: tch, self.CICE_dict["y_dim"]: yck, self.CICE_dict["x_dim"]: xck}
        aice, uvel, vvel = aice.chunk(chunks), uvel.chunk(chunks), vvel.chunk(chunks)

        # Build T-grid speed ONCE over the extended span
        ispd_T = self._compose_T_speed(uvel, vvel, label="base").persist()
        ispd_T = self._with_T_coords(ispd_T, name="ispd_T")
        ispd_out['daily'] = ispd_T

        # DAILY mask on the EXTENDED range, THEN crop
        PI_dly_ext = self._apply_thresholds(ispd_T, aice, label="daily", ice_type="PI", mask_name="PI_mask")
        PI_dly     = PI_dly_ext.sel(time=slice(dt0_str, dtN_str))

        # BINARY-DAYS on the EXTENDED DAILY mask, THEN crop
        self.logger.info("classifying binary-day pack ice mask")
        PI_bin_ext = self.classify_binary_days_pack_ice(PI_dly_ext["PI_mask"],
                                                        bin_win_days = bin_win_days,
                                                        bin_min_days = bin_min_days,
                                                        time_dim     = time_dim,
                                                        centered     = True)
        PI_bin = PI_bin_ext.sel(time=slice(dt0_str, dtN_str))

        # Optional ROLLING on the extended speed, THEN crop
        PI_roll = None
        if enable_rolling_output:
            self.logger.info(f"applying rolling mean on T-grid speed (period={roll_win_days})")
            ispd_T_roll = self._rolling_mean(ispd_T, mean_period=roll_win_days)
            ispd_T_roll = self._with_T_coords(ispd_T_roll, name="ispd_T_roll")
            PI_roll_ext = self._apply_thresholds(ispd_T_roll, aice, label="rolling-mean", ice_type="PI", mask_name="PI_mask")
            PI_roll     = PI_roll_ext.sel(time=slice(dt0_str, dtN_str))
            ispd_out['rolly'] = ispd_T_roll

        return PI_dly, PI_bin, PI_roll, ispd_out

    def classify_sea_ice(self,
                         dt0_str       : str = None,
                         dtN_str       : str = None,
                         time_dim_name : str = None):
        """
        Classify total sea ice (SI) with **no speed thresholding**:

            SI := (aice > icon_thresh)

        Notes
        -----
        - SI is inclusive of both fast ice (FI) and pack ice (PI).
        - There is no binary-day (persistence) product for SI because it is not
          derived from speed-thresholded immobility.
        """
        time_dim = time_dim_name or self.CICE_dict["time_dim"]

        # Load aice only (no need to compute speed)
        self.logger.info(f"loading model data between {dt0_str} and {dtN_str} (aice only)")
        CICE_all = self.load_cice_zarr(slice_hem = False,
                                       variables = ['aice'],
                                       dt0_str   = dt0_str,
                                       dtN_str   = dtN_str)

        aice = CICE_all['aice'].drop_vars(self.CICE_dict["drop_coords"]).astype(np.float32)

        # Chunk in a similar manner to FI/PI products
        yck    = self.CICE_dict.get("y_chunk", 540)
        xck    = self.CICE_dict.get("x_chunk", 1440)
        tlen   = aice.sizes.get(time_dim, 0)
        tch    = max(8, min(31, tlen)) if tlen else 8
        chunks = {self.CICE_dict["time_dim"]: tch, self.CICE_dict["y_dim"]: yck, self.CICE_dict["x_dim"]: xck}
        aice   = aice.chunk(chunks)

        aice_hem = self.slice_hemisphere(aice)
        self.logger.info(f"   thresholding (sea-ice): aice>{self.icon_thresh:0.2f} (no speed constraint)")
        SI_mask  = (aice_hem > self.icon_thresh)

        SI_dly = SI_mask.rename("SI_mask").to_dataset().chunk(self.CICE_dict["FI_chunks"]).compute()

        for dim in self.CICE_dict.get('spatial_dims', []):
            if dim in SI_dly.coords:
                SI_dly = SI_dly.drop_vars(dim)

        # Maintain the same return signature as classify_fast_ice / classify_pack_ice
        return SI_dly, None, None, {}
