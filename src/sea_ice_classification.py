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

    def __init__(self, sim_name=None, logger=None, **kwargs):
        self.sim_name = sim_name
        self.logger = logger or logging.getLogger(__name__)
        for k, v in kwargs.items():
            setattr(self, k, v)

    def _rolling_mean(self, np_mag: np.ndarray) -> np.ndarray:
        from scipy.ndimage import uniform_filter1d
        return uniform_filter1d(np_mag, size=self.mean_period, axis=0, mode="nearest")

    def B_mag_to_Ta_mag_numpy(self, B_mag_np: np.ndarray) -> np.ndarray:
        """
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
        xr_B = self.np3d_to_xr3d(np_B, tgrid=False)
        self.logger.debug(f"input shape to Tx spatial averaging: {np.shape(xr_B)}")
        return self.reG(xr_B)

    def _binary_days(self, mask):
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

    def classify_fast_ice(self, enable_rolling_output=False, dt0_str=None, dtN_str=None):
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
        FI_bin = self.classify_binary_days_fast_ice(FI_dly["FI_mask"]).sel(time=slice(dt0_str, dtN_str))
        FI_dly = FI_dly.sel(time=slice(dt0_str, dtN_str))
        # Optional: prepare rolling classification
        if enable_rolling_output:
            self.logger.info(f"applying rolling mean (period={self.mean_period})")
            ispd_B_roll = self._rolling_mean(ispd_B)
            FI_roll = self._compute_fast_ice_mask(ispd_B_roll, aice, dt0_str, dtN_str, label="rolling-mean").sel(time=slice(dt0_str, dtN_str))
        else:
            FI_roll = None
        return FI_dly, FI_bin, FI_roll