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
- self.BorC2T_type: some combination of {"B","C","Ta","Tb","Tc","Tx"} to average into composite
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
        Compute a centered rolling mean along the model time dimension.

        This helper supports both xarray-backed (optionally Dask-chunked) arrays and
        plain NumPy arrays:

        - For `xarray.DataArray`, uses `DataArray.rolling(...).mean()` which can run
        in parallel under Dask. The time chunk size is adjusted to be no larger
        than the rolling window to improve rolling performance and avoid large
        rechunk overheads.
        - For NumPy input, uses `scipy.ndimage.uniform_filter1d` as a fast fallback.

        Parameters
        ----------
        mag : xarray.DataArray or array-like
            Magnitude-like array with time on axis 0 (NumPy) or on the dimension
            named by `self.CICE_dict["time_dim"]` (xarray). Typically shape
            (time, y, x).
        mean_period : int, optional
            Window length (number of timesteps/days) for the rolling mean. If None,
            defaults to `self.mean_period`.

        Returns
        -------
        out : xarray.DataArray or numpy.ndarray
            Rolling-mean result cast to float32. For xarray input, coordinates and
            dimension names are preserved.

        Notes
        -----
        - The xarray path uses `center=True` and `min_periods=mean_period`, which
        yields NaNs at the edges (no partial-window averaging).
        - The NumPy fallback uses `mode="nearest"`, which does not introduce NaNs at
        edges but will bias edge values relative to the strict full-window xarray
        behaviour.
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
            out = (mag.rolling({time_dim: mean_period},
                                center      = True, 
                                min_periods = mean_period).mean().astype(np.float32))
            return out
        else: # NumPy fallback
            from scipy.ndimage import uniform_filter1d
            return uniform_filter1d(np.asarray(mag), size=mean_period, axis=0, mode="nearest").astype(np.float32)

    def Bavg_methA(self, B_component: np.ndarray, y_len, x_len) -> np.ndarray:
        """
        Average a corner-staggered (B-grid) field onto T-cell centers via a 2x2 mean.

        This method computes the simple arithmetic mean of the four corner values
        surrounding each T-cell:
            avg = (v00 + v01 + v10 + v11) / 4

        It then pads the result back to the target T-grid shape and optionally
        wraps the x-direction by copying the first column into the last column.

        Parameters
        ----------
        B_component : numpy.ndarray
            Corner-staggered field with shape (time, y, x). Typically a velocity
            component (u or v) on B-grid corners.
        y_len : int
            Target T-grid y dimension length. (Note: current implementation reads
            `self.CICE_dict["y_dim_length"]` and ignores the passed value.)
        x_len : int
            Target T-grid x dimension length. (Note: current implementation reads
            `self.CICE_dict["x_dim_length"]` and ignores the passed value.)

        Returns
        -------
        avg_T : numpy.ndarray
            T-centered field with shape (time, y_len, x_len). Padding is filled with
            NaN. dtype is unchanged until returned (typically float).

        Notes
        -----
        - This method does not apply a no-slip treatment: NaNs in the input propagate
        into the 2x2 average.
        - X wrapping is implemented by setting the last x column equal to the first,
        which is appropriate for cyclic/periodic global grids.
        - If the input has insufficient halo/extent for the 2x2 operation, the output
        will contain NaNs introduced by padding.
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
        Average a corner-staggered (B-grid) field onto T-cell centers using a 2x2 mean
        with no-slip handling at land corners (NaNs treated as zero).

        The method forms the mean of the four corner values surrounding each T-cell,
        but first converts NaNs to 0.0 to enforce a no-slip-style boundary at masked
        corners:
            out = 0.25 * (nan_to_num(v00) + nan_to_num(v01) + nan_to_num(v10) + nan_to_num(v11))

        Parameters
        ----------
        B_component : array-like
            Corner-staggered field with shape (time, y, x), e.g., uvel or vvel on a
            CICE B-grid.
        y_len : int
            Target T-grid y dimension length.
        x_len : int
            Target T-grid x dimension length.
        wrap_x : bool
            If True, enforce cyclic boundary in x by copying the first column to the
            last column after padding.

        Returns
        -------
        out_T : numpy.ndarray
            T-centered field with shape (time, y_len, x_len), dtype float32.

        Notes
        -----
        - NaN-to-zero treatment is intentional and corresponds to a no-slip behaviour
        at land/corner masks; this differs from `Bavg_methA`, where NaNs propagate.
        - Padding (if needed) is filled with NaN after averaging.
        - For cyclic grids, the wrap is applied after padding and before trimming.
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
    
    def Cavg_meth(self, C_component, direction: str,
                  y_len       : int  | None = None,
                  x_len       : int  | None = None,
                  wrap_x      : bool | None = None,
                  nan_to_zero : bool = True):
        """
        Average a C-grid edge-centered field onto T-cell centers.

        This is the C-grid analogue of your B-grid corner-to-center averaging, but uses
        a 1D average along the staggered direction:

        - direction="x": average along x (typical for U-grid / east-west faces)
        - direction="y": average along y (typical for V-grid / north-south faces)

        Supports both xarray.DataArray (Dask-friendly) and NumPy arrays.

        Parameters
        ----------
        C_component : xarray.DataArray or array-like
            C-grid edge-centered field. Expected to have time as the first dimension
            and the last two dims as (y, x) in some order consistent with typical CICE
            output. For xarray input, dims are inferred from the last two dims.
        direction : {"x","y"}
            Stagger direction to average over:
            - "x": average adjacent points in x to map U-grid → T-grid
            - "y": average adjacent points in y to map V-grid → T-grid
        y_len, x_len : int, optional
            Target T-grid spatial sizes. Defaults to `self.CICE_dict["y_dim_length"]`
            and `self.CICE_dict["x_dim_length"]`.
        wrap_x : bool, optional
            Whether to enforce cyclic wrap in the x-direction (seam handling). Defaults
            to `self.CICE_dict.get("wrap_x", True)`.
        nan_to_zero : bool, default True
            If True, treat NaNs as 0.0 before averaging (no-slip-like treatment near
            land/masks). If False, NaNs propagate into the average.

        Returns
        -------
        out_T : same type as input
            T-centered field with shape (time, y_len, x_len), float32 for NumPy input.
            For xarray input, coords/dims are preserved as much as possible.

        Notes
        -----
        Size handling:
        - If the staggered dimension has length (target+1), a straightforward adjacent
        slice-average is used and yields the target length.
        - If it has length (target), we average using roll/shift:
            * wrap_x=True uses roll (periodic)
            * wrap_x=False uses shift (non-periodic; last column becomes NaN)
        """
        direction = (direction or "").lower().strip()
        if direction not in ("x", "y"):
            raise ValueError("Cavg_meth: direction must be 'x' or 'y'")
        y_len = int(y_len or self.CICE_dict["y_dim_length"])
        x_len = int(x_len or self.CICE_dict["x_dim_length"])
        if wrap_x is None:
            wrap_x = bool(self.CICE_dict.get("wrap_x", True))
        # -------------------------
        # xarray / Dask-friendly
        # -------------------------
        if isinstance(C_component, xr.DataArray):
            da = C_component
            ydim, xdim = da.dims[-2], da.dims[-1]
            if nan_to_zero:
                da = da.fillna(0.0)
            if direction == "x":
                nx = int(da.sizes.get(xdim, 0))
                if nx == x_len + 1:
                    left  = da.isel({xdim: slice(0, x_len)})
                    right = da.isel({xdim: slice(1, x_len + 1)})
                    out = 0.5 * (left + right)
                else:
                    # length == x_len (or unknown): roll/shift based averaging
                    if wrap_x:
                        shifted = da.roll({xdim: -1}, roll_coords=False)
                    else:
                        shifted = da.shift({xdim: -1})
                    out = 0.5 * (da + shifted)
                    # ensure target length
                    if out.sizes.get(xdim, 0) > x_len:
                        out = out.isel({xdim: slice(0, x_len)})
            else:
                ny = int(da.sizes.get(ydim, 0))
                if ny == y_len + 1:
                    top = da.isel({ydim: slice(0, y_len)})
                    bot = da.isel({ydim: slice(1, y_len + 1)})
                    out = 0.5 * (top + bot)
                else:
                    shifted = da.shift({ydim: -1})  # no wrap in y by default
                    out = 0.5 * (da + shifted)
                    if out.sizes.get(ydim, 0) > y_len:
                        out = out.isel({ydim: slice(0, y_len)})

            # pad to target if needed (rare, but keeps behaviour consistent with Bavg)
            pad_y = max(y_len - int(out.sizes.get(ydim, 0)), 0)
            pad_x = max(x_len - int(out.sizes.get(xdim, 0)), 0)
            if pad_y or pad_x:
                out = out.pad({ydim: (0, pad_y), xdim: (0, pad_x)}, constant_values=np.nan)
            return out
        # -------------------------
        # NumPy fallback
        # -------------------------
        a = np.asarray(C_component, dtype="float64")
        if nan_to_zero:
            a = np.nan_to_num(a, nan=0.0)

        if a.ndim < 3:
            raise ValueError("Cavg_meth (NumPy): expected array with at least 3 dims (time,y,x)")

        T, ny, nx = a.shape[0], a.shape[1], a.shape[2]

        if direction == "x":
            if nx == x_len + 1:
                out = 0.5 * (a[:, :, :-1] + a[:, :, 1:])
            else:
                if wrap_x:
                    out = 0.5 * (a + np.roll(a, -1, axis=2))
                else:
                    core = 0.5 * (a[:, :, :-1] + a[:, :, 1:])
                    out = np.pad(core, ((0, 0), (0, 0), (0, 1)), constant_values=np.nan)
        else:
            if ny == y_len + 1:
                out = 0.5 * (a[:, :-1, :] + a[:, 1:, :])
            else:
                core = 0.5 * (a[:, :-1, :] + a[:, 1:, :])
                out = np.pad(core, ((0, 0), (0, 1), (0, 0)), constant_values=np.nan)
        # pad/crop to target
        pad_y = max(y_len - out.shape[1], 0)
        pad_x = max(x_len - out.shape[2], 0)
        if pad_y or pad_x:
            out = np.pad(out, ((0, 0), (0, pad_y), (0, pad_x)), constant_values=np.nan)
        if wrap_x and out.shape[2] > 1:
            out[:, :, -1] = out[:, :, 0]
        return out[:, :y_len, :x_len].astype(np.float32)

    def B2Ta(self, uB, vB, y_len, x_len):
        """
        Compute T-cell speed magnitude by B->T averaging (method A) of u and v.

        This routine:
        1) Averages corner-staggered u and v to T-cell centers using `Bavg_methA`.
        2) Computes speed magnitude: sqrt(uT^2 + vT^2).

        Parameters
        ----------
        uB, vB : numpy.ndarray
            Corner-staggered velocity components with shape (time, y, x).
        y_len, x_len : int
            Target T-grid spatial dimensions.

        Returns
        -------
        sT : numpy.ndarray
            T-centered speed magnitude with shape (time, y_len, x_len).

        Notes
        -----
        Because `Bavg_methA` propagates NaNs, this speed product will also be NaN
        where any of the contributing corners are NaN.
        """
        uT = self.Bavg_methA(uB, y_len, x_len)
        vT = self.Bavg_methA(vB, y_len, x_len)
        sT = np.sqrt(uT**2 + vT**2)
        return sT

    def B2Tb(self, uB, vB, y_len, x_len, wrap_x):
        """
        Compute T-cell speed magnitude by B->T averaging (method B) of u and v.

        This routine:
        1) Averages corner-staggered u and v to T-cell centers using `Bavg_methB`
            (NaNs treated as 0.0 to enforce no-slip at masked corners).
        2) Computes speed magnitude: sqrt(uT^2 + vT^2).

        Parameters
        ----------
        uB, vB : numpy.ndarray
            Corner-staggered velocity components with shape (time, y, x).
        y_len, x_len : int
            Target T-grid spatial dimensions.
        wrap_x : bool
            Whether to wrap the x-direction cyclically during averaging.

        Returns
        -------
        sT : numpy.ndarray
            T-centered speed magnitude with shape (time, y_len, x_len).

        Notes
        -----
        This method is typically more robust near land masks because NaNs are treated
        as 0.0 prior to averaging.
        """
        uT = self.Bavg_methB(uB, y_len, x_len, wrap_x)
        vT = self.Bavg_methB(vB, y_len, x_len, wrap_x)
        sT = np.sqrt(uT**2 + vT**2)
        return sT#, uT, vT

    def B2Tx(self, uB_da, vB_da):
        """
        Regrid B-grid velocity components to T points via xESMF and compute speed.

        This method uses a pre-defined regridder `self.reG` (constructed by
        `define_reG_weights` / `_ensure_reG_defined`) to map B-grid corner fields onto
        the T-grid. Prior to regridding, NaNs are replaced by 0.0 to enforce a
        no-slip-style boundary condition.

        Parameters
        ----------
        uB_da, vB_da : xarray.DataArray
            Corner-staggered velocity components. Must be compatible with the
            configured xESMF regridder. Typically dims are (time, y, x) or similar.

        Returns
        -------
        ispd_Tx : xarray.DataArray
            T-grid speed magnitude computed as hypot(uT, vT), dtype float32. Uses
            `xr.apply_ufunc(..., dask="parallelized")` for Dask compatibility.

        Notes
        -----
        - This is distinct from the simple 2x2 averaging methods (`B2Ta`, `B2Tb`);
        it uses spatial regridding (often bilinear) and may yield smoother fields.
        - Because NaNs are converted to 0.0, masked/corner land values do not
        propagate into neighbouring T cells in the same way as NaN-preserving
        averages.
        """
        self._ensure_reG_defined()
        uTx = self.reG(uB_da.fillna(0.0))  # no-slip: NaNs->0
        vTx = self.reG(vB_da.fillna(0.0))
        return xr.apply_ufunc(np.hypot, uTx, vTx, dask="parallelized", output_dtypes=[np.float32])
    
    def C2T(self, uvelE, uvelN, vvelE, vvelN,
            y_len: int | None = None,
            x_len: int | None = None,
            wrap_x: bool | None = None,
            combine: str = "mean",
            nan_to_zero: bool = True):
        """
        Reconstruct T-grid east/north velocity components from C-grid edge fields.

        Inputs are assumed to be *physical* components in the local east/north basis:
        - uvelE, uvelN: (E,N) components defined on the U-stagger (x-faces)
        - vvelE, vvelN: (E,N) components defined on the V-stagger (y-faces)

        The mapping is performed by edge-to-center averaging:
        - U-stagger → T via direction="x"
        - V-stagger → T via direction="y"

        Then the two estimates are combined to produce a single (E,N) pair at T.

        Parameters
        ----------
        uvelE, uvelN, vvelE, vvelN : xarray.DataArray or array-like
            C-grid velocity components on the U and V staggers.
        y_len, x_len : int, optional
            Target T-grid spatial sizes. Defaults to `self.CICE_dict["*_dim_length"]`.
        wrap_x : bool, optional
            Seam handling in x. Defaults to `self.CICE_dict.get("wrap_x", True)`.
        combine : {"mean","uv"}, default "mean"
            How to combine U- and V-stagger estimates at T:
            - "mean": E_T = 0.5*(E_from_U + E_from_V), N_T = 0.5*(N_from_U + N_from_V)
            - "uv"  : E_T = E_from_U, N_T = N_from_V (classic C-grid intuition)
        nan_to_zero : bool, default True
            If True, treat NaNs as 0 prior to averaging (no-slip-like).

        Returns
        -------
        velE_T, velN_T : same type as inputs
            Eastward and northward velocity components at T-cell centers.

        Notes
        -----
        If you only need speed, you can compute `hypot(velE_T, velN_T)` downstream
        using `xr.apply_ufunc(np.hypot, ...)` for Dask-parallel execution.
        """
        y_len = int(y_len or self.CICE_dict["y_dim_length"])
        x_len = int(x_len or self.CICE_dict["x_dim_length"])
        if wrap_x is None:
            wrap_x = bool(self.CICE_dict.get("wrap_x", True))
        # Map each stagger to T
        self.logger.info("averaging uvelE to T-grid")
        EU = self.Cavg_meth(uvelE, direction="x", y_len=y_len, x_len=x_len, wrap_x=wrap_x, nan_to_zero=nan_to_zero)
        self.logger.info("averaging uvelN to T-grid")
        NU = self.Cavg_meth(uvelN, direction="x", y_len=y_len, x_len=x_len, wrap_x=wrap_x, nan_to_zero=nan_to_zero)
        self.logger.info("averaging vvelE to T-grid")
        EV = self.Cavg_meth(vvelE, direction="y", y_len=y_len, x_len=x_len, wrap_x=wrap_x, nan_to_zero=nan_to_zero)
        self.logger.info("averaging vvelN to T-grid")
        NV = self.Cavg_meth(vvelN, direction="y", y_len=y_len, x_len=x_len, wrap_x=wrap_x, nan_to_zero=nan_to_zero)
        combine = (combine or "mean").lower().strip()
        if combine == "uv":
            velE_T = EU
            velN_T = NV
        elif combine == "mean":
            velE_T = 0.5 * (EU + EV)
            velN_T = 0.5 * (NU + NV)
        else:
            raise ValueError("C2T: combine must be 'mean' or 'uv'")
        return velE_T, velN_T

    def _binary_days(self, mask, 
                     win     : int, 
                     min_days: int,
                     centered: bool = True):
        """
        Apply a binary-day (persistence) filter to a 0/1 mask using a sliding window.

        For each time index, this computes the number of '1' days within a window of
        length `win`. The output is 1 where the count is at least `min_days`, else 0.

        Parameters
        ----------
        mask : numpy.ndarray
            Input mask of shape (time, nj, ni). Values may be bool or {0,1}. NaNs are
            treated as 0.
        win : int
            Window length in timesteps/days.
        min_days : int
            Minimum number of '1' values required within the window to mark the
            central/trailing time as 1. Values greater than `win` are clipped to `win`.
        centered : bool, default=True
            If True, the window is centered on each time. If False, uses a trailing
            window ending at each time.

        Returns
        -------
        out : numpy.ndarray
            uint8 array of shape (time, nj, ni) with values {0,1}.

        Notes
        -----
        - Uses `numpy.lib.stride_tricks.sliding_window_view` along the time axis,
        producing an intermediate array of shape (T-win+1, nj, ni, win).
        - Edge handling is by zero padding (i.e., edges default to 0 where a full
        window is not available). This differs from xarray rolling with partial
        windows.
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
        Compute a generic binary-day persistence mask using xarray rolling counts.

        In each rolling window of length W, the mask is set to 1 if at least M days
        are True/1 within that window:
            mask_bin(t) = 1 if sum(mask[t-window]) >= M else 0

        Parameters
        ----------
        I_mask : xarray.DataArray
            Input boolean or 0/1 mask with dimensions (time, y, x).
        mask_name : str
            Name of the output variable (e.g., "FI_mask", "PI_mask").
        bin_win_days : int, optional
            Window length W. Defaults to `self.bin_win_days` when None.
        bin_min_days : int, optional
            Minimum days M. Defaults to `self.bin_min_days` when None.
        time_dim : str, optional
            Name of the time dimension. Defaults to `self.CICE_dict["time_dim"]`.
        centered : bool, default=True
            If True, uses a centered window. If False, uses a trailing window.

        Returns
        -------
        ds : xarray.Dataset
            Dataset containing one uint8 variable named `mask_name`.

        Notes
        -----
        - Time chunking is adjusted so that time chunks are not larger than the
        rolling window, which generally improves Dask rolling performance.
        - Rolling uses `min_periods=M`, so partial windows can contribute when at
        least M samples are available; this differs from a strict “full-window only”
        approach.
        - If the input has zero time length, returns an empty dataset (with logging).
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
        Backward-compatible wrapper: binary-day persistence for fast-ice masks.

        Parameters
        ----------
        FI_mask : xarray.DataArray
            Input fast-ice mask (boolean or 0/1) with dimensions (time, y, x).
        bin_win_days, bin_min_days : int, optional
            Window length and minimum count passed to `classify_binary_days_mask`.
        time_dim : str, optional
            Time dimension name passed to `classify_binary_days_mask`.
        centered : bool, default=True
            Center the rolling window.

        Returns
        -------
        xr.Dataset
            Dataset containing a uint8 variable named "FI_mask".
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
        Wrapper: binary-day persistence for pack-ice masks.

        Parameters
        ----------
        PI_mask : xarray.DataArray
            Input pack-ice mask (boolean or 0/1) with dimensions (time, y, x).
        bin_win_days, bin_min_days : int, optional
            Window length and minimum count passed to `classify_binary_days_mask`.
        time_dim : str, optional
            Time dimension name passed to `classify_binary_days_mask`.
        centered : bool, default=True
            Center the rolling window.

        Returns
        -------
        xr.Dataset
            Dataset containing a uint8 variable named "PI_mask".

        Notes
        -----
        “Pack ice” here is the converse classification of fast ice under the same
        concentration threshold and a speed threshold (defined elsewhere).
        """
        return self.classify_binary_days_mask(PI_mask,
                                             mask_name    = "PI_mask",
                                             bin_win_days = bin_win_days,
                                             bin_min_days = bin_min_days,
                                             time_dim     = time_dim,
                                             centered     = centered)
    
    def _wrap_x_last_equals_first(self, da: xr.DataArray) -> xr.DataArray:
        xdim = self.CICE_dict["x_dim"]
        if da.sizes.get(xdim, 0) <= 1:
            return da
        if xdim not in da.coords:
            da = da.assign_coords({xdim: np.arange(da.sizes[xdim], dtype="int32")})
        last_label = da[xdim].isel({xdim: -1})
        first_col  = da.isel({xdim: 0})
        return da.where(da[xdim] != last_label, other=first_col)
    
    def compute_ispdT(self, uvel, vvel=None) -> xr.DataArray:
        """
        Compute T-grid sea-ice speed magnitude using the configured staggering→T strategy.

        Supported conversion modes in `self.BorC2T_type`:
        - "Ta": B-grid corner 2x2 mean (NaNs propagate) via `B2Ta`
        - "Tb": B-grid corner 2x2 mean (NaNs->0, no-slip) via `B2Tb`
        - "Tx": regrid (xESMF) via `B2Tx`
        - "Tc": C-grid edge→center reconstruction from (uvelE/uvelN/vvelE/vvelN) via `C2T`

        IMPORTANT: "Tc" is exclusive. If "Tc" is selected, it must be the only mode
        and the result is NOT averaged with other modes.

        Parameters
        ----------
        uvel, vvel : xarray.DataArray or xarray.Dataset
            For Ta/Tb/Tx:
                uvel and vvel are the usual B-grid corner-staggered components.
            For Tc:
                pass a Dataset (or mapping-like object) containing the four required
                variables: uvelE, uvelN, vvelE, vvelN. `vvel` is ignored in this case.
        label : str, optional
            Optional label for logging; not required.

        Returns
        -------
        ispd_T : xarray.DataArray
            T-grid speed magnitude (float32).

        Raises
        ------
        ValueError
            If "Tc" is combined with other modes, or if required C-grid variables
            are missing.
        """
        sel   = self._b2t_selection_set()
        is_Tc = ("Tc" in sel)
        if is_Tc:
            if isinstance(uvel, xr.Dataset):
                ds = uvel
            elif isinstance(uvel, dict):
                ds = uvel
            else:
                raise ValueError("Tc requires `uvel` to be an xarray.Dataset (or dict-like) containing uvelE, uvelN, vvelE, vvelN.")
            for req in ("uvelE", "uvelN", "vvelE", "vvelN"):
                if req not in ds:
                    raise ValueError(f"Tc requested but missing required C-grid variable '{req}'")
            y_len  = self.CICE_dict["y_dim_length"]
            x_len  = self.CICE_dict["x_dim_length"]
            wrap_x = self.CICE_dict.get("wrap_x", True)
            self.logger.info("Computing ice speeds from C-grid to T-grid")
            velE_T, velN_T = self.C2T(ds["uvelE"], ds["uvelN"], ds["vvelE"], ds["vvelN"],
                                      y_len       = y_len,
                                      x_len       = x_len,
                                      wrap_x      = wrap_x,
                                      combine     = "mean", 
                                      nan_to_zero = True)
            ispd_Tc = xr.apply_ufunc(np.hypot, velE_T, velN_T,
                                     dask          = "parallelized", 
                                     output_dtypes = [np.float32]).astype(np.float32)
            if self.CICE_dict.get("wrap_x", True):
                ispd_Tc = self._wrap_x_last_equals_first(ispd_Tc)
            return ispd_Tc.rename("ispd_Tc")
        # -----------------------------
        # Existing multi-member branch
        # -----------------------------
        y_len   = self.CICE_dict["y_dim_length"]
        x_len   = self.CICE_dict["x_dim_length"]
        wrap_x  = self.CICE_dict.get("wrap_x", True)
        members = []
        if "Ta" in self.BorC2T_type:
            if hasattr(self, "B2Ta"):
                ispd = self.B2Ta(uvel, vvel, y_len, x_len).astype(np.float32)
                ispd_Ta = xr.DataArray(ispd, dims=(self.CICE_dict["three_dims"]), name="ispd_Ta")
                members.append(ispd_Ta)
            else:
                self.logger.warning("Ta requested in BorC2T_type but B2Ta not defined. Skipping Ta.")
        if "Tb" in self.BorC2T_type:
            if hasattr(self, "B2Tb"):
                ispd = self.B2Tb(uvel, vvel, y_len, x_len, wrap_x=wrap_x).astype(np.float32)
                ispd_Tb = xr.DataArray(ispd, dims=(self.CICE_dict["three_dims"]), name="ispd_Tb")
                members.append(ispd_Tb)
            else:
                self.logger.warning("Tb requested in BorC2T_type but B2Tb not defined. Skipping Tb.")
        if "Tx" in self.BorC2T_type:
            if hasattr(self, "B2Tx"):
                ispd_Tx = self.B2Tx(uvel, vvel).astype(np.float32).rename("ispd_Tx")
                members.append(ispd_Tx)
            else:
                self.logger.warning("Tx requested in BorC2T_type but B2Tx not defined. Skipping Tx.")
        if len(members) == 1:
            return members[0].astype(np.float32)
        return xr.concat(members, dim="__concat_dim__").mean("__concat_dim__", skipna=True).astype(np.float32)

    def _apply_thresholds(self,
                          ispd_T   : xr.DataArray,
                          aice     : xr.DataArray,
                          label    : str,
                          ice_type : str = "FI",
                          mask_name: str | None = None) -> xr.Dataset:
        """
        Apply hemisphere slicing and speed/concentration thresholds to build an ice-class mask.

        This method slices both speed and concentration to the configured hemisphere,
        applies class-dependent threshold logic, and returns a compact chunked dataset
        containing a single 0/1 mask variable.

        Parameters
        ----------
        ispd_T : xarray.DataArray
            T-grid ice-speed magnitude with dimensions (time, y, x).
        aice : xarray.DataArray
            Sea-ice concentration on the same grid and dimensions as `ispd_T`.
        label : str
            Descriptive label used in log messages (e.g., "daily", "rolling-mean").
        ice_type : {"FI","PI"}, default="FI"
            Classification type:
            - FI (fast ice): (aice > icon_thresh) AND (0 < speed <= ispd_thresh)
            - PI (pack ice): (aice > icon_thresh) AND (speed > ispd_thresh)
        mask_name : str, optional
            Output variable name. Defaults to "{ice_type}_mask".

        Returns
        -------
        ds : xarray.Dataset
            Dataset with one boolean variable named `mask_name`, chunked according to
            `self.CICE_dict["FI_chunks"]` and computed into memory.

        Raises
        ------
        ValueError
            If `ice_type` is not one of {"FI","PI"}.

        Notes
        -----
        - Hemisphere slicing is performed via `self.slice_hemisphere`.
        - The result is `.compute()`'d to materialize the mask (and avoid deferred
        graphs in downstream workflows).
        - Any coordinate variables that duplicate spatial dims are dropped to keep the
        dataset compact and consistent for writing.
        """
        ice_type_u = (ice_type or "FI").upper()
        mask_name  = mask_name or f"{ice_type_u}_mask"
        ispd_hem   = self.slice_hemisphere(ispd_T)
        aice_hem   = self.slice_hemisphere(aice)
        if ice_type_u == "FI":
            self.logger.info(f"   thresholding ({label}): aice>{self.icon_thresh:0.2f} and 0<speed<={self.ispd_thresh:.3e}")
            I_mask = (aice_hem > self.icon_thresh) & (ispd_hem > 0) & (ispd_hem <= self.ispd_thresh)
        elif ice_type_u == "PI":
            self.logger.info(f"   thresholding ({label}): aice>{self.icon_thresh:0.2f} and speed>{self.ispd_thresh:.3e}" )
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
        Attach standard T-grid auxiliary coordinates (lon, lat, angle, area) to a DataArray.

        This is a convenience utility to ensure diagnostic fields (e.g., speed) carry
        T-grid geolocation and cell metadata. The method assumes the spatial
        dimensions match `self.CICE_dict["spatial_dims"]` and that `self.G_t` has been
        loaded.

        Parameters
        ----------
        da : xarray.DataArray
            DataArray with spatial dims matching the model T-grid, typically with
            overall dims (time, y, x).
        name : str, optional
            New name for the DataArray. If None, keeps `da.name`.

        Returns
        -------
        da_out : xarray.DataArray
            DataArray with added coordinates:
            - lon(y,x), lat(y,x), angle(y,x), area(y,x)

        Notes
        -----
        This attaches coordinates as 2D arrays; it does not modify the underlying data
        or rechunk by itself.
        """
        ydim, xdim = self.CICE_dict["spatial_dims"]
        da         = da.rename(name or da.name)
        da         = da.assign_coords({"lon"  : ((ydim, xdim), self.G_t["lon"].data),
                                       "lat"  : ((ydim, xdim), self.G_t["lat"].data),
                                       "angle": ((ydim, xdim), self.G_t["angle"].data),
                                       "area" : ((ydim, xdim), self.G_t["area"].data)})
        return da
    
    def _b2t_selection_set(self) -> set[str]:
        """
        Normalize self.BorC2T_type into a set of mode tokens, e.g. {"Ta","Tb"} or {"Tc"}.
        Accepts list/tuple/set or a string (optionally comma/space separated).
        """
        b2t = getattr(self, "BorC2T_type", None)

        if isinstance(b2t, (list, tuple, set)):
            sel = {str(t).strip() for t in b2t if str(t).strip()}
            if getattr(self, "logger", None) is not None:
                self.logger.info(f"BorC2T_type raw={b2t!r} (iterable) -> selection={sorted(sel)}")
            return sel

        s = str(b2t).replace(",", " ").strip()
        sel = set(s.split()) if s else set()

        if getattr(self, "logger", None) is not None:
            self.logger.info(f"BorC2T_type raw={b2t!r} (string) -> normalized={s!r} -> selection={sorted(sel)}")

        return sel

    def classify_fast_ice(self,
                        dt0_str               : str  = None,
                        dtN_str               : str  = None,
                        bin_win_days          : int  = None,
                        bin_min_days          : int  = None,
                        enable_rolling_output : bool = False,
                        roll_win_days         : int  = None,
                        time_dim_name         : str  = None):
        """
        Classify fast ice (FI) from model output using daily, persistence, and optional rolling criteria.

        Fast ice is defined where sea-ice concentration exceeds a threshold and ice
        speed is small (immobile) but non-zero:
            FI := (aice > icon_thresh) AND (0 < speed <= ispd_thresh)

        The method returns:
        1) A daily FI mask computed directly from instantaneous (daily) T-grid speed.
        2) A binary-day (persistence) FI mask computed from the *extended-range* daily
            mask using a W-day window requiring at least M FI days.
        3) Optionally, a rolling-mean FI mask computed from a rolling-mean speed field.
        4) Speed diagnostics (daily and optional rolling speed) in `ispd_out`.

        To avoid edge artefacts in persistence and rolling calculations, the input time
        range is extended by half the largest window on each side, masks are computed
        on the extended range, and then cropped back to [dt0_str, dtN_str].

        Parameters
        ----------
        dt0_str, dtN_str : str
            Start and end dates (inclusive) in a pandas-parsable format (e.g., "YYYY-MM-DD").
        bin_win_days : int, optional
            Persistence window length W (days). Defaults to `self.bin_win_days`.
        bin_min_days : int, optional
            Minimum FI days M within the window. Defaults to `self.bin_min_days`.
        enable_rolling_output : bool, default=False
            If True, also compute a rolling-mean speed field and corresponding FI mask.
        roll_win_days : int, optional
            Rolling-mean window length (days). Defaults to `self.mean_period`.
        time_dim_name : str, optional
            Name of the time dimension. Defaults to `self.CICE_dict["time_dim"]`.

        Returns
        -------
        FI_dly : xarray.Dataset
            Daily fast-ice mask dataset with variable "FI_mask" over [dt0_str, dtN_str].
        FI_bin : xarray.Dataset
            Binary-day (persistence) fast-ice mask dataset with variable "FI_mask" over
            [dt0_str, dtN_str].
        FI_roll : xarray.Dataset or None
            Rolling-mean fast-ice mask dataset with variable "FI_mask" over [dt0_str, dtN_str],
            or None if `enable_rolling_output` is False.
        ispd_out : dict
            Dictionary with speed diagnostics:
            - "daily": DataArray "ispd_T"
            - "rolly": DataArray "ispd_T_roll" (if enabled), else None

        Notes
        -----
        - Loads variables ['aice','uvel','vvel'] over an extended time range.
        - Computes T-grid speed once on the extended range and reuses it.
        - Uses Dask-friendly chunking with small time chunks to improve rolling operations.
        - The binary-day mask is computed from the extended daily mask to prevent “all-zero”
        behaviour near the analysis window edges.
        """
        time_dim      = time_dim_name or self.CICE_dict["time_dim"]
        bin_win_days  = int(bin_win_days  or self.bin_win_days)
        bin_min_days  = int(bin_min_days  or self.bin_min_days)
        roll_win_days = int(roll_win_days or self.mean_period)
        sel           = self._b2t_selection_set()
        is_Tc         = ("Tc" in sel)
        if is_Tc and any(t in sel for t in ("Ta", "Tb", "Tx")):
            raise ValueError("BorC2T_type includes 'Tc' but also Ta/Tb/Tx. 'Tc' must be exclusive.")
        self.load_cice_grid()
        if not is_Tc:
            self.define_reG_weights()
        ispd_out = {'daily' : None,
                    'rolly' : None}
        # extend enough time for BOTH methods
        ext_per = max(bin_win_days // 2, roll_win_days // 2 if enable_rolling_output else 0)
        dt0_ext = pd.to_datetime(dt0_str) - pd.Timedelta(days=ext_per)
        dtN_ext = pd.to_datetime(dtN_str) + pd.Timedelta(days=ext_per)
        # Load (extended) data
        if is_Tc:
            variables = ["aice", "uvelE", "uvelN", "vvelE", "vvelN"]
        else:
            variables = ['aice','uvel','vvel']
        self.logger.info(f"loading model data between {dt0_ext:%Y-%m-%d} and {dtN_ext:%Y-%m-%d}")
        CICE_all = self.load_cice_zarr(slice_hem = False,
                                       variables = variables,
                                       dt0_str   = f"{dt0_ext:%Y-%m-%d}",
                                       dtN_str   = f"{dtN_ext:%Y-%m-%d}")
        # Compact coords & set friendly chunks
        aice = CICE_all['aice'].drop_vars(self.CICE_dict["drop_coords"]).astype(np.float32)
        if is_Tc:
            uvelE = CICE_all['uvelE'].drop_vars(self.CICE_dict["drop_coords"]).astype(np.float32)
            uvelN = CICE_all['uvelN'].drop_vars(self.CICE_dict["drop_coords"]).astype(np.float32)
            vvelE = CICE_all['vvelE'].drop_vars(self.CICE_dict["drop_coords"]).astype(np.float32)
            vvelN = CICE_all['vvelN'].drop_vars(self.CICE_dict["drop_coords"]).astype(np.float32)
        else:            
            uvel = CICE_all['uvel'].drop_vars(self.CICE_dict["drop_coords"]).astype(np.float32)
            vvel = CICE_all['vvel'].drop_vars(self.CICE_dict["drop_coords"]).astype(np.float32)
        # chunking is an art form
        yck    = self.CICE_dict.get("y_chunk", 540)
        xck    = self.CICE_dict.get("x_chunk", 1440)
        tch    = max(8, min(bin_win_days, roll_win_days))  # small time chunks help rolling parallelism
        chunks = {self.CICE_dict["time_dim"]: tch, 
                  self.CICE_dict["y_dim"]   : yck, 
                  self.CICE_dict["x_dim"]   : xck}
        aice   = aice.chunk(chunks)
        # Build T-grid speed ONCE over the extended span
        if is_Tc:
            ds_vel = xr.Dataset({"uvelE": uvelE.chunk(chunks),
                                 "uvelN": uvelN.chunk(chunks), 
                                 "vvelE": vvelE.chunk(chunks),
                                 "vvelN": vvelN.chunk(chunks)}).persist()
            ispd_T = self.compute_ispdT(ds_vel).persist()
        else:
            uvel, vvel = uvel.chunk(chunks), vvel.chunk(chunks)
            ispd_T     = self.compute_ispdT(uvel, vvel).persist()
        # DAILY ispd
        ispd_T            = self._with_T_coords(ispd_T, name="ispd_T")
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
        Classify pack ice (PI) from model output using daily, persistence, and optional rolling criteria.

        Pack ice is defined as the complement of fast ice under the same concentration
        threshold, using a speed threshold:
            PI := (aice > icon_thresh) AND (speed > ispd_thresh)

        The method returns:
        1) A daily PI mask computed from instantaneous (daily) T-grid speed.
        2) A binary-day (persistence) PI mask computed from the *extended-range* daily
            mask using a W-day window requiring at least M PI days.
        3) Optionally, a rolling-mean PI mask computed from rolling-mean speed.
        4) Speed diagnostics (daily and optional rolling speed) in `ispd_out`.

        As with fast-ice classification, the input time range is extended to avoid
        edge artefacts in persistence/rolling calculations, and products are cropped
        back to [dt0_str, dtN_str].

        Parameters
        ----------
        dt0_str, dtN_str : str
            Start and end dates (inclusive) in a pandas-parsable format (e.g., "YYYY-MM-DD").
        bin_win_days : int, optional
            Persistence window length W (days). Defaults to `self.bin_win_days`.
        bin_min_days : int, optional
            Minimum PI days M within the window. Defaults to `self.bin_min_days`.
        enable_rolling_output : bool, default=False
            If True, also compute a rolling-mean speed field and corresponding PI mask.
        roll_win_days : int, optional
            Rolling-mean window length (days). Defaults to `self.mean_period`.
        time_dim_name : str, optional
            Name of the time dimension. Defaults to `self.CICE_dict["time_dim"]`.

        Returns
        -------
        PI_dly : xarray.Dataset
            Daily pack-ice mask dataset with variable "PI_mask" over [dt0_str, dtN_str].
        PI_bin : xarray.Dataset
            Binary-day (persistence) pack-ice mask dataset with variable "PI_mask" over
            [dt0_str, dtN_str].
        PI_roll : xarray.Dataset or None
            Rolling-mean pack-ice mask dataset with variable "PI_mask" over [dt0_str, dtN_str],
            or None if `enable_rolling_output` is False.
        ispd_out : dict
            Dictionary with speed diagnostics:
            - "daily": DataArray "ispd_T"
            - "rolly": DataArray "ispd_T_roll" (if enabled), else None

        Notes
        -----
        - Pack ice is only classified where `aice > icon_thresh`.
        - Binary-day persistence is computed on the extended daily mask to avoid
        edge-related bias and truncation artefacts.
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
        ispd_T = self.compute_ispdT(uvel, vvel, label="base").persist()
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
        Classify total sea ice (SI) using only a concentration threshold.

        Total sea ice is defined as:
            SI := (aice > icon_thresh)

        This classification does not use speed, and therefore does not produce
        binary-day (persistence) or rolling-mean products.

        Parameters
        ----------
        dt0_str, dtN_str : str
            Start and end dates (inclusive) in a pandas-parsable format (e.g., "YYYY-MM-DD").
        time_dim_name : str, optional
            Name of the time dimension. Defaults to `self.CICE_dict["time_dim"]`.

        Returns
        -------
        SI_dly : xarray.Dataset
            Daily sea-ice mask dataset with variable "SI_mask" over [dt0_str, dtN_str].
        SI_bin : None
            Always None (no persistence product for SI).
        SI_roll : None
            Always None (no rolling-mean product for SI).
        ispd_out : dict
            Empty dict, returned for API consistency with FI/PI classifiers.

        Notes
        -----
        - Loads only `aice` (no velocity fields required).
        - Applies hemisphere slicing via `self.slice_hemisphere`.
        - Output is chunked using `self.CICE_dict["FI_chunks"]` and computed into memory.
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
