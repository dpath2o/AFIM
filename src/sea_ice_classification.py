import xarray    as xr
import pandas    as pd
import numpy     as np
__all__ = ["SeaIceClassification"]
class SeaIceClassification:
    """
    Classify Antarctic landfast ice, pack ice, and total sea ice from CICE output.

    This class implements a threshold-based classification framework driven by:
      - sea-ice concentration (`aice`), and
      - sea-ice speed magnitude on the model’s analysis grid (nominally the CICE T-grid).

    It supports multiple velocity staggering strategies for constructing a *T-grid*
    speed magnitude (`ispd_T`) from the native model velocity components:

    Staggering → T-grid strategies
    ------------------------------
    The strategy is controlled by `self.BorC2T_type`, which can be a string (e.g. "Tb")
    or an iterable (e.g. ["Ta", "Tb"]). Supported tokens are:

    **B-grid derived modes (legacy / B-grid CICE)**
    - "Ta": Convert corner-staggered B-grid velocities to a T-grid speed using a 2×2
            box mean where NaNs propagate (area-like average; intended to preserve holes).
    - "Tb": Convert corner-staggered B-grid velocities to a T-grid speed using a 2×2
            box mean with NaNs → 0 (no-slip-like fill; tends to remove coastal holes).
    - "Tx": Regrid B-grid velocities to T-grid using precomputed weights (e.g., xESMF).

    **C-grid derived mode (new / C-grid CICE)**
    - "Tc": Reconstruct T-grid (cell-centre) velocities directly from C-grid edge
            components:
              uvelE, uvelN, vvelE, vvelN  →  (velE_T, velN_T)  →  ispd_T
            This mode uses a dedicated C-grid-to-T-grid reconstruction (`C2T`) and is
            conceptually distinct from the B-grid averaging/regridding modes. Because "Tc"
            is a physically different reconstruction (not merely an alternate smoothing),
            it is enforced as **exclusive**: if "Tc" is selected, it must be the only token.

    Classification products
    -----------------------
    For a given date range [dt0, dtN] the class can produce:
      - Daily masks: computed directly from daily `aice` and `ispd_T`.
      - Persistence (binary-day) masks: computed from the extended-range daily mask using
        a sliding window length W (days) and minimum count M (days in window).
      - Optional rolling-mean masks: computed from a rolling-mean speed field and the
        same concentration + speed thresholds.

    Definitions
    -----------
    Using thresholds `icon_thresh` (concentration) and `ispd_thresh` (speed):

      Fast ice (FI):
        (aice > icon_thresh) AND (0 < ispd_T <= ispd_thresh)

      Pack ice (PI):
        (aice > icon_thresh) AND (ispd_T > ispd_thresh)

      Total sea ice (SI):
        (aice > icon_thresh)

    Persistence (“binary-day”) classification:
      A day is classed as persistent if at least `bin_min_days` out of `bin_win_days`
      days in a centred sliding window are classed as FI/PI.

    Required configuration and dependencies
    ---------------------------------------
    The class expects the broader AFIM stack to provide the following attributes:

    Core configuration
    - CICE_dict : dict
        Must include (names may vary across your stack):
          "time_dim", "y_dim", "x_dim", "y_dim_length", "x_dim_length",
          "drop_coords", "FI_chunks", and optionally "wrap_x", "y_chunk", "x_chunk".
    - icon_thresh : float
        Concentration threshold (e.g., 0.15).
    - ispd_thresh : float
        Speed threshold for FI/PI discrimination (e.g., 1e-4 m/s).
    - bin_win_days : int
        Persistence window length W.
    - bin_min_days : int
        Minimum count M within window.
    - mean_period : int
        Rolling mean window length for optional rolling output.
    - BorC2T_type : str | list[str] | set[str]
        Mode tokens controlling `compute_ispdT()`.

    External methods (implemented elsewhere in your toolbox)
    - load_cice_grid()
    - define_reG_weights()
    - load_cice_zarr(...)
    - slice_hemisphere(...)
    - reG(...) or B2Tx(...)  (depending on your architecture)
    - B2Ta(...), B2Tb(...), B2Tx(...)  (if Ta/Tb/Tx are used)
    - C2T(...)  (required for Tc)
    - _apply_thresholds(...), _rolling_mean(...), _with_T_coords(...)
    - classify_binary_days_fast_ice(...), classify_binary_days_pack_ice(...)
    - _wrap_x_last_equals_first(...)  (if wrap_x=True)

    Notes
    -----
    - Many operations are Dask-friendly; however, persistence and rolling calculations
      benefit from careful chunking along time.
    - The high-level classifiers extend the requested time window to avoid edge artefacts
      when applying centred rolling/persistence windows, then crop back to [dt0, dtN].
    - For C-grid runs ("Tc"), velocities are expected as edge components (E/N) and are
      reconstructed rather than averaged/regridded. This typically reduces the systematic
      “west-of-land NaN/zero” artefacts seen in B-grid corner-staggered products.
    """
    
    def __init__(self, sim_name=None, logger=None, **kwargs):
        """
        Initialise the classifier and attach configuration.

        Parameters
        ----------
        sim_name : str, optional
            Simulation identifier used for logging and/or output paths.
        logger : logging.Logger, optional
            Logger instance. If None, a module-level logger is created.
        **kwargs
            Configuration values to be set as attributes on `self`. Typical keys include:
            `CICE_dict`, `icon_thresh`, `ispd_thresh`, `bin_win_days`, `bin_min_days`,
            `mean_period`, `BorC2T_type`, plus grid or path metadata used by your stack.

        Notes
        -----
        - This constructor is intentionally lightweight; grid loading and weight
          construction are deferred to the high-level classifier methods.
        """
        import logging
        self.sim_name = sim_name
        self.logger   = logger or logging.getLogger(__name__)
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
            y_len      : int | None = None,
            x_len      : int | None = None,
            wrap_x     : bool | None = None,
            combine    : str = "mean",
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
        from numpy.lib.stride_tricks import sliding_window_view
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
                                  ice_type    : str | None = None,
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
        ice_type = ice_type or self.ice_type
        time_dim = time_dim or self.CICE_dict["time_dim"]
        W        = int(bin_win_days or self.bin_win_days)
        M        = int(bin_min_days or self.bin_min_days)
        self._check_ice_type(ice_type)
        self.define_ice_mask_name(ice_type=ice_type)
        if I_mask.sizes.get(time_dim, 0) == 0:
            self.logger.warning(f"   {self.mask_name} has zero time length for this slice; returning empty dataset.")
            return I_mask.astype(np.uint8).rename(self.mask_name).to_dataset()
        tlen = I_mask.sizes[time_dim]
        tch  = max(8, min(W, tlen))
        if I_mask.chunks is None or I_mask.chunksizes[time_dim][0] in (0, None) or I_mask.chunksizes[time_dim][0] > tch:
            I_mask = I_mask.chunk({time_dim: tch})
        counts = I_mask.astype("uint8").rolling({time_dim: W}, center=centered, min_periods=M).sum()
        if ice_type=="FI":
            self.logger.info("       fast ice")
            da_bin = (counts >= M).astype("uint8").rename(self.mask_name)
        else:
            self.logger.info("       pack ice")
            da_bin = (counts < M).astype("uint8").rename(self.mask_name)
        return da_bin.to_dataset()
    
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
        Compute T-grid sea-ice speed magnitude (`ispd_T`) using the configured strategy.

        This method converts native velocity components to a speed magnitude on the
        analysis grid (nominally T-grid). The conversion mode(s) are determined by
        `self.BorC2T_type`, which is normalised to a set of tokens via `_b2t_selection_set()`.

        Supported modes
        ---------------
        **B-grid derived modes (can be combined and averaged):**
        - "Ta": Uses `self.B2Ta(uvel, vvel, y_len, x_len)` to produce a T-grid-like speed
          from corner-staggered B-grid components via a 2×2 mean where NaNs propagate.
        - "Tb": Uses `self.B2Tb(uvel, vvel, y_len, x_len, wrap_x=...)` to produce a T-grid-like
          speed via a 2×2 mean where NaNs are converted to 0 (no-slip-like fill).
        - "Tx": Uses `self.B2Tx(uvel, vvel)` to regrid B-grid components to T-grid (e.g., via weights).

        If multiple of {Ta, Tb, Tx} are selected, the method returns the mean of the
        members (skipna=True) to form a composite speed estimate.

        **C-grid derived mode (exclusive):**
        - "Tc": Reconstructs T-grid velocities directly from C-grid edge components using
          `self.C2T(uvelE, uvelN, vvelE, vvelN, ...)`, then computes:
              ispd_T = hypot(velE_T, velN_T)

          "Tc" must be the **only** selected token. It is not averaged with Ta/Tb/Tx,
          because it represents a different staggering and reconstruction pathway.

        Parameters
        ----------
        uvel : xr.DataArray, xr.Dataset, or mapping-like
            Velocity input. Interpretation depends on mode:

            - For Ta/Tb/Tx:
                `uvel` is the B-grid corner-staggered u-component DataArray (typically `uvel`).
                `vvel` must be provided as the matching v-component DataArray.

            - For Tc:
                `uvel` must be a Dataset (or dict-like) containing the four required C-grid
                edge components:
                    "uvelE", "uvelN", "vvelE", "vvelN"
                In Tc mode, the `vvel` argument is ignored.

        vvel : xr.DataArray, optional
            B-grid corner-staggered v-component. Required for Ta/Tb/Tx; ignored for Tc.

        Returns
        -------
        xr.DataArray
            T-grid speed magnitude (float32). Name will reflect the mode:
              - "ispd_Tc" for Tc
              - "ispd_Ta"/"ispd_Tb"/"ispd_Tx" when only one B-grid member is selected
              - composite mean when multiple B-grid members are selected

        Raises
        ------
        ValueError
            If Tc is requested but `uvel` does not provide the required variables.
        ValueError
            If Tc is selected together with any of Ta/Tb/Tx.
        ValueError
            If required C-grid variables ("uvelE", "uvelN", "vvelE", "vvelN") are missing.
        AttributeError
            If a requested conversion method (e.g., B2Ta/B2Tb/B2Tx/C2T) is not defined.

        Notes
        -----
        - In Tc mode, `wrap_x` handling is applied (if `self.CICE_dict["wrap_x"]` is True)
          using `_wrap_x_last_equals_first()` so that the last x-column matches the first.
        - Speed is computed using `np.hypot` via `xr.apply_ufunc(..., dask="parallelized")`.
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
                          ice_type : str = "FI") -> xr.Dataset:
        """
        Apply hemisphere slicing and speed/concentration thresholds to build an ice-class mask.

        FI: (aice > icon_thresh) AND (0 < speed <= ispd_thresh)
        PI: (aice > icon_thresh) AND (speed > ispd_thresh)
        """
        ice_type = ice_type or self.ice_type
        self._check_ice_type(ice_type)
        self.define_ice_mask_name(ice_type=ice_type)
        # Hemisphere slice
        ispd_hem = self.slice_hemisphere(ispd_T)
        aice_hem = self.slice_hemisphere(aice)
        # Align to avoid any latent coordinate mismatch (time or x/y)
        ispd_hem, aice_hem = xr.align(ispd_hem, aice_hem, join="inner")
        # Thresholding
        if ice_type == "FI":
            self.logger.info(f"   thresholding ({label}): aice>{self.icon_thresh:0.2f} and 0<speed<={self.ispd_thresh:.3e}")
            I_mask = (aice_hem > self.icon_thresh) & (ispd_hem > 0) & (ispd_hem <= self.ispd_thresh)
            desc = "Fast-ice mask: (aice>icon_thresh) & (0<ispd_T<=ispd_thresh)"
        else:  # PI
            self.logger.info(f"   thresholding ({label}): aice>{self.icon_thresh:0.2f} and speed>{self.ispd_thresh:.3e}")
            I_mask = (aice_hem > self.icon_thresh) & (ispd_hem > self.ispd_thresh)
            desc = "Pack-ice mask: (aice>icon_thresh) & (ispd_T>ispd_thresh)"
        # Make compact 0/1 mask (uint8) for cheaper IO and stable downstream use
        da = I_mask.rename(self.mask_name).astype("uint8")
        da.attrs.update({"long_name"  : f"{ice_type} ice mask",
                         "description": desc,
                         "ice_type"   : ice_type,
                         "icon_thresh": float(self.icon_thresh),
                         "ispd_thresh": float(self.ispd_thresh),
                         "label"      : str(label)})
        # Chunk + compute
        chunks = self.CICE_dict.get("FI_chunks", None)
        DS = da.to_dataset()
        if chunks is not None:
            DS = DS.chunk(chunks)
        # Drop known “fat” coordinate vars (keep dimension coords unless you *really* want to remove them)
        for c in self.CICE_dict.get("drop_coords", []):
            if c in DS.coords or c in DS.data_vars:
                DS = DS.drop_vars(c)
        # Materialize (your stated design choice)
        DS = DS.compute()
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
        Normalise `self.BorC2T_type` into a set of mode tokens.

        This helper accepts either:
          - an iterable (list/tuple/set) of tokens, e.g. ["Ta", "Tb"], or
          - a string containing tokens separated by commas and/or whitespace,
            e.g. "Ta,Tb" or "Ta Tb".

        Parameters
        ----------
        None

        Returns
        -------
        set[str]
            Set of non-empty token strings (whitespace stripped), e.g. {"Ta", "Tb"} or {"Tc"}.

        Side Effects
        ------------
        Logs the raw value and normalised selection if `self.logger` exists.

        Notes
        -----
        - This method does not validate tokens beyond basic normalisation; validation is
          performed by downstream methods (e.g., `compute_ispdT()` enforcing Tc exclusivity).
        """
        b2t = getattr(self, "BorC2T_type", None)
        if isinstance(b2t, (list, tuple, set)):
            sel = {str(t).strip() for t in b2t if str(t).strip()}
            if getattr(self, "logger", None) is not None:
                self.logger.info(f"BorC2T_type raw={b2t!r} (iterable) -> selection={sorted(sel)}")
            return sel
        s   = str(b2t).replace(",", " ").strip()
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
        Classify fast ice (FI) using concentration + speed thresholds with optional persistence/rolling products.

        Fast ice is defined where concentration exceeds `icon_thresh` and the T-grid
        speed magnitude is small but non-zero:

            FI := (aice > icon_thresh) AND (0 < ispd_T <= ispd_thresh)

        The method computes products over an *extended* time range to avoid edge artefacts
        for centred windows (binary-day persistence and optional rolling-mean speed). It
        then crops the outputs back to the requested interval [dt0_str, dtN_str].

        C-grid "Tc" behaviour (important)
        ---------------------------------
        If `BorC2T_type` includes "Tc", the classification uses C-grid edge velocities:
            uvelE, uvelN, vvelE, vvelN
        and reconstructs T-grid velocities using `C2T`. In this mode:
          - `define_reG_weights()` is not called (no B→T regridding required),
          - the velocity dataset passed to `compute_ispdT()` must contain the four
            C-grid variables, and
          - "Tc" must be exclusive (it cannot be combined with Ta/Tb/Tx).

        Parameters
        ----------
        dt0_str, dtN_str : str
            Start and end dates (inclusive). Must be parseable by pandas, e.g. "YYYY-MM-DD".
        bin_win_days : int, optional
            Persistence window length W (days). Defaults to `self.bin_win_days`.
        bin_min_days : int, optional
            Minimum number of FI days M required within the window. Defaults to `self.bin_min_days`.
        enable_rolling_output : bool, default False
            If True, compute a rolling-mean speed field and corresponding FI mask.
        roll_win_days : int, optional
            Rolling mean window length (days). Defaults to `self.mean_period`.
        time_dim_name : str, optional
            Name of the time dimension. Defaults to `self.CICE_dict["time_dim"]`.

        Returns
        -------
        FI_dly : xr.Dataset
            Daily fast-ice mask dataset with variable "FI_mask" over [dt0_str, dtN_str].
        FI_bin : xr.Dataset
            Binary-day (persistence) fast-ice mask dataset with variable "FI_mask" over
            [dt0_str, dtN_str]. Computed from the *extended-range* daily mask and then cropped.
        FI_roll : xr.Dataset or None
            Rolling-mean fast-ice mask dataset with variable "FI_mask" over [dt0_str, dtN_str],
            or None if `enable_rolling_output` is False.
        ispd_out : dict
            Speed diagnostics:
              - "daily": xr.DataArray "ispd_T" over the *extended* range (with T coords)
              - "rolly": xr.DataArray "ispd_T_roll" over the *extended* range (if enabled), else None

        Raises
        ------
        ValueError
            If "Tc" is selected together with any of {"Ta", "Tb", "Tx"}.
        ValueError
            If required C-grid velocity variables are missing when "Tc" is selected.
        Exception
            Propagates exceptions from dataset loading, coordinate attachment, thresholding,
            or persistence/rolling helpers.

        Notes
        -----
        Implementation outline
        ----------------------
        1) Determine selection set from `BorC2T_type` and enforce Tc exclusivity.
        2) Load grid and, if needed (non-Tc), load/define regridding weights.
        3) Extend requested time range by half the maximum of:
             - persistence window (W/2)
             - rolling window (R/2, if enabled)
        4) Load required variables across the extended range:
             - Tc:  aice + (uvelE,uvelN,vvelE,vvelN)
             - else: aice + (uvel,vvel)
        5) Chunk for rolling-friendly time blocks.
        6) Compute `ispd_T` once on the extended range and persist.
        7) Compute:
             - daily FI mask on extended range, then crop
             - binary-day FI mask on extended daily mask, then crop
             - optional rolling FI mask on rolled speed, then crop
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

    def classify_ice(self,
                     ice_type              : str  = None,
                     dt0_str               : str  = None,
                     dtN_str               : str  = None,
                     bin_win_days          : int  = None,
                     bin_min_days          : int  = None,
                     enable_rolling_output : bool = False,
                     roll_win_days         : int  = None,
                     time_dim_name         : str  = None): 
        """
        Classify ocean ice as fast ice (FI), pack ice (PI) and both of them as sea ice (SI) 
        using concentration (+ speed for fast and pack) thresholds with optional persistence/rolling 
        products.

        Pack ice is defined where concentration exceeds `icon_thresh` and the T-grid
        speed magnitude exceeds the fast-ice speed threshold:

            PI := (aice > icon_thresh) AND (ispd_T > ispd_thresh)
            
        Fast ice is defined where concentration exceeds `icon_thresh` and the T-grid
        speed magnitude is small but non-zero:

            FI := (aice > icon_thresh) AND (0 < ispd_T <= ispd_thresh)

        C-grid "Tc" behaviour
        ---------------------
        If `BorC2T_type` includes "Tc", the classification uses C-grid edge velocities:
            uvelE, uvelN, vvelE, vvelN
        and reconstructs T-grid velocities using `C2T` inside `compute_ispdT()`. In this mode:
        - `define_reG_weights()` is not called,
        - "Tc" must be exclusive (cannot be combined with Ta/Tb/Tx).

        The method computes products over an *extended* time range to avoid edge artefacts
        for centred windows (binary-day persistence and optional rolling-mean speed), then
        crops the outputs back to the requested interval [dt0_str, dtN_str].

        Returns
        -------
        PI_dly : xr.Dataset
            Daily pack-ice mask dataset with variable "PI_mask" over [dt0_str, dtN_str].
        PI_bin : xr.Dataset
            Binary-day (persistence) pack-ice mask dataset with variable "PI_mask" over
            [dt0_str, dtN_str]. Computed from the extended daily mask and then cropped.
        PI_roll : xr.Dataset or None
            Rolling-mean pack-ice mask dataset with variable "I_mask" over [dt0_str, dtN_str],
            or None if `enable_rolling_output` is False.
        ispd_out : dict
            Speed diagnostics:
            - "daily": xr.DataArray "ispd_T" over the *extended* range (with T coords)
            - "rolly": xr.DataArray "ispd_T_roll" over the *extended* range (if enabled), else None

        Raises
        ------
        ValueError
            If "Tc" is selected together with any of {"Ta", "Tb", "Tx"}.
        ValueError
            If required C-grid velocity variables are missing when "Tc" is selected.
        """

        time_dim      = time_dim_name or self.CICE_dict["time_dim"]
        bin_win_days  = int(bin_win_days  or self.bin_win_days)
        bin_min_days  = int(bin_min_days  or self.bin_min_days)
        roll_win_days = int(roll_win_days or self.mean_period)
        ice_type      = ice_type or self.ice_type
        self._check_ice_type(ice_type)
        self.define_ice_mask_name(ice_type)
        # ---- Tc selection rules (mirror classify_fast_ice) ----
        sel   = self._b2t_selection_set()
        is_Tc = ("Tc" in sel)
        if is_Tc and any(t in sel for t in ("Ta", "Tb", "Tx")):
            raise ValueError("BorC2T_type includes 'Tc' but also Ta/Tb/Tx. 'Tc' must be exclusive.")
        # ---- Grid + (optional) B->T regridding weights ----
        self.load_cice_grid()
        if not is_Tc:
            self.define_reG_weights()
        ispd_out = {"daily": None, "rolly": None}
        # ---- extend enough time for BOTH persistence + rolling (if enabled) ----
        ext_per = max(bin_win_days // 2, roll_win_days // 2 if enable_rolling_output else 0)
        dt0_ext = pd.to_datetime(dt0_str) - pd.Timedelta(days=ext_per)
        dtN_ext = pd.to_datetime(dtN_str) + pd.Timedelta(days=ext_per)
        # ---- Load (extended) data ----
        if is_Tc:
            variables = ["aice", "uvelE", "uvelN", "vvelE", "vvelN"]
        else:
            variables = ["aice", "uvel", "vvel"]
        self.logger.info(f"loading model data between {dt0_ext:%Y-%m-%d} and {dtN_ext:%Y-%m-%d}")
        CICE_all = self.load_cice_zarr(slice_hem = False,
                                       variables = variables,
                                       dt0_str   = f"{dt0_ext:%Y-%m-%d}",
                                       dtN_str   = f"{dtN_ext:%Y-%m-%d}")
        # ---- Compact coords & types ----
        aice = CICE_all["aice"].drop_vars(self.CICE_dict["drop_coords"]).astype(np.float32)
        if is_Tc:
            missing = [v for v in ("uvelE", "uvelN", "vvelE", "vvelN") if v not in CICE_all]
            if missing:
                raise ValueError(f"Tc selected but missing required C-grid velocity fields: {missing}")
            uvelE = CICE_all["uvelE"].drop_vars(self.CICE_dict["drop_coords"]).astype(np.float32)
            uvelN = CICE_all["uvelN"].drop_vars(self.CICE_dict["drop_coords"]).astype(np.float32)
            vvelE = CICE_all["vvelE"].drop_vars(self.CICE_dict["drop_coords"]).astype(np.float32)
            vvelN = CICE_all["vvelN"].drop_vars(self.CICE_dict["drop_coords"]).astype(np.float32)
        else:
            uvel  = CICE_all["uvel"].drop_vars(self.CICE_dict["drop_coords"]).astype(np.float32)
            vvel  = CICE_all["vvel"].drop_vars(self.CICE_dict["drop_coords"]).astype(np.float32)
        # ---- Chunking (same “art form” as FI) ----
        yck    = self.CICE_dict.get("y_chunk", 540)
        xck    = self.CICE_dict.get("x_chunk", 1440)
        tch    = max(8, min(bin_win_days, roll_win_days))  # small time chunks help rolling parallelism
        chunks = {self.CICE_dict["time_dim"]: tch,
                  self.CICE_dict["y_dim"]   : yck,
                  self.CICE_dict["x_dim"]   : xck}
        aice = aice.chunk(chunks)
        if ice_type=='SI':
            aice_hem = self.slice_hemisphere(aice)
            self.logger.info(f"   thresholding (sea-ice): aice>{self.icon_thresh:0.2f} (no speed constraint)")
            I_mask = (aice_hem > self.icon_thresh)
            I_dly  = I_mask.rename(self.mask_name).to_dataset().chunk(self.CICE_dict["FI_chunks"]).compute()
            for dim in self.CICE_dict.get('spatial_dims', []):
                if dim in I_dly.coords:
                    I_dly = I_dly.drop_vars(dim)
            # Maintain the same return signature as classify_fast_ice / classify_pack_ice
            return I_dly, None, None, ispd_out
        if is_Tc:
            ds_vel = xr.Dataset({"uvelE": uvelE.chunk(chunks),
                                 "uvelN": uvelN.chunk(chunks), 
                                 "vvelE": vvelE.chunk(chunks),
                                 "vvelN": vvelN.chunk(chunks)}).persist()
            with self._suppress_large_graph_warning():
                ispd_T = self.compute_ispdT(ds_vel).persist()
        else:
            uvel, vvel = uvel.chunk(chunks), vvel.chunk(chunks)
            with self._suppress_large_graph_warning():
                ispd_T = self.compute_ispdT(uvel, vvel).persist()
        # DAILY ispd
        ispd_T            = self._with_T_coords(ispd_T, name="ispd_T")
        ispd_out['daily'] = ispd_T
        # ---- DAILY PI mask on EXTENDED range, THEN crop ----
        I_dly_ext = self._apply_thresholds(ispd_T, aice, label = "daily", ice_type = ice_type )
        I_dly     = I_dly_ext.sel(time=slice(dt0_str, dtN_str))
        # ---- BINARY-DAYS on EXTENDED DAILY mask, THEN crop ----
        self.logger.info("classifying binary-days mask")
        with self._suppress_large_graph_warning():
            I_bin_ext = self.classify_binary_days_mask(I_dly_ext[self.mask_name],     # use extended mask to avoid truncation artefacts
                                                       ice_type     = ice_type,
                                                       bin_win_days = bin_win_days,
                                                       bin_min_days = bin_min_days,
                                                       time_dim     = time_dim,
                                                       centered     = True)
        I_bin = I_bin_ext.sel(time=slice(dt0_str, dtN_str))
        # ---- Optional ROLLING on extended speed, THEN crop ----
        I_roll = None
        if enable_rolling_output:
            self.logger.info(f"applying rolling mean on T-grid speed (period={roll_win_days})")
            ispd_T_roll = self._rolling_mean(ispd_T, mean_period=roll_win_days)
            ispd_T_roll = self._with_T_coords(ispd_T_roll, name="ispd_T_roll")
            I_roll_ext  = self._apply_thresholds(ispd_T_roll, aice,
                                                 label     = "rolling-mean",
                                                 ice_type  = ice_type)
            I_roll            = I_roll_ext.sel(time=slice(dt0_str, dtN_str))
            ispd_out["rolly"] = ispd_T_roll
        return I_dly, I_bin, I_roll, ispd_out

    def classify_sea_ice(self,
                         dt0_str       : str = None,
                         dtN_str       : str = None,
                         time_dim_name : str = None):
        """
        Classify total sea ice (SI) using only a concentration threshold.

        Total sea ice is defined as:

            SI := (aice > icon_thresh)

        This classifier does not require velocity variables and therefore does not
        produce persistence (binary-day) or rolling-mean products.

        Parameters
        ----------
        dt0_str, dtN_str : str
            Start and end dates (inclusive). Must be parseable by pandas, e.g. "YYYY-MM-DD".
        time_dim_name : str, optional
            Name of the time dimension. Defaults to `self.CICE_dict["time_dim"]`.

        Returns
        -------
        SI_dly : xr.Dataset
            Daily sea-ice mask dataset with variable "SI_mask" over [dt0_str, dtN_str].
        SI_bin : None
            Always None (no persistence product for SI).
        SI_roll : None
            Always None (no rolling product for SI).
        ispd_out : dict
            Empty dict, returned for API consistency with FI/PI classifiers.

        Notes
        -----
        - Loads only `aice` across the requested period.
        - Applies hemisphere slicing via `slice_hemisphere()` prior to thresholding.
        - Uses chunking consistent with FI/PI products to keep downstream workflows uniform.
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

    def _miz_mask_from_aice_2d(aice2d: np.ndarray,
                            edge_aice_thresh: float,
                            width_cells: int,
                            aice_min: float | None,
                            aice_max: float | None,
                            exclude_polynyas: bool) -> np.ndarray:
        """
        Compute a 2D MIZ mask (uint8 0/1) from a 2D aice field.

        Strategy:
        1) ice = aice > edge_aice_thresh
        2) define "ocean open water" as open-water connected to the domain boundary (optional)
        3) distance = distance_transform_edt(~ocean_open)  # distance to ocean water (in cells)
        4) MIZ = ice & (distance <= width_cells) & optional(aice_min<=aice<=aice_max)

        Notes:
        - If exclude_polynyas=True, open water inside the pack that is NOT boundary-connected
            does not create local MIZ bands.
        """
        if aice2d is None:
            return None
        # Treat NaNs as open water (conservative for edge distance)
        a = np.asarray(aice2d)
        nanmask = ~np.isfinite(a)
        a = np.where(nanmask, 0.0, a)

        ice = a > edge_aice_thresh
        openw = ~ice

        if exclude_polynyas:
            lbl, n = ndi.label(openw)
            if n == 0:
                ocean_open = openw
            else:
                # Labels touching any boundary are "ocean" (boundary-connected open water)
                b = np.unique(np.concatenate([lbl[0, :], lbl[-1, :], lbl[:, 0], lbl[:, -1]]))
                b = b[b != 0]
                ocean_open = np.isin(lbl, b)
        else:
            ocean_open = openw

        # distance_transform_edt computes distance to nearest zero.
        # We want distance to ocean_open==True, so ocean cells should be zero.
        dist_cells = ndi.distance_transform_edt(~ocean_open)

        miz = ice & (dist_cells <= int(width_cells))

        if aice_min is not None:
            miz &= (a >= float(aice_min))
        if aice_max is not None:
            miz &= (a <= float(aice_max))

        return miz.astype(np.uint8)


    def classify_marginal_ice_zone(self,
                                dt0_str              : str,
                                dtN_str              : str,
                                miz_width_km         : float = 100.0,
                                cell_km              : float | None = None,
                                edge_aice_thresh     : float | None = None,
                                aice_min             : float | None = None,
                                aice_max             : float | None = 0.80,
                                exclude_polynyas     : bool  = True,
                                bin_win_days         : int   | None = None,
                                bin_min_days         : int   | None = None,
                                enable_rolling_output: bool  = False,
                                roll_win_days        : int   | None = None,
                                time_dim_name        : str   | None = None):
        """
        Classify Marginal Ice Zone (MIZ) using an ice-edge distance-band method.

        Default definition (recommended):
        - Ice edge: aice > edge_aice_thresh (defaults to icon_thresh, usually 0.15)
        - MIZ: ice cells within miz_width_km of the open-ocean ice edge
        - Optional aice band constraint: aice_min <= aice <= aice_max (defaults to <=0.80)
        - Optional exclusion of polynyas: open water not boundary-connected does NOT create local MIZ

        Returns
        -------
        MIZ_dly : xr.Dataset
            Daily MIZ mask dataset with variable "MIZ_mask" over [dt0_str, dtN_str].
        MIZ_bin : xr.Dataset or None
            Binary-day (persistence) MIZ mask dataset (same variable name) over [dt0_str, dtN_str],
            or None if bin_win_days/bin_min_days are not provided.
        MIZ_roll : xr.Dataset or None
            Rolling-time-smoothed MIZ mask dataset over [dt0_str, dtN_str], or None if disabled.
            Rolling is applied to aice (not distance), then MIZ is re-computed.
        diag : dict
            Diagnostics: includes "width_cells", "cell_km", "edge_aice_thresh", "aice_min/max".
        """
        time_dim      = time_dim_name or self.CICE_dict["time_dim"]
        roll_win_days = int(roll_win_days or self.mean_period)

        # thresholds/defaults
        edge_thr = float(edge_aice_thresh) if edge_aice_thresh is not None else float(self.icon_thresh)

        # infer cell_km if not provided (approx; recommend passing cell_km explicitly)
        if cell_km is None:
            # try to infer from T-cell area if available
            try:
                self.load_cice_grid()
                if hasattr(self, "G_t") and isinstance(self.G_t, xr.Dataset) and "tarea" in self.G_t:
                    # characteristic cell length = sqrt(area)
                    Lm = np.sqrt(self.G_t["tarea"].median().compute().item())  # meters
                    cell_km = float(Lm / 1000.0)
                else:
                    cell_km = 5.0
            except Exception:
                cell_km = 5.0

        width_cells = max(1, int(np.round(float(miz_width_km) / float(cell_km))))

        # extend for rolling / persistence (edge-distance itself has no temporal window)
        ext_per = 0
        if (bin_win_days is not None) and (bin_min_days is not None):
            ext_per = max(ext_per, int(bin_win_days) // 2)
        if enable_rolling_output:
            ext_per = max(ext_per, roll_win_days // 2)

        dt0_ext = pd.to_datetime(dt0_str) - pd.Timedelta(days=ext_per)
        dtN_ext = pd.to_datetime(dtN_str) + pd.Timedelta(days=ext_per)

        self.logger.info(f"loading model data between {dt0_ext:%Y-%m-%d} and {dtN_ext:%Y-%m-%d}")
        CICE_all = self.load_cice_zarr(slice_hem=False,
                                    variables=["aice"],
                                    dt0_str=f"{dt0_ext:%Y-%m-%d}",
                                    dtN_str=f"{dtN_ext:%Y-%m-%d}")

        aice = CICE_all["aice"].drop_vars(self.CICE_dict["drop_coords"]).astype(np.float32)
        aice = self.slice_hemisphere(aice)

        # IMPORTANT for distance transforms: you want full spatial context.
        # Rechunk so each time chunk sees the full (y,x) slab.
        ydim, xdim = self.CICE_dict["y_dim"], self.CICE_dict["x_dim"]
        tch = max(8, min(int(bin_win_days or 8), int(roll_win_days or 8)))
        aice = aice.chunk({time_dim: tch, ydim: -1, xdim: -1})

        # Compute daily MIZ mask on extended range (vectorized over time)
        miz_da = xr.apply_ufunc(
            _miz_mask_from_aice_2d,
            aice,
            input_core_dims=[[ydim, xdim]],
            output_core_dims=[[ydim, xdim]],
            vectorize=True,
            dask="parallelized",
            output_dtypes=[np.uint8],
            kwargs=dict(edge_aice_thresh=edge_thr,
                        width_cells=width_cells,
                        aice_min=aice_min,
                        aice_max=aice_max,
                        exclude_polynyas=exclude_polynyas)
        ).rename("MIZ_mask")

        MIZ_dly_ext = miz_da.to_dataset()

        # chunk/compute like your FI/PI products
        MIZ_dly_ext = MIZ_dly_ext.chunk(self.CICE_dict["FI_chunks"]).compute()

        # crop to requested interval
        MIZ_dly = MIZ_dly_ext.sel({time_dim: slice(dt0_str, dtN_str)})

        # optional persistence (“binary-days”) on the EXTENDED daily mask, then crop
        MIZ_bin = None
        if (bin_win_days is not None) and (bin_min_days is not None):
            W = int(bin_win_days)
            M = int(bin_min_days)
            self.logger.info(f"classifying binary-day MIZ mask (W={W}, M={M}, centered=True)")
            # rolling sum over full window; require full window (min_periods=W)
            rs = MIZ_dly_ext["MIZ_mask"].rolling({time_dim: W}, center=True, min_periods=W).sum()
            MIZ_bin_ext = (rs >= M).astype(np.uint8).rename("MIZ_mask").to_dataset()
            MIZ_bin_ext = MIZ_bin_ext.chunk(self.CICE_dict["FI_chunks"]).compute()
            MIZ_bin = MIZ_bin_ext.sel({time_dim: slice(dt0_str, dtN_str)})

        # optional rolling: smooth aice in time then recompute MIZ (helps edge jitter)
        MIZ_roll = None
        if enable_rolling_output:
            self.logger.info(f"computing rolling-aice MIZ (period={roll_win_days})")
            aice_roll = self._rolling_mean(aice, mean_period=roll_win_days)

            miz_roll_da = xr.apply_ufunc(
                _miz_mask_from_aice_2d,
                aice_roll,
                input_core_dims=[[ydim, xdim]],
                output_core_dims=[[ydim, xdim]],
                vectorize=True,
                dask="parallelized",
                output_dtypes=[np.uint8],
                kwargs=dict(edge_aice_thresh=edge_thr,
                            width_cells=width_cells,
                            aice_min=aice_min,
                            aice_max=aice_max,
                            exclude_polynyas=exclude_polynyas)
            ).rename("MIZ_mask")

            MIZ_roll_ext = miz_roll_da.to_dataset().chunk(self.CICE_dict["FI_chunks"]).compute()
            MIZ_roll = MIZ_roll_ext.sel({time_dim: slice(dt0_str, dtN_str)})

        diag = dict(width_cells=int(width_cells),
                    miz_width_km=float(miz_width_km),
                    cell_km=float(cell_km),
                    edge_aice_thresh=float(edge_thr),
                    aice_min=None if aice_min is None else float(aice_min),
                    aice_max=None if aice_max is None else float(aice_max),
                    exclude_polynyas=bool(exclude_polynyas))

        return MIZ_dly, MIZ_bin, MIZ_roll, diag