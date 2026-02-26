from __future__ import annotations
import xarray as xr
import numpy  as np
import pandas as pd
from pathlib  import Path
__all__ = ["SeaIceMetrics"]
class SeaIceMetrics:
    """
    Compute diagnostic and skill metrics for processed sea-ice classifications.

    This class aggregates common post-processing calculations used in AFIM workflows,
    including:

      - building standardised input dictionaries for metrics computations
        (fast ice / pack ice / total sea ice),
      - computing hemisphere-integrated time series (area, volume, thickness, rates),
      - computing spatial summary fields (temporal mean thickness; spatial rates),
      - computing inter-comparison skill statistics against observations,
      - computing seasonal statistics (timing, growth/retreat slopes, extrema),
      - computing persistence-distance statistics relative to an Antarctic coastline,
      - loading previously computed metrics from Zarr and padding/reindexing time.

    The principal entry point is `compute_sea_ice_metrics()`, which ingests a dictionary
    of masks and core prognostic/tendency fields (already prepared elsewhere), computes
    a suite of metrics, and optionally writes the results to a consolidated Zarr archive.

    Expected external configuration
    ------------------------------
    This class is designed to operate inside the broader AFIM stack and expects several
    attributes to exist on `self`, typically injected via kwargs or a toolbox manager:

    - logger : logging.Logger
        For progress/debug logging.
    - CICE_dict : dict
        Must contain at least:
          * "time_dim" and "spatial_dims"
          * "drop_coords" (list of coordinate vars to drop before computations)
          * "FI_chunks" (chunking spec for small 1D/2D outputs)
    - dt0_str, dtN_str : str
        Analysis window bounds ("YYYY-MM-DD").
    - D_metrics : str | Path
        Output directory for metrics Zarr archives.
    - metrics_name : str
        Suffix used when constructing default metrics filenames.
    - ice_type : str
        e.g., "FI", "FI_Tb_bin", etc. Cleaned internally to one of {"FI","PI","SI"}.
    - FIC_scale : float
        Scaling factor applied to area (commonly used to convert m² → km² or similar).
    - sim_config : dict
        Metadata inserted into output dataset attributes and/or scalars.
    - AF_FI_dict : dict
        Must include path to AF2020 fast-ice area observations (e.g., "P_AF2020_FIA").
    - BAS_dict : dict
        Must include coastline shapefile path (e.g., "P_Ant_Cstln") for persistence distance.
    - BorC2T_type, ispd_thresh, D_zarr, etc.
        Used by `load_computed_metrics()` to locate metrics on disk.

    Required external methods (defined elsewhere in your stack)
    ----------------------------------------------------------
    - compute_hemisphere_ice_area(...)
    - compute_hemisphere_ice_volume(...)
    - compute_hemisphere_ice_thickness(...)
    - compute_hemisphere_variable_aggregate(...)
    - compute_NSIDC_metrics()  (for PI/SI observational comparisons)
    - persistence_stability_index(...)
    - load_cice_grid(slice_hem=True)
    - _prepare_BAS_coast(...)
    - compute_grounded_iceberg_area() (optional, for GI area accounting)

    Notes
    -----
    - Many methods are written to be safe with Dask-backed inputs. Where appropriate,
      the implementation forces local computation to avoid shipping large task graphs.
    - Output writing uses Zarr v2 (explicitly `zarr_format=2`) for compatibility.
    """

    def __init__(self, **kwargs):
        """
        Initialise the metrics helper.

        Parameters
        ----------
        **kwargs
            Optional configuration injected into `self` by the surrounding workflow.
            Typical keys include `logger`, `CICE_dict`, directory paths, and analysis
            settings (dt0/dtN, thresholds, etc.).

        Notes
        -----
        - The shown implementation is a no-op; in production you typically assign kwargs
          onto `self` (as done in your other classes) or inherit a shared base class.
        """
        return

    @staticmethod
    def _clean_zarr_chunks(ds):
        """
        Remove per-variable chunk encodings and return an unchunked dataset.

        Zarr writes can be sensitive to stale or inconsistent `encoding["chunks"]`
        metadata inherited from prior reads. This helper deletes the `chunks` entry
        from each variable encoding (if present) and returns `ds.chunk(None)` to
        remove active Dask chunking.

        Parameters
        ----------
        ds : xr.Dataset
            Dataset to sanitise prior to writing.

        Returns
        -------
        xr.Dataset
            Dataset with `encoding["chunks"]` removed (where present) and with
            chunking disabled (`chunk(None)`).

        Notes
        -----
        - This does not change data values; it only adjusts encoding and chunk metadata.
        - Use immediately before `.to_zarr(...)` to avoid RTD/runtime confusion and
          Zarr rechunk warnings.
        """
        for var in ds.variables:
            if 'chunks' in ds[var].encoding:
                del ds[var].encoding['chunks']
        return ds.chunk(None)
    
    @staticmethod
    def thickness_tendency_to_mps(dvt: xr.DataArray) -> xr.DataArray:
        u = (dvt.attrs.get("units", "") or "").strip().lower()
        # common cases
        if u in ["m/s", "m s-1", "m s^-1"]:
            return dvt
        if u in ["cm/s", "cm s-1", "cm s^-1"]:
            return dvt * 0.01
        if "cm/day" in u or "cm d-1" in u:
            return dvt * 0.01 / 86400.0
        if "m/day" in u or "m d-1" in u:
            return dvt / 86400.0
        # if unknown, do not silently scale; return as-is and force you to inspect
        raise ValueError(f"Unrecognized DVT units: '{dvt.attrs.get('units','')}'")

    def metrics_data_dict(self, I_mask, I_data, A,
                        ice_type=None, *,
                        required=None,
                        optional=None):
        """
        Build a standardised input dictionary for FI/PI/SI metrics.

        This method:
        - validates required variables exist in `I_data`,
        - adds optional variables if present,
        - always includes the ice mask and T-grid area (`tarea`),
        - and (NEW) auto-adds any referenced cell-area measures (e.g., earea/narea)
            when a variable declares them via CF `cell_measures`, e.g.:
                cell_measures: "area: earea"
            This is needed for lateral-drag stress diagnostics (KuxE/KuxN/KuyE/KuyN).

        Parameters
        ----------
        I_mask : xr.DataArray
            Ice-type mask (e.g., FI/PI/SI), typically with dims (time, nj, ni) or (nj, ni).
        I_data : xr.Dataset or dict-like
            Source dataset/dictionary containing model variables.
        A : xr.DataArray
            T-grid cell area (m^2), stored into the output dict as "tarea".
        ice_type : str, optional
            Ice type key ("FI","PI","SI"). Defaults to `self.ice_type`.
        required : list[str], optional
            Required variable names that must exist in `I_data`.
        optional : list[str], optional
            Optional variable names to include if present in `I_data`.

        Returns
        -------
        dict
            Dictionary containing mask, tarea, and required/optional fields.

        Raises
        ------
        KeyError
            If any required variables are missing.

        Notes
        -----
        - NEW: If a present variable declares an area measure via `cell_measures`,
        and that area variable exists in `I_data`, it is added automatically.
        This supports stress diagnostics with earea/narea weights.
        """
        import re
        ice_type = ice_type or self.ice_type
        self._check_ice_type(ice_type)
        self.define_ice_mask_name(ice_type=ice_type)

        # Defaults (tune as needed)
        required = required or ["aice", "hi", "uvel", "vvel"]
        optional = optional or [
            "strength", "dvidtt", "dvidtd", "daidtt", "daidtd",
            # NEW: lateral-drag stress outputs (if written by CICE)
            "KuxE", "KuxN", "KuyE", "KuyN",
            # NEW: possible area weights (often static 2-D fields)
            "earea", "narea", "uarea",
        ]

        def _has(var):
            if isinstance(I_data, xr.Dataset):
                return var in I_data.data_vars
            try:
                return var in I_data
            except TypeError:
                return hasattr(I_data, var)

        def _get(var):
            if isinstance(I_data, xr.Dataset):
                return I_data[var]
            try:
                return I_data[var]
            except Exception:
                return getattr(I_data, var)

        out = {
            self.mask_name: I_mask,
            "tarea": A,
        }

        missing_required = [v for v in required if not _has(v)]
        if missing_required:
            raise KeyError(f"metrics_data_dict missing required variables in I_data: {missing_required}")

        # Add required
        for v in required:
            out[v] = _get(v)

        # Add optional
        missing_optional = []
        for v in optional:
            if _has(v):
                out[v] = _get(v)
            else:
                missing_optional.append(v)

        if missing_optional:
            try:
                self.logger.info(f"metrics_data_dict: skipping missing optional vars: {missing_optional}")
            except Exception:
                pass

        # ------------------------------------------------------------------
        # NEW: auto-add any referenced area measures from CF cell_measures
        # ------------------------------------------------------------------
        try:
            area_vars_needed = set()
            pat = re.compile(r"area:\s*([A-Za-z0-9_]+)")
            for k, da in list(out.items()):
                if not isinstance(da, xr.DataArray):
                    continue
                cm = da.attrs.get("cell_measures", "") or da.attrs.get("cell_measures:area", "")
                m = pat.search(str(cm))
                if m:
                    area_vars_needed.add(m.group(1))

            for avar in sorted(area_vars_needed):
                if avar in out:
                    continue
                if _has(avar):
                    out[avar] = _get(avar)
                    self.logger.info(f"metrics_data_dict: auto-added area measure '{avar}' referenced by cell_measures")
                else:
                    self.logger.info(f"metrics_data_dict: variable(s) reference area '{avar}' via cell_measures, but '{avar}' not found in I_data")
        except Exception as e:
            self.logger.debug(f"metrics_data_dict: cell_measures parsing skipped ({e})")

        return out

    def compute_sea_ice_metrics(self, da_dict, P_mets_zarr,
                                ice_type       = None,
                                dt0_str        = None,
                                dtN_str        = None,
                                ice_area_scale = None):
        import re
        ice_type = ice_type or self.ice_type
        dt0_str  = dt0_str  or self.dt0_str
        dtN_str  = dtN_str  or self.dtN_str

        self._check_ice_type(ice_type)
        self.define_ice_mask_name(ice_type=ice_type)

        if ice_type == "FI":
            ice_area_scale = self.FIC_scale
        else:
            ice_area_scale = self.SIC_scale

        spatial_dim_names = self.CICE_dict["spatial_dims"]  # NEW: used below for stress aggregation

        self.logger.info(f"¡¡¡ COMPUTING ICE METRICS for {ice_type} !!!")
        self.logger.info(f"    results to {P_mets_zarr}")

        I_mask = da_dict[self.mask_name]

        # --------- Guard helpers ----------
        def has(*keys):
            return all(k in da_dict and da_dict[k] is not None for k in keys)

        def _all_nan(da):
            try:
                return bool(da.isnull().all().compute())
            except Exception:
                try:
                    return bool(np.all(np.isnan(da)))
                except Exception:
                    return False

        def valid(da):
            if da is None:
                return False
            try:
                return not _all_nan(da)
            except Exception:
                return True

        def maybe_compute(tag, fn, req_keys, *, store, out_key, post=None):
            if not has(*req_keys):
                self.logger.info(f"skipping {tag}: missing {set(req_keys) - set(da_dict.keys())}")
                return None
            for k in req_keys:
                if isinstance(da_dict[k], xr.DataArray) and not valid(da_dict[k]):
                    self.logger.info(f"skipping {tag}: input '{k}' is all-NaN")
                    return None
            try:
                out = fn()
                if post is not None:
                    out = post(out)
                store[out_key] = out.load() if isinstance(out, xr.DataArray) else out
                return out
            except Exception as e:
                self.logger.warning(f"{tag} failed: {e}")
                return None

        # --------- Pull what exists (do NOT assume keys) ----------
        I_C   = da_dict.get("aice")
        I_T   = da_dict.get("hi")
        I_S   = da_dict.get("strength")
        I_TVT = da_dict.get("dvidtt")
        I_MVT = da_dict.get("dvidtd")
        I_TAT = da_dict.get("daidtt")
        I_MAT = da_dict.get("daidtd")
        A     = da_dict.get("tarea")

        # NEW: lateral drag stress components (if present)
        KuxE = da_dict.get("KuxE")
        KuxN = da_dict.get("KuxN")
        KuyE = da_dict.get("KuyE")
        KuyN = da_dict.get("KuyN")

        METS = {}

        # --------- Time series metrics (conditional) ----------
        IA = maybe_compute("ice area",
                        lambda: self.compute_hemisphere_ice_area(I_C, A, ice_area_scale=ice_area_scale),
                        req_keys=("aice", "tarea"),
                        store=METS,
                        out_key=f"{ice_type}A")

        IV = maybe_compute("ice volume",
                        lambda: self.compute_hemisphere_ice_volume(I_C, I_T, A),
                        req_keys=("aice", "hi", "tarea"),
                        store=METS,
                        out_key=f"{ice_type}V")

        IT = maybe_compute("ice thickness",
                        lambda: self.compute_hemisphere_ice_thickness(I_C, I_T, A),
                        req_keys=("aice", "hi", "tarea"),
                        store=METS,
                        out_key=f"{ice_type}T")

        # FIXED: strength call should include IS and A
        IS = maybe_compute("ice strength (1D, area-weighted hPa)",
                        lambda: self.compute_area_weighted_strength_hpa(
                            SIC=I_C, HI=I_T, IS=I_S, A=A,
                            spatial_dim_names=spatial_dim_names,
                            sic_threshold=self.icon_thresh,
                            hmin=0.05,
                            assume_IS_units="N/m"),
                        req_keys=("aice", "hi", "strength", "tarea"),
                        store=METS,
                        out_key=f"{ice_type}S")

        ITVR = maybe_compute("thermo volume rate (1D)",
                            lambda: self.compute_hemisphere_ice_volume_rate(DVT=I_TVT, IV=IV, A=A, SIC=I_C,
                                                                            mode="absolute", out_units="per_day"),
                            req_keys=("aice", "dvidtt"),
                            store=METS,
                            out_key=f"{ice_type}TVR")

        IMVR = maybe_compute("dynamic volume rate (1D)",
                            lambda: self.compute_hemisphere_ice_volume_rate(DVT=I_MVT, IV=IV, A=A, SIC=I_C,
                                                                            mode="absolute", out_units="per_day"),
                            req_keys=("aice", "dvidtd"),
                            store=METS,
                            out_key=f"{ice_type}MVR")

        ITAR = maybe_compute("thermo area rate (1D)",
                            lambda: self.compute_hemisphere_ice_area_rate(DAT=I_TAT, IA=IA, A=A,
                                                                          mode="absolute", out_units="per_day"),
                            req_keys=("daidtt", "tarea"),
                            store=METS,
                            out_key=f"{ice_type}TAR",
                            post=lambda x: x) if (I_TAT is not None and IA is not None and A is not None) else None
        if (I_TAT is not None) and (IA is None):
            self.logger.info("skipping thermo area rate: requires IA to be computed first")

        IMAR = maybe_compute("dynamic area rate (1D)",
                            lambda: self.compute_hemisphere_ice_area_rate(DAT=I_MAT, IA=IA, A=A,
                                                                          mode="absolute", out_units="per_day"),
                            req_keys=("daidtd", "tarea"),
                            store=METS,
                            out_key=f"{ice_type}MAR") if (I_MAT is not None and IA is not None and A is not None) else None
        if (I_MAT is not None) and (IA is None):
            self.logger.info("skipping dynamic area rate: requires IA to be computed first")

        IP = maybe_compute("persistence aggregate (2D)",
                        lambda: self.compute_hemisphere_variable_aggregate(I_C),
                        req_keys=("aice",),
                        store=METS,
                        out_key=f"{ice_type}P")

        # ------------------------------------------------------------------
        # NEW PATCH: hemispheric aggregates for lateral-drag stress components
        # ------------------------------------------------------------------
        def _parse_area_name_from_cell_measures(da: xr.DataArray) -> str | None:
            cm = str(da.attrs.get("cell_measures", "") or "")
            m = re.search(r"area:\s*([A-Za-z0-9_]+)", cm)
            return m.group(1) if m else None

        def _get_area_for_tau(da: xr.DataArray) -> tuple[xr.DataArray | None, str | None]:
            # prefer explicit cell_measures area variable (earea/narea/etc)
            aname = _parse_area_name_from_cell_measures(da)
            if aname and (aname in da_dict) and (da_dict[aname] is not None):
                return da_dict[aname], aname
            # fallback to tarea if nothing else available
            if A is not None:
                return A, "tarea"
            return None, None

        def _mask_for_tau(tau: xr.DataArray) -> xr.DataArray | None:
            # Use the ice-type mask if available; align to tau
            if I_mask is None:
                return None
            try:
                m = (I_mask > 0)
                m, tau_al = xr.align(m, tau, join="inner")
                return m.astype(bool)
            except Exception:
                return (I_mask > 0)

        def _add_stress_agg(tau: xr.DataArray, base_name: str):
            """
            Compute and store:
            <base>_mean (signed), <base>_abs_mean, <base>_valid_area_m2
            using appropriate area weights if available.
            """
            if tau is None or not valid(tau):
                self.logger.info(f"skipping {base_name}: missing or all-NaN")
                return

            A_tau, A_name = _get_area_for_tau(tau)
            if A_tau is None:
                self.logger.info(f"skipping {base_name}: no area weights available (expected earea/narea or tarea)")
                return

            m_tau = _mask_for_tau(tau)

            self.logger.info(f"computing **LATERAL-DRAG STRESS AGGREGATE**: {base_name}")
            self.logger.debug(f"  • using area weights : {A_name}")
            self.logger.debug(f"  • output units       : Pa")
            self.logger.debug(f"  • masked by          : {self.mask_name} (>0)")

            try:
                ds_tau = self.compute_hemisphere_area_weighted_stress(
                    tau=tau,
                    A=A_tau,
                    spatial_dim_names=spatial_dim_names,
                    mask=m_tau,
                    out_units="Pa",
                    return_abs=True,
                    min_area_m2=0.0,
                    name=base_name,
                )

                # rename valid area to avoid collisions
                ds_tau = ds_tau.rename({"valid_area_m2": f"{base_name}_valid_area_m2"})

                # store each 1D series
                for v in ds_tau.data_vars:
                    METS[v] = ds_tau[v].load()

            except Exception as e:
                self.logger.warning(f"{base_name} stress aggregate failed: {e}")

        # component aggregates
        _add_stress_agg(KuxE, f"{ice_type}KuxE")
        _add_stress_agg(KuyE, f"{ice_type}KuyE")
        _add_stress_agg(KuxN, f"{ice_type}KuxN")
        _add_stress_agg(KuyN, f"{ice_type}KuyN")

        # magnitude aggregates (more interpretable than signed components)
        if (KuxE is not None) and (KuyE is not None) and valid(KuxE) and valid(KuyE):
            try:
                KuE_mag = np.hypot(KuxE, KuyE)
                KuE_mag.name = "KuE_mag"
                _add_stress_agg(KuE_mag, f"{ice_type}KuE_mag")
            except Exception as e:
                self.logger.warning(f"{ice_type}KuE_mag failed: {e}")

        if (KuxN is not None) and (KuyN is not None) and valid(KuxN) and valid(KuyN):
            try:
                KuN_mag = np.hypot(KuxN, KuyN)
                KuN_mag.name = "KuN_mag"
                _add_stress_agg(KuN_mag, f"{ice_type}KuN_mag")
            except Exception as e:
                self.logger.warning(f"{ice_type}KuN_mag failed: {e}")

        # --------- Spatial summary diagnostics (conditional) ----------
        time_dim = self.CICE_dict["time_dim"]

        if I_T is not None and valid(I_T):
            try:
                self.logger.info("computing **ICE THICKNESS TEMPORAL-MEAN**")
                METS[f"{ice_type}HI"] = I_T.mean(dim=time_dim).load()
            except Exception as e:
                self.logger.warning(f"mean thickness failed: {e}")

        if (I_S is not None) and (I_T is not None) and valid(I_S) and valid(I_T):
            try:
                self.logger.info("computing **ICE STRENGTH TEMPORAL-SUM**; units mPa")
                IST = (I_S / I_T.where(I_T > 0)).sum(dim=time_dim) / 1e6
                METS[f"{ice_type}ST"] = IST.load()
            except Exception as e:
                self.logger.warning(f"strength temporal-sum failed: {e}")

        if I_TVT is not None and valid(I_TVT):
            try:
                self.logger.info("computing **ICE VOLUME TENDENCY (SPATIAL RATE)**; units m/yr")
                METS[f"{ice_type}TVR_YR"] = ((I_TVT * 1e2).mean(dim=time_dim) / 3.65).load()
            except Exception as e:
                self.logger.warning(f"TVR_YR failed: {e}")

        if I_MVT is not None and valid(I_MVT):
            try:
                self.logger.info("computing **ICE VOLUME TENDENCY (SPATIAL RATE)**; units m/yr")
                METS[f"{ice_type}MVR_YR"] = ((I_MVT * 1e2).mean(dim=time_dim) / 3.65).load()
            except Exception as e:
                self.logger.warning(f"MVR_YR failed: {e}")

        if (I_TAT is not None) and (A is not None) and valid(I_TAT) and valid(A):
            try:
                self.logger.info("computing **ICE AREA TENDENCY (SPATIAL RATE)**; units m/yr")
                METS[f"{ice_type}TAR_YR"] = ((I_TAT * A).mean(dim=time_dim) / 31_536_000).load()
            except Exception as e:
                self.logger.warning(f"TAR_YR failed: {e}")

        if (I_MAT is not None) and (A is not None) and valid(I_MAT) and valid(A):
            try:
                self.logger.info("computing **ICE AREA TENDENCY (SPATIAL RATE)**; units m/yr")
                METS[f"{ice_type}MAR_YR"] = ((I_MAT * A).mean(dim=time_dim) / 31_536_000).load()
            except Exception as e:
                self.logger.warning(f"MAR_YR failed: {e}")

        # --------- Skill / seasonal / persistence (conditional) ----------
        IA_skill = {}
        IA_seasonal = {}
        if IA is not None and isinstance(IA, xr.DataArray) and valid(IA):
            try:
                if ice_type == "FI":
                    self.logger.info(f"loading FI obs: {self.AF_FI_dict['P_AF2020_FIA']}")
                    IA_obs = xr.open_dataset(self.AF_FI_dict["P_AF2020_FIA"], engine="netcdf4")["AF2020"]
                elif ice_type in ("PI", "SI"):
                    NSIDC = self.compute_NSIDC_metrics()
                    IA_obs = NSIDC["SIA"]
                else:
                    IA_obs = None

                if IA_obs is not None:
                    IA_obs = IA_obs.load()
                    result = self.compute_skill_statistics(IA, IA_obs)
                    IA_skill = result if isinstance(result, dict) else {}
            except Exception as e:
                self.logger.warning(f"compute_skill_statistics failed: {e}")
                IA_skill = {}

            try:
                result = self.compute_seasonal_statistics(IA, stat_name=f"{ice_type}A")
                if isinstance(result, tuple):
                    IA_seasonal = result[0] if isinstance(result[0], dict) else {}
                elif isinstance(result, dict):
                    IA_seasonal = result
                else:
                    IA_seasonal = {}
            except Exception as e:
                self.logger.warning(f"compute_seasonal_statistics failed: {e}")
                IA_seasonal = {}
        else:
            self.logger.info("skipping skill/seasonal stats: IA not available")

        IP_stab, IP_dist = {}, {}
        if ice_type == "FI":
            if (I_mask is not None) and (A is not None):
                try:
                    result = self.persistence_stability_index(I_mask, A)
                    IP_stab = result if isinstance(result, dict) else {}
                except Exception as e:
                    self.logger.warning(f"persistence_stability_index failed: {e}")
                    IP_stab = {}
            else:
                self.logger.info("skipping persistence_stability_index: missing mask or area")

            if IP is not None:
                try:
                    result = self.persistence_ice_distance_mean_max(IP)
                    IP_dist = result if isinstance(result, dict) else {}
                except Exception as e:
                    self.logger.warning(f"persistence_ice_distance_mean_max failed: {e}")
                    IP_dist = {}
            else:
                self.logger.info("skipping persistence distance: IP not available")

        # --------- Build output dataset ----------
        DS_METS = xr.Dataset()
        for k, v in METS.items():
            if isinstance(v, xr.DataArray):
                DS_METS[k] = v
            else:
                DS_METS[k] = xr.DataArray(v, dims=())

        # --------- Merge metadata / summary dict ----------
        IA_seasonal = IA_seasonal if isinstance(IA_seasonal, dict) else {}
        IP_stab     = IP_stab if isinstance(IP_stab, dict) else {}
        IP_dist     = IP_dist if isinstance(IP_dist, dict) else {}
        IA_skill    = IA_skill if isinstance(IA_skill, dict) else {}

        summary = {**{f"{ice_type}A_{k}": v for k, v in IA_seasonal.items()},
                **IP_stab, **IP_dist, **IA_skill, **self.sim_config}

        for k, v in summary.items():
            if k in self.sim_config:
                DS_METS.attrs[k] = v
            else:
                DS_METS[k] = xr.DataArray(v, dims=())

        # --------- Save to Zarr ----------
        if P_mets_zarr:
            DS_METS = self._clean_zarr_chunks(DS_METS)
            DS_METS.to_zarr(P_mets_zarr, mode="w", consolidated=True, zarr_format=2)
            self.logger.info(f"Metrics written to {P_mets_zarr}")

        return DS_METS

    def _subset_and_pad_time(self, ds: xr.Dataset, dt0_str: str, dtN_str: str, 
                             time_dim: str = "time") -> xr.Dataset:
        """
        Reindex a dataset to a requested inclusive time axis, padding with NaNs as needed.

        This helper is used when loading metrics computed over a slightly different time span
        than the current analysis window. It creates a desired time coordinate from `dt0_str`
        to `dtN_str` (inclusive) using the dataset’s native time type and step, then reindexes
        onto that coordinate. Missing days are introduced as NaNs.

        Parameters
        ----------
        ds : xr.Dataset
            Dataset to subset/reindex. If `time_dim` is absent, `ds` is returned unchanged.
        dt0_str, dtN_str : str
            Start/end dates ("YYYY-MM-DD") for the desired inclusive time axis.
        time_dim : str, default "time"
            Name of the time coordinate/dimension.

        Returns
        -------
        xr.Dataset
            Dataset reindexed to the desired time axis. Values are preserved where the
            original dataset overlaps, and NaNs are inserted where the dataset did not
            originally contain those times.

        Raises
        ------
        ValueError
            If `dt0_str` is later than `dtN_str`.
        RuntimeError
            If the time coordinate appears to be cftime/object and `cftime` is unavailable.

        Notes
        -----
        - Supports both numpy datetime64 time coordinates and cftime calendars.
        - The step is inferred from the first two time points when available; otherwise
          defaults to one day.
        """
        if time_dim not in ds.dims and time_dim not in ds.coords:
            # Nothing to do
            return ds
        if dt0_str is None or dtN_str is None:
            return ds
        t = ds[time_dim]
        if t.size == 0:
            return ds
        # Determine whether we're working with numpy datetime64 or cftime objects
        t0 = t.values[0]
        is_datetime64 = np.issubdtype(t.values.dtype, np.datetime64)
        # Build start/end in matching type
        if is_datetime64:
            start = np.datetime64(dt0_str)
            end   = np.datetime64(dtN_str)
        else:
            # cftime (dtype typically object)
            try:
                import cftime  # noqa: F401
            except Exception as e:
                raise RuntimeError("Dataset time coordinate appears to be cftime/object, but cftime is not available.") from e
            # Try to preserve calendar if present; default to 'standard'
            cal = t.encoding.get("calendar", None) or t.attrs.get("calendar", None) or "standard"
            start = xr.cftime_range(start=dt0_str, periods=1, calendar=cal)[0]
            end   = xr.cftime_range(start=dtN_str, periods=1, calendar=cal)[0]
        # Sanity
        if start > end:
            raise ValueError(f"Requested start > end: {dt0_str} > {dtN_str}")
        # Infer step from dataset time (fallback: 1 day)
        if t.size >= 2:
            step = t.values[1] - t.values[0]
        else:
            step = np.timedelta64(1, "D") if is_datetime64 else __import__("datetime").timedelta(days=1)
        # Build the requested time axis (inclusive)
        if is_datetime64:
            # ensure inclusive end (np.arange is end-exclusive)
            desired = np.arange(start, end + step, step)
        else:
            desired = []
            cur = start
            # robust loop for cftime + timedelta
            while cur <= end:
                desired.append(cur)
                cur = cur + step
            desired = np.array(desired, dtype=object)
        # Reindex onto the desired axis -> pads with NaNs where outside original coverage
        # (keeps original values where they overlap)
        ds_out = ds.reindex({time_dim: desired})
        return ds_out
    
    def load_computed_metrics(self,
                              class_method   : str  = "binary-days",  # "raw", "rolling-mean", "binary-days"
                              BorC2T_type    : str  = None,
                              ice_type       : str  = None,
                              ispd_thresh    : str  = None,
                              zarr_directory : str  = None,
                              clip_to_self   : bool = True,
                              time_dim       : str  = "time"):
        """
        Load a previously computed metrics Zarr store and optionally clip/pad to the current analysis window.

        This method resolves the expected metrics path based on:
          - the requested ice type (FI/PI/SI),
          - the classification method (raw / rolling-mean / binary-days),
          - the staggering→T strategy (BorC2T_type),
          - the speed threshold label (ispd_thresh),
          - the configured directory layout (via `define_classification_dir()` and related helpers).

        Parameters
        ----------
        class_method : {"raw","rolling-mean","binary-days"}, default "binary-days"
            Classification product to load for FI/PI. Used to build `self.FI_class`.
        BorC2T_type : str, optional
            Velocity staggering token(s) used for classification, e.g. "Tb", "Tx", "Tc".
            Defaults to `self.BorC2T_type`.
        ice_type : str, optional
            Ice type to load ("FI", "PI", or "SI"). Defaults to `self.ice_type`.
        ispd_thresh : str, optional
            Speed threshold label (often embedded in directory naming). Defaults to `self.ispd_thresh`.
        zarr_directory : str, optional
            Root Zarr directory. Defaults to `self.D_zarr`.
        clip_to_self : bool, default True
            If True, reindex to `[self.dt0_str, self.dtN_str]` using `_subset_and_pad_time()`.
        time_dim : str, default "time"
            Name of the time dimension in the stored metrics.

        Returns
        -------
        xr.Dataset
            Metrics dataset loaded from disk, optionally reindexed/padded to the current window.

        Raises
        ------
        FileNotFoundError
            If the resolved metrics store does not exist.
        Exception
            Propagates errors from path-resolution helpers or xarray open operations.

        Notes
        -----
        - This method assumes helper methods exist:
            * define_classification_dir(...)
            * define_fast_ice_class_name(...)
          and that these set `self.D_class` and `self.FI_class` consistently with your on-disk layout.
        """
        BorC2T_type = BorC2T_type    or self.BorC2T_type
        ice_type    = ice_type       or self.ice_type
        ispd_thresh = ispd_thresh    or self.ispd_thresh
        D_zarr      = zarr_directory or self.D_zarr
        P_mets_zarr = self.define_metrics_zarr(D_zarr       = D_zarr,
                                               ice_type     = ice_type,
                                               ispd_thresh  = ispd_thresh, 
                                               BorC2T_type  = BorC2T_type,
                                               class_method = class_method)
        ds = xr.open_dataset(P_mets_zarr)
        if clip_to_self:
            ds = self._subset_and_pad_time(ds, self.dt0_str, self.dtN_str, time_dim=time_dim)
        return ds

    def compute_hemisphere_area_weighted_stress(
        self,
        tau: xr.DataArray,
        A: xr.DataArray,
        spatial_dim_names: list | None = None,
        mask: xr.DataArray | None = None,
        out_units: str = "Pa",                  # "Pa" | "kPa" | "hPa"
        return_abs: bool = True,                # if True, also return |tau|
        min_area_m2: float = 0.0,               # if >0, mask timesteps with too little valid area
        name: str | None = None,
    ) -> xr.Dataset:
        """
        Compute a hemispheric, area-weighted aggregate of a stress component (or magnitude).

        Designed for stress-like diagnostics such as lateral-drag stresses (e.g., KuxE, KuyN)
        whose native units are N/m^2 (= Pa). This yields a hemispheric mean stress time series
        with **sensible physical units** (Pa, hPa, or kPa) and a robust treatment of masks.

        The aggregate is an area-weighted mean over the valid region:

            <tau>_A(t) = sum( tau(t,...) * A(...) ) / sum( A(...) )

        If `return_abs=True`, this also returns the area-weighted mean of |tau|, which is often
        more interpretable for signed stress components that alternate in sign over space.

        Parameters
        ----------
        tau : xr.DataArray
            Stress component field in N/m^2 (Pa). Must include a time dimension and spatial
            dimensions in `spatial_dim_names`.
        A : xr.DataArray
            Grid-cell area (m^2) corresponding to `tau`'s grid (e.g., earea for E-grid,
            narea for N-grid). Broadcastable to tau's spatial dims.
        spatial_dim_names : list[str], optional
            Spatial dimensions to reduce over (e.g., ["nj", "ni"]). Defaults to
            `self.CICE_dict["spatial_dims"]`.
        mask : xr.DataArray, optional
            Boolean mask (True=keep) broadcastable to tau. Use this to restrict aggregation
            to a region (e.g., fast-ice mask, SIC threshold mask, coastal mask, etc.).
            If None, uses all finite tau cells.
        out_units : {"Pa","hPa","kPa"}, default "Pa"
            Unit for the returned time series.
        return_abs : bool, default True
            If True, return both signed mean stress and mean absolute stress.
        min_area_m2 : float, default 0.0
            If >0, timesteps where the valid area sum is <= min_area_m2 are masked (set to NaN).
            Helpful to prevent noisy values when the valid mask becomes tiny.
        name : str, optional
            Name root used for output variables. Defaults to `tau.name` if present, else "tau".

        Returns
        -------
        xr.Dataset
            Dataset containing:
            - f"{name}_mean"     : area-weighted mean stress (signed)
            - f"{name}_abs_mean" : area-weighted mean |stress| (if return_abs=True)
            - "valid_area_m2"    : total valid area used in weighting

            Units of *_mean variables are Pa/hPa/kPa as requested.

        Notes
        -----
        - N/m^2 is exactly Pascal. So converting to hPa uses /100, and to kPa uses /1000.
        - Prefer |tau| (absolute mean) when comparing “strength” of drag across experiments,
        because signed components can cancel spatially.
        - If you want a vector magnitude across components (e.g., sqrt(Kux^2+Kuy^2)), compute
        that first and pass as `tau`.
        """
        spatial_dim_names = spatial_dim_names if spatial_dim_names is not None else self.CICE_dict["spatial_dims"]

        nm = name or getattr(tau, "name", None) or "tau"

        self.logger.info(f"Computing **AREA-WEIGHTED STRESS AGGREGATE** for '{nm}'")
        self.logger.debug(f"  • out_units            : {out_units}")
        self.logger.debug(f"  • return_abs           : {return_abs}")
        self.logger.debug(f"  • min_area_m2          : {min_area_m2:.3e} m^2")
        self.logger.debug(f"  • Spatial dimensions   : {spatial_dim_names}")
        self.logger.debug(f"  • tau units (attrs)    : {tau.attrs.get('units', 'unknown')}")
        self.logger.debug(f"  • tau long_name        : {tau.attrs.get('long_name', 'unknown')}")

        self.logger.debug("\n[Stress Aggregation Steps]\n"
                        "  1. Define valid mask: finite(tau) & optional user mask\n"
                        "  2. Compute area-weighted mean: sum(tau*A)/sum(A)\n"
                        "  3. Optionally compute area-weighted mean absolute stress: sum(|tau|*A)/sum(A)\n"
                        "  4. Convert Pa -> requested units (Pa/hPa/kPa)\n"
                        "  5. Optionally mask timesteps where valid area is too small")

        # base validity mask
        valid = np.isfinite(tau)
        if mask is not None:
            valid = valid & mask

        w = A.where(valid)
        valid_area = w.sum(dim=spatial_dim_names)

        # signed mean
        tau_mean = (tau.where(valid) * w).sum(dim=spatial_dim_names) / valid_area

        # absolute mean
        if return_abs:
            tau_abs_mean = (np.abs(tau.where(valid)) * w).sum(dim=spatial_dim_names) / valid_area

        # unit conversion from Pa
        if out_units not in ("Pa", "hPa", "kPa"):
            raise ValueError(f"Unknown out_units='{out_units}'. Expected 'Pa', 'hPa', or 'kPa'.")

        if out_units == "hPa":
            scale = 1.0 / 100.0
        elif out_units == "kPa":
            scale = 1.0 / 1000.0
        else:
            scale = 1.0

        tau_mean = tau_mean * scale
        tau_mean.attrs["units"] = out_units
        tau_mean.attrs["long_name"] = f"Area-weighted mean stress ({nm})"
        tau_mean.attrs["aggregation"] = f"area_weighted_mean_over={spatial_dim_names}"

        out = {
            f"{nm}_mean": tau_mean,
            "valid_area_m2": valid_area,
        }

        if return_abs:
            tau_abs_mean = tau_abs_mean * scale
            tau_abs_mean.attrs["units"] = out_units
            tau_abs_mean.attrs["long_name"] = f"Area-weighted mean absolute stress ({nm})"
            tau_abs_mean.attrs["aggregation"] = f"area_weighted_mean_over={spatial_dim_names}"
            out[f"{nm}_abs_mean"] = tau_abs_mean

        ds_out = xr.Dataset(out)

        if min_area_m2 and float(min_area_m2) > 0.0:
            self.logger.debug("Applying min_area_m2 gating to output time series")
            gate = ds_out["valid_area_m2"] > min_area_m2
            for v in ds_out.data_vars:
                if v != "valid_area_m2":
                    ds_out[v] = ds_out[v].where(gate)

        return ds_out


    def compute_hemisphere_lateral_drag_stress(
        self,
        KuxE: xr.DataArray | None = None,
        KuxN: xr.DataArray | None = None,
        KuyE: xr.DataArray | None = None,
        KuyN: xr.DataArray | None = None,
        earea: xr.DataArray | None = None,
        narea: xr.DataArray | None = None,
        spatial_dim_names: list | None = None,
        maskE: xr.DataArray | None = None,
        maskN: xr.DataArray | None = None,
        out_units: str = "Pa",
        min_area_m2: float = 0.0,
    ) -> xr.Dataset:
        """
        Convenience wrapper to compute hemispheric aggregates for lateral-drag stress components.

        This computes area-weighted mean and mean-absolute stresses for any subset of:
        KuxE, KuxN, KuyE, KuyN

        using their appropriate grid-cell area weights:
        - E-grid stresses weighted by `earea`
        - N-grid stresses weighted by `narea`

        Parameters
        ----------
        KuxE, KuxN, KuyE, KuyN : xr.DataArray, optional
            Stress component fields (Pa) on E or N grids.
        earea, narea : xr.DataArray, optional
            Grid-cell areas (m^2) for E- and N-grids.
        spatial_dim_names : list[str], optional
            Spatial reduction dims. Defaults to `self.CICE_dict["spatial_dims"]`.
        maskE, maskN : xr.DataArray, optional
            Optional boolean masks for E-grid and N-grid fields.
        out_units : {"Pa","hPa","kPa"}, default "Pa"
            Units for output time series.
        min_area_m2 : float, default 0.0
            Gating threshold for valid area per timestep.

        Returns
        -------
        xr.Dataset
            Dataset containing time series for each provided field:
            - <name>_mean
            - <name>_abs_mean
            - <name>_valid_area_m2

            Plus, if both x and y are present on the same grid, a magnitude aggregate:
            - K_uE_mag_mean / abs_mean (from KuxE,KuyE)
            - K_uN_mag_mean / abs_mean (from KuxN,KuyN)

        Notes
        -----
        - For interpretation across experiments, the **absolute mean** and the **magnitude**
        are usually the most stable.
        """
        spatial_dim_names = spatial_dim_names if spatial_dim_names is not None else self.CICE_dict["spatial_dims"]
        self.logger.info("Computing **HEMISPHERE LATERAL-DRAG STRESS AGGREGATES**")
        self.logger.debug(f"  • out_units          : {out_units}")
        self.logger.debug(f"  • min_area_m2        : {min_area_m2:.3e} m^2")
        self.logger.debug(f"  • Spatial dimensions : {spatial_dim_names}")

        out = []

        # --- E grid components ---
        if (KuxE is not None) and (earea is None):
            raise ValueError("KuxE provided but earea is None.")
        if (KuyE is not None) and (earea is None):
            raise ValueError("KuyE provided but earea is None.")

        if KuxE is not None:
            ds = self.compute_hemisphere_area_weighted_stress(
                tau=KuxE, A=earea, spatial_dim_names=spatial_dim_names,
                mask=maskE, out_units=out_units, return_abs=True,
                min_area_m2=min_area_m2, name="KuxE"
            )
            ds = ds.rename({"valid_area_m2": "KuxE_valid_area_m2"})
            out.append(ds)

        if KuyE is not None:
            ds = self.compute_hemisphere_area_weighted_stress(
                tau=KuyE, A=earea, spatial_dim_names=spatial_dim_names,
                mask=maskE, out_units=out_units, return_abs=True,
                min_area_m2=min_area_m2, name="KuyE"
            )
            ds = ds.rename({"valid_area_m2": "KuyE_valid_area_m2"})
            out.append(ds)

        if (KuxE is not None) and (KuyE is not None):
            KmagE = np.hypot(KuxE, KuyE)
            KmagE.name = "K_uE_mag"
            ds = self.compute_hemisphere_area_weighted_stress(
                tau=KmagE, A=earea, spatial_dim_names=spatial_dim_names,
                mask=maskE, out_units=out_units, return_abs=True,
                min_area_m2=min_area_m2, name="K_uE_mag"
            )
            ds = ds.rename({"valid_area_m2": "K_uE_mag_valid_area_m2"})
            out.append(ds)

        # --- N grid components ---
        if (KuxN is not None) and (narea is None):
            raise ValueError("KuxN provided but narea is None.")
        if (KuyN is not None) and (narea is None):
            raise ValueError("KuyN provided but narea is None.")

        if KuxN is not None:
            ds = self.compute_hemisphere_area_weighted_stress(
                tau=KuxN, A=narea, spatial_dim_names=spatial_dim_names,
                mask=maskN, out_units=out_units, return_abs=True,
                min_area_m2=min_area_m2, name="KuxN"
            )
            ds = ds.rename({"valid_area_m2": "KuxN_valid_area_m2"})
            out.append(ds)

        if KuyN is not None:
            ds = self.compute_hemisphere_area_weighted_stress(
                tau=KuyN, A=narea, spatial_dim_names=spatial_dim_names,
                mask=maskN, out_units=out_units, return_abs=True,
                min_area_m2=min_area_m2, name="KuyN"
            )
            ds = ds.rename({"valid_area_m2": "KuyN_valid_area_m2"})
            out.append(ds)

        if (KuxN is not None) and (KuyN is not None):
            KmagN = np.hypot(KuxN, KuyN)
            KmagN.name = "K_uN_mag"
            ds = self.compute_hemisphere_area_weighted_stress(
                tau=KmagN, A=narea, spatial_dim_names=spatial_dim_names,
                mask=maskN, out_units=out_units, return_abs=True,
                min_area_m2=min_area_m2, name="K_uN_mag"
            )
            ds = ds.rename({"valid_area_m2": "K_uN_mag_valid_area_m2"})
            out.append(ds)

        if not out:
            raise ValueError("No stress fields provided. Supply at least one of KuxE/KuxN/KuyE/KuyN.")

        return xr.merge(out)

    def compute_hemisphere_ice_area_rate(self, 
                                         DAT: xr.DataArray,
                                         IA: xr.DataArray,
                                         A: xr.DataArray,
                                         spatial_dim_names: list | None = None,
                                         mode: str = "fractional",               # "fractional" (1/s or 1/day) OR "absolute" (km^2/s or km^2/day)
                                         out_units: str = "per_day",             # "per_day" | "per_second"
                                         IA_min_m2: float = 5e10,                # denominator floor (m^2) for fractional mode
                                         mask_invalid: bool = True) -> xr.DataArray:
        """
        Compute a hemisphere-aggregated ice-area tendency rate with optional denominator gating.

        This function aggregates a gridded area-tendency field over the hemisphere and returns either:

        (A) an **absolute** integrated tendency:
                d(Area)/dt  [km^2/s] or [km^2/day]

        (B) a **fractional** tendency (normalised by instantaneous integrated ice area):
                (1/Area) d(Area)/dt  [1/s] or [1/day]

        The fractional form is often the most comparable across experiments, but it can become
        numerically unstable when the integrated ice area is small. To mitigate this, a minimum
        area floor `IA_min_m2` is applied in fractional mode (values are masked where IA <= IA_min_m2).

        Parameters
        ----------
        DAT : xr.DataArray
            Local ice-area tendency field, typically in [1/s] (e.g., dynamic area tendency),
            with dimensions including time and the spatial dimensions in `spatial_dim_names`.
            NOTE: `DAT` is assumed to act on the same area basis as `A` (grid-cell area).
        IA : xr.DataArray
            Hemisphere-integrated ice area time series in [m^2], used as the denominator in
            fractional mode. Must share the same time coordinate as `DAT`.
        A : xr.DataArray
            Grid-cell area in [m^2], broadcastable to the spatial dimensions of `DAT`.
        spatial_dim_names : list[str], optional
            Spatial dimensions to sum over (e.g., ["nj", "ni"]). Defaults to `self.CICE_dict["spatial_dims"]`.
        mode : {"fractional","absolute"}, default "fractional"
            - "fractional": return (DAT*A).sum / IA  -> [1/s] or [1/day]
            - "absolute"  : return (DAT*A).sum       -> [m^2/s] converted to [km^2/s] or [km^2/day]
        out_units : {"per_day","per_second"}, default "per_day"
            Output time unit scaling. For "per_day", multiplies by 86400.
        IA_min_m2 : float, default 5e10
            Minimum integrated ice area (m^2) required to compute fractional rates. Values where
            IA <= IA_min_m2 are masked (set to NaN) in fractional mode.
            A typical order-of-magnitude choice is 1e10–1e11 m^2 (~10^4–10^5 km^2).
        mask_invalid : bool, default True
            If True, masks non-finite values after calculation (NaN/inf).

        Returns
        -------
        xr.DataArray
            Hemisphere-aggregated area tendency rate.
            Units:
            - mode="absolute"  : km^2/s or km^2/day
            - mode="fractional": 1/s   or 1/day

        Notes
        -----
        - This function replaces the previous hard-coded `/1e9` scaling and instead returns
        physically interpretable units (km^2/time for absolute, 1/time for fractional).
        - If your `IA` is not in m^2, convert it before calling this function.
        - If `DAT` is not in 1/s, units will differ accordingly; the calculations here assume
        the common CICE convention for area tendencies.
        """
        spatial_dim_names = spatial_dim_names if spatial_dim_names is not None else self.CICE_dict["spatial_dims"]
        self.logger.info("Computing **HEMISPHERE ICE-AREA TENDENCY RATE**")
        self.logger.debug(f"  • mode                          : {mode}")
        self.logger.debug(f"  • out_units                     : {out_units}")
        self.logger.debug(f"  • IA_min_m2 (fractional gating) : {IA_min_m2:.3e} m^2")
        self.logger.debug(f"  • Spatial dimensions            : {spatial_dim_names}")

        # 1) integrate local tendency over area: (1/s)*m^2 -> m^2/s
        self.logger.debug("\n[Area Tendency Aggregation Steps]\n"
                        "  1. Multiply DAT by grid-cell area A and sum over spatial dims -> dA/dt in m^2/s\n"
                        "  2. If mode='fractional', divide by IA (with IA floor) -> fractional rate in 1/s\n"
                        "  3. Convert to requested output units (per_day vs per_second) and (if absolute) km^2 scaling")

        dA_dt_m2_s = (DAT * A).sum(dim=spatial_dim_names)

        if mode not in ("fractional", "absolute"):
            raise ValueError(f"Unknown mode='{mode}'. Expected 'fractional' or 'absolute'.")

        if out_units not in ("per_day", "per_second"):
            raise ValueError(f"Unknown out_units='{out_units}'. Expected 'per_day' or 'per_second'.")

        if mode == "fractional":
            # denominator gating to avoid blow-ups
            IA_safe = IA.where(IA > IA_min_m2)
            rate = dA_dt_m2_s / IA_safe  # 1/s (if DAT is 1/s and IA is m^2)
            rate.attrs["long_name"] = "Hemisphere fractional ice-area tendency rate"
            rate.attrs["units"] = "1/s"
        else:
            # absolute: convert m^2/s -> km^2/s
            rate = dA_dt_m2_s / 1e6
            rate.attrs["long_name"] = "Hemisphere absolute ice-area tendency rate"
            rate.attrs["units"] = "km^2/s"

        if out_units == "per_day":
            rate = rate * 86400.0
            # adjust units string
            if mode == "fractional":
                rate.attrs["units"] = "1/day"
            else:
                rate.attrs["units"] = "km^2/day"

        if mask_invalid:
            rate = rate.where(np.isfinite(rate))

        # keep a breadcrumb for provenance
        rate.attrs["aggregation"] = f"sum_over={spatial_dim_names}"
        if mode == "fractional":
            rate.attrs["IA_min_m2"] = float(IA_min_m2)

        return rate

    def compute_hemisphere_ice_volume_rate(self,
                                        DVT: xr.DataArray,
                                        IV: xr.DataArray | None,
                                        A: xr.DataArray,
                                        SIC: xr.DataArray | None = None,
                                        spatial_dim_names: list | None = None,
                                        mode: str = "fractional",           # "fractional" | "absolute"
                                        out_units: str = "per_day",         # "per_day" | "per_second"
                                        IV_min_m3: float = 5e10,            # denominator floor (m^3) for fractional mode
                                        sic_threshold: float | None = None, # optional mask (SIC > thresh)
                                        mask_invalid: bool = True,
                                        assume_units: str = "auto"          # "auto" | "cm/day" | "m/day" | "m/s"
                                        ) -> xr.DataArray:
        """
        Compute a hemisphere-aggregated ice-volume tendency rate with sensible units and optional
        denominator gating.

        This function aggregates a gridded "volume tendency" field (CICE history typically stores
        `dvidtt` / `dvidtd` in **cm/day**, representing an equivalent thickness tendency) into:

        (A) an **absolute** integrated volume tendency:
                dV/dt  [km^3/s] or [km^3/day]

        (B) a **fractional** (normalised) tendency:
                (1/V) dV/dt  [1/s] or [1/day]

        The fractional form is often the most comparable across experiments, but it can become
        numerically unstable when the integrated ice volume is small. To mitigate this, a minimum
        volume floor `IV_min_m3` is applied in fractional mode (values are masked where IV <= IV_min_m3).

        Parameters
        ----------
        DVT : xr.DataArray
            Local ice "volume tendency" field. For standard CICE diagnostics:
            - dvidtt / dvidtd typically have units "cm/day" (equivalent thickness tendency).
            Must have dims including time and the spatial dims in `spatial_dim_names`.
        IV : xr.DataArray or None
            Hemisphere-integrated ice volume time series in [m^3], used as the denominator in
            fractional mode. Must share the same time coordinate as `DVT`.
            Required if mode="fractional"; can be None if mode="absolute".
        A : xr.DataArray
            Grid-cell area in [m^2], broadcastable to the spatial dimensions of `DVT`.
            (For CICE dvidtt/dvidtd: thickness tendency * area -> volume tendency.)
        SIC : xr.DataArray, optional
            Sea-ice concentration used to mask contributions. If provided, cells with
            SIC <= sic_threshold are excluded. If None, no SIC masking is applied.
        spatial_dim_names : list[str], optional
            Spatial dimensions to sum over (e.g., ["nj","ni"]). Defaults to `self.CICE_dict["spatial_dims"]`.
        mode : {"fractional","absolute"}, default "fractional"
            - "fractional": return (dV/dt)/IV   -> [1/s] or [1/day]
            - "absolute"  : return dV/dt        -> [km^3/s] or [km^3/day]
        out_units : {"per_day","per_second"}, default "per_day"
            Output time unit. "per_day" scales by 86400 relative to per-second base.
        IV_min_m3 : float, default 5e10
            Minimum integrated ice volume (m^3) required to compute fractional rates.
            Values where IV <= IV_min_m3 are masked (set to NaN) in fractional mode.
            Order-of-magnitude guidance:
            - fast-ice V ~ 1e10–1e11 m^3 (10–100 km^3)
            - sea-ice V  ~ 1e12–1e13 m^3 (1000–10000 km^3)
        sic_threshold : float, optional
            SIC threshold for masking. Defaults to `self.icon_thresh` when SIC is provided.
        mask_invalid : bool, default True
            If True, mask non-finite results (NaN/inf) after calculation.
        assume_units : {"auto","cm/day","m/day","m/s"}, default "auto"
            How to interpret `DVT` units. If "auto", tries to parse `DVT.attrs["units"]`.
            If parsing fails and you have `self.thickness_tendency_to_mps`, it will be used as fallback.

        Returns
        -------
        xr.DataArray
            Hemisphere-aggregated volume tendency rate time series.
            Units:
            - mode="absolute"  : km^3/s or km^3/day
            - mode="fractional": 1/s   or 1/day

        Notes
        -----
        - For CICE dvidtt/dvidtd (cm/day): the physically consistent absolute tendency is:
            dV/dt [m^3/day] = (DVT * 1e-2) * A  summed over space
        then converted to km^3/day by dividing by 1e9.
        - This function does NOT assume any hard-coded scaling like /1e9*86400 unless justified by units.
        """

        spatial_dim_names = spatial_dim_names if spatial_dim_names is not None else self.CICE_dict["spatial_dims"]

        if mode not in ("fractional", "absolute"):
            raise ValueError(f"Unknown mode='{mode}'. Expected 'fractional' or 'absolute'.")
        if out_units not in ("per_day", "per_second"):
            raise ValueError(f"Unknown out_units='{out_units}'. Expected 'per_day' or 'per_second'.")
        if (mode == "fractional") and (IV is None):
            raise ValueError("mode='fractional' requires IV (hemisphere-integrated volume time series in m^3).")

        sic_threshold = sic_threshold if sic_threshold is not None else getattr(self, "icon_thresh", 0.15)

        self.logger.info("Computing **HEMISPHERE ICE-VOLUME TENDENCY RATE**")
        self.logger.debug(f"  • mode                           : {mode}")
        self.logger.debug(f"  • out_units                      : {out_units}")
        self.logger.debug(f"  • IV_min_m3 (fractional gating)  : {IV_min_m3:.3e} m^3")
        self.logger.debug(f"  • SIC mask applied?              : {'yes' if SIC is not None else 'no'}")
        if SIC is not None:
            self.logger.debug(f"  • sic_threshold                  : {sic_threshold:.2f}")
        self.logger.debug(f"  • Spatial dimensions             : {spatial_dim_names}")

        self.logger.debug(
            "\n[Volume Tendency Aggregation Steps]\n"
            "  1. Convert DVT to thickness tendency in m/s (unit-aware)\n"
            "  2. Optionally apply SIC mask (SIC > threshold)\n"
            "  3. Multiply by grid-cell area A and sum over spatial dims -> dV/dt in m^3/s\n"
            "  4. If mode='fractional', divide by IV (with IV floor) -> fractional rate in 1/s\n"
            "  5. Convert to requested output units (per_day vs per_second) and (if absolute) km^3 scaling"
        )

        # ---- unit-aware conversion to m/s ----
        def _to_m_per_s(da: xr.DataArray) -> xr.DataArray:
            if assume_units != "auto":
                u = assume_units.lower()
            else:
                u = str(da.attrs.get("units", "")).strip().lower()

            # Common CICE cases
            if ("cm" in u) and ("/day" in u or "d-1" in u or "day-1" in u):
                # cm/day -> m/s
                return da * 1e-2 / 86400.0
            if ("m" in u) and ("/day" in u or "d-1" in u or "day-1" in u):
                # m/day -> m/s
                return da / 86400.0
            if ("cm" in u) and ("/s" in u or "s-1" in u):
                # cm/s -> m/s
                return da * 1e-2
            if ("m" in u) and ("/s" in u or "s-1" in u):
                # m/s -> m/s
                return da

            # Fallback to your existing helper if available
            if hasattr(self, "thickness_tendency_to_mps"):
                self.logger.debug("  • units not recognised; falling back to self.thickness_tendency_to_mps(DVT)")
                return self.thickness_tendency_to_mps(da)

            raise ValueError(f"Could not infer DVT units='{da.attrs.get('units','')}'. "
                            f"Set assume_units='cm/day'|'m/day'|'m/s' or implement thickness_tendency_to_mps().")

        DVT_mps = _to_m_per_s(DVT)

        # ---- optional SIC mask ----
        if SIC is not None:
            mask = (SIC > float(sic_threshold))
            DVT_mps = DVT_mps.where(mask)

        # ---- integrate: (m/s)*m^2 -> m^3/s ----
        dV_dt_m3_s = (DVT_mps * A).sum(dim=spatial_dim_names)

        if mode == "fractional":
            IV_safe = IV.where(IV > float(IV_min_m3))
            rate = dV_dt_m3_s / IV_safe  # 1/s
            rate.attrs["long_name"] = "Hemisphere fractional ice-volume tendency rate"
            rate.attrs["units"] = "1/s"
            rate.attrs["IV_min_m3"] = float(IV_min_m3)
        else:
            # absolute: m^3/s -> km^3/s
            rate = dV_dt_m3_s / 1e9
            rate.attrs["long_name"] = "Hemisphere absolute ice-volume tendency rate"
            rate.attrs["units"] = "km^3/s"

        if out_units == "per_day":
            rate = rate * 86400.0
            if mode == "fractional":
                rate.attrs["units"] = "1/day"
            else:
                rate.attrs["units"] = "km^3/day"

        if mask_invalid:
            rate = rate.where(np.isfinite(rate))

        # provenance breadcrumbs
        rate.attrs["aggregation"] = f"sum_over={spatial_dim_names}"
        if SIC is not None:
            rate.attrs["sic_threshold"] = float(sic_threshold)

        return rate

    def compute_area_weighted_strength_hpa(self,
                                           SIC: xr.DataArray,
                                           HI: xr.DataArray,
                                           IS: xr.DataArray,
                                           A: xr.DataArray,
                                           spatial_dim_names: list | None = None,
                                           sic_threshold: float | None = None,
                                           hmin: float = 0.05,
                                           assume_IS_units: str = "N/m",          # "N/m" -> divide by HI to get Pa; "Pa" -> use IS directly
                                           mask_invalid: bool = True) -> xr.DataArray:
        """
        Compute an area-weighted mean ice "pressure-equivalent" in hectopascals (hPa).

        This is a robust hemispheric aggregate designed to avoid two common pathologies:
        1) **summing** a pressure-like field instead of averaging (domain-size dependence),
        2) **division-by-small-thickness** spikes (handled via `hmin` and masking).

        The calculation proceeds as:
        - Define a valid-ice mask: SIC > sic_threshold AND HI > hmin AND finite(IS)
        - Convert the internal strength variable to a pressure-like field in Pa:
            * if IS is in N/m, then P = IS / HI  (Pa)
            * if IS is already in Pa, then P = IS
        - Compute an area-weighted mean:
            P_mean = sum(P * A) / sum(A)   over the masked region
        - Convert Pa -> hPa (divide by 100)

        Parameters
        ----------
        SIC : xr.DataArray
            Sea ice concentration [0–1].
        HI : xr.DataArray
            Sea ice thickness [m].
        IS : xr.DataArray
            Internal ice strength field. Expected units depend on `assume_IS_units`:
            - "N/m": common for CICE strength diagnostics; converted to Pa via division by HI.
            - "Pa" : if your diagnostic is already pressure-like.
        A : xr.DataArray
            Grid-cell area [m^2], broadcastable to SIC/HI/IS.
        spatial_dim_names : list[str], optional
            Spatial dimensions to reduce over. Defaults to `self.CICE_dict["spatial_dims"]`.
        sic_threshold : float, optional
            Concentration threshold for including cells. Defaults to `self.icon_thresh`.
        hmin : float, default 0.05
            Minimum thickness (m) required to include a cell. Prevents spikes from dividing by
            very small HI and excludes near-zero ice.
        assume_IS_units : {"N/m","Pa"}, default "N/m"
            Controls conversion to Pa.
        mask_invalid : bool, default True
            If True, masks non-finite values after calculation.

        Returns
        -------
        xr.DataArray
            Area-weighted mean pressure-equivalent in hPa.

        Notes
        -----
        - This is an *analysis* diagnostic: it translates a CICE internal strength-like quantity
        into a pressure-equivalent for hemispheric comparison. Interpret with care.
        - If you want the same mask to be applied across multiple experiments, pass a common mask
        or use a common SIC/HI selection upstream.
        """
        spatial_dim_names = spatial_dim_names if spatial_dim_names is not None else self.CICE_dict["spatial_dims"]
        sic_threshold = sic_threshold if sic_threshold is not None else self.icon_thresh

        self.logger.info("Computing **AREA-WEIGHTED INTERNAL ICE PRESSURE (hPa)**")
        self.logger.debug(f"  • SIC threshold          : {sic_threshold:.2f}")
        self.logger.debug(f"  • HI minimum (hmin)      : {hmin:.3f} m")
        self.logger.debug(f"  • assume_IS_units        : {assume_IS_units}")
        self.logger.debug(f"  • Spatial dimensions     : {spatial_dim_names}")

        self.logger.debug("\n[Area-Weighted Strength Steps]\n"
                        f"  1. Build mask: (SIC > {sic_threshold:.2f}) & (HI > {hmin:.3f}) & finite(IS)\n"
                        "  2. Convert to Pa: P = IS/HI if IS in N/m, else P = IS if already Pa\n"
                        "  3. Area-weighted mean over mask: sum(P*A)/sum(A)\n"
                        "  4. Convert Pa -> hPa by dividing by 100")

        mask = (SIC > sic_threshold) & (HI > hmin) & np.isfinite(IS)

        if assume_IS_units not in ("N/m", "Pa"):
            raise ValueError(f"Unknown assume_IS_units='{assume_IS_units}'. Expected 'N/m' or 'Pa'.")

        if assume_IS_units == "N/m":
            P_pa = (IS / HI).where(mask)
        else:
            P_pa = IS.where(mask)

        w = A.where(mask)

        # avoid 0/0 if mask is empty
        denom = w.sum(dim=spatial_dim_names)
        P_pa_mean = (P_pa * w).sum(dim=spatial_dim_names) / denom

        P_hpa = P_pa_mean / 100.0
        if mask_invalid:
            P_hpa = P_hpa.where(np.isfinite(P_hpa))

        P_hpa.attrs["long_name"] = "Area-weighted internal ice pressure-equivalent"
        P_hpa.attrs["units"] = "hPa"
        P_hpa.attrs["sic_threshold"] = float(sic_threshold)
        P_hpa.attrs["hmin_m"] = float(hmin)
        P_hpa.attrs["assume_IS_units"] = assume_IS_units
        P_hpa.attrs["aggregation"] = f"area_weighted_mean_over={spatial_dim_names}"

        return P_hpa

    def compute_hemisphere_ice_strength(self, SIC, HI, IS,
                                        spatial_dim_names  : list  = None,
                                        sic_threshold      : float = None,
                                        ice_strength_scale : float = 100):
        """
        Compute hemispheric mean sea ice pressure in hectopascals (hPa), based on masked concentration,
        thickness, strength, and area fields.

        This function calculates an area-weighted mean ice thickness (m) by dividing total
        ice volume by total ice-covered area. It then converts internal ice strength (N/m) 
        into an equivalent pressure field (Pa = N/m²), and finally into hectopascals (hPa). 
        Grid cells below a threshold sea ice concentration are masked from the analysis.

        INPUTS:
        -------------------
        SIC               : xr.DataArray; Sea ice concentration (unitless, typically 0–1), per grid cell.
        HI                : xr.DataArray; Sea ice thickness in meters, per grid cell.
        IS                : xr.DataArray; Internal ice strength in N/m, per grid cell.
        A                 : xr.DataArray; Grid cell area in m².
        spatial_dim_names : list, optional; Names of the spatial dimensions to reduce over (e.g., ["nj", "ni"]).
                            If not provided, defaults to values in self.CICE_dict['spatial_dims'].
        sic_threshold     : float, optional; Threshold for masking sea ice concentration (default: self.icon_thresh).

        OUTPUTS:
        -------
        P : xr.DataArray; Sea ice pressure in hectopascals (hPa), with masked cells excluded
            and computed as (IS / mean_thickness) / 100.

        Notes
        -----
        - The output pressure is an approximation derived from CICE's internal strength field and is not a direct model output.
        - Use with caution where sea ice thickness is very small or near-zero.
        """
        sic_threshold     = sic_threshold     if sic_threshold     is not None else self.icon_thresh
        spatial_dim_names = spatial_dim_names if spatial_dim_names is not None else self.CICE_dict['spatial_dims']
        self.logger.info("Computing **INTERNAL ICE PRESSURE** for hemisphere")
        self.logger.debug(f"  • Using SIC threshold               : {sic_threshold:.2f}\n"
                          f"  • Spatial dimension names           : {spatial_dim_names}")
        self.logger.debug("\n[Sea Ice Pressure Computation Steps]\n"
                          f"  1. Apply SIC mask (SIC > {sic_threshold:.2f})\n"
                          f"  2. Divide internal ice strength (N/m = kg*m/s^2) by thickness (m) to get ice pressure (kg/m*s^2 = Pa) then divide by 100 to get hPa")
        mask = SIC > sic_threshold
        HI   = HI.where(mask)
        IS   = IS.where(mask)
        return (IS / HI).sum(dim=spatial_dim_names) / ice_strength_scale

    def compute_sector_ice_area(self, I_mask, area_grid, sector_defs, 
                                GI_area = False):
        """
        Compute ice area per geographic sector and a domain total from an arbitrary ice mask.

        This method generalises sector-based area computations for any ice type (FI/PI/SI)
        by integrating a boolean mask over an area grid within sector bounding boxes.

        Parameters
        ----------
        I_mask : xr.DataArray
            Ice mask (boolean or 0/1) on the same grid as `area_grid`.
        area_grid : xr.DataArray
            Grid-cell area field (typically m²). Must have coordinate variables `lon` and `lat`
            (either as DataArray attributes or attached coords) compatible with the sector
            definitions.
        sector_defs : dict
            Sector definition dictionary. For each sector key, expects:
              sector_defs[sec_name]["geo_region"] == (lon_min, lon_max, lat_min, lat_max)
        GI_area : bool, default False
            If True, include grounded iceberg area in the domain total. This uses
            `compute_grounded_iceberg_area()` (called internally) and adds the result
            (converted to km²) to the total.

        Returns
        -------
        IA_da : xr.DataArray
            Sector ice areas (same units as `area_grid` integration), indexed by "sector".
        IA_tot : float
            Domain total ice area (sum of sectors) plus optional grounded iceberg contribution.

        Notes
        -----
        - Current implementation assumes sectors do not cross the dateline in a way that
          requires wrap-aware longitude masking. If any sector bounds cross the dateline,
          adjust the masking logic accordingly.
        - If `I_mask` or `area_grid` are Dask-backed, each sector sum is computed independently;
          this is typically acceptable because the number of sectors is small.
        """
        if GI_area:
            self.compute_grounded_iceberg_area()
            GI_ttl_area = self.compute_grounded_iceberg_area() / 1e6
            self.logger.info(f"adding {GI_ttl_area} to ice area computation")
        else:
            GI_ttl_area = 0
        sector_names = list(sector_defs.keys())
        ia_values    = []
        for sec_name in sector_names:
            lon_min, lon_max, lat_min, lat_max = sector_defs[sec_name]["geo_region"]
            sec_mask = I_mask.where(
                (area_grid.lon >= lon_min) & (area_grid.lon <= lon_max) &
                (area_grid.lat >= lat_min) & (area_grid.lat <= lat_max)
            )
            tmp = (sec_mask * area_grid).sum(skipna=True)
            if hasattr(tmp, "compute"):
                ia_val = float(tmp.compute())
            else:
                ia_val = tmp.item() if hasattr(tmp, "item") else float(tmp)
            ia_values.append(ia_val)
        IA_da  = xr.DataArray(ia_values, dims=["sector"], coords={"sector": sector_names})
        IA_tot = sum(ia_values) + GI_ttl_area
        return IA_da, IA_tot

    def compute_sector_FIA(self, FI_mask, area_grid, sector_defs, GI_area=False):
        """
        Backwards-compatible wrapper for fast-ice area (FIA) sector computation.

        This method preserves the original FI-specific API by delegating to the
        ice-type-agnostic `compute_sector_ice_area()`.

        Parameters
        ----------
        FI_mask : xr.DataArray
            Fast-ice mask on the same grid as `area_grid`.
        area_grid : xr.DataArray
            Grid-cell area field with lon/lat coordinates.
        sector_defs : dict
            Sector definition dictionary (see `compute_sector_ice_area()`).
        GI_area : bool, default False
            If True, add grounded iceberg area to the domain total.

        Returns
        -------
        (xr.DataArray, float)
            Sector areas and domain total, as returned by `compute_sector_ice_area()`.
        """
        return self.compute_sector_ice_area(FI_mask, area_grid, sector_defs, GI_area=GI_area)

    def compute_hemisphere_ice_area(self, SIC, A,
                                    ice_area_scale            = None,
                                    spatial_dim_names         = None,
                                    sic_threshold             = None,
                                    add_grounded_iceberg_area = None,
                                    grounded_iceberg_area     = None,
                                    region                    = None):
        """
        Compute the total sea ice area (IA) by integrating sea ice concentration (SIC) over the grid area, and optionally 
        including the grounded iceberg area (GI_total_area).

        This method computes the sea ice area by multiplying the sea ice concentration (SIC) by the grid cell area (GC_area),
        summing the result over the specified spatial dimensions, and optionally adding the grounded iceberg area. The final
        result is then scaled by an optional ice area scale factor.

        Parameters:
        -----------
        SIC : xarray.DataArray
            Sea ice concentration field with dimensions (nj, ni) representing the concentration of sea ice at each grid point.
            
        A : xarray.DataArray
            Grid cell area with the same dimensions (nj, ni) as SIC, representing the area of each grid cell in square meters.
            
        ice_area_scale : float, optional
            A scale factor for the sea ice area calculation. If not provided, the default scale is used from the `FIC_scale` attribute.
            
        spatial_dim_names : list of str, optional
            The dimension names over which to sum the sea ice area (typically latitude and longitude). If not provided, the default spatial dimensions 
            are used from the `CICE_dict` attribute.

        sic_threshold : float, optional
            Minimum SIC value to be included in the area calculation. Points with SIC below this threshold are excluded.
            Defaults to `self.icon_thresh`.
            
        add_grounded_iceberg_area : bool, optional
            A flag indicating whether to include the grounded iceberg area in the ice area calculation. Defaults to the class attribute `use_gi`.
            
        grounded_iceberg_area : float, optional
            The grounded iceberg area in square meters to be added to the sea ice area. If not provided, the grounded iceberg area is computed using 
            `compute_grounded_iceberg_area()`.

        Returns:
        --------
        IA : xarray.DataArray
            The total sea ice area, including the optional grounded iceberg area, with the same spatial dimensions as SIC and GC_area.
            The value is returned after summing over the spatial dimensions and scaling by the provided `ice_area_scale`.

        Notes:
        ------
        - The grounded iceberg area is added to the sea ice area calculation if `add_grounded_iceberg_area` is `True`. If no grounded iceberg area 
        is provided, the method will compute it using the `compute_grounded_iceberg_area()` method.
        - The sea ice area is computed as the sum of SIC * GC_area across the specified spatial dimensions, and then scaled by `ice_area_scale`.
        """
        add_grounded_iceberg_area = add_grounded_iceberg_area if add_grounded_iceberg_area is not None else self.use_gi
        sic_threshold             = sic_threshold             if sic_threshold             is not None else self.icon_thresh
        spatial_dim_names         = spatial_dim_names         if spatial_dim_names         is not None else self.CICE_dict['spatial_dims']
        ice_area_scale            = ice_area_scale            if ice_area_scale            is not None else self.FIC_scale
        if region is not None:
            self.logger.info(f"Computing Ice **AREA** for {region} ")
        else:
            self.logger.info("Computing Ice **AREA** ")
        self.logger.debug(f"  • Using SIC threshold               : {sic_threshold:.2f}\n"
                          f"  • Spatial dimension names           : {spatial_dim_names}\n"
                          f"  • Ice area scaling factor           : {ice_area_scale:.2e} (to convert m² → km² or other units)\n"
                          f"  • Include grounded iceberg area     : {add_grounded_iceberg_area}")
        if add_grounded_iceberg_area:
            GI_total_area = (grounded_iceberg_area if grounded_iceberg_area is not None else self.compute_grounded_iceberg_area() )
        else:
            GI_total_area = 0
        self.logger.debug("\n[Sea Ice Area Computation Steps]\n"
                        f"  1. Apply SIC mask (SIC > {sic_threshold:.2f})\n"
                        f"  2. Multiply SIC × area and sum over {spatial_dim_names}\n"
                        f"  3. Add grounded iceberg area (GIA) = {GI_total_area:.2f} m²\n"
                        f"  4. Divide by scale factor {ice_area_scale:.2e} to convert units")
        mask = SIC > sic_threshold
        SIC  = SIC.where(mask)
        A    = A.where(mask)
        IA   = ((SIC * A).sum(dim=spatial_dim_names))
        return (IA + GI_total_area) / ice_area_scale 

    def compute_hemisphere_ice_volume(self, SIC, HI, A,
                                      add_grounded_iceberg_area = None,
                                      ice_volume_scale          = 1e12, 
                                      spatial_dim_names         = None,
                                      sic_threshold             = None,
                                      grounded_iceberg_area     = None):
        """
        Compute total sea ice volume by integrating sea ice concentration, thickness, and grid cell area.

        This method calculates the total sea ice volume as the sum of the product of sea ice concentration (SIC),
        thickness (HI), and grid cell area (GC_area) across the model domain. A SIC threshold can be applied to exclude
        grid cells with low concentrations. Optionally includes grounded iceberg volume and applies a scaling factor 
        for unit conversion.

        Parameters
        ----------
        SIC : xarray.DataArray
            Sea ice concentration (unitless, typically between 0 and 1).
            
        HI : xarray.DataArray
            Sea ice thickness in meters.
            
        A : xarray.DataArray
            Grid cell area in square meters (m²).

        ice_volume_scale : float, optional
            Scale factor for the output volume. Default is `1e12`, converting m³ to 1000 km³.

        spatial_dim_names : list of str, optional
            Names of spatial dimensions to sum over (e.g., ['nj', 'ni']). Defaults to `self.CICE_dict['spatial_dims']`.

        sic_threshold : float, optional
            Minimum SIC value to be included in the volume calculation. Grid cells with SIC ≤ `sic_threshold` are masked.
            Defaults to `self.icon_thresh`.

        Returns
        -------
        float
            Total sea ice volume, optionally including grounded ice and scaled (e.g., in 1000 km³ if default scale is used).

        Notes
        -----
        - SIC, HI, and GC_area must share the same grid and be broadcastable.
        - Grid cells below `sic_threshold` are excluded from the volume calculation.
        - If `self.use_gi` is True, a grounded iceberg area is converted to volume and included.
        """
        add_grounded_iceberg_area = add_grounded_iceberg_area if add_grounded_iceberg_area is not None else self.use_gi
        spatial_dim_names         = spatial_dim_names         if spatial_dim_names         is not None else self.CICE_dict['spatial_dims']
        sic_threshold             = sic_threshold             if sic_threshold             is not None else self.icon_thresh
        self.logger.info("Computing Ice **VOLUME** ")
        self.logger.debug(f"  • Using SIC threshold               : {sic_threshold:.2f}\n"
                          f"  • Spatial dimension names           : {spatial_dim_names}\n"
                          f"  • Ice volume scaling factor         : {ice_volume_scale:.2e} (to convert {1e6:.1e}-km^3)\n"
                          f"  • Include grounded iceberg area     : {add_grounded_iceberg_area}")
        if add_grounded_iceberg_area:
            GI_total_area = (grounded_iceberg_area if grounded_iceberg_area is not None else self.compute_grounded_iceberg_area() )
        else:
            GI_total_area = 0
        self.logger.debug("\n[Sea Ice Volume Computation Steps]\n"
                         f"  1. Apply SIC mask (SIC > {sic_threshold:.2f})\n"
                         f"  2. Multiply SIC × thickness × area and sum over {spatial_dim_names}\n"
                         f"  3. Add grounded iceberg area (GIA) = {GI_total_area:.2f} m²\n"
                         f"  4. Divide by scale factor {ice_volume_scale:.2e} to convert units")
        mask = SIC > sic_threshold
        HI   = HI.where(mask)
        SIC  = SIC.where(mask)
        A    = A.where(mask)
        IV   = (SIC * HI * A).sum(dim=spatial_dim_names)#.chunk({'time': -1}).compute()
        return (IV + GI_total_area) / ice_volume_scale 

    def compute_hemisphere_ice_thickness(self, SIC, HI, A, 
                                         spatial_dim_names = None,
                                         sic_threshold     = None):
        """
        Compute average sea ice thickness weighted by grid cell area and sea ice concentration.

        This method calculates the domain-averaged sea ice thickness as the ratio of the total sea ice volume to 
        the total sea ice area. A SIC threshold is used to exclude low-concentration grid cells from both the numerator 
        and denominator.

        Parameters
        ----------
        HI : xarray.DataArray
            Sea ice thickness in meters.

        SIC : xarray.DataArray
            Sea ice concentration (unitless, typically between 0 and 1).

        GC_area : xarray.DataArray
            Grid cell area in square meters (m²).

        spatial_dim_names : list of str, optional
            Names of spatial dimensions over which to compute the sums. Defaults to `self.CICE_dict['spatial_dims']`.

        sic_threshold : float, optional
            Minimum SIC value to include in the average. Grid cells with SIC ≤ `sic_threshold` are excluded.
            Defaults to `self.icon_thresh`.

        Returns
        -------
        float or xarray.DataArray
            Domain-averaged sea ice thickness, computed as the total volume divided by total area.

        Notes
        -----
        - Computed as: ∑(HI × area) / ∑(SIC × area), with masking by `sic_threshold`.
        - Excludes cells with SIC below threshold from both numerator and denominator.
        - Units are in meters.
        """
        spatial_dim_names = spatial_dim_names if spatial_dim_names is not None else self.CICE_dict['spatial_dims']
        sic_threshold     = sic_threshold     if sic_threshold     is not None else self.icon_thresh
        self.logger.info("Computing Ice **THICKNESS** ")
        self.logger.debug(f"  • Using SIC threshold               : {sic_threshold:.2f}\n"
                          f"  • Spatial dimension names           : {spatial_dim_names}")
        self.logger.debug("\n[Sea Ice Thickness Computation Steps]\n"
                        f"  1. Apply SIC mask (SIC > {sic_threshold:.2f})\n"
                        f"  2. Compute ice area: Multiply SIC × area and sum over {spatial_dim_names}\n"
                        f"  3. Compute ice volume: Multiply SIC × thickness × area and sum over {spatial_dim_names}\n"
                        f"  4. Divide by area / volume = thickness in metres")
        mask = SIC > sic_threshold
        HI   = HI.where(mask)
        SIC  = SIC.where(mask)
        A    = A.where(mask)
        IA   = (SIC * A).sum(dim=spatial_dim_names)
        IV   = (HI * A).sum(dim=spatial_dim_names)
        return (IV / IA)

    def compute_hemisphere_variable_aggregate(self, da, time_coord_name='time', flat_mean=False, da2=None, scale=None):
        """
        Compute a time-aggregated field from a DataArray across the specified time dimension.

        This method calculates the mean of the input variable over time, with an optional
        normalisation step based on the maximum value across the time dimension.
        Designed for Dask-backed DataArrays to support scalable computation.

        Parameters
        ----------
        da : xarray.DataArray
            Input data array to aggregate over time. Typically a geophysical variable such as
            sea ice thickness (m), strength (N/m), or velocity (m/s), with time as one of the dimensions.
        time_coord_name : str, optional
            Name of the time dimension to aggregate over. Default is 'time'.
        flat_mean : bool, optional
            If True, the data has simple mean over time, Default is False.

        Returns
        -------
        xarray.DataArray
            A time-aggregated version of the input array. If `normalise_max` is True,
            the result represents the normalised mean occurrence across the time dimension;
            otherwise, it is the simple time mean.

        Notes
        -----
        - If the specified `time_coord_name` is not found in the input DataArray,
        the original array is returned with a warning.
        - If `normalise_max` is enabled, the output represents a relative frequency-like
        measure, especially useful for thresholded or binary input fields (e.g., fast ice presence).
        - Assumes that all non-time dimensions are to be preserved in the output.
        """
        self.logger.info("Computing Ice **AGGREGATE OVER TIME** (percentage of occurence)")
        if time_coord_name not in da.dims:
            self.logger.warning(f"No time dimension '{time_coord_name}' found. Returning original array.")
            return da
        time_len = da[time_coord_name].sizes.get(time_coord_name, 1)
        if flat_mean:
            return da.mean(dim=time_coord_name)
        elif da2 and scale:
            return (da/da2).sum(dim=time_coord_name) / scale
        else:
            return da.sum(dim=time_coord_name) / time_len

    # SEASONAL METRICS
    @staticmethod
    def _ensure_odd(n: int) -> int:
        return n if n % 2 == 1 else n + 1

    @staticmethod
    def _safe_slope(y: np.ndarray, t_days: np.ndarray, scaling: float) -> float | None:
        """Least-squares slope without polyfit overhead; returns None if flat/NaN."""
        mask = np.isfinite(y) & np.isfinite(t_days)
        if mask.sum() < 2:
            return None
        y = y[mask]
        t = t_days[mask]
        t = t - t[0]  # improves conditioning
        t_mean = t.mean()
        y_mean = y.mean()
        denom = np.sum((t - t_mean) ** 2)
        if denom == 0:
            return None
        slope = np.sum((t - t_mean) * (y - y_mean)) / denom
        return float(slope * scaling)

    @staticmethod
    def _year_index_slices(times: np.ndarray) -> dict[int, np.ndarray]:
        years = pd.DatetimeIndex(times).year.values
        uniq = np.unique(years)
        return {int(y): np.where(years == y)[0] for y in uniq}

    @staticmethod
    def _drop_leap_day(times: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        # Remove Feb 29 for DOY-consistent filters
        di = pd.DatetimeIndex(times)
        keep = ~((di.month == 2) & (di.day == 29))
        return times[keep], y[keep]

    @staticmethod
    def _summary(x):
        x = [v for v in x if v is not None]
        return (float(np.mean(x)), float(np.std(x))) if len(x) else (None, None)

    def compute_seasonal_statistics(self, da,
                                    stat_name           : str = 'FIA',
                                    window              : int = 15,
                                    polyorder           : int = 2,
                                    min_onset_doy       : int = 50,
                                    min_retreat_doy     : int = 240,
                                    growth_range        : tuple[int, int] = (71, 273),
                                    retreat_early_range : tuple[int, int] = (273, 330),
                                    retreat_late_range  : tuple[tuple[int, int], tuple[int, int]] = ((331, 365), (1, 71)),
                                    scaling_factor      : float = 1e6,
                                    drop_leap_day       : bool = True,
                                    return_per_year     : bool = False):
        """
        Compute seasonal timing and rate statistics from a daily 1D time series.

        This method is optimised for long daily time series (multi-year) by:
          - loading the full 1D series once,
          - computing Savitzky–Golay derivatives once per year,
          - using NumPy slicing for seasonal windows (minimising xarray overhead).

        It returns summary statistics (mean and standard deviation across years) for:
          - annual maxima and minima,
          - day-of-year (DOY) of maxima and minima,
          - onset timing (first positive derivative after `min_onset_doy`),
          - retreat timing (last negative curvature after `min_retreat_doy`),
          - duration (onset → retreat, wrap-aware),
          - growth slope (within `growth_range`),
          - retreat slope(s): early, late, and combined.

        Parameters
        ----------
        da : xr.DataArray
            Daily 1D time series with dimension "time" (e.g., ice area in km²).
        stat_name : str, default "FIA"
            Label used for logging and downstream reporting.
        window : int, default 15
            Savitzky–Golay filter window length (days). Will be coerced to an odd integer
            and adjusted upward if too small for `polyorder`.
        polyorder : int, default 2
            Polynomial order for Savitzky–Golay derivatives.
        min_onset_doy : int, default 50
            Earliest DOY at which onset is permitted (limits spurious early detections).
        min_retreat_doy : int, default 240
            Earliest DOY at which retreat is permitted.
        growth_range : tuple[int, int], default (71, 273)
            DOY window used to estimate growth slope.
        retreat_early_range : tuple[int, int], default (273, 330)
            DOY window used to estimate early retreat slope within the same year.
        retreat_late_range : tuple[tuple[int,int], tuple[int,int]], default ((331,365), (1,71))
            Late retreat window that may span the year boundary:
              - first tuple: late-year DOY window in current year
              - second tuple: early-year DOY window in next year
        scaling_factor : float, default 1e6
            Multiplier applied to slope estimates to yield convenient units (e.g., km²/day).
        drop_leap_day : bool, default True
            If True, removes Feb 29 to simplify DOY-consistent filtering.
        return_per_year : bool, default False
            If True, also return a per-year table (pandas.DataFrame) with raw annual values.

        Returns
        -------
        summary : dict[str, float | None]
            Dictionary of mean/std summary metrics across years.
        per_year : pandas.DataFrame, optional
            Returned only if `return_per_year=True`. Indexed by year, containing annual
            extrema, DOYs, and slope estimates.

        Raises
        ------
        ValueError
            If `da` does not have a "time" dimension.

        Notes
        -----
        - The routine iterates over years except the final year (by design) because the
          late-retreat window can require samples from the subsequent year.
        - Slope estimates use a lightweight least-squares implementation (`_safe_slope`)
          instead of `np.polyfit` to reduce overhead.
        """
        from scipy.signal import savgol_filter
        from tqdm         import tqdm
        if 'time' not in da.dims:
            raise ValueError("Input DataArray must have a 'time' dimension")
        self.logger.info(f"computing SEASONAL STATISTICS for {stat_name} ... ")
        # 1) Load once (small 1D array) and optionally drop leap days
        self.logger.debug("   1. loading array ")
        da = da.load()  # safe for 1D daily series
        times = da.time.values
        y_full = da.values.astype(np.float32)
        self.logger.debug("   2. drop leap-days")
        if drop_leap_day:
            times, y_full = self._drop_leap_day(times, y_full)
        self.logger.debug("   3. compute DOY indeces once")
        doy_full = pd.DatetimeIndex(times).dayofyear.values
        idx_by_year = self._year_index_slices(times)
        years_sorted = sorted(idx_by_year.keys())
        self.logger.debug("   4. prepare statistical lists")
        per_year_rows = []
        max_list, min_list = [], []
        grow_list, ret_list = [], []
        ret_early_list, ret_late_list = [], []
        min_doy_list, max_doy_list = [], []
        onset_list = [] 
        duration_list = []
        self.logger.debug("   5. Ensure Savitzky–Golay (savgol) settings are valid")
        window = self._ensure_odd(window)
        if window < polyorder + 2:
            window = polyorder + 3 if (polyorder + 3) % 2 == 1 else polyorder + 4
        self.logger.debug("   6. compute derivatives once per-year")
        for i, year in enumerate(years_sorted[:-1]):
            self.logger.debug(f"     ....{year}")
            idx_this = idx_by_year[year]
            idx_next = idx_by_year[years_sorted[i + 1]]
            t_this = times[idx_this]
            y_this = y_full[idx_this]
            doy_this = doy_full[idx_this]
            # Skip tiny years
            if y_this.size < max(window, polyorder + 2):
                continue
            # Derivatives once per year
            try:
                dy_dt = savgol_filter(y_this, window_length=window, polyorder=polyorder, deriv=1, mode='interp')
                d2y_dt2 = savgol_filter(y_this, window_length=window, polyorder=polyorder, deriv=2, mode='interp')
            except ValueError:
                # window too long for this year, skip
                continue
            # --- Max/Min & their DOY ---
            if np.all(~np.isfinite(y_this)):
                # skip year with all NaN
                continue
            with np.errstate(invalid='ignore'):
                max_idx = int(np.nanargmax(y_this))
                min_idx = int(np.nanargmin(y_this))
            max_val = float(y_this[max_idx])
            min_val = float(y_this[min_idx])
            max_doy = int(doy_this[max_idx])
            min_doy = int(doy_this[min_idx])
            max_list.append(max_val); min_list.append(min_val)
            max_doy_list.append(max_doy); min_doy_list.append(min_doy)
            # --- Onset via derivative (first positive slope after min_onset_doy) ---
            onset_candidates = np.flatnonzero((dy_dt > 0) & (doy_this >= min_onset_doy))
            onset_doy = int(doy_this[onset_candidates[0]]) if onset_candidates.size else None
            # --- Retreat via inflection (last negative curvature after min_retreat_doy) ---
            retreat_candidates = np.flatnonzero((d2y_dt2 < 0) & (doy_this >= min_retreat_doy))
            retreat_doy = int(doy_this[retreat_candidates[-1]]) if retreat_candidates.size else None
            # --- Duration (onset → retreat, handle wrap) ---
            if onset_doy is not None and retreat_doy is not None:
                duration = (retreat_doy - onset_doy) if (retreat_doy > onset_doy) else (365 - onset_doy + retreat_doy)
                duration_list.append(int(duration))
                onset_list.append(int(onset_doy))
            # --- Slopes (growth / retreat early / retreat late / overall retreat) ---
            # Build masks once (NumPy, no xarray .where)
            doy_mask = doy_this  # alias
            g0, g1 = growth_range
            re0, re1 = retreat_early_range
            rl0a, rl1a = retreat_late_range[0]
            rl0b, rl1b = retreat_late_range[1]
            # growth season (this year)
            g_mask = (doy_mask >= g0) & (doy_mask <= g1)
            y_grow = y_this[g_mask]
            t_grow = (t_this[g_mask] - t_this[g_mask][0]) / np.timedelta64(1, 'D') if y_grow.size else np.array([])
            slope_grow = self._safe_slope(y_grow, t_grow, scaling_factor)
            if slope_grow is not None:
                grow_list.append(slope_grow)
            # retreat early (this year)
            e_mask = (doy_mask >= re0) & (doy_mask <= re1)
            y_early = y_this[e_mask]
            t_early = (t_this[e_mask] - t_this[e_mask][0]) / np.timedelta64(1, 'D') if y_early.size else np.array([])
            slope_early = self._safe_slope(y_early, t_early, scaling_factor)
            if slope_early is not None:
                ret_early_list.append(slope_early)
            # retreat late (next year)
            t_next = times[idx_next]
            y_next = y_full[idx_next]
            doy_next = doy_full[idx_next]
            l_mask_a = (doy_mask >= rl0a) & (doy_mask <= rl1a)
            l_mask_b = (doy_next >= rl0b) & (doy_next <= rl1b)
            y_late = np.concatenate([y_this[l_mask_a], y_next[l_mask_b]])
            t_late = np.concatenate([(t_this[l_mask_a] - t_this[l_mask_a][0]) / np.timedelta64(1, 'D') if l_mask_a.any() else np.array([]),
                                      # Continue time axis across year boundary
                                     ((t_next[l_mask_b] - t_this[l_mask_a][0]) / np.timedelta64(1, 'D')) if l_mask_b.any() and l_mask_a.any() else
                                     ((t_next[l_mask_b] - t_next[l_mask_b][0]) / np.timedelta64(1, 'D')) if l_mask_b.any() else np.array([])]) if (l_mask_a.any() or l_mask_b.any()) else np.array([])
            slope_late = self._safe_slope(y_late, t_late, scaling_factor)
            if slope_late is not None:
                ret_late_list.append(slope_late)
            # overall retreat (early + late concatenated)
            if y_early.size or y_late.size:
                y_ret = np.concatenate([y_early, y_late])
                # stitch time so it's strictly increasing
                t_ret = np.concatenate([t_early if y_early.size else np.array([]), (t_late + (t_early[-1] + 1) if y_early.size and y_late.size else t_late)])
                slope_ret = self._safe_slope(y_ret, t_ret, scaling_factor)
                if slope_ret is not None:
                    # Keep sign convention negative (ablation)
                    ret_list.append(-abs(slope_ret))
            if return_per_year:
                per_year_rows.append({'year'                : year,
                                      'max'                 : max_val,   'min'         : min_val,
                                      'doy_max'             : max_doy,   'doy_min'     : min_doy,
                                      'onset_doy'           : onset_doy, 'retreat_doy' : retreat_doy,
                                      'duration_days'       : duration if (onset_doy is not None and retreat_doy is not None) else None,
                                      'growth_slope'        : slope_grow,
                                      'retreat_early_slope' : slope_early,
                                      'retreat_late_slope'  : slope_late,
                                      'retreat_slope'       : slope_ret if slope_ret is not None else None,})
        self.logger.debug("   7. compile summary and return")
        summary = {"Maximum Mean"       : self._summary(max_list)[0],       "Maximum Std"       : self._summary(max_list)[1],
                   "Minimum Mean"       : self._summary(min_list)[0],       "Minimum Std"       : self._summary(min_list)[1],
                   "Growth Mean"        : self._summary(grow_list)[0],      "Growth Std"        : self._summary(grow_list)[1],
                   "Retreat Mean"       : self._summary(ret_list)[0],       "Retreat Std"       : self._summary(ret_list)[1],
                   "Retreat Early Mean" : self._summary(ret_early_list)[0], "Retreat Early Std" : self._summary(ret_early_list)[1],
                   "Retreat Late Mean"  : self._summary(ret_late_list)[0],  "Retreat Late Std"  : self._summary(ret_late_list)[1],
                   "Duration-days Mean" : self._summary(duration_list)[0],  "Duration-days Std" : self._summary(duration_list)[1],
                   "DOY Min Mean"       : self._summary(min_doy_list)[0],   "DOY Min Std"       : self._summary(min_doy_list)[1],
                   "DOY Max Mean"       : self._summary(max_doy_list)[0],   "DOY Max Std"       : self._summary(max_doy_list)[1],
                   "DOY Onset Mean"     : self._summary(onset_list)[0],     "DOY Onset Std"     : self._summary(onset_list)[1],}
        if return_per_year:
            return summary, pd.DataFrame(per_year_rows).set_index('year')
        return summary

    # PERSISTENCE METRICS
    def persistence_stability_index(self, I_mask, A,
                                    persistence_threshold: float = 0.8,
                                    winter_months        : tuple = (5, 6, 7, 8, 9, 10),
                                    area_scale           : float = 1e9) -> dict:
        """
        Compute FIPSI (Fast-Ice Persistence Stability Index) from a binary fast-ice mask.

        Definition
        ----------
        FIPSI = (area of grid cells whose WINTER fast-ice persistence ≥ threshold)
                / (area of grid cells that have fast ice at least once in WINTER),
        where persistence is the fraction of winter days (e.g., May–Oct) that a cell is fast.

        Parameters
        ----------
        I_mask : xarray.DataArray or xarray.Dataset
            Binary fast-ice mask (0/1) with a 'time' dimension; if a Dataset is given,
            a variable named 'FI_mask' is used.
        A : xarray.DataArray or xarray.Dataset
            Grid-cell area on the same grid; if a Dataset is given, the first of
            {'tarea','area','TAREA'} is used. If A has a 'time' dim, the first slice is used.
        persistence_threshold : float, optional
            Minimum winter persistence to classify a cell as persistent (default 0.8).
        winter_months : tuple[int], optional
            Months to include for the winter persistence (default May–Oct).

        Returns
        -------
        dict:
            {'persistence_stability_index': float,   # FIPSI in [0,1] or NaN if undefined
             'area_persistent_winter': float,        # area of persistent winter FI (units of A)
             'area_ever_winter': float               # area with any winter FI (units of A)}

        Notes
        -----
        • Uses the class's `spatial_dims` to sum areas.
        • Assumes `I_mask` is truly binary (0/1). For classification outputs that are boolean, they will be cast to float internally.
        • Designed to be Dask-friendly: only two area reductions are computed eagerly.
        """
        import dask
        self.logger.info("Computing FIPSI: persistent-winter area / ever-winter area")
        I_mask = self._as_da_mask(I_mask).astype("float32")
        A      = self._as_da_area(A)
        # Ensure A has only spatial dims
        if "time" in A.dims:
            A = A.isel(time=0)
            if "time" in A.coords:
                # safer than drop_vars which can produce a Dataset
                A = A.drop_vars("time", errors="ignore")
            A = A.squeeze(drop=True)
        spat_dims = self.CICE_dict["spatial_dims"]
        # restrict to winter months and form persistence (fraction in [0,1])
        in_winter = I_mask["time"].dt.month.isin(list(winter_months))
        FI_winter = I_mask.where(in_winter)
        prstnc    = FI_winter.mean(dim="time", skipna=True)
        # masks for footprints
        persistent_mask = prstnc >= float(persistence_threshold)       # persistent winter FI
        ever_mask       = FI_winter.max(dim="time", skipna=True) > 0   # ever-winter FI
        # # area-weighted sums over spatial dims
        prstnc_A = (A.where(persistent_mask)).sum(dim=spat_dims, skipna=True)
        ttl_A    = (A.where(ever_mask)).sum(dim=spat_dims, skipna=True)
        with dask.config.set(scheduler="threads"):   # or "synchronous"
            prstnc_A = prstnc_A.compute()
            ttl_A    = ttl_A.compute()
        prstnc_A = self._to_float_scalar(prstnc_A)
        ttl_A    = self._to_float_scalar(ttl_A)
        psi      = (prstnc_A / ttl_A) if (np.isfinite(ttl_A) and ttl_A > 0.0) else np.nan
        return {"persistence_stability_index": float(psi),
                "area_persistent_winter"     : prstnc_A/area_scale,
                "area_ever_winter"           : ttl_A/area_scale}

    def _prepare_BAS_coast(self,
                        path_coast_shape: str | None = None,
                        crs_out: str = "EPSG:3031",
                        target_spacing_km: float = 1.0,
                        cache_dir: str | None = None,
                        rebuild_cache: bool = False):
        """
        Load BAS high-res coastline and build a KD-tree of densified coastline points in a metric CRS.

        This implementation is optimised for repeated calls:
        1) It supports an on-disk cache of densified coastline points (NPZ), keyed by:
            (shapefile path, shapefile mtime, crs_out, target_spacing_km).
        2) It uses Shapely 2.x `segmentize()` + `get_coordinates()` to densify in C (fast).
            If Shapely 2.x isn't available, it falls back to a slower Python interpolation loop.

        Parameters
        ----------
        path_coast_shape : str, optional
            Path to the BAS polygon coastline shapefile. Defaults to `self.BAS_dict["P_Ant_Cstln"]`.
        crs_out : str, default "EPSG:3031"
            Target projected CRS for metric distances.
        target_spacing_km : float, default 1.0
            Approx spacing along coastline for densification (km). Clipped to >= 0.1 km (100 m).
        cache_dir : str, optional
            Directory to store NPZ cache files. If None, uses:
            - `self.D_cache/coast_kdtree` if `self.D_cache` exists
            - else `~/.cache/AFIM/coast_kdtree`
        rebuild_cache : bool, default False
            If True, ignore any existing cache and rebuild from shapefile.

        Sets
        ----
        self._coast_kdtree  : scipy.spatial.cKDTree
        self._coast_xy_proj : tuple(ndarray, ndarray) of coastline points (meters)
        self._coast_crs     : str (crs_out)
        self._coast_cache   : str (path to cache file used, if any)
        """
        import os
        import hashlib
        from pathlib import Path

        import numpy as np
        import geopandas as gpd
        from scipy.spatial import cKDTree
        from shapely.ops import unary_union

        # --- inputs / paths ---
        P_shp = Path(path_coast_shape) if path_coast_shape is not None else Path(self.BAS_dict["P_Ant_Cstln"])
        if not P_shp.exists():
            raise FileNotFoundError(f"Coastline shapefile not found: {P_shp}")

        spacing_m = max(100.0, float(target_spacing_km) * 1000.0)  # >= 100 m

        # --- determine cache location ---
        if cache_dir is not None:
            P_cache_dir = Path(cache_dir)
        else:
            if hasattr(self, "D_cache") and self.D_cache is not None:
                P_cache_dir = Path(self.D_cache) / "coast_kdtree"
            else:
                P_cache_dir = Path.home() / ".cache" / "AFIM" / "coast_kdtree"
        P_cache_dir.mkdir(parents=True, exist_ok=True)

        # --- build cache key (path + mtime + crs + spacing) ---
        mtime = os.path.getmtime(P_shp)
        sig = f"{str(P_shp.resolve())}|{int(mtime)}|{crs_out}|{spacing_m:.3f}"
        h = hashlib.md5(sig.encode("utf-8")).hexdigest()[:12]
        P_npz = P_cache_dir / f"bas_coast_{h}.npz"

        self.logger.info("Preparing BAS coastline KD-tree")
        self.logger.debug(f"  • shapefile         : {P_shp}")
        self.logger.debug(f"  • crs_out           : {crs_out}")
        self.logger.debug(f"  • target_spacing_km : {target_spacing_km}")
        self.logger.debug(f"  • spacing_m         : {spacing_m:.1f}")
        self.logger.debug(f"  • cache_file        : {P_npz}")
        self.logger.debug(f"  • rebuild_cache     : {rebuild_cache}")

        # --- fast path: load cached points ---
        if P_npz.exists() and not rebuild_cache:
            self.logger.info(f"Coast cache hit: {P_npz.name}")
            z = np.load(P_npz)
            x_coast = z["x"].astype(float, copy=False)
            y_coast = z["y"].astype(float, copy=False)
            kdt = cKDTree(np.column_stack([x_coast, y_coast]))

            self._coast_kdtree  = kdt
            self._coast_xy_proj = (x_coast, y_coast)
            self._coast_crs     = crs_out
            self._coast_cache   = str(P_npz)
            self.logger.info(f"Loaded coastline points: {x_coast.size:,} (cached)")
            return

        # --- slow path: read + densify + cache ---
        self.logger.info(f"Coast cache miss: building from shapefile {P_shp.name}")
        gdf = gpd.read_file(P_shp)

        # CRS handling
        if gdf.crs is None:
            self.logger.warning("Input CRS is None; assuming EPSG:4326 (lon/lat).")
            gdf = gdf.set_crs("EPSG:4326", allow_override=True)
        gdf = gdf.to_crs(crs_out)

        # Prefer union of BOUNDARIES (linework) instead of polygon union (usually faster)
        self.logger.info("Extracting coastline boundary linework...")
        boundary_series = gdf.geometry.boundary

        # Merge boundaries to remove duplicates (still can be expensive, but less than polygon union)
        self.logger.info("Merging boundary linework...")
        boundary_geom = unary_union(boundary_series.values)

        # Densify
        # Shapely 2 path: segmentize + get_coordinates (fast)
        try:
            import shapely  # shapely 2.x provides segmentize/get_coordinates at top-level

            have_segmentize = hasattr(shapely, "segmentize") and hasattr(shapely, "get_coordinates")
            self.logger.debug(f"  • shapely version   : {getattr(shapely, '__version__', 'unknown')}")
            self.logger.debug(f"  • segmentize avail? : {have_segmentize}")

            if have_segmentize:
                self.logger.info("Densifying coastline using shapely.segmentize (fast path)...")
                densified = shapely.segmentize(boundary_geom, spacing_m)
                coords = shapely.get_coordinates(densified)  # (N,2)
                if coords.size == 0:
                    raise RuntimeError("Shapely densification produced zero coordinates.")
                x_coast = coords[:, 0].astype(float, copy=False)
                y_coast = coords[:, 1].astype(float, copy=False)
            else:
                raise ImportError("segmentize/get_coordinates not available")

        except Exception as e:
            # Fallback: Python interpolation along linework (slower)
            self.logger.warning(f"Falling back to Python densification (shapely2 not available or failed): {e}")

            from shapely.geometry import LineString, MultiLineString, GeometryCollection

            def _iter_lines(g):
                if g is None:
                    return
                gt = getattr(g, "geom_type", "").lower()
                if gt == "linestring":
                    yield g
                elif gt == "multilinestring":
                    for gg in g.geoms:
                        yield gg
                elif gt == "geometrycollection":
                    for gg in g.geoms:
                        yield from _iter_lines(gg)
                else:
                    # last resort: try boundary
                    try:
                        yield from _iter_lines(g.boundary)
                    except Exception:
                        return

            xs, ys = [], []
            n_lines = 0
            for ln in _iter_lines(boundary_geom):
                try:
                    if ln.length <= 0 or not np.isfinite(ln.length):
                        continue
                    n_lines += 1
                    n = max(2, int(np.ceil(ln.length / spacing_m)) + 1)
                    dists = np.linspace(0.0, float(ln.length), n)
                    for d in dists:
                        p = ln.interpolate(d)
                        x, y = p.coords[0]
                        xs.append(x)
                        ys.append(y)
                except Exception:
                    continue

            if len(xs) == 0:
                raise RuntimeError("Fallback densification produced zero points.")
            x_coast = np.asarray(xs, dtype=float)
            y_coast = np.asarray(ys, dtype=float)
            self.logger.info(f"Fallback densified {n_lines} lines -> {x_coast.size:,} points")

        # Optional de-duplication (helps reduce KD-tree size if union missed duplicates)
        # Keep it cheap: round to 1 m and unique.
        self.logger.debug("De-duplicating coastline points (rounded to 1 m)...")
        coords_i = np.column_stack([np.round(x_coast, 0), np.round(y_coast, 0)]).astype(np.int64, copy=False)
        # Unique rows (stable-ish)
        _, idx = np.unique(coords_i, axis=0, return_index=True)
        idx.sort()
        x_coast = x_coast[idx]
        y_coast = y_coast[idx]

        if x_coast.size == 0:
            raise RuntimeError("Coastline point set is empty after processing.")

        self.logger.info(f"Built coastline point set: {x_coast.size:,} pts @ ~{target_spacing_km} km spacing")

        # Build KD-tree
        kdt = cKDTree(np.column_stack([x_coast, y_coast]))

        # Cache to disk
        try:
            np.savez_compressed(P_npz, x=x_coast.astype(np.float64), y=y_coast.astype(np.float64))
            self.logger.info(f"Wrote coast cache: {P_npz}")
        except Exception as e:
            self.logger.warning(f"Failed to write coast cache '{P_npz}': {e}")

        # Cache on self
        self._coast_kdtree  = kdt
        self._coast_xy_proj = (x_coast, y_coast)
        self._coast_crs     = crs_out
        self._coast_cache   = str(P_npz)

    def persistence_ice_distance_mean_max(self, ice_prstnc,
                                          persistence_threshold : float = 0.8,
                                          path_coast_shape      : str = None,
                                          crs_out               : str = "EPSG:3031"):
        """
        Compute mean and maximum distance from the Antarctic coastline for persistent ice.

        This method measures how far persistent ice (e.g., persistent fast ice) extends
        from the coastline by:
          1) thresholding a persistence field `ice_prstnc` in [0, 1],
          2) projecting grid coordinates into a metric CRS (default EPSG:3031),
          3) querying nearest coastline points via a pre-built KD-tree, and
          4) returning mean and maximum distances (km) of persistent grid cells.

        Parameters
        ----------
        ice_prstnc : xr.DataArray
            Persistence field in [0, 1], typically computed from a binary mask over a
            seasonal window (e.g., May–Oct). Must be on the same spatial grid as `self.G_t`.
        persistence_threshold : float, default 0.8
            Cells with persistence >= threshold are treated as “persistent” for distance
            computations.
        path_coast_shape : str, optional
            Path to a coastline shapefile. Defaults to `self.BAS_dict["P_Ant_Cstln"]`.
        crs_out : str, default "EPSG:3031"
            Output CRS used for distance computations. Must be a metric CRS.

        Returns
        -------
        dict
            Dictionary with keys:
              - "persistence_mean_distance": mean distance (km)
              - "persistence_max_distance" : max distance (km)

        Raises
        ------
        KeyError
            If longitude/latitude fields cannot be found in `self.G_t`.
        ValueError
            If `ice_prstnc` shape does not match the loaded grid coordinates.
        Exception
            Propagates errors from coordinate transforms or KD-tree construction.

        Notes
        -----
        - The coastline KD-tree is cached on `self` as `_coast_kdtree` and rebuilt only if
          `crs_out` changes.
        - The method forces local computation of `ice_prstnc` (via `.compute()`) using the
          threaded scheduler to avoid distributing a large Dask graph.
        """
        import dask
        from pyproj import Transformer
        P_shp = path_coast_shape if path_coast_shape is not None else self.BAS_dict["P_Ant_Cstln"]
        spat_dims = self.CICE_dict["spatial_dims"]
        self.logger.info("Computing persistence distance (USNIC coastline KD-Tree)")
        with dask.config.set(scheduler="threads"):
            ice_prstnc = ice_prstnc.compute()  # instead of .load()
        if not getattr(self, "bgrid_loaded", False):
            self.load_cice_grid(slice_hem=True)
        Gt       = self.G_t
        lon_name = 'TLON' if 'TLON' in Gt else ('lon' if 'lon' in Gt else None)
        lat_name = 'TLAT' if 'TLAT' in Gt else ('lat' if 'lat' in Gt else None)
        if lon_name is None or lat_name is None:
            raise KeyError("Could not find lon/lat fields (TLON/TLAT or lon/lat) in self.G_t.")
        grid_shape = tuple(ice_prstnc.sizes[d] for d in spat_dims)
        if grid_shape != Gt[lon_name].shape:
            self.logger.warning("Grid mismatch detected; reloading b-grid with slice_hem=True")
            self.load_cice_grid(slice_hem=True)
            Gt = self.G_t
            if tuple(ice_prstnc.sizes[d] for d in spat_dims) != Gt[lon_name].shape:
                raise ValueError("ice_prstnc and grid coordinate shapes still mismatch after reload.")
        # Prepare coastline KD-tree once (unchanged)
        need_kdtree = (not hasattr(self, "_coast_kdtree")) or (getattr(self, "_coast_crs", None) != crs_out)
        if need_kdtree:
            self._prepare_BAS_coast(P_shp, crs_out=crs_out, target_spacing_km=1.0)
        # Build persistent mask from fully-realized array
        prst = (ice_prstnc >= float(persistence_threshold)).values
        if not np.any(prst):
            self.logger.warning("No persistent fast-ice cells found at this threshold.")
            return {"persistence_mean_distance": float('nan'),
                    "persistence_max_distance":  float('nan')}
        # cache dict on self
        if not hasattr(self, "_grid_xy_cache"):
            self._grid_xy_cache = {}
        # key includes CRS + grid shape + hemisphere slicing state (shape is usually enough)
        proj_key = (crs_out, lon_name, lat_name, tuple(Gt[lon_name].shape))
        if proj_key in self._grid_xy_cache:
            self.logger.debug(f"Using cached projected grid coords for {crs_out}")
            x_all, y_all = self._grid_xy_cache[proj_key]
        else:
            self.logger.info(f"Projecting grid lon/lat -> {crs_out} (cached thereafter)")
            lon = np.asarray(Gt[lon_name].values)
            lat = np.asarray(Gt[lat_name].values)
            transformer = Transformer.from_crs("EPSG:4326", crs_out, always_xy=True)
            x_all, y_all = transformer.transform(lon, lat)
            self._grid_xy_cache[proj_key] = (x_all, y_all)
        # lon = np.asarray(Gt[lon_name].values)
        # lat = np.asarray(Gt[lat_name].values)
        # transformer = Transformer.from_crs("EPSG:4326", crs_out, always_xy=True)
        x_all, y_all = transformer.transform(lon, lat)
        iy, ix = np.where(prst)  # spat_dims order (nj, ni)
        x_p = x_all[iy, ix]
        y_p = y_all[iy, ix]
        if x_p.size == 0:
            return {"persistence_mean_distance": float('nan'),
                    "persistence_max_distance":  float('nan')}
        d_m, _ = self._coast_kdtree.query(np.column_stack([x_p, y_p]), k=1)
        d_km = d_m / 1000.0
        return {"persistence_mean_distance": float(np.nanmean(d_km)),
                "persistence_max_distance":  float(np.nanmax(d_km))}

    # INTER-COMPARISONS
    @staticmethod
    def _skill_stats(x, y):
        from sklearn.metrics import mean_squared_error, mean_absolute_error
        from scipy.stats     import pearsonr
        valid = (~np.isnan(x)) & (~np.isnan(y))
        x, y  = x[valid], y[valid]
        return {"Bias"     : np.mean(x - y),
                "RMSE"     : np.sqrt(mean_squared_error(y, x)),
                "MAE"      : mean_absolute_error(y, x),
                "Corr"     : pearsonr(x, y)[0],
                "SD_Model" : np.std(x),
                "SD_Obs"   : np.std(y)}

    @staticmethod
    def _safe_load_array_to_memory(x):
        # --- Force local computation to avoid serializing Dask HLGs to the cluster ---
        # Keep memory light and avoid shipping graphs to distributed workers
        try:
            return x.astype("float32").compute(scheduler="single-threaded")
        except Exception:
            return x.compute(scheduler="single-threaded") if hasattr(x, "compute") else x
        
    def compute_skill_statistics(self, mod, obs, min_points=100):
        """
        Compute model-vs-observation skill statistics for a time series.

        The method attempts a direct time-series comparison by intersecting model and
        observation timestamps. If there are fewer than `min_points` intersecting times,
        it falls back to a climatological comparison using day-of-year means.

        Parameters
        ----------
        mod : xr.DataArray
            Modelled time series (must have a "time" coordinate for direct alignment).
        obs : xr.DataArray
            Observational time series. Ideally shares the same time basis as `mod`.
        min_points : int, default 100
            Minimum number of intersecting timestamps required to compute direct time-series
            statistics. Otherwise, a climatology comparison is used.

        Returns
        -------
        dict
            Dictionary of skill metrics:
              - Bias
              - RMSE
              - MAE
              - Corr
              - SD_Model
              - SD_Obs

            Returns NaNs for all metrics if time alignment fails and climatology cannot be formed.

        Notes
        -----
        - Uses `_safe_load_array_to_memory()` to compute Dask-backed series locally and
          minimise graph shipping.
        - Climatology fallback uses `groupby("time.dayofyear").mean("time")` for both series.
        """
        mod = self._safe_load_array_to_memory(mod)
        obs = self._safe_load_array_to_memory(obs)
        if hasattr(mod, 'time') and hasattr(obs, 'time'):
            dt_mod  = mod.time.values
            dt_obs  = obs.time.values
            dt_xsct = np.intersect1d(dt_mod, dt_obs)
            if len(dt_xsct) > min_points:
                self.logger.warning(f"more than {min_points} date-times intersected, therefore do stats on time series")
                mod_data = mod.sel(time=dt_xsct).values
                obs_data = obs.sel(time=dt_xsct).values
                return self._skill_stats(mod_data, obs_data)
            else:
                self.logger.warning("date-time mismatch detected — falling back to climatology comparison.")
                mod_clim = mod.groupby('time.dayofyear').mean('time')
                obs_clim = obs.groupby('time.dayofyear').mean('time')
                doy_xsct = np.intersect1d(mod_clim.dayofyear.values, obs_clim.dayofyear.values)
                mod_data = mod_clim.sel(dayofyear=doy_xsct).values
                obs_data = obs_clim.sel(dayofyear=doy_xsct).values
                return self._skill_stats(mod_data, obs_data)
        self.logger.warning("Time dimension not found or could not align. Returning empty stats.")
        return {"Bias": np.nan, "RMSE": np.nan, "MAE": np.nan, "Corr": np.nan,
                "SD_Model": np.nan, "SD_Obs": np.nan}

    def antarctic_year_labels(self, time_da: xr.DataArray) -> xr.DataArray:
        """
        Construct Antarctic-year labels for a time coordinate (Jul→Jun “years”).

        Antarctic-year (AY) grouping is commonly defined such that AY YYYY runs from
        July YYYY through June YYYY+1. This helper returns an integer label "ayear"
        for each timestamp:
          - months Jan–Jun are assigned to the previous calendar year,
          - months Jul–Dec are assigned to the current calendar year.

        Parameters
        ----------
        time_da : xr.DataArray
            Time coordinate DataArray supporting `.dt.year` and `.dt.month`.

        Returns
        -------
        xr.DataArray
            Integer DataArray named "ayear" aligned with `time_da`.

        Examples
        --------
        - 2000-06-15 -> ayear = 1999
        - 2000-07-01 -> ayear = 2000
        """
        year  = time_da.dt.year
        month = time_da.dt.month
        ayear = xr.where(month <= 6, year - 1, year)
        return ayear.rename("ayear")

    def extrema_means_antarctic_year(self, da_ts: xr.DataArray):
        """
        Compute mean maxima/minima across Antarctic years for a 1D time series.

        Parameters
        ----------
        da_ts : xr.DataArray
            1D time series with dimension "time" (e.g., fast-ice area in km²).

        Returns
        -------
        dict
            Dictionary containing:
              - "max_mean": mean of per-AY maxima (xr.DataArray scalar)
              - "min_mean": mean of per-AY minima (xr.DataArray scalar)
              - "max_ts"  : per-AY maxima time series (dims: ayear)
              - "min_ts"  : per-AY minima time series (dims: ayear)

        Notes
        -----
        - Uses `antarctic_year_labels()` and xarray groupby reductions.
        - For short 1D series, it is acceptable to fully load/compute in memory.
        """
        # for tiny 1D series, it's fine to fully chunk or even compute up front
        # da_ts = da_ts.chunk({"time": -1})
        ay = self.antarctic_year_labels(da_ts["time"])
        g  = da_ts.groupby(ay)
        max_ts = g.max(skipna=True)                  # dims: ayear
        min_ts = g.min(skipna=True)
        return {
            "max_mean": max_ts.mean("ayear", skipna=True),
            "min_mean": min_ts.mean("ayear", skipna=True),
            "max_ts": max_ts,
            "min_ts": min_ts,
        }

    #######################################################################
