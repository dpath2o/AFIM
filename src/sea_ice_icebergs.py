import os, shutil
import xarray      as xr
import numpy       as np
import pandas      as pd
from pathlib       import Path
from datetime      import datetime


class SeaIceIcebergs:

    def __init__(self, **kwargs):
        return

    ######################################################################
    # small helpers
    ######################################################################
    @staticmethod
    def _first_existing_column(columns, candidates):
        for c in candidates:
            if c in columns:
                return c
        return None

    def _grid_da(self, values, name, dtype=None, attrs=None):
        arr = np.asarray(values, dtype=dtype) if dtype is not None else np.asarray(values)
        da = xr.DataArray(
            data=arr,
            dims=self.spatial_dims,
            coords={
                "lat": (self.spatial_dims, self.G_t["lat"].values),
                "lon": (self.spatial_dims, self.G_t["lon"].values),
            },
            name=name,
        )
        if attrs:
            da.attrs.update(attrs)
        return da

    def _dedup_gi_features(self, df, uid_col="uid"):
        """
        Collapse repeated observations of the same iceberg UID into one feature row.
        """
        if uid_col not in df.columns:
            return df.reset_index(drop=True)

        uid = df[uid_col]
        if uid.isna().all():
            return df.reset_index(drop=True)

        agg = {}
        for c in ["x_m", "y_m", "lon", "lat"]:
            if c in df.columns:
                agg[c] = "mean"
        for c in ["area_m2", "perim_m", "bed_depth"]:
            if c in df.columns:
                agg[c] = "median"
        for c in ["timestamp", "date_range", "orbit", "swath_mode",
                  "source_type", "geometry_type", "source_path"]:
            if c in df.columns:
                agg[c] = "first"

        out = df.groupby(uid_col, as_index=False).agg(agg)
        out["obs_count"] = (
            df.groupby(uid_col).size().reindex(out[uid_col]).to_numpy(dtype=np.int32)
        )
        return out.reset_index(drop=True)

    ######################################################################
    # canonical source loader
    ######################################################################
    def load_grounded_iceberg_source(
        self,
        GI_P_raw=None,
        GI_layer="grounded_icebergs",
        csv_lon_col=None,
        csv_lat_col=None,
        dedup_uid=True,
        normalise_lon_to="0-360",
        use_attr_area=True,
    ):
        """
        Read CSV or GPKG grounded iceberg input and return a canonical feature table.

        Returns columns as available:
          uid, lon, lat, x_m, y_m, area_m2, perim_m, bed_depth, obs_count,
          source_type, geometry_type, source_path
        """
        GI_P_raw = GI_P_raw if GI_P_raw is not None else self.GI_dict["P_raw"]
        P_raw = Path(GI_P_raw)
        ext = P_raw.suffix.lower()

        if ext == ".csv":
            raw = pd.read_csv(P_raw)

            if csv_lon_col is None:
                csv_lon_col = self._first_existing_column(
                    raw.columns,
                    ["X", "x", "lon", "Lon", "LON", "longitude", "LONGITUDE"],
                )
            if csv_lat_col is None:
                csv_lat_col = self._first_existing_column(
                    raw.columns,
                    ["Y", "y", "lat", "Lat", "LAT", "latitude", "LATITUDE"],
                )

            if csv_lon_col is None or csv_lat_col is None:
                raise ValueError(
                    f"Could not identify lon/lat columns in CSV: {P_raw}\n"
                    f"Available columns: {list(raw.columns)}"
                )

            df = self.read_grounded_iceberg_csv(
                P_csv=P_raw,
                lon_col=csv_lon_col,
                lat_col=csv_lat_col,
                normalise_lon_to=normalise_lon_to,
            ).copy()

            if "uid" not in df.columns:
                df["uid"] = np.arange(len(df), dtype=np.int64)

            df["x_m"] = np.nan
            df["y_m"] = np.nan
            if "area_m2" not in df.columns:
                df["area_m2"] = np.nan
            if "perim_m" not in df.columns:
                df["perim_m"] = np.nan
            if "bed_depth" not in df.columns:
                df["bed_depth"] = np.nan

            df["obs_count"] = 1
            df["source_type"] = "csv_point"
            df["geometry_type"] = "Point"
            df["source_path"] = str(P_raw)

        elif ext == ".gpkg":
            df = self.read_grounded_iceberg_gpkg(
                P_gpkg=P_raw,
                layer=GI_layer,
                dedup_uid=False,          # dedup ourselves so we can preserve obs_count
                use_attr_area=use_attr_area,
                normalise_lon_to=normalise_lon_to,
            ).copy()

            df["obs_count"] = 1
            df["source_type"] = "gpkg_polygon"
            df["geometry_type"] = "Polygon"
            df["source_path"] = str(P_raw)

            if dedup_uid:
                df = self._dedup_gi_features(df, uid_col="uid")

        else:
            raise ValueError(
                f"Unsupported grounded iceberg input: {P_raw}\n"
                f"Expected .csv or .gpkg"
            )

        self.GI_features = df.reset_index(drop=True)
        self.logger.info(
            f"Loaded grounded iceberg source: {P_raw} "
            f"(n_features={len(self.GI_features)})"
        )
        return self.GI_features

    ######################################################################
    # mapping to model grid
    ######################################################################
    def assign_grounded_icebergs_to_tgrid(
        self,
        GI_features=None,
        proj_crs="EPSG:3031",
        max_dist_km=None,
    ):
        """
        Map canonical GI features to nearest T-grid cell.

        Notes
        -----
        map_xy_to_tgrid / map_points_to_tgrid return:
            ii -> ni index
            jj -> nj index
        """
        if GI_features is None:
            if not hasattr(self, "GI_features"):
                raise ValueError("No GI_features found. Run load_grounded_iceberg_source first.")
            GI_features = self.GI_features

        self.load_cice_grid(slice_hem=False, build_faces=True)

        df = GI_features.copy()

        has_xy = (
            ("x_m" in df.columns) and ("y_m" in df.columns)
            and np.isfinite(df["x_m"]).any() and np.isfinite(df["y_m"]).any()
        )

        if has_xy:
            ni_i, nj_i, dist_m, valid = self.map_xy_to_tgrid(
                x_m=df["x_m"].to_numpy(dtype=float),
                y_m=df["y_m"].to_numpy(dtype=float),
                proj_crs=proj_crs,
                max_dist_km=max_dist_km,
            )
        else:
            ni_i, nj_i, dist_m, valid = self.map_points_to_tgrid(
                lon_deg=df["lon"].to_numpy(dtype=float),
                lat_deg=df["lat"].to_numpy(dtype=float),
                proj_crs=proj_crs,
                max_dist_km=max_dist_km,
            )

        df["ni_i"] = ni_i.astype(np.int32)
        df["nj_i"] = nj_i.astype(np.int32)
        df["map_dist_m"] = dist_m.astype(np.float64)
        df["map_valid"] = valid.astype(bool)

        df = df[df["map_valid"]].reset_index(drop=True)

        self.GI_assignments = df
        self.GI_clean = df.copy()  # backward compatibility
        self.logger.info(
            f"Mapped grounded icebergs to T-grid: "
            f"n_valid={len(df)}"
        )
        return self.GI_assignments

    ######################################################################
    # cellwise aggregation
    ######################################################################
    def aggregate_grounded_icebergs_on_tgrid(self, GI_assignments=None):
        """
        Aggregate mapped features to model-grid cell diagnostics.
        """
        if GI_assignments is None:
            if not hasattr(self, "GI_assignments"):
                raise ValueError("No GI_assignments found. Run assign_grounded_icebergs_to_tgrid first.")
            GI_assignments = self.GI_assignments

        self.load_cice_grid(slice_hem=False, build_faces=True)

        df = GI_assignments.copy()
        nj, ni = self.G_t["lon"].shape

        if "obs_count" not in df.columns:
            df["obs_count"] = 1
        if "area_m2" not in df.columns:
            df["area_m2"] = np.nan
        if "perim_m" not in df.columns:
            df["perim_m"] = np.nan
        if "bed_depth" not in df.columns:
            df["bed_depth"] = np.nan

        grouped = (
            df.groupby(["nj_i", "ni_i"], sort=False)
              .agg(
                  GI_count=("nj_i", "size"),
                  GI_obs_count=("obs_count", "sum"),
                  GI_unique_uid=("uid", "nunique") if "uid" in df.columns else ("nj_i", "size"),
                  GI_area_sum_m2=("area_m2", "sum"),
                  GI_area_max_m2=("area_m2", "max"),
                  GI_perim_sum_m=("perim_m", "sum"),
                  GI_bed_depth_mean=("bed_depth", "mean"),
                  GI_map_dist_mean_m=("map_dist_m", "mean"),
              )
              .reset_index()
        )

        def zeros_int():
            return np.zeros((nj, ni), dtype=np.int32)

        def zeros_float():
            return np.full((nj, ni), np.nan, dtype=np.float32)

        GI_count = zeros_int()
        GI_obs_count = zeros_int()
        GI_unique_uid = zeros_int()
        GI_area_sum_m2 = zeros_float()
        GI_area_max_m2 = zeros_float()
        GI_perim_sum_m = zeros_float()
        GI_bed_depth_mean = zeros_float()
        GI_map_dist_mean_m = zeros_float()

        jj = grouped["nj_i"].to_numpy(dtype=np.int32)
        ii = grouped["ni_i"].to_numpy(dtype=np.int32)

        GI_count[jj, ii] = grouped["GI_count"].to_numpy(dtype=np.int32)
        GI_obs_count[jj, ii] = grouped["GI_obs_count"].to_numpy(dtype=np.int32)
        GI_unique_uid[jj, ii] = grouped["GI_unique_uid"].to_numpy(dtype=np.int32)
        GI_area_sum_m2[jj, ii] = grouped["GI_area_sum_m2"].to_numpy(dtype=np.float32)
        GI_area_max_m2[jj, ii] = grouped["GI_area_max_m2"].to_numpy(dtype=np.float32)
        GI_perim_sum_m[jj, ii] = grouped["GI_perim_sum_m"].to_numpy(dtype=np.float32)
        GI_bed_depth_mean[jj, ii] = grouped["GI_bed_depth_mean"].to_numpy(dtype=np.float32)
        GI_map_dist_mean_m[jj, ii] = grouped["GI_map_dist_mean_m"].to_numpy(dtype=np.float32)

        ds = xr.Dataset(
            data_vars={
                "GI_count": self._grid_da(GI_count, "GI_count", dtype=np.int32),
                "GI_obs_count": self._grid_da(GI_obs_count, "GI_obs_count", dtype=np.int32),
                "GI_unique_uid": self._grid_da(GI_unique_uid, "GI_unique_uid", dtype=np.int32),
                "GI_area_sum_m2": self._grid_da(GI_area_sum_m2, "GI_area_sum_m2", dtype=np.float32),
                "GI_area_max_m2": self._grid_da(GI_area_max_m2, "GI_area_max_m2", dtype=np.float32),
                "GI_perim_sum_m": self._grid_da(GI_perim_sum_m, "GI_perim_sum_m", dtype=np.float32),
                "GI_bed_depth_mean": self._grid_da(GI_bed_depth_mean, "GI_bed_depth_mean", dtype=np.float32),
                "GI_map_dist_mean_m": self._grid_da(GI_map_dist_mean_m, "GI_map_dist_mean_m", dtype=np.float32),
            }
        )

        ds["GI_occupied"] = (ds["GI_count"] >= 1).astype(np.int8)
        ds["GI_occupied"].attrs["long_name"] = "binary occupied GI cell"

        self.GI_cell_stats = ds
        self.GI_cnts = ds["GI_count"].rename("GI_counts")  # backward compatibility
        self.logger.info(
            f"Aggregated grounded icebergs to grid: "
            f"occupied_cells={int((ds['GI_count'].values > 0).sum())}"
        )
        return self.GI_cell_stats

    ######################################################################
    # old public entrypoint, now upgraded
    ######################################################################
    def process_raw_grounded_iceberg_database(
        self,
        GI_P_raw=None,
        GI_layer="grounded_icebergs",
        csv_lon_col=None,
        csv_lat_col=None,
        dedup_uid=True,
        normalise_lon_to="0-360",
        proj_crs="EPSG:3031",
        max_dist_km=None,
        use_attr_area=True,
    ):
        """
        Load, map, and aggregate grounded iceberg source data.

        Backward-compatible public entry point.
        """
        self.load_grounded_iceberg_source(
            GI_P_raw=GI_P_raw,
            GI_layer=GI_layer,
            csv_lon_col=csv_lon_col,
            csv_lat_col=csv_lat_col,
            dedup_uid=dedup_uid,
            normalise_lon_to=normalise_lon_to,
            use_attr_area=use_attr_area,
        )
        self.assign_grounded_icebergs_to_tgrid(
            proj_crs=proj_crs,
            max_dist_km=max_dist_km,
        )
        self.aggregate_grounded_icebergs_on_tgrid()
        return self.GI_cell_stats

    ######################################################################
    # thinning
    ######################################################################
    def _is_isolated(self, y, x, mask):
        ny, nx = mask.shape
        adjacent_coords = [
            (y-1, x-1), (y-1, x), (y-1, x+1),
            (y,   x-1),           (y,   x+1),
            (y+1, x-1), (y+1, x), (y+1, x+1),
        ]
        for adj_y, adj_x in adjacent_coords:
            if 0 <= adj_y < ny and 0 <= adj_x < nx:
                if mask[adj_y, adj_x]:
                    return False
        return True

    def _thin_mask(
        self,
        mask,
        counts,
        p_min=0.2,
        p_max=0.8,
        scaling_exponent=2.0,
        random_state=0,
        protect_isolated=True,
    ):
        """
        Legacy probabilistic thinning, now reproducible via random_state.
        """
        rng = np.random.default_rng(random_state)

        ny, nx = mask.shape
        thinned_mask = np.zeros_like(mask, dtype=bool)

        count_min = np.nanmin(counts)
        count_max = np.nanmax(counts)

        if count_max > count_min:
            norm_counts = (counts - count_min) / (count_max - count_min)
            norm_counts = norm_counts ** scaling_exponent
        else:
            norm_counts = np.ones_like(counts, dtype=float)

        retain_probs = p_min + norm_counts * (p_max - p_min)

        retained_count = 0
        total_count = int(np.count_nonzero(mask))

        for y in range(ny):
            for x in range(nx):
                if not mask[y, x]:
                    continue

                keep = False
                if protect_isolated and self._is_isolated(y, x, mask):
                    keep = True
                elif rng.random() < retain_probs[y, x]:
                    keep = True

                if keep:
                    thinned_mask[y, x] = True
                    retained_count += 1

        GI_thinning_factor = 1.0 - (retained_count / total_count) if total_count > 0 else 0.0
        return thinned_mask, GI_thinning_factor

    def _ranked_thin_mask(
        self,
        occupied,
        score,
        keep_fraction=0.5,
        protect_isolated=True,
    ):
        """
        Deterministic thinning by retaining the highest-scoring occupied cells.
        """
        occupied = np.asarray(occupied, dtype=bool)
        score = np.asarray(score, dtype=float)

        total = int(np.count_nonzero(occupied))
        if total == 0:
            return np.zeros_like(occupied, dtype=bool), 0.0

        keep_fraction = float(np.clip(keep_fraction, 0.0, 1.0))
        n_keep_target = int(np.ceil(keep_fraction * total))

        keep = np.zeros_like(occupied, dtype=bool)

        if protect_isolated:
            ys, xs = np.where(occupied)
            for y, x in zip(ys, xs):
                if self._is_isolated(y, x, occupied):
                    keep[y, x] = True

        remaining_need = max(n_keep_target - int(keep.sum()), 0)
        eligible = occupied & (~keep)

        if remaining_need > 0 and np.any(eligible):
            flat = np.flatnonzero(eligible)
            flat_scores = score.ravel()[flat]
            flat_scores = np.where(np.isfinite(flat_scores), flat_scores, -np.inf)

            order = flat[np.argsort(flat_scores)[::-1]]
            chosen = order[:remaining_need]
            keep.ravel()[chosen] = True

        thinning_factor = 1.0 - (int(keep.sum()) / total)
        return keep, thinning_factor

    def thin_grounded_icebergs(
        self,
        strategy="legacy_count_prob",
        keep_fraction=None,
        rank_var="GI_area_sum_m2",
        p_min=0.1,
        p_max=0.9,
        scaling_exponent=1.5,
        random_state=0,
        protect_isolated=True,
    ):
        """
        Thin occupied GI cells using either the legacy stochastic method or a new
        deterministic ranked method.
        """
        if not hasattr(self, "GI_cell_stats"):
            raise ValueError("GI_cell_stats not found. Run process_raw_grounded_iceberg_database first.")

        counts = self.GI_cell_stats["GI_count"].values
        occupied = counts >= 1

        if strategy == "legacy_count_prob":
            thin_mask, thin_fact = self._thin_mask(
                mask=occupied,
                counts=counts,
                p_min=p_min,
                p_max=p_max,
                scaling_exponent=scaling_exponent,
                random_state=random_state,
                protect_isolated=protect_isolated,
            )
            score_used = counts.astype(float)

        elif strategy in {"ranked_area", "ranked_count", "ranked"}:
            if keep_fraction is None:
                keep_fraction = 0.5

            if strategy == "ranked_count":
                score_name = "GI_count"
            elif strategy == "ranked_area":
                score_name = "GI_area_sum_m2"
            else:
                score_name = rank_var

            if score_name not in self.GI_cell_stats:
                raise ValueError(
                    f"Requested rank_var='{score_name}' not found in GI_cell_stats. "
                    f"Available: {list(self.GI_cell_stats.data_vars)}"
                )

            score_used = self.GI_cell_stats[score_name].values
            thin_mask, thin_fact = self._ranked_thin_mask(
                occupied=occupied,
                score=score_used,
                keep_fraction=keep_fraction,
                protect_isolated=protect_isolated,
            )

        else:
            raise ValueError(
                f"Unknown thinning strategy '{strategy}'. "
                f"Supported: legacy_count_prob, ranked_area, ranked_count, ranked"
            )

        self.GI_keep_mask = thin_mask
        self.GI_thin_mask = thin_mask
        self.GI_thin_fact = thin_fact
        self.GI_thin_strategy = strategy
        self.GI_keep_score = self._grid_da(score_used, "GI_keep_score", dtype=np.float32)
        return thin_mask, thin_fact

    ######################################################################
    # modified landmask
    ######################################################################
    def modify_landmask_with_grounded_icebergs(
        self,
        GI_P_raw=None,
        GI_layer="grounded_icebergs",
        csv_lon_col=None,
        csv_lat_col=None,
        dedup_uid=True,
        normalise_lon_to="0-360",
        proj_crs="EPSG:3031",
        max_dist_km=None,
        use_attr_area=True,
        strategy="legacy_count_prob",
        keep_fraction=None,
        rank_var="GI_area_sum_m2",
        p_min=0.1,
        p_max=0.9,
        scaling_exponent=1.5,
        random_state=0,
        protect_isolated=True,
    ):
        """
        Full GI-to-KMT workflow.

        This remains backward-compatible with the old usage, but now supports:
          - CSV or GPKG input
          - legacy stochastic thinning
          - deterministic ranked thinning
        """
        if (GI_P_raw is not None) or (not hasattr(self, "GI_cell_stats")):
            self.process_raw_grounded_iceberg_database(
                GI_P_raw=GI_P_raw,
                GI_layer=GI_layer,
                csv_lon_col=csv_lon_col,
                csv_lat_col=csv_lat_col,
                dedup_uid=dedup_uid,
                normalise_lon_to=normalise_lon_to,
                proj_crs=proj_crs,
                max_dist_km=max_dist_km,
                use_attr_area=use_attr_area,
            )

        thin_mask, thin_fact = self.thin_grounded_icebergs(
            strategy=strategy,
            keep_fraction=keep_fraction,
            rank_var=rank_var,
            p_min=p_min,
            p_max=p_max,
            scaling_exponent=scaling_exponent,
            random_state=random_state,
            protect_isolated=protect_isolated,
        )

        counts = self.GI_cell_stats["GI_count"].values
        GI_thinned = counts * thin_mask

        self.GI_thin_da = self._grid_da(
            GI_thinned,
            "GI_counts",
            dtype=np.int32,
            attrs={
                "thinning_strategy": strategy,
                "thinning_factor": float(thin_fact),
            },
        )

        kmt = xr.open_dataset(self.P_KMT_org).kmt.values.copy()
        thin_nj, thin_ni = np.where(thin_mask)
        kmt[thin_nj, thin_ni] = 0

        self.GI_KMT = self._grid_da(
            kmt,
            "kmt",
            dtype=kmt.dtype,
            attrs={
                "modified_with_grounded_icebergs": 1,
                "thinning_strategy": strategy,
                "thinning_factor": float(thin_fact),
            },
        )

        self.logger.info(
            f"Modified landmask with grounded icebergs: "
            f"strategy={strategy}, thinning_factor={thin_fact:.3f}"
        )
        return self.GI_KMT

    def compute_grounded_iceberg_cell_count(self, hemisphere_only=True, return_area=False):
        """
        Count modified landmask cells introduced by grounded icebergs.

        A grounded-iceberg cell is defined here as:
            (kmt_org > 0) AND (kmt_mod == 0)

        Parameters
        ----------
        hemisphere_only : bool, default True
            If True, apply the configured hemisphere slice before counting.
        return_area : bool, default False
            If True, also return the summed T-cell area (m^2) of those cells.

        Returns
        -------
        int
            Number of grounded-iceberg-modified T-grid cells.

        or

        tuple
            (n_cells, area_m2) if return_area=True
        """
        self.load_cice_grid(slice_hem=hemisphere_only, build_faces=False)
        if not hasattr(self, "kmt_org"):
            raise ValueError("Original landmask (kmt_org) not loaded.")
        if not hasattr(self, "kmt_mod"):
            raise ValueError("Modified landmask (kmt_mod) not loaded. Set use_gi=True and/or define P_KMT_mod.")
        kmt_org = self.kmt_org["kmt_org"].values
        kmt_mod = self.kmt_mod["kmt_mod"].values
        gi_mask = (kmt_org > 0) & (kmt_mod == 0)
        n_cells = int(np.count_nonzero(gi_mask))
        if not return_area:
            return n_cells
        area_m2 = float(np.nansum(self.G_t["area"].values[gi_mask]))
        return n_cells, area_m2

    def write_modified_landmask_and_counts_to_disk(self,
                                                   GI_thin_fact   : float|int = None,
                                                   GI_version_str : str       = None,
                                                   iteration      : int       = 0,
                                                   overwrite      : bool      = False):
        """
        Save the thinned grounded iceberg counts and modified KMT landmask to NetCDF files.
        File paths are generated from the config template using GI thinning factor and version.
        """
        GI_thin_fact     = GI_thin_fact   if GI_thin_fact   is not None else self.GI_thin_fact
        GI_version_str   = GI_version_str if GI_version_str is not None else self.GI_version_str
        GI_thin_fact_str = f"{GI_thin_fact:.2f}".replace('.', 'p')
        F_thin           = self.GI_dict['GI_thin_fmt'].format(GI_thin=GI_thin_fact_str, version=GI_version_str, iteration=iteration)
        F_KMT_mod        = self.GI_dict['KMT_mod_fmt'].format(GI_thin=GI_thin_fact_str, version=GI_version_str, iteration=iteration)
        self.GI_P_counts = os.path.join(self.GI_dict['D_GI_thin'], F_thin)
        self.P_KMT_mod   = os.path.join(self.GI_dict['D_GI_thin'], F_KMT_mod)
        print(f"modified KMT     : {self.GI_P_counts}")
        print(f"thinned GI count : {self.P_KMT_mod}")
        if os.path.exists(self.GI_P_counts) or os.path.exists(self.P_KMT_mod):
            if not overwrite:
                print(f"GI counts file and/or KMT file already exists:\n  - {self.GI_P_counts}\n  - {self.P_KMT_mod}\n"
                      "These files may have been used in model simulations. Set `overwrite=True` in the constructor to allow overwriting." )
            else:
                date_tag   = datetime.now().strftime("%Y-%m-%d")
                base_gi    = os.path.splitext(self.GI_P_counts)[0]
                base_kmt   = os.path.splitext(self.P_KMT_mod)[0]
                backup_gi  = f"{base_gi}_{date_tag}.nc"
                backup_kmt = f"{base_kmt}_{date_tag}.nc"
                shutil.copy2(self.GI_P_counts, backup_gi)
                shutil.copy2(self.P_KMT_mod, backup_kmt)
                print("*** WARNING OVER-WRITING EXISTING FILES ***")
                print(f"Backups created:\n - {backup_gi}\n - {backup_kmt}")
                print("Saving thinned GI count and modified KMT files")
                self.GI_thin_da.to_netcdf(self.GI_P_counts, mode='w')
                self.GI_KMT.to_netcdf(self.P_KMT_mod, mode='w')
        else:
            print("Modified KMT and GI-count-file do not exist. Writing out the above files")
            self.GI_thin_da.to_netcdf(self.GI_P_counts, mode='w')
            self.GI_KMT.to_netcdf(self.P_KMT_mod, mode='w')

# import os, shutil
# import xarray      as xr
# import numpy       as np
# import pandas      as pd
# import geopandas   as gpd
# from scipy.spatial import cKDTree
# from pathlib       import Path
# from datetime      import datetime
# from scipy.spatial import cKDTree

# class SeaIceIcebergs:

#     def __init__(self, **kwargs):
#         return

#     def load_GI_lon_lats(self):
#         self.load_cice_grid(slice_hem=True)
#         lon  = self.G_t["lon"].reset_coords(drop=True)
#         lat  = self.G_t["lat"].reset_coords(drop=True)
#         mask = self.G_GI["mask"].reset_coords(drop=True).astype(bool)
#         return {"lon": lon.data[mask.data].ravel(),
#                 "lat": lat.data[mask.data].ravel()}

#     def _region_lon_mask(self, lon2d, lon_min, lon_max):
#         """
#         Build a longitude mask that works for either:
#         - non-wrapping bounds (lon_min <= lon_max)
#         - seam-crossing bounds (lon_min > lon_max) e.g. 330..10 on a 0..360 grid
#         """
#         if lon_min <= lon_max:
#             return (lon2d >= lon_min) & (lon2d <= lon_max)
#         else: # crosses seam
#             return (lon2d >= lon_min) | (lon2d <= lon_max)

#     def _lonlat_to_unit_xyz(lon, lat):
#         """
#         Convert lon/lat in degrees to 3D unit-sphere Cartesian coordinates.
#         This makes nearest-neighbour search more robust near the pole/dateline.
#         """
#         lon     = np.asarray(lon, dtype=float)
#         lat     = np.asarray(lat, dtype=float)
#         lon_rad = np.deg2rad(((lon + 180.0) % 360.0) - 180.0)
#         lat_rad = np.deg2rad(lat)
#         clat    = np.cos(lat_rad)
#         x       = clat * np.cos(lon_rad)
#         y       = clat * np.sin(lon_rad)
#         z       = np.sin(lat_rad)
#         return np.column_stack((x, y, z))

#     def _first_existing_column(columns, candidates):
#         for c in candidates:
#             if c in columns:
#                 return c
#         return None

#     def _warn(logger, msg):
#         if logger is not None:
#             logger.warning(msg)
#         else:
#             print(f"WARNING: {msg}")

#     def compute_grounded_iceberg_area(self, region=None, scale=1e6):
#         """
#         Returns:
#         - dict of km^2 by region if region is provided
#         - total GI area (m^2) if region is None (or km^2 if you divide by scale)
#         """
#         self.load_cice_grid(slice_hem=True)
#         # IMPORTANT: use FULL-grid lon/lat for region tests (no NaNs)
#         lon2d = self.G_t["lon"].values
#         lat2d = self.G_t["lat"].values
#         # GI area field already has NaNs outside GI in your loader
#         area2d = self.G_GI["area"].values
#         # IMPORTANT: convert mask back to boolean for indexing
#         gi_mask = self.G_GI["mask"].values.astype(bool)
#         # Logging (nan-safe)
#         self.logger.info(f"GI-area_calc: grid lon min/max: {np.nanmin(lon2d):.2f} / {np.nanmax(lon2d):.2f}")
#         self.logger.info(f"GI-area_calc: grid lat min/max: {np.nanmin(lat2d):.2f} / {np.nanmax(lat2d):.2f}")
#         if region is None:
#             # total GI area (m^2)
#             return np.nansum(area2d[gi_mask])
#         # Determine grid lon convention
#         grid_is_0360 = (np.nanmin(lon2d) >= 0.0) and (np.nanmax(lon2d) > 180.0)
#         area_dict = {}
#         for reg_name, reg_cfg in region.items():
#             lon_min, lon_max, lat_min, lat_max = reg_cfg["plot_region"]
#             # Convert region bounds if grid is 0..360 but bounds are negative / -180..180 style
#             if grid_is_0360:
#                 lon_min_c = self.normalise_longitudes(lon_min)
#                 lon_max_c = self.normalise_longitudes(lon_max)
#             else:
#                 lon_min_c, lon_max_c = lon_min, lon_max
#             lon_mask = self._region_lon_mask(lon2d, lon_min_c, lon_max_c)
#             lat_mask = (lat2d >= lat_min) & (lat2d <= lat_max)
#             region_mask = lon_mask & lat_mask
#             # Sum GI areas within region
#             area_km2 = np.nansum(area2d[gi_mask & region_mask]) / scale
#             area_dict[reg_name] = area_km2
#             self.logger.info(f"GI-area_calc: region {reg_name} GI area = {area_km2:0.2f} km^2")
#         return area_dict

#     def check_GI_coverage(self, da, varname="GI_counts", lon_name="lon", lat_name="lat"):
#         """
#         Assess grounded iceberg (GI) or landmask coverage across Antarctica.

#         Parameters
#         ----------
#         da : xr.DataArray
#             The GI or KMT data.
#         varname : str
#             Determines masking logic: 'GI_counts' uses >0, anything else uses >1.
#         lon_name : str, default 'lon'
#             Name of the longitude coordinate.
#         lat_name : str, default 'lat'
#             Name of the latitude coordinate.

#         Returns
#         -------
#         report : dict
#             Dictionary of summary statistics and Antarctic sector coverage.
#         """
#         report = {}
#         lon = da[lon_name]
#         lat = da[lat_name]
#         lon_vals = lon.values
#         lat_vals = lat.values
#         # Handle longitude wrapping
#         if lon_vals.max() <= 180:
#             lon_wrapped = np.mod(lon_vals + 360, 360)
#             if lon_vals.ndim == 1:
#                 da = da.assign_coords({lon_name: (da[lon_name].dims, lon_wrapped)})
#                 da = da.sortby(lon_name)
#             elif lon_vals.ndim == 2:
#                 da = da.assign_coords({lon_name: (da[lon_name].dims, lon_wrapped)})
#                 self.logger.warning("Skipping sortby(): longitude is 2D.")
#             else:
#                 raise ValueError("Longitude data must be 1D or 2D.")
#         # Construct binary coverage mask
#         if varname == "GI_counts":
#             mask = da.values > 0
#         else:
#             mask = da.values > 1
#         # Basic spatial stats
#         report['lon_range'] = (float(np.nanmin(lon_vals)), float(np.nanmax(lon_vals)))
#         report['lat_range'] = (float(np.nanmin(lat_vals)), float(np.nanmax(lat_vals)))
#         report['n_valid_cells'] = int(np.count_nonzero(mask))
#         report['n_total_cells'] = int(np.prod(mask.shape))
#         report['coverage_fraction'] = report['n_valid_cells'] / report['n_total_cells']
#         report['nonfinite_count'] = int(np.sum(~np.isfinite(da.values)))
#         # Value stats
#         report['mean_value'] = float(np.nanmean(da.values)) if np.any(mask) else float('nan')
#         report['max_value'] = float(np.nanmax(da.values)) if np.any(mask) else float('nan')
#         report['min_nonzero'] = float(np.nanmin(da.values[da.values > 0])) if np.any(da.values > 0) else float('nan')
#         # Sector-based coverage check
#         if lon_vals.ndim == 2 and lat_vals.ndim == 2:
#             try:
#                 sector_edges = np.linspace(0, 360, 9)  # 8 sectors
#                 presence_by_sector = []
#                 for i in range(len(sector_edges) - 1):
#                     lon_min, lon_max = sector_edges[i], sector_edges[i + 1]
#                     sector_mask = (
#                         (lon_vals >= lon_min) & (lon_vals < lon_max) &
#                         (lat_vals < -60)  # Limit to Antarctic region
#                     )
#                     if np.any(sector_mask):
#                         sector_data = da.values[sector_mask]
#                         presence_by_sector.append(np.any(sector_data > 0))
#                     else:
#                         presence_by_sector.append(False)
#                 report['sector_gaps'] = [i for i, present in enumerate(presence_by_sector) if not present]
#                 report['sector_gap_count'] = len(report['sector_gaps'])
#             except Exception as e:
#                 report['sector_check_failed'] = str(e)
#         return report

#     def load_existing_thinned(self):
#         """
#         Load previously saved thinned grounded iceberg data and modified KMT file.

#         Sets the attributes `self.GI_thin_da`, `self.GI_thin_mask`, and `self.GI_KMT`
#         if their respective NetCDF files exist.
#         """
#         if os.path.exists(self.P_GI_thin):
#             self.GI_thin_da = xr.open_dataarray(self.P_GI_thin)
#             self.logger.info(f"GROUNDED ICEBERG COUNTS PER GRID CELL REPORT {self.P_GI_thin}")
#             self.logger.info(self.check_GI_coverage(self.GI_thin_da))
#             self.GI_thin_mask = self.GI_thin_da.values > 0
#         if os.path.exists(self.P_KMT_mod):
#             self.GI_KMT = xr.open_dataarray(self.P_KMT_mod)
#             self.logger.info(f"MODIFIED LANDMASK REPORT {self.P_KMT_mod}")
#             self.logger.info(self.check_GI_coverage(self.GI_KMT, varname="kmt"))

#     def _read_grounded_iceberg_points(self,
#         P_raw,
#         layer=None,
#         x_col=None,
#         y_col=None,
#         assume_gpkg_crs="EPSG:3031",
#         logger=None,
#     ):
#         """
#         Read grounded iceberg point data from CSV or GPKG and return a pandas DataFrame
#         with lon/lat stored in columns 'X' and 'Y'.

#         Parameters
#         ----------
#         P_raw : str or Path
#             Input CSV or GPKG path.
#         layer : str, optional
#             GeoPackage layer name. If None, geopandas will use its default behaviour.
#         x_col, y_col : str, optional
#             Explicit coordinate column names for CSV input.
#         assume_gpkg_crs : str, default 'EPSG:4326'
#             CRS to assign if a GPKG has no CRS metadata.
#         logger : logging.Logger, optional

#         Returns
#         -------
#         pandas.DataFrame
#             DataFrame with at least columns ['X', 'Y'].
#         """
#         P_raw = Path(P_raw)
#         suffix = P_raw.suffix.lower()
#         if suffix == ".csv":
#             df = pd.read_csv(P_raw)
#             if x_col is None:
#                 x_col = _first_existing_column(
#                     df.columns,
#                     ["X", "x", "lon", "Lon", "LON", "longitude", "LONGITUDE"],
#                 )
#             if y_col is None:
#                 y_col = _first_existing_column(
#                     df.columns,
#                     ["Y", "y", "lat", "Lat", "LAT", "latitude", "LATITUDE"],
#                 )

#             if x_col is None or y_col is None:
#                 raise ValueError(
#                     f"Could not identify longitude/latitude columns in CSV: {P_raw}\n"
#                     f"Available columns: {list(df.columns)}\n"
#                     f"Pass x_col=... and y_col=... explicitly."
#                 )

#             out = df.copy()
#             out["X"] = pd.to_numeric(out[x_col], errors="coerce")
#             out["Y"] = pd.to_numeric(out[y_col], errors="coerce")
#             out = out.dropna(subset=["X", "Y"]).copy()
#             return out

#         elif suffix == ".gpkg":
#             try:
#                 import geopandas as gpd
#             except ImportError as e:
#                 raise ImportError(
#                     "Reading GeoPackage files requires geopandas. "
#                     "Install geopandas (plus pyogrio/fiona) in your environment."
#                 ) from e

#             gdf = gpd.read_file(P_raw, layer=layer)

#             if gdf.empty:
#                 raise ValueError(f"No features found in GeoPackage: {P_raw}")

#             if gdf.geometry is None:
#                 raise ValueError(f"No geometry column found in GeoPackage: {P_raw}")

#             gdf = gdf[gdf.geometry.notna()].copy()

#             # Expand MultiPoint to individual Point features where relevant
#             try:
#                 gdf = gdf.explode(index_parts=False, ignore_index=True)
#             except TypeError:
#                 # compatibility with older geopandas
#                 gdf = gdf.explode()

#             if gdf.crs is None:
#                 _warn(
#                     logger,
#                     f"{P_raw} has no CRS metadata; assuming {assume_gpkg_crs}.",
#                 )
#                 gdf = gdf.set_crs(assume_gpkg_crs)

#             # Convert to lon/lat for consistency with model grid lookup
#             gdf = gdf.to_crs("EPSG:4326")

#             non_point = ~gdf.geometry.geom_type.isin(["Point"])
#             if np.any(non_point):
#                 _warn(
#                     logger,
#                     f"{P_raw} contains non-Point geometries; using representative points.",
#                 )
#                 geom = gdf.geometry.copy()
#                 geom.loc[non_point] = gdf.loc[non_point].representative_point()
#                 gdf = gdf.set_geometry(geom)

#             out = pd.DataFrame(gdf.drop(columns=[gdf.geometry.name]))
#             out["X"] = gdf.geometry.x.to_numpy()
#             out["Y"] = gdf.geometry.y.to_numpy()
#             out = out.dropna(subset=["X", "Y"]).copy()
#             return out
#         else:
#             raise ValueError(f"Unsupported grounded iceberg file type: {P_raw.suffix}. "
#                              f"Supported types are .csv and .gpkg")

#     def process_raw_grounded_iceberg_database(
#         self,
#         GI_P_raw=None,
#         GI_layer=None,
#         x_col=None,
#         y_col=None,
#         assume_gpkg_crs="EPSG:4326",
#     ):
#         """
#         Load and spatially map raw grounded iceberg locations from CSV or GPKG to
#         the model grid.

#         Uses a KD-tree nearest-neighbour search on 3D unit-sphere coordinates to
#         assign iceberg coordinates to grid cells.

#         Sets
#         ----
#         self.GI_clean : pandas.DataFrame
#             Original iceberg table with attached model-grid indices (nj_i, ni_i)
#             and standardised lon/lat columns (X, Y).
#         self.GI_cnts : xarray.DataArray
#             2D grounded iceberg counts on the model grid.

#         Parameters
#         ----------
#         GI_P_raw : str or Path, optional
#             Path to raw grounded iceberg database (.csv or .gpkg). If omitted,
#             uses self.GI_dict['P_raw'].
#         GI_layer : str, optional
#             GeoPackage layer name, only used for .gpkg input.
#         x_col, y_col : str, optional
#             Explicit lon/lat column names for CSV input.
#         assume_gpkg_crs : str, default 'EPSG:4326'
#             CRS to assume if a GPKG is missing CRS metadata.
#         """
#         self.load_cice_grid()
#         GI_P_raw = GI_P_raw if GI_P_raw is not None else self.GI_dict["P_KJ"]
#         logger = getattr(self, "logger", None)
#         GI_clean = self._read_grounded_iceberg_points(
#             P_raw=GI_P_raw,
#             layer=GI_layer,
#             x_col=x_col,
#             y_col=y_col,
#             assume_gpkg_crs=assume_gpkg_crs,
#             logger=logger,
#         )
#         # Model grid lon/lat
#         G_lon = np.asarray(self.G_t["lon"].values)
#         G_lat = np.asarray(self.G_t["lat"].values)
#         # KD-tree in 3D Cartesian coordinates on the unit sphere
#         G_pts = _lonlat_to_unit_xyz(G_lon.ravel(), G_lat.ravel())
#         tree = cKDTree(G_pts)
#         GI_locs = _lonlat_to_unit_xyz(
#             GI_clean["X"].to_numpy(),
#             GI_clean["Y"].to_numpy(),
#         )
#         _, nrst_i = tree.query(GI_locs)
#         nj, ni = G_lon.shape
#         nj_i, ni_i = np.unravel_index(nrst_i, (nj, ni))
#         GI_clean = GI_clean.copy()
#         GI_clean["nj_i"] = nj_i.astype(np.int32)
#         GI_clean["ni_i"] = ni_i.astype(np.int32)
#         self.GI_clean = GI_clean
#         GI_cnts = np.zeros((nj, ni), dtype=np.int32)
#         GI_cnts_df = (
#             self.GI_clean
#             .groupby(["nj_i", "ni_i"], sort=False)
#             .size()
#             .reset_index(name="count")
#         )

#         GI_cnts[
#             GI_cnts_df["nj_i"].to_numpy(),
#             GI_cnts_df["ni_i"].to_numpy()
#         ] = GI_cnts_df["count"].to_numpy()

#         GI_cnts_da = xr.DataArray(
#             GI_cnts,
#             dims=self.spatial_dims,
#             name="GI_counts",
#         )

#         GI_cnts_da.attrs["source_file"] = str(GI_P_raw)
#         if GI_layer is not None:
#             GI_cnts_da.attrs["source_layer"] = GI_layer

#         self.GI_cnts = GI_cnts_da

#     def modify_landmask_with_grounded_icebergs(self, p_min=0.1, p_max=0.9, scaling_exponent=1.5):
#         """
#         Apply probabilistic thinning to the grounded iceberg mask based on iceberg counts.

#         Parameters:
#         -----------
#         p_min : float
#             Minimum probability of retaining a cell with low counts.
#         p_max : float
#             Maximum probability of retaining a cell with high counts.
#         scaling_exponent : float
#             Controls steepness of probability scaling with count.

#         Sets:
#         ------
#         self.GI_thin_da : xarray.DataArray
#         self.GI_thin_mask : ndarray (boolean)
#         self.GI_thin_fact : float (fraction thinned out)
#         """
#         kmt                   = xr.open_dataset(self.P_KMT_org).kmt.values
#         mask                  = self.GI_cnts.values >= 1
#         counts                = self.GI_cnts.values
#         thin_mask, thin_fact  = self._thin_mask(mask, counts, p_min, p_max, scaling_exponent)
#         GI_thinned            = counts * thin_mask
#         self.GI_thin_da       = xr.DataArray(data   = GI_thinned,
#                                              dims   = self.spatial_dims,
#                                              coords = {'lat': (self.spatial_dims, self.G_t['lat'].values),
#                                                        'lon': (self.spatial_dims, self.G_t['lon'].values)},
#                                              name   = 'GI_counts')
#         self.GI_thin_mask     = thin_mask
#         self.GI_thin_fact     = thin_fact
#         thin_nj, thin_ni      = np.where(thin_mask)
#         kmt[thin_nj, thin_ni] = 0
#         self.GI_KMT           = xr.DataArray(data   = kmt,
#                                              dims   = self.spatial_dims,
#                                              coords = {'lat': (self.spatial_dims, self.G_t['lat'].values),
#                                                        'lon': (self.spatial_dims, self.G_t['lon'].values)},
#                                              name   = 'kmt' )

#     def write_modified_landmask_and_counts_to_disk(self,
#                                                    GI_thin_fact   : float|int = None,
#                                                    GI_version_str : str       = None,
#                                                    overwrite      : bool      = None):
#         """
#         Save the thinned grounded iceberg counts and modified KMT landmask to NetCDF files.
#         File paths are generated from the config template using GI thinning factor and version.
#         """
#         GI_thin_fact     = GI_thin_fact   if GI_thin_fact   is not None else self.GI_thin_fact
#         GI_version_str   = GI_version_str if GI_version_str is not None else self.GI_version_str
#         overwrite        = overwrite      if overwrite      is not None else self.overwrite
#         GI_thin_fact_str = f"{GI_thin_fact:.2f}".replace('.', 'p')
#         F_thin           = self.GI_dict['GI_thin_fmt'].format(GI_thin_fact=GI_thin_fact_str, version=GI_version_str)
#         F_KMT_mod        = self.GI_dict['KMT_mod_fmt'].format(GI_thin_fact=GI_thin_fact_str, version=GI_version_str)
#         self.GI_P_counts = os.path.join(self.GI_dict['D_GI_thin'], F_thin)
#         self.P_KMT_mod   = os.path.join(self.GI_dict['D_GI_thin'], F_KMT_mod)
#         print(f"modified KMT     : {self.GI_P_counts}")
#         print(f"thinned GI count : {self.P_KMT_mod}")
#         if os.path.exists(self.GI_P_counts) or os.path.exists(self.P_KMT_mod):
#             if not overwrite:
#                 print(f"GI counts file and/or KMT file already exists:\n  - {self.GI_P_counts}\n  - {self.P_KMT_mod}\n"
#                       "These files may have been used in model simulations. Set `overwrite=True` in the constructor to allow overwriting." )
#             else:
#                 date_tag   = datetime.now().strftime("%Y-%m-%d")
#                 base_gi    = os.path.splitext(self.GI_P_counts)[0]
#                 base_kmt   = os.path.splitext(self.P_KMT_mod)[0]
#                 backup_gi  = f"{base_gi}_{date_tag}.nc"
#                 backup_kmt = f"{base_kmt}_{date_tag}.nc"
#                 shutil.copy2(self.GI_P_counts, backup_gi)
#                 shutil.copy2(self.P_KMT_mod, backup_kmt)
#                 print("*** WARNING OVER-WRITING EXISTING FILES ***")
#                 print(f"Backups created:\n - {backup_gi}\n - {backup_kmt}")
#                 print("Saving thinned GI count and modified KMT files")
#                 self.GI_thin_da.to_netcdf(self.GI_P_counts, mode='w')
#                 self.GI_KMT.to_netcdf(self.P_KMT_mod, mode='w')
#         else:
#             print("Modified KMT and GI-count-file do not exist. Writing out the above files")
#             self.GI_thin_da.to_netcdf(self.GI_P_counts, mode='w')
#             self.GI_KMT.to_netcdf(self.P_KMT_mod, mode='w')

#     def _is_isolated(self, y, x, mask):
#         """
#         Determine if a given grid cell is isolated (no grounded icebergs in 8-connected neighbors).

#         Parameters:
#         -----------
#         y, x : int
#             Grid coordinates to test.
#         mask : ndarray
#             Boolean mask where True = grounded iceberg.

#         Returns:
#         --------
#         bool
#             True if isolated, False otherwise.
#         """
#         ny, nx = mask.shape
#         adjacent_coords = [(y-1, x-1), (y-1, x)  , (y-1, x+1),
#                            (y,   x-1), (y,   x+1),
#                            (y+1, x-1), (y+1, x)  , (y+1, x+1) ]
#         for adj_y, adj_x in adjacent_coords:
#             if 0 <= adj_y < ny and 0 <= adj_x < nx:
#                 if mask[adj_y, adj_x]:
#                     return False
#         return True

#     def _thin_mask(self, mask, counts, p_min=0.2, p_max=0.8, scaling_exponent=2.0):
#         """
#         Perform probabilistic thinning of grounded iceberg cells based on their counts.

#         Parameters:
#         -----------
#         mask : ndarray
#             Boolean array of original grounded iceberg presence.
#         counts : ndarray
#             Corresponding iceberg counts.
#         p_min : float
#             Minimum probability of retaining a cell.
#         p_max : float
#             Maximum probability of retaining a cell.
#         scaling_exponent : float
#             Controls how sharply retention probability scales with count.

#         Returns:
#         --------
#         thinned_mask : ndarray (bool)
#             Mask after thinning.
#         GI_thinning_factor : float
#             Fraction of cells that were removed.
#         """
#         ny, nx               = mask.shape
#         thinned_mask         = np.zeros_like(mask, dtype=bool)
#         count_min, count_max = np.nanmin(counts), np.nanmax(counts)
#         if count_max > count_min:
#             norm_counts = (counts - count_min) / (count_max - count_min)
#             norm_counts = norm_counts ** scaling_exponent
#         else:
#             norm_counts = np.ones_like(counts)
#         retain_probs   = p_min + norm_counts * (p_max - p_min)
#         retained_count = 0
#         total_count    = np.count_nonzero(mask)
#         for y in range(ny):
#             for x in range(nx):
#                 if mask[y, x]:
#                     retain_prob = retain_probs[y, x]
#                     if self._is_isolated(y, x, mask) or np.random.rand() < retain_prob:
#                         thinned_mask[y, x] =  True
#                         retained_count     += 1
#         GI_thinning_factor = 1 - (retained_count / total_count) if total_count > 0 else 0
#         return thinned_mask, GI_thinning_factor

#     def regrid_gi_to_fip_bool(self, gi_mask_da: xr.DataArray, gi_lon2d: xr.DataArray, gi_lat2d: xr.DataArray, fip_lon2d: xr.DataArray, fip_lat2d: xr.DataArray,
#                               dilate_pixels: int = 1):   # set 0 to disable adjacency expansion
#         """
#         Regrids a boolean GI mask from the model grid onto the FIP grid via nearest-neighbour
#         on great-circle distance (xesmf 'nearest_s2d'). Returns a boolean mask on the FIP grid.
#         Optionally dilates by N pixels on the target grid to include coast-adjacent cells.
#         """
#         from scipy.ndimage import binary_dilation
#         import xesmf as xe
#         # xESMF expects float fields; cast True->1, False->0
#         gi_float = gi_mask_da.astype(np.float32)
#         # Build source/target 'grid' dictionaries from 2D lon/lat
#         src = {'lon': gi_lon2d.values, 'lat': gi_lat2d.values}
#         dst = {'lon': fip_lon2d.values, 'lat': fip_lat2d.values}
#         # Regridder (nearest source-to-destination). reuse_weights=True if you call repeatedly.
#         regridder       = xe.Regridder(src, dst, method='nearest_s2d', locstream_in=False, locstream_out=False)
#         gi_on_fip_float = regridder(gi_float.values)  # ndarray on target shape
#         gi_on_fip       = xr.DataArray(gi_on_fip_float > 0.5, coords=fip_lon2d.coords, dims=fip_lon2d.dims)
#         # Optional: 1-pixel dilation to catch coast-adjacent cells next to GI
#         if dilate_pixels and dilate_pixels > 0:
#             # ndimage expects numpy array; preserve NaNs as False for mask ops
#             mask_np   = np.asarray(gi_on_fip.fillna(False).values, dtype=bool)
#             structure = np.ones((1 + 2*dilate_pixels, 1 + 2*dilate_pixels), dtype=bool)
#             mask_np   = binary_dilation(mask_np, structure=structure)
#             gi_on_fip = xr.DataArray(mask_np, coords=fip_lon2d.coords, dims=fip_lon2d.dims)
#         gi_on_fip.name = "GI_mask_on_FIP"
#         gi_on_fip.attrs.update(dict(long_name="Grounded iceberg mask (regridded to FIP)", comment="nearest_s2d + optional dilation"))
#         return gi_on_fip
