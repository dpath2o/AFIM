# AFIM Sea Ice Toolbox: `sea_ice_config.json` reference

This page documents the **global configuration JSON** consumed by `SeaIceToolbox` and its mixin sub-classes:

- `SeaIceClassification`
- `SeaIceMetrics`
- `SeaIcePlotter`
- `SeaIceIcebergs`
- `SeaIceObservations`
- `SeaIceACCESS`
- `SeaIceGridWork`
- `SeaIceRegridder`

`SeaIceToolbox` loads this JSON at construction time and uses it to:
1. Define **filesystem locations** (inputs/outputs, logs, grids, observational archives).
2. Define **canonical names** for key variables and dimensions (e.g., `aice`, `hi`, `time`, `nj`, `ni`).
3. Provide **standard thresholds** and analysis defaults (speed threshold, concentration threshold, persistence window, etc.).
4. Provide **plotting metadata** (labels, colormaps, panel defaults, sector extents).
5. Provide **dataset-specific settings** (AF2020 bins, NSIDC formats, ORAS/ERA5 access patterns, etc.).

> Important: this is the **global** config. During runtime, `SeaIceToolbox.parse_simulation_metadata()` may also write a **per-simulation** JSON (e.g., `ice_in_AFIM_subset_<sim>.json`) derived from `ice_diag.d`. That per-simulation metadata file is **separate** from `sea_ice_config.json` and is not described in full here.

---

## Conventions used in this document

- **Required** means the toolbox will raise a `KeyError` or otherwise fail without the field (based on current code paths).
- **Optional** means the toolbox has a fallback default or the field is only needed for certain workflows.
- Types are given in Python terms: `str`, `int`, `float`, `bool`, `list[...]`, `dict[...]`.

---

## Top-level (non-dictionary) keys

These keys control the default analysis window and common thresholds.

| Key | Type | Required | Meaning |
|---|---:|:---:|---|
| `dt0_str` | `str` | Optional | Default analysis start date (`YYYY-MM-DD`). |
| `dtN_str` | `str` | Optional | Default analysis end date (`YYYY-MM-DD`). |
| `leap_year` | `int` | Optional | A leap year used in certain climatology/time utilities. |
| `hemisphere` | `str` | Optional | Default hemisphere (`"south"` or `"north"`). |
| `iceh_freq` | `str` | Optional | Default history cadence (`"daily"`, `"monthly"`, etc.). |
| `cm2m_fact` | `float` | Optional | Unit conversion factor (centimetres → metres), used where needed. |
| `mean_period` | `int` | Optional | Rolling mean window (days) used in speed-based methods. |
| `bin_win_days` | `int` | Optional | Binary-days window length (days). |
| `bin_min_days` | `int` | Optional | Minimum number of “fast” days within `bin_win_days`. |
| `valid_BorC2T_types` | `list[str]` | Optional | Allowed B→T / regrid modes (e.g., `["B","Ta","Tb","Tc","Tx"]`). |
| `BorC2T_type` | `list[str]` or `str` | Optional | Default mode(s) used to build speed fields. |
| `valid_ice_types` | `list[str]` | Optional | Allowed classification types (e.g., `["FI","PI","SI"]`). |
| `ice_type` | `str` or `list[str]` | Optional | Default classification type(s). |
| `ice_speed_thresh_hi` | `float` | Optional | Default “strict” speed threshold (m/s), e.g. `5e-4`. |
| `ice_speed_thresh_lo` | `float` | Optional | Default “lenient” speed threshold (m/s). |
| `ice_conc_thresh` | `float` | Optional | Ice concentration threshold for masking (fraction, e.g. `0.15`). |
| `FIC_scale` | `float` | Optional | Scale factor used in fast-ice metrics (project-specific). |
| `FI_thick_max` | `float` | Optional | Upper cap for fast-ice thickness diagnostics (m). |
| `SIC_scale` | `float` | Optional | Scale factor used in sea-ice volume/area metrics (project-specific). |
| `metrics_name` | `str` | Optional | Tag used in metrics directory/file naming. |

---

## Dictionary keys (modules and their config blocks)

### 1) `FI_class_types` (classification naming)

**Purpose:** Maps classification methods to suffixes used in output naming.

**Schema:**
- **Required keys:** whatever methods you use in code (commonly `binary-days`, `rolling`).
- **Values:** short string tags used in constructed names.

Example:
```json
"FI_class_types": {
  "binary-days": "bin",
  "rolling": "roll"
}
```

---

### 2) `D_dict` (top-level directory registry)

**Purpose:** Central location for all major directory roots used throughout AFIM.

**Common requirements (practical minimum):**
- `AFIM_out` (model output root)
- `logs` (log root)
- `tmp` (temporary working directory)
- `graph` (figure output root)

**Schema (recommended):**
| Key | Type | Required | Meaning |
|---|---:|:---:|---|
| `AFIM_out` | `str` | **Yes** | Root directory containing simulations (`<AFIM_out>/<sim_name>/...`). |
| `AFIM_in` | `str` | Usually | Root directory for input archives (obs, forcing, etc.). |
| `logs` | `str` | **Yes** | Directory for toolbox log files. |
| `tmp` | `str` | **Yes** | Temporary workspace (Dask spill, intermediate files). |
| `graph` | `str` | **Yes** | Root directory for saved figures. |
| `grids` | `str` | Optional | Root for grid files if you keep them external to the repo. |

---

### 3) `CICE_dict` (CICE conventions, dimensions, and required variables)

**Purpose:** Defines how to interpret CICE grid/history data and how to chunk/select variables.

**Fields typically required by core workflows:**
- Grid files/paths: `P_G` (CICE grid NetCDF)
- Dimension names: `time_dim`, `x_dim`, `y_dim`, `spatial_dims`, `three_dims`
- Coordinate names: `lon_coord_name`, `lat_coord_name`, and/or `tcoord_names`, `bcoord_names`
- Variable name conventions: `SIC_name`, `SIT_name`, `SI_vel_names`, `cice_vars_reqd`
- Regridding weights path: `P_reG_u2t_weights` (if using xESMF regridding)
- Chunking presets: `FI_chunks` (if using Zarr + Dask patterns)

**Schema:**
| Key | Type | Required | Meaning |
|---|---:|:---:|---|
| `D_input` | `str` | Optional | Root for CICE input artifacts (optional if paths elsewhere). |
| `P_G` | `str` | **Yes** | Path to the CICE grid file (contains `tlat/tlon/ulat/ulon/angle/...`). |
| `P_B` | `str` | Optional | Path to bathymetry/topography file used by some pipelines. |
| `P_KMT` | `str` | Optional | Path to a landmask/bathymetry mask file (if used directly). |
| `P_reG_u2t_weights` | `str` | Required for xESMF | Filename for saved xESMF weights (U→T). |
| `G_res` | `float` | Optional | Grid resolution in degrees (metadata/plotting). |
| `SIC_name` | `str` | Optional | Standard sea-ice concentration variable name (e.g. `aice`). |
| `SIT_name` | `str` | Optional | Standard sea-ice thickness variable name (e.g. `hi`). |
| `SI_vel_names` | `list[str]` | Optional | Velocity component names (e.g. `["uvel","vvel"]`). |
| `area_name` | `str` | Optional | Cell area variable name on grid (e.g. `tarea`). |
| `time_dim` | `str` | **Yes** | Time dimension name (often `"time"`). |
| `x_dim`, `y_dim` | `str` | **Yes** | Horizontal dimension names (often `"ni"`, `"nj"`). |
| `x_dim_length`, `y_dim_length` | `int` | Optional | Expected full-grid sizes (used in padding/consistency checks). |
| `wrap_x` | `bool` | Optional | Whether the x-dimension is cyclic (global wrap). |
| `spatial_dims` | `list[str]` or `tuple[str,str]` | **Yes** | Horizontal dims, e.g. `["nj","ni"]`. |
| `three_dims` | `list[str]` | Optional | Spatiotemporal dims, e.g. `["time","nj","ni"]`. |
| `lon_coord_name`, `lat_coord_name` | `str` | Optional | Default coordinate names in CICE datasets (e.g. `TLON/TLAT`). |
| `tcoord_names` | `list[str]` | Optional | T-grid lon/lat coordinate names used in some routines. |
| `bcoord_names` | `list[str]` | Required for B→T averaging | B/U-grid lon/lat coordinate names (e.g. `["ULON","ULAT"]`). |
| `drop_coords` | `list[str]` | Optional | Coordinates to drop when cleaning datasets for merge/concat. |
| `drop_vars` | `list[str]` | Optional | Variables to drop when cleaning datasets. |
| `cice_vars_reqd` | `list[str]` | **Yes** | Minimal variable set required for AFIM workflows (masking/metrics). |
| `cice_vars_ext` | `list[str]` | Optional | Extra variables to include when `extra_cice_vars=True`. |
| `FI_chunks` | `dict[str,int]` | Optional | Recommended Dask chunking for fast-ice/classification products. |

---

### 4) `GI_dict` (grounded iceberg dataset + modified landmask formats)

**Purpose:** Supports workflows that modify the landmask to represent grounded icebergs (GI), and stores naming conventions for GI products.

**Schema:**
| Key | Type | Required | Meaning |
|---|---:|:---:|---|
| `D_GI` | `str` | Optional | Root directory holding grounded iceberg source products. |
| `D_GI_thin` | `str` | Often | Directory used to store thinned GI products and KMT masks. |
| `GI_thin_fmt` | `str` | Optional | Filename template for GI thinning products. |
| `KMT_org_fmt` | `str` | Optional | Filename for the “original” KMT file (no GI). |
| `KMT_mod_fmt` | `str` | Optional | Filename template for “modified” KMT including GI. |

---

### 5) `AF_FI_dict` (Fraser et al. 2020 fast-ice observational conventions)

**Purpose:** Configures access and interpretation of the AF2020 fast-ice dataset (binning, variables, file layout).

**Common fields:**
| Key | Type | Required | Meaning |
|---|---:|:---:|---|
| `D_AF2020` | `str` | Usually | Root directory for AF2020 dataset. |
| `F_AF2020_fmt` | `str` | Usually | Filename pattern for AF2020 files (often includes year). |
| `AF2020_var_name` | `str` | Usually | Variable name for fast-ice mask/field in AF2020. |
| `DOY_vals` | `list[int]` | **Yes** | Start day-of-year values defining AF2020 temporal bins. |

---

### 6) `NSIDC_dict`, `Sea_Ice_Obs_dict`, `AOM2_dict`, `MOM_dict`, `ORAS_dict`, `ERA5_dict`

These dictionaries serve as *registries* for dataset-specific paths, variable mappings, and conventions. Their exact fields are project-specific, but the intent is consistent:
- provide **directory roots** and **file naming patterns**
- define **variable names** and **dimension/coordinate names**
- define **time coverage** and dataset identifiers

---

### 7) `hemispheres_dict`, `Ant_8sectors`, `Ant_2sectors`, `specific_regions`

These control **spatial slicing** and **plotting extents/projections**, ensuring consistent regional figures across workflows.

---

### 8) `pygmt_dict`, `pygmt_FIA_dict`, `pygmt_FI_panel`, `plot_var_dict`

These control **plot styling** (labels, legend entries, axis presets, colormaps, and panel defaults).

---

## Practical validation checklist

Before running large workflows, validate that:

1. **Paths exist**:
   - `D_dict.AFIM_out`, `D_dict.logs`, `D_dict.tmp`, `D_dict.graph`
   - `CICE_dict.P_G`

2. **Dimension conventions are consistent**:
   - `CICE_dict.time_dim`, `CICE_dict.x_dim`, `CICE_dict.y_dim`
   - `CICE_dict.spatial_dims` matches your datasets (typically `("nj","ni")`)

3. **Required CICE variables are present**:
   - `CICE_dict.cice_vars_reqd`

4. **If using xESMF**:
   - `CICE_dict.P_reG_u2t_weights` points to a writable location on first run.
