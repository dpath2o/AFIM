# Lateral Drag Form Factor Creation

This workflow constructs **unitless lateral drag form factors** $(\mathrm{F2}_{x}, \mathrm{F2}_{y})$ on a **T-grid** (intended for use with [CICE](https://cice-consortium-cice.readthedocs.io/en/main/index.html) model) using:

1. A **high-resolution polygon coastline** (Antarctic land + ice-shelf surfaces) to compute a *coast-only* form factors $(\mathrm{F2cst}_{x}, \mathrm{F2cst}_{y})$; and
2. A **grounded-iceberg polygon dataset** to add sub-grid obstacle form drag, producing a *combined* (high-resolution coastline + grounded-iceberg locations).

The implementation follows the cell-based formulation of [**Liu et al. (2022)**](https://onlinelibrary.wiley.com/doi/abs/10.1029/2022JC018413) for coastal geometry, then appends grounded-iceberg contributions via a separate methodolgy.

---

## 1. Definitions and target quantities

For each T-grid-cell $(i,j)$, two **cell-based** form factors are defined:

$$
F_{2x}(i,j) \;=\; \frac{1}{\Delta x(i,j)} \sum_{n \in \mathcal{S}(i,j)} \left| \ell_n \cos\theta_n \right|,
\qquad
F_{2y}(i,j) \;=\; \frac{1}{\Delta y(i,j)} \sum_{n \in \mathcal{S}(i,j)} \left| \ell_n \sin\theta_n \right|.
$$

Where:

- $\ell_n$ is the **geodesic length** (m) of coastline segment $n$,
- $\Delta x(i,j)$, $\Delta y(i,j)$ are **local grid metric lengths** (m) on the T-grid,
- $\theta_n$ is the **segment orientation** expressed in the **local model coordinate frame**,
- $\mathcal{S}(i,j)$ is the set of segments assigned to cell $(i,j)$ by nearest-neighbour mapping in a projected coordinate reference system (CRS).

A convenient plotted diagnostic is the magnitude:

$$
\left|F_2\right|(i,j) \;=\; \sqrt{F_{2x}(i,j)^2 + F_{2y}(i,j)^2}.
$$

---

## 2. Required grid geometry (CICE C-grid file)

The grid file supplies, per T-cell:

- T-cell centres: $\lambda(i,j)$ (lon), $\phi(i,j)$ (lat),
- local rotation angle: $\alpha(i,j)$ (radians), and
- metric lengths: $\Delta x(i,j)$, $\Delta y(i,j)$ (m).

**Unit handling** is robust:
- lon/lat are inferred as radians vs degrees by magnitude and converted to degrees if needed;
- $\Delta x,\Delta y$ are converted to meters using metadata when available (cm → m is supported).

Longitudes are normalised to $[-180^\circ, 180^\circ]$ for stable Antarctic projection transforms.

---

## 3. Coastline ingestion and conditioning (high-res polygon shapefile)

### 3.1 Feature filtering
The coastline layer is filtered by `surface` classes to retain only polygonal coastal boundaries representing:
- `land`, `ice shelf`, `ice tongue`, `rumple`.

### 3.2 Geometry repair and dissolve
To avoid artefacts from internal polygon boundaries (e.g., grounding-line edges between land and ice shelf polygons), geometries are optionally:
- repaired (e.g., `make_valid` or `buffer(0)` fallback), then
- **dissolved / unioned** in the native projected CRS prior to reprojection.

### 3.3 Exterior rings as lon/lat
The dissolved geometry is reprojected to EPSG:4326 and each polygon’s **exterior ring** is emitted as a vertex sequence:
$$
(\lambda_k, \phi_k)_{k=1}^K,
$$
with consecutive vertices forming line segments.

---

## 4. Segment geometry on the sphere (WGS84 geodesic)

For each coastline segment connecting $(\lambda_0,\phi_0)$ → $(\lambda_1,\phi_1)$:

- compute geodesic distance $\ell$ (m), and
- compute forward azimuth ($\mathrm{az}$) (degrees, clockwise from north).

The implementation converts azimuth to a segment angle in **east-north** convention:
$$
\gamma \;=\; \mathrm{rad}\left(90^\circ - \mathrm{az}\right),
$$
so $\gamma=0$ corresponds to +east.

---

## 5. Assigning segments to model cells (projected KDTree)

### 5.1 Project to Antarctic CRS
Segment endpoints are projected to an Antarctic planar CRS (default `EPSG:3031`). Segment midpoints:
$$
(x_m, y_m) \;=\; \tfrac{1}{2}(x_0+x_1,\, y_0+y_1)
$$
are used for mapping.

### 5.2 KDTree over selected T-cells
A KDTree is built from projected T-cell centres $(x_{ij}, y_{ij})$, with optional restrictions:

- latitude subset: $\phi(i,j) \le \phi_{\max}$ (default $-30^\circ$ for Antarctic relevance/performance;
- optionally restrict the KDTree to **coastal-ocean cells** only (ocean band within a small dilation of land), to reduce erroneous assignment to interior grounding-line/ice-shelf regions.

Each segment midpoint is assigned to the nearest KDTree T-cell; assignments beyond a maximum distance are rejected:
$$
d_{\min} \le d_{\max} \quad \text{(default } d_{\max}=50\text{ km)}.
$$

---

## 6. Local-angle transform and accumulation

For a segment assigned to cell $(i,j)$, the **local** segment angle is:
$$
\theta \;=\; \gamma \;-\; \alpha(i,j),
$$
where $\alpha(i,j)$ is the grid rotation angle.

The segment contributes projected lengths:
$$
s_x \;=\; \left|\ell \cos\theta\right|,
\qquad
s_y \;=\; \left|\ell \sin\theta\right|.
$$

Accumulators sum these within each cell:
$$
S_x(i,j)=\sum s_x, \qquad S_y(i,j)=\sum s_y,
$$
and final coast-only form factors are:
$$
F_{2x}^{\mathrm{coast}}(i,j)=\frac{S_x(i,j)}{\Delta x(i,j)},
\qquad
F_{2y}^{\mathrm{coast}}(i,j)=\frac{S_y(i,j)}{\Delta y(i,j)}.
$$

The coast-only output NetCDF stores $(F2cst_{2}, F2cst_{y})$, plus (optionally thinned) coastline vertices for provenance.

---

## 7. Grounded-iceberg (GI) contributions from a polygon GeoPackage

Kaihong’s Jiao grounded iceberg dataset is provided as polygons in `EPSG:3031`. For each polygon $p$:

- geometry area $\mathrm{A}_p$ (m$^2$) and perimeter $\mathrm{P}_p$ (m) are computed directly from the geometry;
- an attribute area (e.g., `Area_Mean_km2`) may be used when present, falling back to geometry area if missing;
- features can be **deduplicated by unique ID** (e.g., `Global_UID`) to avoid multiple detections of the same iceberg.

Each grounded iceberg is mapped to the nearest T-grid-cell using the same projected KDTree strategy, with a maximum mapping distance constraint.

### 7.1 Default "simple-geometry" GI parameterisation
Each mapped iceberg contributes an **isotropic projected length scale** $\mathrm{L}_p$ to both directions, scaled by a tunable coefficient $\mathrm{C}_{\mathrm{gi}}$:
$$
\mathrm{F2gi}_{x}(i,j) \;{+}{=}\; \mathrm{C}_{\mathrm{gi}} \frac{\mathrm{L}_p}{\Delta x(i,j)},
\qquad
\mathrm{F2gi}_{y}(i,j) \;{+}{=}\; \mathrm{C}_{\mathrm{gi}} \frac{\mathrm{L}_p}{\Delta y(i,j)},
$$

The default $\mathrm{L}_p$ is based on perimeter where available:
$$
\mathrm{L}_p \;=\; \frac{2}{\pi} \mathrm{P}_p,
$$
which corresponds to an isotropic mean absolute projection assumption (using $\mathbb{E}[|\cos\Theta|]=2/\pi$ for uniformly distributed $\Theta$).

If perimeter is missing or *unusable*, a circular-equivalent scale derived from area is used:
$$
\mathrm{L}_p \;=\; 4\sqrt{\frac{\mathrm{A}_p}{\pi}}.
$$

This method produces $\mathrm{F2gi}_{x}$ and $\mathrm{F2gi}_{y}$ plus diagnostics (mapped cell indices, mapping distances, area/perimeter vectors, and per-cell GI counts).

### 7.2 Optional “cluster-axis” GI parameterisation (for dense GI fields)

A second option clusters nearby grounded icebergs (DBSCAN-like) in projected space, estimates principal axes via a principal component analysis (PCA), and distributes an oriented cluster-scale contribution back to grounded icebergs. This is intended for treating densely packed iceberg “fields” as *anisotropic* obstacles. In the current workflow, **simple-geometry** is used.

---

## 8. Combined coast + GI form factors and NetCDF outputs

The final combined form factors are:

$$
\mathrm{F2}_{x} \;=\; \mathrm{F2cst}_{x} + \mathrm{F2gi}_{x},
\qquad
\mathrm{F2}_{y} \;=\; \mathrm{F2cst}_{y} + \mathrm{F2gi}_{y}.
$$

The combined NetCDF includes:
- totals: $\mathrm{F2}_{x}, \mathrm{F2}_{y}$,
- $x$-components: $\mathrm{F2cst}_{x}, \mathrm{F2gi}_{x}$
- $y$-components: $\mathrm{F2cst}_{y}, \mathrm{F2gi}_{y}$
- coastline vertex vectors (if stored), and
- GI diagnostic vectors (if enabled).

Key tunables include:
- coastline-to-cell mapping CRS and max distance,
- coastal-ocean band restriction (buffer cells),
- GI scaling $\mathrm{C}_{\mathrm{gi}}$,
- GI length-scale definition (perimeter vs area),
- optional clustering controls (if cluster-axis is enabled).

---

## 9. Practical interpretation

- Large $\mathrm{F2}$ values typically occur where the coastline (or dense GI) is geometrically complex within/near a cell, implying stronger parameterised form drag.
- The result is **cell-based** (T-grid) by construction; conversion to velocity points ($u$/$v$) can be deferred to the dynamical core as needed.
