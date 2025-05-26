# Methodology: Antarctic (Land)Fast (Sea) Ice Modelling and Metrics

This document describes the core methodologies implemented in the AFIM (Antarctic Fast Ice Modelling) framework, corresponding to functions and constructs in the repository [`src/`](https://github.com/dpath2o/AFIM/tree/main/src).

---

When classifying fast ice from numerical model sea ice simulation results, it is common practice to use sea ice concentration (`aice`) and sea ice speed (`ispd`). The former is a tracer variable derived on the spatial grid centre (commonly referred to as the /T-grid/). In contrast, `ispd` is derived from the momentum components (`u`, `v`), which are defined at displaced locations relative to the /T-grid/, forming what is referred to as the Arakawa /B-grid/. [This figure](https://raw.githubusercontent.com/dpath2o/AFIM/main/docs/figures/bgrid.png) shows a depiction of the computational [Arakawa B-grid](https://doi.org/10.1016/B978-0-12-460817-7.50009-4) in the vicinity of an island. Note that in this depiction, only one grid cell has been masked out as land. This is important, as it motivates the classification strategy: due to computational constraints, a no-slip boundary is imposed on any velocity cell adjacent to land (i.e., touching land). These velocity components are explicitly set to $0$ by the model. Velocity is defined on staggered /B-grid/ locations, this results in cells near the coast having artificially zero speed regardless of the physical ice state. This introduces ambiguity: a zero speed may result either from a valid physical condition (e.g., grounded ice) or from an imposed boundary condition.

In [CICE v6](https://github.com/CICE-Consortium/CICE), the momentum equation evolves horizontal velocity components ($u$, $v$) under various forces:

$$
\rho h \frac{\partial \vec{u}}{\partial t} = \vec{\tau}_a + \vec{\tau}_o - mf\hat{z} \times \vec{u} - mg\nabla H + \nabla \cdot \boldsymbol{\sigma} - C_d \vec{u}
$$

Where:
- $\vec{u}$ is ice velocity
- $\vec{\tau}_a$, $\vec{\tau}_o$ are wind and ocean stresses
- $f$ is the Coriolis parameter, $g$ is gravity
- $H$ is surface height
- $\boldsymbol{\sigma}$ is the internal ice stress tensor
- $C_d$ is the drag coefficient

In mathematical terms, the ice speed is computed as:

$$
|\vec{u}| = \sqrt{u^2 + v^2}
$$

But near land, either $u = 0$, $v = 0$, or both are imposed, and so $|\vec{u}| = 0$ even though ice may not be landfast. Because of this, AFIM implements multiple strategies to mitigate false classifications near land by using interpolated speeds or spatially averaged approaches (see Sections 2.1 and 2.2).

---

## 1. Primary Classification

A grid cell is classified as fast ice if:

$$
a \geq a_\text{thresh} \quad \text{and} \quad |\vec{u}| \leq u_\text{thresh}
$$

where:
* **sea ice concentration threshold**, $a_\text{thresh} = 0.15$
* **sea ice speed threshold**, $u_\text{thresh} \in{\{10^{-3}, 5 \times 10^{-4}\, 2.5 \times 10^{-4}}~\text{m/s}}$

Fast ice is identified from sea ice concentration ($a$) and speed ($|\vec{u}|$) using multiple thresholding methods.

Since we are masking (thresholding) values that are two grids ($a$ is on the /T-grid/ and $\vec{u}$ is on the /B-grid/) we need to decide on either re-gridding $a$ to the /B-grid/ or $\vec{u}$ to the /T-grid/. Since the underlying landmask file that I have been using for simulation is based on the /T-grid/ I choose to re-grid $u$ and $v$ sea ice velocity components to the /T-grid/. Other established results ([Lemieux et.al](https://doi.org/10.1002/2016JC012006) and [VanAchter et.al](https://doi.org/10.1016/j.ocemod.2021.101920)) have used this same approach. However, to be thorough in my reporting and analysis I have chosen to preserve the non-re-gridded ice speeds and hence have created the following four sea ice speed categories.

### 1.1 `ispd` categories:

1. `ispd_B`  :: no re-gridding
2. `ispd_Ta` :: spatially-averaged
3. `ispd_Tx` :: spatially-weighted-average
4. `ispd_BT` :: composite-mean of `ispd_B`, `ispd_Ta` and `ispd_Tx`

#### 1.1.1 `ispd_B`

From [`ispd_B`](https://github.com/dpath2o/AFIM/blob/main/src/sea_ice_processor.py#L228)

$$
|\vec{u}|_B = \sqrt{u^2 + v^2}
$$

#### 1.1.2 `ispd_Ta`

From [`ispd_Ta`](https://github.com/dpath2o/AFIM/blob/main/src/sea_ice_processor.py#L238)

To approximate the speed on the T-grid ($|\vec{u}|_T$), AFIM applies a spatial average of the B-grid speed $|\vec{u}|_B$ from the four surrounding corners:

$$
|\vec{u}|_T(i,j) = \frac{1}{4} \Big[ |\vec{u}|_B(i,j) + |\vec{u}|_B(i+1,j) + |\vec{u}|_B(i,j+1) + |\vec{u}|_B(i+1,j+1) \Big]
$$

This is equivalent to averaging the velocity magnitudes from the four B-grid corners around a T-grid center, mitigating the impact of single-point no-slip anomalies.

#### 1.1.3 `ispd_Tx`

From [`ispd_Tx`](https://github.com/dpath2o/AFIM/blob/273090c740618e4db7a5d835e614fa855a9fc793/src/sea_ice_processor.py#L190)

More generally, this type of spatial interpolation is an instance of **bilinear interpolation**, a common spatial re-gridding method in numerical modeling. The bilinear interpolation at a point $(x, y)$ within a cell bounded by $(x_1, y_1)$, $(x_2, y_2)$, with values $Q_{11}$, $Q_{21}$, $Q_{12}$, $Q_{22}$ at each corner, is given by:

$$
\begin{aligned}
f(x, y) &= \frac{1}{(x_2 - x_1)(y_2 - y_1)} \Big[ \\
&\quad Q_{11}(x_2 - x)(y_2 - y) + Q_{21}(x - x_1)(y_2 - y) \\
&\quad + Q_{12}(x_2 - x)(y - y_1) + Q_{22}(x - x_1)(y - y_1) \Big]
\end{aligned}
$$

This expression performs linear interpolation in $x$ followed by $y$, and provides smooth, continuous estimates across grid cells while not necessarily conserving quantities like mass or energy. More information on this specific method that I used can be obtained from [ESMF Regridding Documentation](https://earthsystemmodeling.org/regrid/) and [xESMF](https://xesmf.readthedocs.io/en/latest/notebooks/Compare_algorithms.html).

### 1.2 masking methods:

From the above four `ispd` categories there are then three ways in which to apply the masking/thresholding:
1. on daily-averaged $a$ and $\vec{u}$, **or** (see 1.2.1 below)
2. N-day-average $\bar{a}$ and $\bar{\vec{u}}$ (see 1.2.2 below)
3. persistence method, which uses daily-averaged $a$ and $\vec{u}$ (see 1.2.3 below)

Lastly, I then apply additional mainly geo-spatial criteria (see section 1.2.4 below) to ensure that resulting gridded dataset is southern hemisphere and that landmask file is correctly applied. **Then classify/mask for fast ice using***:

$$
FImask_{ispd-cat} = \bar{a} \geq a_\text{thresh} \quad \text{and} \quad \bar{u} \leq u_\text{thresh}
$$
where $u_{\text{thresh}}$ is one of the four sea ice speed categories: `ispd_B`, `ispd_Ta`, `ispd_Tx`, or `ispd_BT`.

This results in eight different fast ice classifications (or conversely **pack ice** classifications). Hence my naming scheme I have chosen to use the following nomenclature for brevity and remaining consistent with the underlying sea ice speed categories:
1. `B`, `Ta`, `Tx`, `BT`                     :: **no** temporal-averaging of $a$ and $\vec{u}$
2. `B_roll`, `Ta_roll`, `Tx_roll`, `BT_roll` :: rolling-averaging of $N$-days on \bar{a}$ and $\bar{\vec{u}}$

#### 1.2.1 daily-average method:
[Full method](https://github.com/dpath2o/AFIM/blob/273090c740618e4db7a5d835e614fa855a9fc793/src/sea_ice_processor.py#L434)

#### 1.2.2 rolling-average method:
More generally the $N$-day-average that [I employed](https://github.com/dpath2o/AFIM/blob/273090c740618e4db7a5d835e614fa855a9fc793/src/sea_ice_processor.py#L540) can be expressed as a rolling average over $N$-days (default $N=15$):

$$
\bar{a}(t) = \frac{1}{N} \sum_{\tau = t - N/2}^{t + N/2} a(\tau)
$$

#### 1.2.3 persistence method:

To identify fast ice based not only on instantaneous conditions but also on temporal persistence, AFIM applies a rolling boolean mask (or binary-days mask) as an alternate to the rolling average method described above. The conceptual approach is given in  grid cell is flagged as fast ice if it satisfies the fast ice condition (e.g., $a \geq 0.15$, $|\vec{u}| \leq \varepsilon$) for at least $M$ days within a centred window of $N$ days. This additional method for classifying fast ice was proposed by @adfraser, and that is to use the daily-averaged $a$ and $\vec{u}$ to determine if fast ice is /present/ in a grid cell. After $N$-days count the precedding $N$-days and if $M$ of those days were classified as fast ice then mark that cell as fast ice for those $N$-days.

Mathematically, define the daily binary mask for fast ice presence as:

$$
F(t, i, j) = \begin{cases} 1, & \text{if fast ice condition is met at } (i,j) \text{ on day } t \\ 0, & \text{otherwise} \end{cases}
$$

Then the boolean fast ice mask is:

$$
F_{\text{bool}}(t, i, j) = \begin{cases} 1, & \displaystyle \sum_{\tau = t - N/2}^{t + N/2} F(\tau, i, j) \geq M \\ 0, & \text{otherwise} \end{cases}
$$

Where:
- $N$ is the rolling window length (e.g., 7 days)
- $M$ is the minimum count threshold (e.g., 6 days)

This persistence filter is applied using xarray’s `rolling(...).construct(...).sum(...)` approach, centred in time, [here](https://github.com/dpath2o/AFIM/blob/273090c740618e4db7a5d835e614fa855a9fc793/src/sea_ice_processor.py#L620C1-L626C23)

#### 1.2.4 Additional criteria imposed:

#### 1.3.4.1 re-apply landmask:
After either doing (or not doing the temporal averaging) the dataset is then sub-set for particular hemisphere: `north` or `south` (default `south`).

#### 1.3.4.2 hemisphere masking:
After either doing (or not doing the temporal averaging) the dataset is then sub-set for particular hemisphere: `north` or `south` (default `south`).

---

## 2. Fast Ice Metrics

Metrics are computed via [`compute_fast_ice_metrics.py`](https://github.com/dpath2o/AFIM/blob/main/scripts/sea_ice_metrics/compute_fast_ice_metrics.py), with methods:

### 2.1 Fast Ice Area (`FIA`)  
From [`compute_fast_ice_area`](https://github.com/dpath2o/AFIM/blob/main/src/sea_ice_processor.py#L620):

$$
\text{FIA}(t) = \sum_{i,j} A_{i,j} M_{i,j}(t)
$$
Where:
- $M_{i,j}(t)$ is the boolean fast ice mask
- $A_{i,j}$ is the grid cell area

### 2.2 Fast Ice Persistence (`FIP`)
From [`compute_variable_aggregate`](https://github.com/dpath2o/AFIM/blob/273090c740618e4db7a5d835e614fa855a9fc793/src/sea_ice_processor.py#L636)

AFIM computes the temporal mean of fast ice concentration or other variables across a defined time window using a simple arithmetic mean. Let $a(t)_{ispd-cat}_$ represent a variable (e.g., sea ice concentration masked for fast ice category) at time $t$ defined on a fixed spatial grid. Then the temporal mean over $T$ time steps is defined by:

$$
\bar{a_{ispd-cat}} = \frac{1}{T} \sum_{t=1}^{T} a(t)_{ispd-cat}_
$$

This operation is used within AFIM to collapse a time series into a single climatological estimate, such as the 1993–1999 mean FIA.


---

## 3. Threshold Sensitivity
AFIM supports experiments with multiple $u_\text{thresh}$ values:
* $10^{-3}~\text{m/s}$ (default)
* $5 \times 10^{-4}~\text{m/s}$ (for sensitivity testing)

These are applied to both boolean and rolling classifications.

---

## 📁 Source Files
All methods above are implemented in:
* [`src/sea_ice_processor.py`](https://github.com/dpath2o/AFIM/blob/main/src/sea_ice_processor.py)
* [`scripts/process_fast_ice/process_fast_ice.py`](https://github.com/dpath2o/AFIM/blob/main/scripts/process_fast_ice/process_fast_ice.py)
* [`scripts/sea_ice_metrics/compute_fast_ice_metrics.py`](https://github.com/dpath2o/AFIM/blob/main/scripts/sea_ice_metrics/compute_fast_ice_metrics.py)

See also: `config.json` files for applied thresholds and flags.

