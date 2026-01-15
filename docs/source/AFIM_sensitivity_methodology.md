# AFIM Sensitivity Study

This document describes the simulation setup and core methodologies implemented for an Antarctic Fast Ice Modelling (AFIM) sensitivity study and published here: [link to be provided once published]

## 0. Model Configuration & Setup

+ **Simulation period** : [1993-01-01 Fri] to [1999-12-31 Tue]
+ **Grid** : $1/4^{\circ}$ global, *tripole* Arakawa-B (`B-grid`)
+ **Modified Landmask** : a grid cell is determined to contain a (or a number of) sub-grid scale grounded icebergs and a land cell is created where normally an ocean cell would exist [see sea_ice_toolbox.sea_ice_icebergs.modify_landmask_with_grounded_icebergs()](https://github.com/dpath2o/AFIM/blob/1f8881284c82781579cabf898f6776cde7bc18df/src/sea_ice_icebergs.py#L208)
+ `dt`   : 1800 seconds
+ `ndte` : 240
+ `kdyn` : 1 $\textrightarrow$ EVP solver
+ Initial conditions : None

### 0.1 forcing:

#### ocean

+ [ECMWF Ocean Re-analysis version 5 (ORAS)](https://www.cen.uni-hamburg.de/en/icdc/data/ocean/easy-init-ocean/ecmwf-oras5.html)
+ regridded to the grid file above for 30-year period [1993-01-01 Fri] to [2023-12-31 Sun]

#### atmosphere

+ `ERA5`
+ regridded to the grid file above for 30-year period [1993-01-01 Fri] to [2023-12-31 Sun]

---

## Fast-ice classification (`sea_ice_classification.py`)

This module classifies **fast ice (FI)**, **pack ice (PI)**, and **total sea ice (SI)** from CICE output using:
1) a concentration threshold (`aice`), and  
2) a speed threshold applied to a **T-grid** speed magnitude `ispd_T`.

Because concentration is defined on the **T-grid** (cell centres) but model velocities may be staggered (legacy **B-grid** corners, or **C-grid** edges), the module first constructs `ispd_T` using one of several reconstruction strategies. It then forms a **daily candidate** mask and optionally applies (i) a centred rolling-mean diagnostic and (ii) a binary-day persistence filter.

### 1. Notation

Let:
- `a(t,i,j)` be sea-ice concentration (`aice`) at day `t` on a T-cell `(i,j)`;
- `|\vec{u}|_T(t,i,j)` be the speed magnitude on the analysis grid (T-grid);
- `a_th = icon_thresh` (typically 0.15);
- `u_th = ispd_thresh` (e.g. 1e-3, 5e-4, 2.5e-4 m/s).

The speed is always computed as:
$$
|\vec{u}| = \sqrt{u^2+v^2}
$$

### 2. Constructing T-grid speed (`ispd_T`)

The selection is controlled by `BorC2T_type`. Supported tokens:

#### 2.1 B-grid derived options (can be combined and averaged)

These options start from corner-staggered B-grid components `u_B, v_B`.

**(Ta) 2×2 corner mean (NaNs propagate)**  
$$
u_{Ta}(t,i,j)=\frac14\sum_{\delta i,\delta j\in\{0,1\}}u_B(t,i+\delta i,j+\delta j),
\quad
v_{Ta}(t,i,j)=\frac14\sum_{\delta i,\delta j\in\{0,1\}}v_B(t,i+\delta i,j+\delta j),
$$
$$
|\vec{u}|_{Ta}=\sqrt{(u_{Ta})^2+(v_{Ta})^2}.
$$

**(Tb) 2×2 corner mean with no-slip fill (NaNs → 0)**  
Same as Ta, but replace NaNs with 0 before averaging:
$$
\tilde{u}_B=\mathrm{nan2zero}(u_B),\quad \tilde{v}_B=\mathrm{nan2zero}(v_B).
$$

**(Tx) Regridding (weights, e.g. xESMF)**
Corner-staggered `u_B, v_B` are mapped to T-grid using a precomputed regridder, after NaNs are replaced by 0:
$$
u_{Tx}=\mathcal{R}(\tilde{u}_B),\quad v_T^{Tx}=\mathcal{R}(\tilde{v}_B),\quad
|\vec{u}|_{Tx}=\sqrt{(u_{Tx})^2+(v_{Tx})^2}.
$$

**Composite (if multiple B-grid modes are selected)**  
If more than one of `{Ta, Tb, Tx}` is requested, the module forms a composite speed:
$$
|\vec{u}|_T = \mathrm{mean}\left(\,|\vec{u}|_{Ta},|\vec{u}|_{Tb},|\vec{u}|_{Tx}\,\right)
$$
with missing values skipped.

#### 2.2 C-grid derived option (Tc; exclusive)

For C-grid output, velocities are provided as east/north components on **edge staggers**:
- U-stagger: `uvelE, uvelN`
- V-stagger: `vvelE, vvelN`

Edges are mapped to centres by adjacent averaging along the staggered direction, then combined (default: mean of U→T and V→T estimates):
$$
E_T=\tfrac12(E_{U\to T}+E_{V\to T}),\qquad
N_T=\tfrac12(N_{U\to T}+N_{V\to T}),
\qquad
|\vec{u}|_{Tc}=\sqrt{E_T^2+N_T^2}.
$$

**Important:** `Tc` is enforced as **exclusive** and is not mixed with `Ta/Tb/Tx`, because it is a different reconstruction pathway.

### 3. Daily candidate fast‑ice mask

For a chosen T‑grid speed product $s_T$, AFIM defines a *daily candidate* mask

$$
\mathcal{M}_{\mathrm{FI,day}}(t,i,j) =
\begin{cases}
1, & a(t,i,j) > a_{\text{thresh}} \;\wedge\; 0 < s_T(t,i,j) \le s_{\text{thresh}},\\
0, & \text{otherwise}.
\end{cases}
$$

The thresholds are:

- $a_{\text{thresh}} = 0.15$ (15% concentration)
- $s_{\text{thresh}} = \varepsilon$ (default $\varepsilon = 5\times 10^{-4}\,\mathrm{m\,s^{-1}}$, with sensitivity tests at other values)

**Important:** the strict inequality $s_T>0$ is deliberate: it reduces false positives from land‑adjacent *no‑slip* zeros.

#### 3.1 Pack ice (PI):
$$
\mathcal{M}_{PI,day}(t,i,j)=\mathbb{I}[a(t,i,j)>a_{th}]\ \mathbb{I}[|\vec{u}|_T(t,i,j)>u_{th}].
$$

#### 3.2 Total sea ice (SI):
$$
\mathcal{M}_{SI,day}(t,i,j)=\mathbb{I}[a(t,i,j)>a_{th}].
$$

### 4. Rolling‑mean fast‑ice mask (thresholding the smoothed speed)

AFIM also implements a rolling‑mean approach that mirrors common practice in the literature. Let $P$ be the rolling‑mean period (default $P=15$ days). Define the centred rolling mean speed:

$$
\overline{s}_T(t,i,j)=\frac{1}{P}\sum_{\tau\in\mathcal{P}(t)} s_T(\tau,i,j),
$$

and then apply the same daily threshold logic to $\overline{s}_T$:

$$
\mathcal{M}_{\mathrm{FI,roll}}(t,i,j) =
\mathbf{1}\left(a(t,i,j)>a_{\text{thresh}}\right)\,
\mathbf{1}\left(0<\overline{s}_T(t,i,j)\le s_{\text{thresh}}\right).
$$

AFIM returns (when requested) the triplet of masks:
${\mathcal{M}_{\mathrm{FI,day}}, \mathcal{M}_{\mathrm{FI,bin}}, \mathcal{M}_{\mathrm{FI,roll}}}$.

### 5. Binary‑days persistence fast‑ice mask (primary diagnostic)

The binary‑days method imposes persistence using a centred rolling count. Let $W$ be the window length (days) and $N$ the required number of “fast‑ice days” in the window. Define:

$$
\mathrm{C}(t,i,j) = \sum_{\tau \in \mathcal{W}(t)} \mathcal{M}_{\mathrm{FI,day}}(\tau,i,j),
$$

where $\mathcal{W}(t)$ is the centred $W$-day window around $t$. The persistent fast‑ice mask is then:

$$
\mathcal{M}_{\mathrm{FI,bin}}(t,i,j) =
\begin{cases}
1, & \mathrm{C}(t,i,j) \ge N,\\
0, & \text{otherwise}.
\end{cases}
$$

Default AFIM values are $W=11$ and $N=9$, allowing up to $W-N=2$ “mobile” days inside the window while still classifying the cell as persistently fast.

Implementation details:

- AFIM uses `rolling(time=W, center=True, min_periods=N)` so edge times can still be classified when at least $N$ days are present.
- Internally, the time range is extended by $\max(\lfloor W/2\rfloor, \lfloor P/2\rfloor)$ days (where $P$ is the rolling‑mean period) to minimise edge effects, then cropped back to the requested $[t_0,t_N]$.

Operational detail: the workflow evaluates daily masks over an **extended time span** (padding by approximately $⌊W/2⌋$ days, and by $⌊R/2⌋$ when rolling output is enabled), then **crops back** to the requested analysis interval. This reduces edge artefacts from centred windows while retaining daily timing information.

Unless otherwise stated, AFIM’s primary fast-ice diagnostic uses $\mathcal{M}_{FI,bin}$ rather than the instantaneous candidate mask $\mathcal{M}_{FI,day}$.

---

## 2. Fast‑ice metrics (SeaIceMetrics)

AFIM computes time‑series metrics (area/volume/thickness) and gridded diagnostics (persistence, stability index, persistence distance) using generic “hemispheric” operators applied to fast‑ice‑masked fields.

### 2.1 Notation and masking convention

Let:

- \(A_{ij}\) be the grid‑cell area (`tarea`; m\(^2\)).
- \(a(t,i,j)\) be sea‑ice concentration (`aice`; 0–1).
- \(h(t,i,j)\) be sea‑ice thickness (`hi`; m).
- \(\mathcal{M}(t,i,j)\in\{0,1\}\) be a chosen fast‑ice mask (typically \(\mathcal{M}_{\mathrm{FI,bin}}\)).

Define masked fields:

\[
a_{\mathrm{FI}}(t,i,j)=a(t,i,j)\,\mathcal{M}(t,i,j), \qquad
h_{\mathrm{FI}}(t,i,j)=h(t,i,j)\,\mathcal{M}(t,i,j).
\]

AFIM then computes metrics from \(a_{\mathrm{FI}}\) and \(h_{\mathrm{FI}}\), optionally applying a concentration threshold (default \(a>0.15\)) within the metric operators.

---

### 2.2 Fast‑ice area (FIA)

The concentration‑weighted fast‑ice area time series is:

$$
\mathrm{FIA}(t)=\sum_{i,j} A_{ij}\,a_{\mathrm{FI}}(t,i,j).
$$

Many studies report a purely Boolean area, $\sum \mathrm{A}_{ij}\mathcal{M}$; the concentration‑weighted form above is what AFIM’s area operator computes when given a concentration field already masked by $\mathcal{M}$.

**Grounded icebergs (optional):** in experiments with grounded‑iceberg (GI) masking, AFIM reports:

$$
\mathrm{FIA}_{\mathrm{tot}}(t)=\mathrm{FIA}(t)+A_{\mathrm{GI}},
$$

where $A_{\mathrm{GI}}$ is the (time‑invariant) grounded‑iceberg footprint area on the model grid.

---

### 2.3 Fast‑ice volume (FIV)

Fast‑ice volume is computed as the area integral of concentration‑weighted thickness:

$$
\mathrm{FIV}(t)=\sum_{i,j} A_{ij}\,a_{\mathrm{FI}}(t,i,j)\,h(t,i,j).
$$

(Equivalently, $\sum \mathrm{A}_{ij}\,h_{\mathrm{FI}}\,a$; both are identical given the masking convention.)

---

### 2.4 Mean fast‑ice thickness (FIT)

The mean thickness of fast ice is defined as the volume‑to‑area ratio:

$$
\mathrm{FIT}(t)=\frac{\mathrm{FIV}(t)}{\mathrm{FIA}(t)}
=\frac{\sum_{i,j} A_{ij}\,a_{\mathrm{FI}}(t,i,j)\,h(t,i,j)}
{\sum_{i,j} A_{ij}\,a_{\mathrm{FI}}(t,i,j)}.
$$

This is an area‑weighted thickness over the fast‑ice footprint.

---

### 2.5 Fast‑ice persistence (FIP)

Given a binary mask $\mathcal{M}(t,i,j)$, the persistence over an interval $\mathcal{T}$ of length $T$ days is:

$$
\mathrm{FIP}(i,j)=\frac{1}{T}\sum_{t\in\mathcal{T}} \mathcal{M}(t,i,j).
$$

$\mathrm{FIP}\in[0,1]$ is a *fraction of days* classified as fast ice. It is often reported as:

- **Percent persistence:** $100\,\mathrm{FIP}$ (%)
- **Days of persistence:** $T\,\mathrm{FIP}$ (days)

For climatological persistence maps, $\mathcal{T}$ is typically the full analysis period or a season (e.g., austral winter).

---

### 2.6 Persistence Stability Index (FIPSI / PSI)

AFIM implements a “persistence stability index” (PSI) over winter months (default May–October). Let $\mathcal{T}_w$ be the set of winter days and define winter persistence:

$$
P_w(i,j)=\frac{1}{|\mathcal{T}_w|}\sum_{t\in\mathcal{T}_w}\mathcal{M}(t,i,j).
$$

Define:

- A *persistent* winter mask: $\mathcal{P}(i,j)=\mathbf{1}\left(P_w(i,j)\ge p_{\text{th}}\right)$, with default $p_{\text{th}}=0.8$.
- An *ever‑fast* winter mask: $\mathcal{E}(i,j)=\mathbf{1}\left(\max_{t\in\mathcal{T}_w}\mathcal{M}(t,i,j)>0\right)$.

Convert these to areas:

$$
A_{\text{pers}}=\sum_{i,j} \mathrm{A}_{ij}\,\mathcal{P}(i,j), \qquad
A_{\text{ever}}=\sum_{i,j} \mathrm{A}_{ij}\,\mathcal{E}(i,j),
$$

and define the PSI:

$$
\mathrm{PSI} = \frac{\mathrm{A}_{\text{pers}}}{\mathrm{A}_{\text{ever}}}\in[0,1].
$$

Interpretation: PSI quantifies how much of the *ever‑fast* winter footprint is *highly persistent*.

---

### 2.7 Persistence distance metrics (mean and max)

AFIM also computes distance‑based persistence statistics relative to the coastline:

1. Identify persistent fast‑ice grid cells (e.g., $P_w \ge p_{\text{th}}$).
2. Compute the minimum Euclidean distance from each persistent cell centre to the nearest coastline point (using Antarctic polar stereographic coordinates, `EPSG:3031`).
3. Report:
   - Mean persistence distance (km): $\overline{d}$
   - Maximum persistence distance (km): $d_{\max}$

These diagnostics summarise how far the most persistent fast ice extends offshore.

---

## 3. Recommended reporting for publications

For reproducibility, AFIM results should report:

- Speed product (`Ta`, `Tb`, `Tx`, `BT`, or `Tc`)
- Mask method `day`, `bin` (with $W,N$), or `roll` (with $P$)
- Thresholds ($a_{\text{thresh}}, s_{\text{thresh}}$)
- Whether GI footprint area is included in FIA totals
- The analysis window and (if relevant) the season used for persistence/PSI

---

## 4. Thresholds

### 4.1 sea ice speed threshold ( $u_\text{thresh}$ )

Choosing the appropriate $u_\text{thresh}$ has a significant effect on the classification of fast ice. These are values that have been used thus far, and their physical representation. 

+ $10^{-3}~\text{m/s}$, which translates to roughly 86 meters of **distributed** ice movement with*in** an Antarctic coastal grid cell (see **Grid** under model configuration above)
+ $5 \times 10^{-4}~\text{m/s}$, which translates to roughly 43 meters of **distributed** ice movement with*in** an Antarctic coastal grid cell (see **Grid** under model configuration above)
+ $2.5 \times 10^{-4}~\text{m/s}$, which translates to roughly 22 meters of **distributed** ice movement with*in** an Antarctic coastal grid cell (see **Grid** under model configuration above)

### 4.2 sea ice concentration threshold ( $a_\text{thresh}$  )

The selection of $a$ (sea ice concentration) has been kept at 15% of a grid cell.

