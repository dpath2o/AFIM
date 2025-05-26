# Methodology: Antarctic (Land)Fast (Sea) Ice Modelling and Metrics

This document describes the core methodologies implemented in the AFIM (Antarctic Fast Ice Modelling) framework, corresponding to functions and constructs in the repository [`src/`](https://github.com/dpath2o/AFIM/tree/main/src).

---

When classifying fast ice from numerical model sea ice simulation results, it is common practice to use sea ice concentration (`aice`) and sea ice speed (`ispd`). The former is a tracer variable derived on the spatial grid centre (commonly referred to as the T-grid). In contrast, `ispd` is derived from the momentum components (`u`, `v`), which are defined at displaced locations relative to the T-grid, forming what is referred to as the Arakawa B-grid.

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

☝️ [This figure](https://raw.githubusercontent.com/dpath2o/AFIM/main/docs/figures/bgrid.png) shows a depiction of the computational [Arakawa B-grid](https://doi.org/10.1016/B978-0-12-460817-7.50009-4) in the vicinity of an island.

Note that in this depiction, only one grid cell has been masked out as land. This is important, as it motivates the classification strategy: due to computational constraints, a no-slip boundary is imposed on any velocity cell adjacent to land (i.e., touching land). These velocity components are explicitly set to $0$ by the model.

Because velocity is defined on staggered B-grid locations, this results in many cells near the coast having artificially zero speed regardless of the physical ice state. This introduces ambiguity: a zero speed may result either from a valid physical condition (e.g., grounded ice) or from an imposed boundary condition.

In mathematical terms, the ice speed is computed as:

$$
|\vec{u}| = \sqrt{u^2 + v^2}
$$

But near land, either $u = 0$, $v = 0$, or both are imposed, and so $|\vec{u}| = 0$ even though ice may not be landfast.

Because of this, AFIM implements multiple strategies to mitigate false classifications near land by using interpolated speeds or spatially averaged approaches (see Sections 2.1 and 2.2).

---


## 1. Fast Ice Masking Criteria

Fast ice is identified from sea ice concentration ($a$) and speed ($|\vec{u}|$) using multiple thresholding methods.

### 1.1 Sea ice speed calculation methods

### 1.1 Native B-grid speed

From [`ispd_B`](https://github.com/dpath2o/AFIM/blob/main/src/sea_ice_processor.py#L228):

$$
|\vec{u}|_B = \sqrt{u^2 + v^2}
$$

### 1.2 `Ta`: spatially-averaged B-grid to T-grid
From [`ispd_Ta`](https://github.com/dpath2o/AFIM/blob/main/src/sea_ice_processor.py#L238):

$$
|\vec{u}|_Ta = \text{mean of 4 surrounding B-grid corners}
$$



### 1.3 `Tx`: ESMF spatially weighted average using a /bilinear/ approach to map B-grid to T-grid
From [`ispd_Tx`](https://github.com/dpath2o/AFIM/blob/main/src/sea_ice_processor.py#L315):

$$
|\vec{u}|_Ta = \text{mean of 4 surrounding B-grid corners}
$$

### 1.1 Sea ice speed mask 
From [`compute_fast_ice_speeds`](https://github.com/dpath2o/AFIM/blob/main/src/sea_ice_processor.py#L224):

A grid cell is classified as fast ice if:

$$
a \geq a_\text{thresh} \quad \text{and} \quad |\vec{u}| \leq u_\text{thresh}
$$

where:
- $a_\text{thresh} = 0.15$ (sea ice concentration threshold)
- $u_\text{thresh} \in \{10^{-3}, 5 \times 10^{-4}\}~\text{m/s}$

---

### 1.2 Rolling Mean Then Masking  
From [`compute_fast_ice_rolling`](https://github.com/dpath2o/AFIM/blob/main/src/sea_ice_processor.py#L287):

Apply a rolling average over $N$ days (default $N=15$):

$$
\bar{a}(t) = \frac{1}{N} \sum_{\tau = t - N/2}^{t + N/2} a(\tau)
$$

Then classify as fast ice using:

$$
\bar{a} \geq a_\text{thresh} \quad \text{and} \quad \bar{u} \leq u_\text{thresh}
$$

---



---

## 3. Daily vs. Rolling-Averaged Workflow

### 3.1 Daily Processing  
From [`process_daily_cice`](https://github.com/dpath2o/AFIM/blob/main/src/sea_ice_processor.py#L434):
Masking is applied directly to daily fields.

### 3.2 Rolling Averaged Processing  
From [`process_rolling_cice`](https://github.com/dpath2o/AFIM/blob/main/src/sea_ice_processor.py#L481):
A rolling average is applied before thresholding.

---

## 4. Fast Ice Metrics

Metrics are computed via [`compute_fast_ice_metrics.py`](https://github.com/dpath2o/AFIM/blob/main/scripts/sea_ice_metrics/compute_fast_ice_metrics.py), with methods:

### 4.1 Fast Ice Area (FIA)  
From [`compute_fast_ice_area`](https://github.com/dpath2o/AFIM/blob/main/src/sea_ice_processor.py#L620):

$$
\text{FIA}(t) = \sum_{i,j} A_{i,j} M_{i,j}(t)
$$
Where:
- $M_{i,j}(t)$ is the boolean fast ice mask
- $A_{i,j}$ is the grid cell area

### 4.2 Fast Ice Concentration (FIC)  
From [`compute_fast_ice_concentration`](https://github.com/dpath2o/AFIM/blob/main/src/sea_ice_processor.py#L628):

$$
\text{FIC}_{i,j}(t) = M_{i,j}(t) \cdot a_{i,j}(t)
$$

### 4.3 Fast Ice Thickness (FIH)  
From [`compute_fast_ice_thickness`](https://github.com/dpath2o/AFIM/blob/main/src/sea_ice_processor.py#L636):

$$
\text{FIH}_{i,j}(t) = M_{i,j}(t) \cdot h_{i,j}(t)
$$

---

## 5. Threshold Sensitivity

AFIM supports experiments with multiple $u_\text{thresh}$ values:
- $10^{-3}~\text{m/s}$ (default)
- $5 \times 10^{-4}~\text{m/s}$ (for sensitivity testing)

These are applied to both boolean and rolling classifications.

---

## 📁 Source Files
All methods above are implemented in:
- [`src/sea_ice_processor.py`](https://github.com/dpath2o/AFIM/blob/main/src/sea_ice_processor.py)
- [`scripts/sea_ice_metrics/compute_fast_ice_metrics.py`](https://github.com/dpath2o/AFIM/blob/main/scripts/sea_ice_metrics/compute_fast_ice_metrics.py)

See also: `config.json` files for applied thresholds and flags.

