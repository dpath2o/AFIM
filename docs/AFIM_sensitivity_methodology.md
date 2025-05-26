# Methodology: Antarctic Fast Ice Modelling and Metrics

This document describes the core methodologies implemented in the AFIM (Antarctic Fast Ice Modelling) framework, corresponding to functions and constructs in the repository [`src/`](https://github.com/dpath2o/AFIM/tree/main/src).

---

## 1. Fast Ice Masking Criteria

Fast ice is identified from sea ice concentration ($a$) and speed ($|\vec{u}|$) using multiple thresholding methods.

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

## 2. Ice Speed Calculation Methods

### 2.1 Native B-grid Speed  
From [`compute_speed_B`](https://github.com/dpath2o/AFIM/blob/main/src/sea_ice_processor.py#L305):

$$
|\vec{u}|_B = \sqrt{u^2 + v^2}
$$

### 2.2 Interpolated to T-grid (Averaged)  
From [`compute_speed_Ta`](https://github.com/dpath2o/AFIM/blob/main/src/sea_ice_processor.py#L315):

$$
|\vec{u}|_T = \text{mean of 4 surrounding B-grid corners}
$$

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

