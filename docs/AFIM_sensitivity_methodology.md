# Methodology: Antarctic Fast Ice Modelling and Metrics

This document describes the core methodologies implemented in the AFIM (Antarctic Fast Ice Modelling) framework, corresponding to functions and constructs in the repository [src directory](../src). Code references are linked directly to the GitHub repository.

---

## 1. Fast Ice Masking Criteria

Fast ice is identified from sea ice concentration (`aice`) and speed (`|\vec{u}|`) using multiple thresholding methods:

### 1.1 Boolean Masking
From [`compute_fast_ice_boolean`](https://github.com/dpath2o/AFIM/blob/main/src/sea_ice_processor.py#L224):

A grid cell is classified as fast ice if:

\[ a_{i,j} \geq a_\text{thresh} \quad \text{and} \quad |\vec{u}_{i,j}| \leq u_\text{thresh} \]

where:
- \( a_\text{thresh} = 0.85 \) (sea ice concentration threshold)
- \( u_\text{thresh} \in \{10^{-3}, 5 \times 10^{-4}\} \, \text{m/s} \)

---

### 1.2 Rolling Mean First, Then Mask
From [`compute_fast_ice_rolling`](https://github.com/dpath2o/AFIM/blob/main/src/sea_ice_processor.py#L287):

A rolling mean is applied over a time window (default \( N=15 \) days):

\[ \bar{a}_{i,j}(t) = \frac{1}{N} \sum_{\tau = t-N/2}^{t+N/2} a_{i,j}(\tau) \]

Masking is then applied to the smoothed fields using the same boolean rule:

\[ \bar{a}_{i,j} \geq a_\text{thresh} \quad \text{and} \quad \bar{u}_{i,j} \leq u_\text{thresh} \]

---

## 2. Ice Speed Calculation Methods

Fast ice speed \( |\vec{u}| \) is derived using various techniques:

### 2.1 Native B-grid Speed
From [`compute_speed_B`](https://github.com/dpath2o/AFIM/blob/main/src/sea_ice_processor.py#L305):

\[ |\vec{u}|_\text{B} = \sqrt{u^2 + v^2} \]

### 2.2 Interpolated to T-grid (Average)
From [`compute_speed_Ta`](https://github.com/dpath2o/AFIM/blob/main/src/sea_ice_processor.py#L315):

\[ |\vec{u}|_\text{T} = \text{mean}(|\vec{u}|_\text{B} \text{ from surrounding corners}) \]

---

## 3. Daily vs. Averaged Processing

### 3.1 Daily Classification
[`process_daily_cice`](https://github.com/dpath2o/AFIM/blob/main/src/sea_ice_processor.py#L434):
Masking is applied directly to each daily field without any temporal smoothing.

### 3.2 Rolling-Averaged Classification
[`process_rolling_cice`](https://github.com/dpath2o/AFIM/blob/main/src/sea_ice_processor.py#L481):
A rolling average is applied before masking.

---

## 4. Fast Ice Metrics
Post-processed masks are analyzed in:
[`compute_fast_ice_metrics.py`](https://github.com/dpath2o/AFIM/blob/main/scripts/sea_ice_metrics/compute_fast_ice_metrics.py)

Key metrics include:

### 4.1 Fast Ice Area (FIA)
[`compute_fast_ice_area`](https://github.com/dpath2o/AFIM/blob/main/src/sea_ice_processor.py#L620):

\[ \text{FIA}(t) = \sum_{i,j} A_{i,j} M_{i,j}(t) \]
Where:
- \( M_{i,j}(t) \) is the fast ice mask
- \( A_{i,j} \) is grid cell area

### 4.2 Fast Ice Concentration (FIC)
[`compute_fast_ice_concentration`](https://github.com/dpath2o/AFIM/blob/main/src/sea_ice_processor.py#L628):

\[ \text{FIC}_{i,j}(t) = M_{i,j}(t) \cdot a_{i,j}(t) \]

### 4.3 Fast Ice Thickness (FIH)
[`compute_fast_ice_thickness`](https://github.com/dpath2o/AFIM/blob/main/src/sea_ice_processor.py#L636):

\[ \text{FIH}_{i,j}(t) = M_{i,j}(t) \cdot h_{i,j}(t) \]

---

## 5. Threshold Sensitivity
Simulations were run with two primary thresholds for ice speed masking:
- \( u_\text{thresh} = 10^{-3} \, \text{m/s} \)
- \( u_\text{thresh} = 5 \times 10^{-4} \, \text{m/s} \)

These were used in both boolean and rolling classifications.

---

## 📁 Source Reference
All methodology is implemented under [`src/`](https://github.com/dpath2o/AFIM/tree/main/src) and [`scripts/`](https://github.com/dpath2o/AFIM/tree/main/scripts).

For reproducibility, see `SeaIceProcessor`, `compute_fast_ice_metrics.py`, and `config.json` for input configuration.
