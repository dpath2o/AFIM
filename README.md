# Antarctic Fast Ice Modelling (AFIM)

AFIM is a Python library for analysing, processing, and visualising sea ice model output, particularly focused on **fast ice** and **pack ice** derived from the [CICE](https://github.com/CICE-Consortium/CICE) sea ice model.

This package is actively developed as part of a PhD research project at the University of Tasmania, in partnership with the Australian Antarctic Program Partnership (AAPP), and integrates a number of post-processing utilities for Antarctic sea ice metrics derived from CICE.

There are currently no other open-source tools dedicated to post-processing **CICE** model outputs for fast and pack ice with this level of structure, specificity, and automation.

+ [CICE documentation (ReadTheDocs)](https://cice-consortium-cice.readthedocs.io/en/main/)
+ [AFIM documentation (ReadTheDocs)](https://AFIM.readthedocs.io/en/main/)

[![Documentation Status](https://readthedocs.org/projects/afim/badge/?version=latest)](https://afim.readthedocs.io/en/latest/?badge=latest)
[![PyPI version](https://badge.fury.io/py/afim.svg)](https://badge.fury.io/py/afim)

### AFIM Sensitivity Study

+ [AFIM Sensitivity Experiment Design and Methodology](https://dpath2o.github.io/AFIM/AFIM_sensitivity_methodology.html)

+ [AFIM table describing the simulation names provided in the results presented below](https://dpath2o.github.io/AFIM/afim_simulation_configs.html)  

Explore key outputs from AFIM simulations:

| Gallery Type                      | Description                                                                                                             | Link                                                                           |
|-----------------------------------|-------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------|
| **FIP Simulation Tree Gallery**   | Simulation > Region > threshold > figures                                                                               | [View Gallery](https://dpath2o.github.io/AFIM/fip_tree_gallery.html)           |
| **Intra-model FIP Comparison** | Compare different FIP group types (`FI_BT_bool`, `FI_Ta_roll`, etc.) for a single simulation, region, and `ispd_thresh` | [View Gallery »](https://dpath2o.github.io/AFIM/fip_intra_model_gallery.html)  |
| **Inter-model FIP Comparison** | Compare multiple simulations for the same region, FIP group, and `ispd_thresh`                                          | [View Gallery »](https://dpath2o.github.io/AFIM/fip_inter_model_gallery.html)  |
| **FIA Timeseries (1993–1999)** | Smoothed fast ice area (FIA) timeseries across simulations                                                              | [View Timeseries »](https://dpath2o.github.io/AFIM/timeseries_gallery.html)    |
| **AFIM Archive Status**        | Table of all simulations and their processed outputs (FI/PI/SO, metrics, etc.)                                          | [View Status Table »](https://dpath2o.github.io/AFIM/AFIM_archive_status.html) |

---

## Highlights

### `SeaIceToolbox` — Unified API

- **All-in-one class** for Antarctic sea ice analysis:
  - Model I/O, fast/pack ice masking, Zarr output
  - Flexible plotting (maps, time series, faceted views)
  - Grounded iceberg masking and landmask modification
  - Observational data integration (Fraser et al. 2020, NSIDC)
- Boolean and rolling fast ice classification, spatial/temporal averaging, and regional metrics
- Regridding, derived fields, and observational overlays built-in
- Highly configurable via JSON

### Batch Processing & HPC-Ready

- Command-line and scriptable workflows for large-scale or HPC batch processing
- PBS job scripts and wrappers for automated monthly/seasonal runs

### Advanced Plotting

- Generate spatial maps, regional/faceted views, and time series
- Overlay observational and grounded iceberg data

### Grounded Iceberg Handling

- Mask grounded iceberg regions and modify landmasks for improved coastal fast ice simulation

### Regridding and Derived Fields

- Supports B-grid, T-grid (average and xESMF) interpolation
- Computes vector magnitudes for stress and velocity fields

### Interactive Analysis

- Jupyter notebook [`fi_anal.ipynb`](https://github.com/dpath2o/AFIM/blob/main/notebooks/fi_anal.ipynb) demonstrates full workflows

---

## Installation

```bash
git clone https://github.com/dpath2o/AFIM.git
cd AFIM
pip install -e ./src
```

Or using Conda:

```bash
conda env create -f src/environment.yml
conda activate afim
```

---

## Example: Interactive Python Usage

```python
from sea_ice_toolbox import SeaIceToolbox

SI_tools = SeaIceToolbox(
    sim_name='gi-mid',
    dt0_str='1998-07-01',
    dtN_str='1999-07-01',
    ice_speed_threshold=1e-3
)
FI = SI_tools.process_daily_cice()
SI_tools.plot_ice_area(FI)
```

---

## Project Layout

```
├── notebooks/                  # Jupyter notebooks
│   └── fi_anal.ipynb           # Main interactive analysis
├── scripts/
│   └── process_fast_ice/
│       ├── process_fast_ice.py           # Legacy CLI script (use SeaIceToolbox for new workflows)
│       ├── process_fast_ice.pbs          # PBS job script
│       └── process_fast_ice_pbs_wrapper.sh # Wrapper to launch monthly jobs
├── src/
│   ├── sea_ice_toolbox.py      # Unified toolbox class (SeaIceToolbox)
│   ├── sea_ice_classification.py       # Model I/O and fast/pack ice logic (subclassed)
│   ├── sea_ice_plotter.py      # Plotting (subclassed)
│   ├── sea_ice_icebergs.py     # Grounded iceberg logic (subclassed)
│   ├── sea_ice_observations.py # Observational data integration (subclassed)
│   └── grounded_iceberg_processor.py
└── docs/                       # Documentation (Sphinx)
```

---

## Documentation

- Full API docs and method descriptions under `docs/`
- [fi_anal.ipynb](https://github.com/dpath2o/AFIM/blob/main/notebooks/fi_anal.ipynb) is the best place to start interactively

---

## Contributions

Contributions and feature requests welcome. Feel free to open an [Issue](https://github.com/dpath2o/AFIM/issues) or submit a pull request!

---

## Background

In early 2020, I chose to embark on a PhD in oceanography. This decision was influenced by my background and interest in [ocean modelling](http://www.cmar.csiro.au/staff/oke/pubs/England_and_Oke_2001.pdf), the [Southern Ocean](https://tos.org/oceanography/issue/volume-25-issue-03), and [Antarctica](https://www.scar.org). This interest stemmed in large part from my [previous](https://www.cencoos.org) [professional life](http://imos.org.au) as a [coastal oceanographer](https://scripps.ucsd.edu/research/topics/coastal-oceanography). That work focused mainly on a remote sensing technology called [high frequency radar](https://tos.org/oceanography/assets/docs/10-2_paduan1.pdf), which continues to be widely used to understand upper ocean dynamics via the digital signal processing of Doppler-shifted Bragg frequencies.

After a decade in that field, I stepped away to diversify my career skillset while also pursuing service-oriented goals. I spent four years [learning to drive Navy ships](https://www.navy.gov.au/sites/default/files/documents/Warfare_Officers_Career_Handbook.pdf), and then specialised as a [Meteorological and Oceanographic Officer](https://www.defencejobs.gov.au/jobs/reserves/navy/meteorologist-and-oceanographer) in the [Royal Australian Navy](https://www.navy.gov.au), applying my [scientific background](https://oceansci.ucsc.edu/academics/graduate/ms.html) to real-time operational contexts. During this time, I became fascinated by Antarctica, both for its [climatic importance](https://tos.org/oceanography/article/southern-ocean-warming) and its [strategic relevance](https://defence.gov.au/adc/Publications/AJDSS/documents/volume3-number2/Where-to-from-here-The-Australian-Defence-Forces-pursuit-of-national-security-and-the-2020-Defence-Strategic-update.pdf) ([see also](https://www.antarctica.gov.au/about-us/antarctic-strategy-and-action-plan/)). This led me to begin searching for a PhD project that would allow me to pursue these interests more deeply.

It wasn’t long before I was introduced to [Alex Fraser](https://tasmanian.com.au/stories/alex-fraser/) and a [project](./ResearchPlan/project_proposal/PROJECT_PROPOSAL.pdf) he had been holding onto that aligned perfectly with my interests. I can't recall exactly when I first learned about the different types of [sea ice](https://en.wikipedia.org/wiki/Sea_ice) or their role in [polar oceanography](https://tos.org/oceanography/issue/volume-24-issue-03), but my early understanding was definitely [Arctic-centric](http://nsidc.org/arcticseaicenews/). I had only a vague awareness that [landfast sea ice (fast ice)](https://arctic.noaa.gov/Report-Card/Report-Card-2018/ArtMID/7878/ArticleID/788/Landfast-Sea-Ice-in-a-Changing-Arctic) played a crucial role in the Arctic system.

When Alex outlined his idea for a fast ice modelling project — and noted that fast ice had been [largely neglected](https://dipot.ulb.ac.be/dspace/bitstream/2013/336850/1/doi_320494.pdf) in circumpolar sea ice modelling efforts — I was immediately intrigued. It was clear that this was a problem with depth, and one that hadn’t yet received the attention it deserved.

In mid-2021, I enrolled as a part-time PhD student at the [University of Tasmania](https://www.utas.edu.au), within the [Institute for Marine and Antarctic Studies (IMAS)](https://www.imas.utas.edu.au), through the [Australian Antarctic Program Partnership](https://aappartnership.org.au). I spent the second half of that year [reviewing the literature](./references) and drafting an initial [research plan](./ResearchPlan/doc/researchplan.pdf). As the project evolved, it became increasingly clear that I should align my work with a well-supported Australian sea ice modelling framework. This would not only benefit my development as a modeller but would also create a stronger pathway for integrating fast ice representation into a nationally supported climate model.

Following a brief pause in early 2022, I resumed my PhD project with a sharpened focus: to explore and implement fast ice modelling in the context of the [COSIMA](http://cosima.org.au) model framework.

In 2025, I was awarded the Australian Defence Force Chief of Defence Force Fellowship in recognition of my research into Antarctic fast ice and its relevance to climate and strategic studies.

