# Antarctic Fast Ice Modelling (AFIM)

AFIM is a Python library for analysing, processing, and visualising sea ice model output, particularly focused on **fast ice** and **pack ice** derived from the [CICE](https://github.com/CICE-Consortium/CICE) sea ice model.

This package is actively developed as part of a PhD research project at the University of Tasmania, in partnership with the Australian Antarctic Program Partnership (AAPP), and integrates a number of post-processing utilities for Antarctic sea ice metrics derived from CICE.

There are currently no other open-source tools dedicated to post-processing **CICE** model outputs for fast and pack ice with this level of structure, specificity, and automation.

📖 [CICE documentation (ReadTheDocs)](https://cice-consortium-cice.readthedocs.io/en/main/)
📘 [AFIM documentation (ReadTheDocs)](https://AFIM.readthedocs.io/en/main/)

[![Documentation Status](https://readthedocs.org/projects/afim/badge/?version=latest)](https://afim.readthedocs.io/en/latest/?badge=latest)
[![PyPI version](https://badge.fury.io/py/afim.svg)](https://badge.fury.io/py/afim)

---

## Features

- 📦 `SeaIceProcessor`: Compute either landfast sea ice or pack ice or sea ice metrics area.
- 📊 `SeaIcePlotter`: Plot spatial maps and time series of fast/pack ice, generate region-faceted maps and animations.
- 🌍 `GroundedIcebergProcessor` grounded iceberg masking and landmask modification (via `GroundedIcebergProcessor`).
- 🧊 Support for regridding, speed calculation, climatological masking, and multiple hemispheres.

---

## Installation

Clone the repository and install in editable mode:

```bash
git clone https://github.com/dpath2o/AFIM.git
cd AFIM
pip install -e ./src
```

Or create a full environment:

```bash
conda env create -f src/environment.yml
conda activate afim
```

---

## Usage

### Simple one-off use:

```python
from sea_ice_processor import SeaIceProcessor
dt0_str  = "1998-08-01"
dtN_str  = "1999-03-31"
sim_name = 'baseline'
SI_proc  = SeaIceProcessor(sim_name            = sim_name, 
                           ice_speed_threshold = 1e-4)
FI_lo_spd = SI_proc.process_window(dt0_str    = dt0_str,
                                    dtN_str    = dtN_str, 
                                    write_zarr = False,
                                    ow_zarrs   = True)
SI_proc2  = SeaIceProcessor(sim_name            = sim_name, 
                           ice_speed_threshold = 1e-3)
FI_hi_spd = SI_proc2.process_window(dt0_str    = dt0_str,
                                    dtN_str    = dtN_str, 
                                    write_zarr = False,
                                    ow_zarrs   = True)
```

### Batch Processing for Full Simulation:

See [`compute_sea_ice.py`](./scripts/compute_sea_ice.py) for a CLI-driven loop script:

```bash
python compute_sea_ice.py Rothrock --sea_ice
```

---

## Project Structure

```
AFIM/
├── src/
│   ├── JSONs/                        # JSON config files for variable metadata and plotting
│   ├── sea_ice_processor.py          # sea ice processor 
│   ├── sea_ice_plotter.py            # Visualisation for fast and pack ice data
│   ├── grounded_iceberg_processor.py # Handles grounded iceberg landmask integration
│   ├── __init__.py                   # Package initializer
│   ├── requirements.txt              # Python dependencies
│   └── environment.yml               # Conda environment file
├── scripts/
│   ├── compute_sea_ice.py           # Looping script for computing fast ice
├── notebooks/
│   └── fi_anal.ipynb                 # Jupyter notebook example
├── README.md                         # Project overview and usage instructions
├── setup.py                          # Installable package configuration
└── docs/                             # Sphinx documentation source (optional)
```

---

## Background
In early 2020 I choose to embark on a PhD in oceanography. This decision was heavily influenced by my background and my interest in [ocean modelling](http://www.cmar.csiro.au/staff/oke/pubs/England_and_Oke_2001.pdf), the [Southern Ocean](https://tos.org/oceanography/issue/volume-25-issue-03) and [Antarctica](https://www.scar.org). This was spurned largely from my
[previous](https://www.cencoos.org) [professional life](http://imos.org.au) as a [coastal oceanographer](https://scripps.ucsd.edu/research/topics/coastal-oceanography). That previous pursuit centred mainly around a remote sensing technology called [high frequency radar](https://tos.org/oceanography/assets/docs/10-2_paduan1.pdf) which is, and continues to be, used in a broad number of ways to understand the dynamics (the motion of the of the upper ocean) through the digital application of signal processing of the physical representaiton of a Doppler shifted Bragg frequency. 

After 10 years in that field I left to diversify my career skillset while ticking some service-related philosophies about being [contributing to a new nation](https://en.wikipedia.org/wiki/National_service). After four years of [learning how-to drive Navy ships](https://www.navy.gov.au/sites/default/files/documents/Warfare_Officers_Career_Handbook.pdf) I specialised as a [Meteorologic and Oceanographic Officer](https://www.defencejobs.gov.au/jobs/reserves/navy/meteorologist-and-oceanographer) in that [organisation](https://www.navy.gov.au) and [operational-ised](https://www.youtube.com/watch?v=_1roFUwV7ss) my [scientific skillset](https://oceansci.ucsc.edu/academics/graduate/ms.html). At about the same time I became fascinated with Antarctica and its criticality both [climatically](https://tos.org/oceanography/article/southern-ocean-warming) and [strategically](https://defence.gov.au/adc/Publications/AJDSS/documents/volume3-number2/Where-to-from-here-The-Australian-Defence-Forces-pursuit-of-national-security-and-the-2020-Defence-Strategic-update.pdf) (and [here](https://www.antarctica.gov.au/about-us/antarctic-strategy-and-action-plan/)). Hence I began searching for a PhD project that would fulfil my interest.

Fortunately, it did not take too long before I was introduced to [Alex Fraser](https://tasmanian.com.au/stories/alex-fraser/) and [a project](./ResearchPlan/project_proposal/PROJECT_PROPOSAL.pdf) that he had on his digital shelf that from the outset ticked all my interest boxes. I can\'t recall exactly where or when in I was introduced to the various forms of [sea ice](https://en.wikipedia.org/wiki/Sea_ice) and [its role in polar oceanography](https://tos.org/oceanography/issue/volume-24-issue-03) but I do remember it being very [Arctic-centric](http://nsidc.org/arcticseaicenews/) with this vague concept that [landfast sea ice (fast ice)](https://arctic.noaa.gov/Report-Card/Report-Card-2018/ArtMID/7878/ArticleID/788/Landfast-Sea-Ice-in-a-Changing-Arctic)played an important role in Arctic sea ice. So when Alex provided me with a rough outline of the project he wanted to pursue modelling fast ice and it being [largely neglected thus far](https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&ved=2ahUKEwiWo5bfuPv2AhWTjeYKHSlPCycQFnoECCcQAQ&url=https%3A%2F%2Fdipot.ulb.ac.be%2Fdspace%2Fbitstream%2F2013%2F336850%2F1%2Fdoi_320494.pdf&usg=AOvVaw3fsCaFuWB9oz-dzVqWsQUW) in circumpolar modelling efforts, I knew the meat on this bone was likely well marbled.

In mid-2021 I was accepted and enrolled as a part-time PhD student at the [University of Tasmania](https://www.utas.edu.au) [Institute of Marine Science](https://www.imas.utas.edu.au) through the
[Australian Antarctic Partnership Program](https://aappartnership.org.au). I spent the second-half of that year [reading](./references) and constructing an initial [Reasearch Plan](./ResearchPlan/doc/researchplan.pdf). However, as all good projects evolve, it became apparent that I should be aligning my project with a well-supported Australian sea ice modelling effort, not only for my own support, but also to aim for a more impactful goal of incorporating fast ice into a nationally recognised/supported model.

After a brief interlude in the first quarter of 2022 I\'m now pursuing fast ice modelling with an eye towards incorporating it into [COSIMA](http://cosima.org.au).
