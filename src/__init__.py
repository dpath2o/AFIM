"""
AFIM (Antarctic Fast Ice and Pack Ice Metrics)

This package provides tools to process and visualise sea ice simulations from the CICE model,
including fast ice and pack ice diagnostics, spatial regridding, and PyGMT-based visualisation.

Submodules
----------
fast_ice_processor         : Compute fast ice diagnostics over time windows.
pack_ice_processor         : Compute pack ice metrics from rolling window CICE outputs.
sea_ice_plotter            : Visualise sea ice concentration and metrics.
grounded_iceberg_processor : Load and process grounded iceberg masks and land configurations.
"""
__version__ = '0.1.0'
from .fast_ice_processor import FastIceProcessor
from .pack_ice_processor import PackIceProcessor
from .sea_ice_plotter import SeaIcePlotter
from .grounded_iceberg_processor import GroundedIcebergProcessor
