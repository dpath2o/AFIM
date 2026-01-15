# Configuration file for the Sphinx documentation builder.
import os
import sys
sys.path.insert(0, os.path.abspath('../../src'))
project = 'AFIM'
copyright = '2025, Daniel Atwater'
author = 'Daniel Atwater'
release = '0.1.0'
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.autosummary',
    'myst_parser',
    'nbsphinx',  # if using Jupyter notebooks
]
myst_enable_extensions = ["dollarmath", "amsmath", "colon_fence"]
nbsphinx_allow_errors = True
autosummary_generate = True
templates_path = ['_templates']
exclude_patterns = []
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
autodoc_mock_imports = [
    "pygmt",
    "geopandas",
    "pyogrio",
    "shapely",
    "pyproj",
    "xesmf",
    "ESMF",
    "netCDF4",
    "dask",
    "distributed",
    "cartopy",
    "cmocean",
    "rasterio",
]

