# Configuration file for the Sphinx documentation builder.
import os
import sys
from unittest.mock import MagicMock
from pathlib import Path
# Ensure AFIM sources are importable
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
# Only on RTD: mock heavy/optional deps so autodoc can import AFIM modules
if os.environ.get("READTHEDOCS") == "True":
    MOCK_MODULES = [
        # PyGMT / GMT
        "pygmt", "pygmt.helpers", "pygmt.clib",
        # Cartopy stack
        "cartopy", "cartopy.crs", "cartopy.feature", "cartopy.io", "cartopy.io.shapereader",
        # Geo stack
        "geopandas", "pyogrio", "shapely", "shapely.geometry", "shapely.ops", "pyproj",
        "rasterio",
        # Regridding / ESMF
        "xesmf", "ESMF", "esmpy",
        # Odds and ends sometimes pulled in by plotting
        "cmocean",
    ]
    for name in MOCK_MODULES:
        sys.modules.setdefault(name, MagicMock())
# sys.path.insert(0, os.path.abspath('../../src'))
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
# -- Autodoc imports ---------------------------------------------------
autodoc_mock_imports = [ "dask.distributed",
    # plotting stack (not available on RTD)
    "pygmt",
    "geopandas",
    "PIL",
    "PIL.Image",
    "imageio",
    # regridding stack (xesmf often not installable on RTD)
    "xesmf",
    "ESMF",
    "esmpy",
    "pyresample",
    "pyresample.geometry",
    "pyresample.kd_tree",
    "pyproj",
    # ACCESS intake catalog
    "intake",
]