from setuptools import setup, find_packages

setup(
    name="AFIM",
    version="0.2.0",
    description="AFIM: Antarctic Fast Ice and Pack Ice Metrics toolkit for CICE model outputs",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    author="Daniel Patrick Atwater",
    author_email="daniel.atwater@utas.edu.au",
    url="https://github.com/dpath2o/AFIM",
    project_urls={
        "Documentation": "https://afim.readthedocs.io/en/latest/",
        "Source": "https://github.com/dpath2o/AFIM",
    },
    license="MIT",

    # You have module files in src/, not a package directory
    py_modules=[
        "sea_ice_toolbox",
        "sea_ice_gridwork",
        "sea_ice_classification",
        "sea_ice_observations",
        "sea_ice_metrics",
        "sea_ice_plotter",
        "sea_ice_regridder",
        "sea_ice_icebergs",
        "sea_ice_ACCESS",
    ],
    package_dir={"": "src"},
    include_package_data=True,

    # Keep this minimal and RTD-friendly
    install_requires=[
        "numpy",
        "pandas",
        "xarray",
        "scipy",
        "tqdm",
    ],

    extras_require={
        "analysis": [
            "dask",
            "distributed",
            "netCDF4",
        ],
        "geo": [
            "pyproj",
            "shapely",
            "geopandas",
            "pyogrio",
            "rasterio",
            "cartopy",
        ],
        "regrid": [
            "xesmf",
        ],
        "plot": [
            "matplotlib",
            "cmocean",
            "seaborn",
            "pygmt",
        ],
        "docs": [
            "sphinx>=7",
            "sphinx-rtd-theme",
            "myst-parser",
            "nbsphinx",
        ],
        "all": [
            "dask",
            "distributed",
            "netCDF4",
            "pyproj",
            "shapely",
            "geopandas",
            "pyogrio",
            "rasterio",
            "cartopy",
            "xesmf",
            "matplotlib",
            "cmocean",
            "seaborn",
            "pygmt",
        ],
    },

    python_requires=">=3.9",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Programming Language :: Python :: 3",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Atmospheric Science",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
