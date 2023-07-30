# Documentation
## [CICE5](https://cesmcice.readthedocs.io/en/latest/index.html)
## CICE6
* [latest version](https://cice-consortium-cice.readthedocs.io/en/main/)
* [CICE_6.2 PDF](https://cice-consortium-cice.readthedocs.io/_/downloads/en/cice6.2.0/pdf/)
# [CICE Consortium (NCAR-CESM)](https://www.cesm.ucar.edu/models/cesm2/sea-ice)

## [CICE Github Repository](https://github.com/CICE-Consortium)
### [Resources Index](https://github.com/CICE-Consortium/About-Us/wiki/Resource-Index)
### [CICE Consortium Test Data](https://github.com/CICE-Consortium/CICE/wiki/CICE-Input-Dat)
#### [Standalone Forcing](https://cice-consortium-cice.readthedocs.io/en/master/developer_guide/dg_forcing.html)
Describes details and helpful when pulling in more forcing data than the *default*. The following two lists show the fields that are required for forcing. The documentation states that all the
variables have reasonable static defaults that will be used in *default* mode. To advance the forcing, the subroutines `get_forcing_atmo` and `get_forcing_ocn` are called each timestep from the step loop. That subroutine computes the forcing year (`fyear`), calls the appropriate forcing data method, and then calls `prepare_forcing` which converts the input data fields to model forcing fields.
## GX3 Description
Global Three-Degree Cartesian

### GX3 Forcing
Three datasets are used for forcing [*GX3 JRA55*]
1.  GX3 JRA55
2.  GX3 NCAR
   The contents of the files can be dedecued from the file names:
      filename   my guess
      ---------- ------------------------------
      cldf       cloud
      prec       precipitation
      swdn       net downward shortwave
      dn10       dew-point temperature @ 10 m
      q10        specific humidity @ 10 m
      t10        air temperature @ 10 m
      u10        zonal wind speed @ 10 m
      v10        meridional wind speed @ 10 m

3.  GX3 WW3
    One wave spectral file is given with a total of four time steps and 25 frequencies. 

### GX3 Initial Conditions
There are a total of 12 initial condition files provided. With each file containing NIL time step across the This is interesting to me as I thought only one would be required. Nonetheless, the header of the first file is given below.

## GX1 Description
Global One-Degree Cartesian

### GX1 Forcing

Four datasets are used for forcing *GX1 CESM*

1.  GX1 CESM
    One file provides monthly summarised oceanographic data in the netCDF format. Note that the time dimension is only 12 increments in length.

2.  GX1 COREII
    These files are given as binary with two sub-directories, 4XDAILY and MONTHLY. The file names indicate the representative variable contained within.

3.  GX1 JRA55
    The forcing files are nearly identical to *GX3 JRA55*

4.  GX1 WOA
    Monthly statistics for silicate and nitrate are provided across twofiles. 

### GX1 Initial Conditions
These files are nearly identical to the *GX3 Initial Conditions*


## TX1 Description

Global One-Degree Tri-polar

### TX1 Grid
Here this grid file is given in both NetCDF and binary format. The header contents indicates are fair amount more information than the *GX3 Grid* file.

### TX1 Forcing
The only dataset utilised for the tri-polar forcing is JRA55. It is equivalent in it\'s structure and content to the other two JRA55 (*GX3 JRA55* and *GX1 JRA55*), but just re-gridded onto the *TX1 Grid*.
### TX1 Initial Conditions
There is no version 6 files provided *CICE Consortium*

## Initial Conditions Discussion

All three (GX3, GX1 and TX1) cases of *CICE Consortium Test Data* have the same structure of the initial condition files. Each case provides twelve files named for each month. The contents of those files show 56 gridded parameters over two and three dimensions (one categorical dimension and two spatial dimensions). There is not a time dimension. Hence time-step is implied as monthly with 12 files named for each month -- it is not known wether this is a snap snot or an average over the month. The table below is made in order to consolidate the initial condition parameters should at some point it be determined that initial conditions should be provided. The short name and description (the first two columns) are obtained from correlating the *CICE Consortium Test Data* information and the information found in the documentation, online [here](https://cice-consortium-cice.readthedocs.io/en/main/cice_index.html).

  | CICE6 short       |  CICE6 description                                 |    Dims      |   Units|
  |-------------------|----------------------------------------------------|--------------|--------|
  |uvel               | sea ice velocity, eastward                         | nj, ni       | m/s    |
  |vvel               | sea ice velocity, northward                        | nj, ni       | m/s    |
  |scale factor       | timestep offsets in shortwave radiation components | nj, ni       |  \-    |
  |swvdr               |incoming shortwave radiation, visible direct       | nj, ni       | W/m^2^ |
  |swvdf               |incoming shortwave radiation, visible diffuse      | nj, ni       | W/m^2^ |
  |swidr               |incoming shortwave radiation, infrared direct      | nj, ni       |  W/m^2^|
  |swidf               |incoming shortwave radiation, infrared diffuse     | nj, ni       |  W/m^2^|
  |strocnxT            |ice-ocean stress                                   | nj, ni       |  N/m   |
  |strocnyT            |  ice-ocean stress                                 | nj, ni       |  N/m   |
  |stressp\_\[1-4\]    |internal strain tensor                             |  nj, ni      |   N/m  |
  |stressm\_\[1-4\]    |  internal strain tensor                           |  nj, ni      |   N/m  |
  |stress12\_\[1-4\]   |internal strain tensor                             |  nj, ni      |  N/m   |
  |iceumask            |the mask where there is ice in a gridcell          |  nj, ni      |  \-    |
  |sst                 |sea surface temperature                            |  nj, ni      |  C     |
  |frzmlt              |ice ocean heat exhange                             |  nj, ni      |  W/m^2^|
  |frzonset            |Julian day in which the onset of freezing occurs   |  nj, ni      |   \-   |
  |ulat                |u-grid latitudes                                   |  nj, ni      |   \-   |
  |ulon                |u-grid longitudes                                  |  nj, ni      |   \-   |
  |tlat                |t-grid latitudes                                   |  nj, ni      |   \-   |
  |tlon                |t-grid latitudes                                   |  nj, ni      |   \-   |
  |aicen               |ice fraction per category                          |  ncat, nj, ni|   \-   |
  |vicen               |ice volume per unit area per category              |  ncat, nj, ni|   \-   |
  |vsnon               |snow volume per unit area per category             |  ncat, nj, ni|   \-   |
  |Tsfcn               |surface snow/ice temperature per category          |  ncat, nj, ni|   C    |
  |iage                |ice age                                            |  ncat, nj, ni|   \-   |
  |FY                  |first year ice area                                |  ncat, nj, ni|   \-   |
  |alvl                |fraction of level ice area                         |  ncat, nj, ni|   \-   |
  |vlvl                |volume of the level ice area                       |  ncat, nj, ni|   \-   |
  |apnd                |fraction of ponds per cell                         |  ncat, nj, ni|   \-   |
  |hpnd                |depth of ponds per cell                            |  ncat, nj, ni|   \-   |
  |ipnd                |other pond characteristics                         |  ncat, nj, ni|   \-   |
  |dhs                 |\-                                                 |  ncat, nj, ni|   \-   |
  |ffrac               |\-                                                 |  ncat, nj, ni|   \-   |
  |fbrs                |\-                                                 |  ncat, nj, ni|   \-   |
  |sice00\[1-5\]       |\-                                                 |  ncat, nj, ni|   \-   |
  |qice00\[1-5\]       |\-                                                 |  ncat, nj, ni|   \-   |
  |qsno00\[1-5\]       |\-                                                 |  ncat, nj, ni|   \-   |

# Climatological Datasets
The following external datasets are relevant to stand-alone *CICE6*. Below is given information that is helpful to reference for this purpose.

## COREv2
[COREv2](https://data1.gfdl.noaa.gov/nomads/forms/core/COREv2/CIAF_v2.html) is a high quality, but now aging, global one degree and one-hour atmospheric climate dataset with from years 1948-2009.
## C-GLORS
The CMCC Global Ocean Physical Reanalysis System ([C-GLORS](http://c-glors.cmcc.it/index/index.html)) is used to simulate the state of the ocean in the last decades. It consists of a variational data assimilation system (OceanVar), capable of assimilating all in-situ observations along with altimetry data, and a forecast step performed by the ocean model NEMO coupled with the LIM2 sea-ice model. 
## ERA5
$1/4^{\circ}$ spatial and 1-hour temporal resolutions. A thorough description of all the fields in ERA5 is given [here](https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview).
[This is also](https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation#ERA5:datadocumentation-Table3), another good location to resource for information on the ERA5 fields.

[COSIMA GitHub Working Group/Issue](https://github.com/COSIMA/access-om2/issues/242). For a comprehensive listing on gadi of mapped ERA5 fields [see Andrew Kiss\'s](https://github.com/COSIMA/access-om2/issues/242#issuecomment-913908910) helpful table -- that thread also contains other helpful for information re. ERA5 setup for CICE.

-   [Latitudes and Longitudes (grid) in text (ascii) format](https://rda.ucar.edu/datasets/ds633.0/docs/ERA5.025deg_lats_lons.txt)
-   [Pressure Levels (levels) in text (ascii) format](https://rda.ucar.edu/datasets/ds633.0/docs/ERA5.025deg_std_pres_levels.txt)
-   [ERA5 Atmospheric Surface Analysis](https://rda.ucar.edu/datasets/ds633.0/docs/ds633.0.e5.oper.an.sfc.grib1.table.web.txt)
-   [ERA5 Atmospheric Pressure Level Analysis](https://rda.ucar.edu/datasets/ds633.0/docs/ds633.0.e5.oper.an.pl.grib1.table.web.txt)
-   [ERA5 Accummulated Fields](https://rda.ucar.edu/datasets/ds633.0/docs/ds633.0.e5.oper.fc.sfc.accumu.grib1.table.web.txt)
-   [ERA5 Instantaneous Fields](https://rda.ucar.edu/datasets/ds633.0/docs/ds633.0.e5.oper.fc.sfc.instan.grib1.table.web.txt)

### The following ERA5 fields have been regridded:

#### 10u: Ten-metre zonal surface winds
#### 10v: Ten-metre meridional surface winds
#### 2d: Two-metre dew-point temperature
#### 2t: Two-metre air temperature
#### metss: Mean eastward turbulent surface stress
#### mntss: Mean northward turbulent surface stress
#### mror: Mean runoff rate
#### msdrswrf: Mean surface direct short-wave radiation flux
#### msdrswrfcs: Mean surface direct short-wave radiation flux, clear sky
#### msdwlwrf: Mean surface downward long-wave radiation flux
#### msdwlwrfcs: Mean surface downward long-wave radiation flux, clear sky
#### msdwswrf: Mean surface downward short-wave radiation flux
#### msdwswrfcs: Mean surface downward short-wave radiation flux, clear sky
#### msl: Mean sea level pressure
#### msnlwrf: Mean surface net long-wave radiation flux
#### msnlwrfcs: Mean surface net long-wave radiation flux, clear sky
#### msnswrf: Mean surface net short-wave radiation flux
#### msnswrfcs: Mean surface net short-wave radiation flux, clear sky
#### msr: Mean snowfall rate
#### msror: Mean surface runoff rate
#### mtdwswrf: Mean top downward short-wave radiation flux
#### mtnlwrf: Mean top net long-wave radiation flux
#### mtnlwrfcs: Mean top net long-wave radiation flux, clear sky
#### mtnswrf: Mean top net short-wave radiation flux
#### mtnswrfcs: Mean top net short-wave radiation flux, clear sky
#### mtpr: Mean total precipitation rate
#### sp: Surface pressure
#### z: Geopotential
## JRA55do
[JRA55do](https://climate.mri-jma.go.jp/pub/ocean/JRA55-do/) is a $1/2^{\circ}$ spatial and 3-hour temporal resolution atmospheric climate dataset. Version 1.5 has been regridded on the following variables:
### huss: Near-surface specific humidity
### prra: Rainfall flux
### prsn: Snowfall flux
### psl: Mean surface pressure
### rlds: Surface downwelling longwave flux
### rsds: Surface downwelling shortwave flux
### tas: Near-surface air temperature
### ts: Surface temperature
### uas: Eastward near-surface wind
### vas: Northward near-surface wind
## MERRAv2
[MERRAv2](https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/) is a $1/2^{\circ}$ spatial and 1-hour temporal resolution atmospheric climate dataset maintained by NASA. gadi has a repository described
[here](http://climate-cms.wikis.unsw.edu.au/MERRA2) and an [excel spreadsheet](file:///Users/dpath2o/PHD/references/data/MERRA2_variables.xlsx) has also been provided to further describe this repository. 2010-2019 are retained locally.

### Surface Flux Diagnostics (\'flx\')
### Radiation Diagnostics (\'rad\')
### Single-level Diagnostics (\'slv\')
### Ocean Surface Diagnostics (\'ocn\')
Not presently on any projects associated with gadi. I have inquired with the one of the project leaders at UTas if there is any scope to \'host\' (for the lack of a better term) this portion of MERRAv2 in gadi project \'ua8\'

## SOSE
The following files were downloaded from SOSE [Iteration 122](http://sose.ucsd.edu/BSOSE6_iter122_solution.html) (2013-2017).

### Bottom pressure
### Fresh water from atmosphere and land
### Mixed layer depth
### Sea surface height anomaly
### Sea ice area
### Sea ice covered part of SIqnet
### Ocean surface freshwater flux
### Sea ice effective ice thickness
### Sea ice surface temperature over sea ice
### Sea ice zonal transport effective thickness
### Sea ice zonal velocity
### Sea ice meridional transport effective thickness
### Sea ice meridional velocity
### Sea ice effective snow thickness
### Sea ice snow precipitation over sea ice
### Zonal wind stress at surface
### Meridional wind stress at surface
### Total salt flux
### Total heat flux
### Atmospheric surface pressure
## BRAN2020
[CSIRO maintained dataset](https://geonetwork.nci.org.au/geonetwork/srv/eng/catalog.search#/metadata/f9372_7752_2015_3718), but available on gadi via the [BRAN data](https://dapds00.nci.org.au/thredds/catalog/gb6/BRAN/BRAN2020/catalog.html) repository.

## CAWCR
The [CICE GitHub online setup guide](https://github.com/CICE-Consortium/CICE/wiki/CICE-Input-Data#descriptions-of-tar-files) explains that ocean surface gravity wave information is necessary for
floe size distribution (FSD). Explanation is given that illustrates from just one day with the only variable present from [WaveWatch III,](https://github.com/NOAA-EMC/WW3) `efreq` ($m^2/s$) the importance of including wave information in sea ice modelling. The [CICE Documentation](https://cice-consortium-cice.readthedocs.io/en/main/#) states the basic problem in sea ice modelling, the one CICE aims to model, is ice thickness distribution (ITD). However, CICE can optionally attempt to model FSD. It does this via a probability distribution function that charcacterises the variability of the horizontal floe size through vertical and lateral growth, new freezing, breaking waves and welding of floes. This work is based on \\[@horvatPrognosticModelSeaice2015]
-   Available to download from [THREDDS](https://data-cbr.csiro.au/thredds/catalog/catch_all/CMAR_CAWCR-Wave_archive/CAWCR_Wave_Hindcast_aggregate/gridded/catalog.html)
# Pre-processing Forcing and Initial Conditions
## Grids
There are two approaches to providing a stand-alone grid to CICE6: using someone elses or creating one from scratch. Methods accounting both for my research with CICE6 stand-alone are given below.
### ACCESS-OM2
[ACCESS-OM2](http://cosima.org.au/index.php/models/access-om2/) has $1^{\circ}$, $1/4^{\circ}$, and $1/10^{\circ}$ grid files that are Arakawa-B grid types. A good description of the varying types of grids is given in [@nGridsNumericalWeather2013]. [CICE6 allows for both Arakawa-B and Arakawa-C](https://cice-consortium-cice.readthedocs.io/en/main/user_guide/ug_implementation.html#grid-boundary-conditions-and-masks) grids.

Regridding of *ERA5* and *BRAN* datasets is required to map the essential forcing variables onto the same GRID that will be used for CICE6 studies herein.
[afim.py](file:///Users/dpath2o/PHD/src/python/afim.py) was written for this purpose. 

1.  ACCESS-OM2 1/4-Degree CICE5

    I downloaded the ACCESS-OM2 grids from gadi and stored them locally.
    There are a lot of grids on gadi and I simply choose the most recent
    grid files after a brief discussion with Paul Spence on the topic.
    The directory listing of the ACCESS-OM2 CICE5 1/4-degree grid is:

    > (afim) ➜ \~ ll
    > /Volumes/ioa03/model~input~/ACCESS-OM2/input~20200530~/cice~025deg~
    > total 510416 -rw-r-xr-- 1 dpath2o staff 119M 26 May 2020 grid.nc
    > lrwxrwxrwx 1 dpath2o staff 39B 26 May 2020 i2o.nc -\>
    > ../../input~20200422~/cice~025deg~/i2o.nc -rw-r-xr-- 1 dpath2o
    > staff 12M 26 May 2020 kmt.nc -rw-r-xr-- 1 dpath2o staff 24M 26 May
    > 2020 monthly~sstsss~.nc -rw-r-xr-- 1 dpath2o staff 83M 26 May 2020
    > o2i.nc -rw-r-xr-- 1 dpath2o staff 12M 26 May 2020 u~star~.nc

    1.  grid description

        The grid file is structured like this:

        > gridtype = generic gridsize = 1555200 xsize = 1440 ysize =
        > 1080 cdo griddes: Processed 10 variables \[0.01s 24MB\].
        >
        > (afim) ➜ \~ cdo sinfon
        > /Volumes/ioa03/model~input~/ACCESS-OM2/input~20200530~/cice~025deg~/grid.nc
        > File format : NetCDF4 -1 : Institut Source T Steptype Levels
        > Num Points Num Dtype : Parameter name 1 : unknown unknown c
        > instant 1 1 1555200 1 F64 : ulat 2 : unknown unknown c instant
        > 1 1 1555200 1 F64 : ulon 3 : unknown unknown c instant 1 1
        > 1555200 1 F64 : tlat 4 : unknown unknown c instant 1 1 1555200
        > 1 F64 : tlon 5 : unknown unknown c instant 1 1 1555200 1 F64 :
        > htn 6 : unknown unknown c instant 1 1 1555200 1 F64 : hte 7 :
        > unknown unknown c instant 1 1 1555200 1 F64 : angle 8 :
        > unknown unknown c instant 1 1 1555200 1 F64 : angleT 9 :
        > unknown unknown c instant 1 1 1555200 1 F64 : tarea 10 :
        > unknown unknown c instant 1 1 1555200 1 F64 : uarea Grid
        > coordinates : 1 : generic : points=1555200 (1440x1080)
        > Vertical coordinates : 1 : surface : levels=1 cdo sinfon:
        > Processed 10 variables \[0.04s 23MB\].
        >
        > netcdf grid { dimensions: nx = 1440 ; ny = 1080 ; nc = 4 ;
        > variables: double ulat(ny, nx) ; ulat:units = \"radians\" ;
        > ulat:title = \"Latitude of U points\" ; double ulon(ny, nx) ;
        > ulon:units = \"radians\" ; ulon:title = \"Longitude of U
        > points\" ; double tlat(ny, nx) ; tlat:units = \"radians\" ;
        > tlat:title = \"Latitude of T points\" ; double tlon(ny, nx) ;
        > tlon:units = \"radians\" ; tlon:title = \"Longitude of T
        > points\" ; double htn(ny, nx) ; htn:units = \"cm\" ; htn:title
        > = \"Width of T cells on North side.\" ; double hte(ny, nx) ;
        > hte:units = \"cm\" ; hte:title = \"Width of T cells on East
        > side.\" ; double angle(ny, nx) ; angle:units = \"radians\" ;
        > angle:title = \"Rotation angle of U cells.\" ; double
        > angleT(ny, nx) ; angleT:units = \"radians\" ; angleT:title =
        > \"Rotation angle of T cells.\" ; double tarea(ny, nx) ;
        > tarea:units = \"m^2^\" ; tarea:title = \"Area of T cells.\" ;
        > double uarea(ny, nx) ; uarea:units = \"m^2^\" ; uarea:title =
        > \"Area of U cells.\" ; }

    2.  re-gridding

        The first step in re-gridding is to pull out the required
        information into two separate files. I used the following python
        code to do this job.

        ``` python
        import xarray as xr
        import numpy as np
        import afim

        G0p25  = xr.open_dataset('/Volumes/ioa03/cice-dirs/input/AFIM/grid/0p25/grid.nc')
        Gt0p25 = xr.Dataset({"tlat"   : G0p25.tlat,
                             "tlon"   : G0p25.tlon,
                             "htn"    : G0p25.htn,
                             "hte"    : G0p25.hte,
                             "angleT" : G0p25.angleT,
                             "tarea"  : G0p25.tarea})
        Gu0p25 = xr.Dataset({"ulat"   : G0p25.ulat,
                             "tlon"   : G0p25.ulon,
                             "angle"  : G0p25.angle,
                             "uarea"  : G0p25.uarea})

        ```

2.  ACCESS-OM2 1/10-degree CICE5

    There is no discernable directory/path difference in the 1/10-degree
    grid from the 1/4-degree grid:

    > -rw-r-xr-- 1 dpath2o staff 1.9G 26 May 2020 grid.nc lrwxrwxrwx 1
    > dpath2o staff 38B 26 May 2020 i2o.nc -\>
    > ../../input~20200422~/cice~01deg~/i2o.nc -rw-r-xr-- 1 dpath2o
    > staff 74M 26 May 2020 kmt.nc -rw-r-xr-- 1 dpath2o staff 148M 26
    > May 2020 monthly~sstsss~.nc -rw-r-xr-- 1 dpath2o staff 519M 26 May
    > 2020 o2i.nc -rw-r-xr-- 1 dpath2o staff 74M 26 May 2020 u~star~.nc

    1.  grid description

        The grid file is slightly different than the 1/4-degree grid in
        that it contains information on the corner of each grid cell and
        hence appears to be an Arakawa-C grid. This may or may not be
        significant.

        > gridtype = generic gridsize = 9720000 xsize = 3600 ysize =
        > 2700 cdo griddes: Processed 14 variables \[8.19s 23MB\].
        >
        > (afim) ➜ \~ cdo sinfon
        > /Volumes/ioa03/model~input~/ACCESS-OM2/input~20200530~/cice~01deg~/grid.nc
        > File format : NetCDF4 -1 : Institut Source T Steptype Levels
        > Num Points Num Dtype : Parameter name 1 : unknown unknown c
        > instant 1 1 9720000 1 F64 : ulat 2 : unknown unknown c instant
        > 1 1 9720000 1 F64 : ulon 3 : unknown unknown c instant 1 1
        > 9720000 1 F64 : tlat 4 : unknown unknown c instant 1 1 9720000
        > 1 F64 : tlon 5 : unknown unknown c instant 4 2 9720000 1 F64 :
        > clon~t~ 6 : unknown unknown c instant 4 2 9720000 1 F64 :
        > clat~t~ 7 : unknown unknown c instant 4 2 9720000 1 F64 :
        > clon~u~ 8 : unknown unknown c instant 4 2 9720000 1 F64 :
        > clat~u~ 9 : unknown unknown c instant 1 1 9720000 1 F64 : htn
        > 10 : unknown unknown c instant 1 1 9720000 1 F64 : hte 11 :
        > unknown unknown c instant 1 1 9720000 1 F64 : angle 12 :
        > unknown unknown c instant 1 1 9720000 1 F64 : angleT 13 :
        > unknown unknown c instant 1 1 9720000 1 F64 : tarea 14 :
        > unknown unknown c instant 1 1 9720000 1 F64 : uarea Grid
        > coordinates : 1 : generic : points=9720000 (3600x2700)
        > Vertical coordinates : 1 : surface : levels=1 2 : generic :
        > levels=4 cdo sinfon: Processed 14 variables \[0.01s 24MB\].
        >
        > netcdf grid { dimensions: nx = 3600 ; ny = 2700 ; nc = 4 ;
        > variables: double ulat(ny, nx) ; ulat:units = \"radians\" ;
        > ulat:title = \"Latitude of U points\" ; double ulon(ny, nx) ;
        > ulon:units = \"radians\" ; ulon:title = \"Longitude of U
        > points\" ; double tlat(ny, nx) ; tlat:units = \"radians\" ;
        > tlat:title = \"Latitude of T points\" ; double tlon(ny, nx) ;
        > tlon:units = \"radians\" ; tlon:title = \"Longitude of T
        > points\" ; double clon~t~(nc, ny, nx) ; clon~t~:units =
        > \"radians\" ; clon~t~:title = \"Longitude of T cell corners\"
        > ; double clat~t~(nc, ny, nx) ; clat~t~:units = \"radians\" ;
        > clat~t~:title = \"Latitude of T cell corners\" ; double
        > clon~u~(nc, ny, nx) ; clon~u~:units = \"radians\" ;
        > clon~u~:title = \"Longitude of U cell corners\" ; double
        > clat~u~(nc, ny, nx) ; clat~u~:units = \"radians\" ;
        > clat~u~:title = \"Latitude of U cell corners\" ; double
        > htn(ny, nx) ; htn:units = \"cm\" ; htn:title = \"Width of T
        > cells on North side.\" ; double hte(ny, nx) ; hte:units =
        > \"cm\" ; hte:title = \"Width of T cells on East side.\" ;
        > double angle(ny, nx) ; angle:units = \"radians\" ; angle:title
        > = \"Rotation angle of U cells.\" ; double angleT(ny, nx) ;
        > angleT:units = \"radians\" ; angleT:title = \"Rotation angle
        > of T cells.\" ; double tarea(ny, nx) ; tarea:units = \"m^2^\"
        > ; tarea:title = \"Area of T cells.\" ; double uarea(ny, nx) ;
        > uarea:units = \"m^2^\" ; uarea:title = \"Area of U cells.\" ;
        > }

    2.  

### Grid Creation

The US [National Snow and Ice Data Center](https://nsidc.org/home) is a
wealth of resources and information and one thing that they have
developed is a [Polar Stereographic
Projection](https://nsidc.org/data/user-resources/help-center/guide-nsidcs-polar-stereographic-projection#anchor-5)
that represents the polar regions with minimal distoration. Thanks to
[Till Soya Rasmussen](mailto:tar@dmi.dk) for turning me on to this
projection. Fortunately there is also a python resource for transforming
\`square\' (Cartesian) latitudes and longitudes into [NSIDC Polar
Stereographic Southern Hemisphere Projection (EPSG
3976)](https://epsg.io/3976).

## Mimicing CICE\'s current options for forcing

In the [CICE6 documentation for the forcing namelist
fields](https://cice-consortium-cice.readthedocs.io/en/main/user_guide/ug_case_settings.html#forcing-nml)
the `atm_data_type` and the `ocn_data_type` appear to be the only
location from which to control the type of atmosphere and ocean forcing
datasets, respectively, that one can provide CICE6. It\'s noted that
neither [*ERA5*]{.spurious-link target="ERA5"} or
[*BRAN*]{.spurious-link target="BRAN"} are offered as supported options.
Hence if I am to get CICE to be forced with [*ERA5*]{.spurious-link
target="ERA5"} and [*BRAN*]{.spurious-link target="BRAN"} then I\'ve got
essentially two options:

1.  write my own Fortran sub-routine that can be used directly by CICE6,
    or
2.  re-configure [*ERA5*]{.spurious-link target="ERA5"} and
    [*BRAN*]{.spurious-link target="BRAN"} to \`look\' like an existing
    supported CICE6 forcing option.

Do to my lack of skill in writing Fortran and my initial thought that
getting [*ERA5*]{.spurious-link target="ERA5"} and
[*BRAN*]{.spurious-link target="BRAN"} to \`look\' like an existing
supported CICE6 forcing option, I chose the second option. This turned
out to be less than trivial.

### Finding which dataset to mimic

Fortunately all the Fortran sub-routines that perform the reading-in of
physical forcing datasets are found in the
[ice~forcing~.F90](file:///Users/dpath2o/src/CICE/cicecore/cicedynB/general/ice_forcing.F90)
file. The sub-routine `JRA55_tx1_files` performs the initial reading-in
of JRA55 files. This seemed like the most logical place to start
mimicing the atmosphere. Hence using the [*TX1 Forcing*]{.spurious-link
target="TX1 Forcing"} (from above) I began in earnest to re-grid two
years worth of [*ERA5*]{.spurious-link target="ERA5"} dataset.
Similarily, for the ocean forcing, `ocn_data_ncar_init_3D` performs the
reading-in of 3D ocean data. The beginning comments in that sub-routine
provide some useful information and just a few lines below those
comments the variable names are given. Using this information I began in
earnest to re-grid two years worth of [*BRAN*]{.spurious-link
target="BRAN"} dataset.

### Creating Python module [afim.py](file:///Users/dpath2o/PHD/src/python/afim.py)

After beginning to do the above re-gridding/mimicing in a Python
notebook it became apparent that the resources used by iPython/Jupyter
were too taxing and hence a script-based (command-line) interface is
best to utilise [afim.py](file:///Users/dpath2o/PHD/src/python/afim.py).
At present, the class `cice_prep` performs the re-gridding of
[*ERA5*]{.spurious-link target="ERA5"} and [*BRAN*]{.spurious-link
target="BRAN"} via a required JSON, here\'s the [current
example](file:///Users/dpath2o/PHD/src/python/afim.json). So performing
a re-gridding is simply done via the following commands, which are [here
in this
example](file:///Users/dpath2o/PHD/src/python/CICE/prepare_cice_forcing_script.py):

``` python
import os
import afim
cp = afim.cice_prep('/Users/dpath2o/PHD/src/python/afim.json')
# do ERA5 regridding and derivations for forcing
cp.regrid_wrapper(ds_name='ERA5')
# do BRAN regridding and derivations for forcing
cp.regrid_wrapper(ds_name='BRAN')
```

1.  Re-gridding with [NetCDF Operators
    (NCO)](~/PHD/admin/computing.org::*NetCDF Operators (NCO))

    [A good
    tutorial](https://acme-climate.atlassian.net/wiki/spaces/DOC/pages/754286611/Regridding+E3SM+Data+with+ncremap)
    explaining how to use
    [NCO](~/PHD/admin/computing.org::*NetCDF Operators (NCO)), and a bit
    closer to home \[\[

    1.  ERA5 re-gridding 1/10-degree

## One-to-One Forcing Discusion

[*As discussed earlier*]{.spurious-link target="CICE Consortium"},
CICE6, \`out-of-the-box\' comes with three example or test runs based on
the underlying spatial grid: a [*three-degree global
cartesian*]{.spurious-link target="GX3 Description"}, a [*one-degree
global cartesian*]{.spurious-link target="GX1 Description"}, and a
[*one-degree global tri-polar*]{.spurious-link
target="TX1 Description"}. [The CICE6
documentation](https://cice-consortium-cice.readthedocs.io/en/main/user_guide/ug_case_settings.html#forcing-nml)
indicates a very narrow range of supported forcing datasets, however,
the documentation does describe what variables are nice-to-haves for
running CICE6 in stand-alone. However, users are cautioned that
stand-alone is not the intended \`scientific\' mode of CICE6, and that
coupling with an ocean and/or atmosphere is preferred for obtaining more
scientific results. From that [discussion in the
documentation](https://cice-consortium-cice.readthedocs.io/en/main/developer_guide/dg_forcing.html),
I created a table for future use in mapping [*ERA5*]{.spurious-link
target="ERA5"} and [*BRAN*]{.spurious-link target="BRAN"} data to in
case it becomes useful/beneficial to force CICE6 writing my own Fortran
sub-routine.

### Table of Stand-alone CICE6 Forcing Fields

  CICE6 short    CICE6 description                                Units         Dataset   Short Name   Long Name                                         Derivation   Notes
  -------------- ------------------------------------------------ ------------- --------- ------------ ------------------------------------------------- ------------ ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  zlvl           atmosphere level height                          m             \-        \-           \-                                                \-           
  uatm           model grid i-direction wind velocity component   m/s           ERA5      10u          10 metre U wind component                         no           
  vatm           model grid j-direction wind velocity component   m/s           ERA5      10v          10 metre V wind component                         no           
  strax          model grid i-direction wind stress               N/m^2^        \-        \-           Mean eastward turbulent surface stress            no           
  stray          model grid j-direction wind stress               N/m^2^        \-        \-           Mean northward turbulent surface stress           no           
  potT           air potential temperature                        K             \-        \-           \-                                                yes          [derivation](https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.potential_temperature.html)
  Tair           air termperature                                 K             ERA5      2t           2 metre temperature                               no           
  Qa             specific humidity                                kg/kg         ERA5      \-           \-                                                yes          [derivation](https://confluence.ecmwf.int/pages/viewpage.action?pageId=171411214)
  rhoa           air density                                      kg/m^3^       ERA5      \-           \-                                                yes          compute [RH](https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.relative_humidity_from_specific_humidity.html) then [mixing ratio](https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.mixing_ratio_from_relative_humidity.html?highlight=mixing%20ratio#metpy.calc.mixing_ratio_from_relative_humidity) then [rho](https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.density.html?highlight=density)
  flw            incoming longwave radiation                      W/m^2^        ERA5      msdwlwrf     Mean surface downward long-wave radiation flux    no           
  fsw            incoming shortwave radiation                     W/m^2^        ERA5      msdwswrf     Mean surface downward short-wave radiation flux   no           
  swvdr          sw down, visible, direct                         W/m^2^        \-        \-           \-                                                no           
  swvdf          sw down, visible, diffuse                        W/m^2^        \-        \-           \-                                                yes          
  swidr          sw down, near IR, direct                         W/m^2^        \-        \-           \-                                                yes          
  swidf          sw down, near IR, diffuse                        W/m^2^        \-        \-           \-                                                yes          
  frain          rainfall rate                                    kg/(m^2^ s)   ERA5      mtpr         Mean total precipitation rate                     no           
  fsnow          snowfall rate                                    kg/(m^2^ s)   \-        \-           \-                                                no           
  uocn           ocean current, x-direction                       m/s           BRAN      u            eastward velocity                                 no           
  vocn           ocean current, y-direction                       m/s           BRAN      v            northward velocity                                no           
  ss~tltx~       sea surface slope, x-direction                   m             BRAN      \-           \-                                                yes          see [eqn:seaslope](eqn:seaslope) ; use `eta`
  ss~tlty~       sea surface slope, y-direction                   m             BRAN      \-           \-                                                yes          see [eqn:seaslope](eqn:seaslope) ; use `eta`
  sss            sea surface salinity                             ppt           BRAN      salt         salinity                                          no           
  sst            sea surface temperature                          C             BRAN      temp         temperature                                       no           
  frzmlt         freezing/melting potential                       W/m^2^        \-        \-           \-                                                             
  frzmlt~init~   frzmlt used in current time step                 W/m^2^        \-        \-           \-                                                             
  Tf             freezing temperature                             C             \-        \-           \-                                                             
  qdp            deep ocean heat flux                             W/m^2^        \-        \-           \-                                                             [eqn:qdp](eqn:qdp)
  hmix           mixed layer depth                                m             BRAN      mld          mixed layer depth                                 no           

1.  Computing Short-wave downward visibile diffuse

    ```{=org}
    #+NAME: eqn:swvdf
    ```
    ```{=latex}
    \begin{equation}
    SW_{diff} = 0.952 - 1.041\exp{-\exp{(2.3-4.702k_t)}}
    \end{equation}
    ```
    where,

    ```{=latex}
    \begin{equation}

    k_t = SW / SW_{TOA}
    \end{equation}
    ```
    and $SW_{TOA}$ is the short-wave at the top-of-the-atmosphere

2.  Computing two-metre wind speed components

    ```{=org}
    #+NAME: eqn:u10to2
    ```
    ```{=latex}
    \begin{equation}
    u2 = (u10 * 4.87) / ln ((67.8 * 10) − 5.42)
    \end{equation}
    ```

3.  Computing sea surface slope

    ```{=org}
    #+NAME: eqn:seaslope
    ```
    ```{=latex}
    \begin{equation}
    ss_tltx = eta/dx
    ss_tlty = eta/dy
    \end{equation}
    ```
    where $dx$ and $dy$ is the grid spacing in meters in the
    longitudinal and latitudinal, respectively

4.  Computing deep ocean heat flux

    To compute the ocean heat flux at the mixed layer depth the
    following reference was used
    [@fofonoffAlgorithmsComputationFundamental1983]

    ```{=org}
    #+NAME: eqn:qdp
    ```
    ```{=latex}
    \begin{equation}
        c0     =  4.2174e3
        c1     = -3.720283
        c2     =   .1412855
        c3     = -2.654387e-3
        c4     = 2.093236e-5
        a0     = -7.64357
        a1     =   .1072763
        a2     = -1.38385e-3
        b0     =   .1770383
        b1     = -4.07718e-3
        b2     =  5.148e-5
        Cpst0  = c0 + c1*T_{mld} + c2*T_{mld}^2 + c3*T_{mld}^3 + c4*T_{mld}^4 + (a0 + a1*T_{mld} + a2*T_{mld}^2)*S_{mld} + (b0 + b1*T_{mld} + b2*T_{mld}^2)*S_{mld}*\sqrt{S_{mld}}
        a0     = -4.9592e-1
        a1     =  1.45747e-2
        a2     = -3.13885e-4
        a3     =  2.0357e-6
        a4     =  1.7168e-8
        b0     =  2.4931e-4
        b1     = -1.08645e-5
        b2     =  2.87533e-7
        b3     = -4.0027e-9
        b4     =  2.2956e-11
        c0     = -5.422e-8
        c1     =  2.6380e-9
        c2     = -6.5637e-11
        c3     =  6.136e-13
        dCp0t0 = (a0 + a1*T_{mld} + a2*T_{mld}^2 + a3*T_{mld}^3 + a4*T_{mld}^4)*P_{mld} + (b0 + b1*T_{mld} + b2*T_{mld}^2 + b3*T_{mld}^3 + b4*T_{mld}^4)*P_{mld}^2 + (c0 + c1*T_{mld} + c2*T_{mld}^2 + c3*T_{mld}^3)*P_{mld}^3
        d0     =  4.9247e-3
        d1     = -1.28315e-4
        d2     =  9.802e-7
        d3     =  2.5941e-8
        d4     = -2.9179e-10
        e0     = -1.2331e-4
        e1     = -1.517e-6
        e2     =  3.122e-8
        f0     = -2.9558e-6
        f1     =  1.17054e-7
        f2     = -2.3905e-9
        f3     =  1.8448e-11
        g0     =  9.971e-8
        h0     =  5.540e-10
        h1     = -1.7682e-11
        h2     =  3.513e-13
        j1     = -1.4300e-12
        S3_2   = S_{mld}*\sqrt(S_{mld})
        dCpstp = ((d0 + d1*T_{mld} + d2*T_{mld}^2 + d3*T_{mld}^3 + d4*T_{mld}^4)*S_{mld} + (e0 + e1*T_{mld} + e2*T_{mld}^2)*S3_2)*P_{mld} + ((f0 + f1*T_{mld} + f2*T_{mld}^2 + f3*T_{mld}^3)*S_{mld} + g0*S3_2)*P_{mld}^2 + ((h0 + h1*T_{mld} + h2*T_{mld}^2)*S_{mld} + j1*T_{mld}*S3_2)*P_{mld}^3
        cp     = Cpst0 + dCp0t0 + dCpstp
        qdp    = cp * (sst-T_{mld})
    \end{equation}
    ```
    where $T_{mld}$, $S_{mld}$ and $P_{mld}$ are the temperature (C),
    salinity (psu) and pressure (dbar) at the mixed layer depth, and
    $sst$ is the sea surface temperature (C).

# External Drives

## ioa01 (5000GB capacity -- 4940GB used)

### Datasets

1.  BRAN (3000GB)

2.  JRA55do (140GB)

3.  MERRAv2 (1800GB -- incomplete)

    -   currently only contains surface fluxes (\'flx\') and radiative
        (\'rad\') components of the dataset
    -   only partial \'rad\' upto mid-2015 (add another 400GB)
    -   also \'ocn\' and \'slv\' (see [*MERRAv2*]{.spurious-link
        target="MERRAv2"} above) need to be downloaded, with estimated
        600GB and 1400GB space required, respectively
    -   estimated total storage space for one decade is 4300GB, hence
        will require own 5TB disk

## ioa02 (5000GB capacity -- 4170GB used)

### Datasets

1.  ERA5 (2300GB)

2.  CAWCR (925GB)

3.  SOSE (675GB -- incomplete)

    missing U and V (both roughly 500GB each)

4.  C-GLORS (7GB)

5.  COREv2 (30GB)

## ioa03 (5000GB capacity -- varying used)

Fast Ice (8GB)

model~input~

model~output~

observations

orography

satellite

## Tracking Data Downloads

  ------------ ----------- -------- ------------ ------------ ------ --------------------------------------------------- ----------------------------------------------------------------------------------------------------------------------------------------- -------------------------------------------------------------------- -------------
  Downloaded   Mandatory   Might    Reanalysis                                                                           Download                                                                                                                                                                                                       
               CICE        Be       or           Field        Time   Time                                                Location                                                                                                                                  Field                                                                
               Field       Useful   Gridded      Short        Step   Coveage                                                                                                                                                                                       Long                                                                 Units
                                    Model        Name                                                                                                                                                                                                              Name                                                                 
  X            x                    ERA5         2d           1h     \[2010-01-01 Fri 00:00\]-\[2019-12-31 Tue 23:00\]   gadi                                                                                                                                      2 metre dewpoint temperature                                         K
  X            x                    ERA5         2t           1h     \[2010-01-01 Fri 00:00\]-\[2019-12-31 Tue 23:00\]   gadi                                                                                                                                      2 metre temperature                                                  K
  x            x                    ERA5         10u          1h     \[2010-01-01 Fri 00:00\]-\[2019-12-31 Tue 23:00\]   gadi                                                                                                                                      10 metre u wind component                                            m/s
  x            x                    ERA5         10v          1h     \[2010-01-01 Fri 00:00\]-\[2019-12-31 Tue 23:00\]   gadi                                                                                                                                      10 metre v wind component                                            m/s
  x            x                    ERA5         metss        1h     \[2010-01-01 Fri 00:00\]-\[2019-12-31 Tue 23:00\]   gadi                                                                                                                                      Mean eastward turbulent surface stress                               N/m^2^
  x            x                    ERA5         mntss        1h     \[2010-01-01 Fri 00:00\]-\[2019-12-31 Tue 23:00\]   gadi                                                                                                                                      Mean northward turbulent surface stress                              N/m^2^
  x                        x        ERA5         mror         1h     \[2010-01-01 Fri 00:00\]-\[2019-12-31 Tue 23:00\]   gadi                                                                                                                                      Mean runoff rate                                                     kg/(m^2^ s)
  x            x                    ERA5         msdrswrf     1h     \[2010-01-01 Fri 00:00\]-\[2019-12-31 Tue 23:00\]   gadi                                                                                                                                      Mean surface direct short-wave radiation flux                        W/m^2^
  x                        x        ERA5         msdrswrfcs   1h     \[2010-01-01 Fri 00:00\]-\[2019-12-31 Tue 23:00\]   gadi                                                                                                                                      Mean surface direct short-wave radiation flux, clear sky             W/m^2^
  x                        x        ERA5         msdwlwrf     1h     \[2010-01-01 Fri 00:00\]-\[2019-12-31 Tue 23:00\]   gadi                                                                                                                                      Mean surface downward long-wave radiation flux                       W/m^2^
  x                        x        ERA5         msdwlwrfcs   1h     \[2010-01-01 Fri 00:00\]-\[2019-12-31 Tue 23:00\]   gadi                                                                                                                                      Mean surface downward long-wave radiation flux, clear sky            W/m^2^
  x                        x        ERA5         msdwswrf     1h     \[2010-01-01 Fri 00:00\]-\[2019-12-31 Tue 23:00\]   gadi                                                                                                                                      Mean surface downward short-wave radiation flux                      W/m^2^
  x                        x        ERA5         msdwswrfcs   1h     \[2010-01-01 Fri 00:00\]-\[2019-12-31 Tue 23:00\]   gadi                                                                                                                                      Mean surface downward short-wave radiation flux, clear sky           W/m^2^
  x            x                    ERA5         msl          1h     \[2010-01-01 Fri 00:00\]-\[2019-12-31 Tue 23:00\]   gadi                                                                                                                                      Mean sea level pressure                                              Pa
  x            x                    ERA5         msnlwrf      1h     \[2010-01-01 Fri 00:00\]-\[2019-12-31 Tue 23:00\]   gadi                                                                                                                                      Mean surface net long-wave radiation flux                            W/m^2^
  x                        x        ERA5         msnlwrfcs    1h     \[2010-01-01 Fri 00:00\]-\[2019-12-31 Tue 23:00\]   gadi                                                                                                                                      Mean surface net long-wave radiation flux, clear sky                 W/m^2^
  x            x                    ERA5         msnswrf      1h     \[2010-01-01 Fri 00:00\]-\[2019-12-31 Tue 23:00\]   gadi                                                                                                                                      Mean surface net short-wave radiation flux                           W/m^2^
  x                        x        ERA5         msnswrfcs    1h     \[2010-01-01 Fri 00:00\]-\[2019-12-31 Tue 23:00\]   gadi                                                                                                                                      Mean surface net short-wave radiation flux, clear sky                W/m^2^
  x            x                    ERA5         msr          1h     \[2010-01-01 Fri 00:00\]-\[2019-12-31 Tue 23:00\]   gadi                                                                                                                                      Mean snowfall rate                                                   kg/(m^2^ s)
  x                        x        ERA5         msror        1h     \[2010-01-01 Fri 00:00\]-\[2019-12-31 Tue 23:00\]   gadi                                                                                                                                      Mean surface runoff rate                                             kg/(m^2^ s)
                           x        ERA5         mtdwswrf     1h     \[2010-01-01 Fri 00:00\]-\[2019-12-31 Tue 23:00\]   gadi                                                                                                                                      Mean top downward short-wave radiation flux                          W/m^2^
  x                        x        ERA5         mtnlwrf      1h     \[2010-01-01 Fri 00:00\]-\[2019-12-31 Tue 23:00\]   gadi                                                                                                                                      Mean top net long-wave radiation flux                                W/m^2^
  x                        x        ERA5         mtnlwrfcs    1h     \[2010-01-01 Fri 00:00\]-\[2019-12-31 Tue 23:00\]   gadi                                                                                                                                      Mean top net long-wave radiation flux, clear sky                     W/m^2^
  x                        x        ERA5         mtnswrf      1h     \[2010-01-01 Fri 00:00\]-\[2019-12-31 Tue 23:00\]   gadi                                                                                                                                      Mean top net short-wave radiation flux                               W/m^2^
  x                        x        ERA5         mtnswrfcs    1h     \[2010-01-01 Fri 00:00\]-\[2019-12-31 Tue 23:00\]   gadi                                                                                                                                      Mean top net short-wave radiation flux, clear sky                    W/m^2^
  x            x                    ERA5         mtpr         1h     \[2010-01-01 Fri 00:00\]-\[2019-12-31 Tue 23:00\]   gadi                                                                                                                                      Mean total precipitation rate                                        kg/(m^2^ s)
  x            x                    ERA5         sp           1h     \[2010-01-01 Fri 00:00\]-\[2019-12-31 Tue 23:00\]   gadi                                                                                                                                      Surface pressure                                                     Pa
  x                        x        CAWCR        dataset      1h     \[2010-01-01 Fri 00:00\]-\[2021-12-31 Fri 23:00\]   [CAWCR](https://data-cbr.csiro.au/thredds/catalog/catch_all/CMAR_CAWCR-Wave_archive/CAWCR_Wave_Hindcast_aggregate/gridded/catalog.html)   complete dataset                                                     various
  x            x                    BRAN         u            1h     \[2010-01-01 Fri 00:00\]-\[2019-12-31 Tue 23:00\]   gadi                                                                                                                                      eastward velocity                                                    m/s
  x            x                    BRAN         v            1h     \[2010-01-01 Fri 00:00\]-\[2019-12-31 Tue 23:00\]   gadi                                                                                                                                      northward velocity                                                   m/s
  x            x                    BRAN         eta          1h     \[2010-01-01 Fri 00:00\]-\[2019-12-31 Tue 23:00\]   gadi                                                                                                                                      surface height on T cells \[Boussinesq (volume conserving) model\]   m
  x            x                    BRAN         salt         1h     \[2010-01-01 Fri 00:00\]-\[2019-12-31 Tue 23:00\]                                                                                                                                                                                                                  
  x            x                    BRAN         temp         1h                                                                                                                                                                                                                                                                        
  x            x                    BRAN         mld          1h                                                                                                                                                                                                   mixed layer depth determined by density criteria                     m
  x                        x        BRAN         ice~force~   1h                                                                                                                                                                                                   contains two variables of short-wave that might be useful            W/m^2^
  ------------ ----------- -------- ------------ ------------ ------ --------------------------------------------------- ----------------------------------------------------------------------------------------------------------------------------------------- -------------------------------------------------------------------- -------------

# Documenting Runs

## Initial Trials

### CICE Installation

[*CICE Github Repository*]{.spurious-link
target="CICE Github Repository"} was cloned onto local machine

### Grid, Forcing and Initial Conditions

Downloaded [*CICE Consortium*]{.spurious-link target="CICE Consortium"}
forcing data, unzipped and copied data into directory structure that
matched [CICE github online setup
guide](https://github.com/CICE-Consortium/CICE/wiki/CICE-Input-Data#directory-structure)

### Test Run \`\`Vanilla\'\'

A successful test run of a \`\`vanilla\'\' setup of CICE. Test directory
setup is preserved in
[case1](file:///Users/dpath2o/PHD/MODELS/src/CICE/case1) and output is
preserved in [case1](file:///Users/dpath2o/cice-dirs/runs/case1). What I
mean by \`\`vanilla\'\' is that I made no changes to any internal files
and simply ran the following (as per the instructions):

``` bash
cd ~/PHD/MODELS/src/CICE/
./cice.setup --case case1 --mach conda --env macos
cd case1
./cice.build
./cice.run
```

## Sandbox

### GX1 fail \[2022-05-23 Mon 06:00\]

Edited `ice_in` file replacing `gx3` enteries with `gx1` . While this
ran through to competion the results were garbarge

### GX1 success \[2022-05-23 Mon 14:00\]

``` {.bash org-language="sh"}
./cice.setup -m conda -e macos -c sandbox_gx1 -g gx1
cd sandbox_gx1
./cice.build
./cice.run
```

![This is the \'sice\' field for GX1
vanilla](./CICE_runs/sandbox_gx1_bin/cice6-gx1_sice.png){#fig:gx1_vanilla}

### GX3 success \[2022-05-23 Mon 14:30\]

``` {.bash org-language="sh"}
./cice.setup -m conda -e macos -c sandbox_gx3 -g gx3
cd sandbox_gx3
./cice.build
./cice.run
```

![This is the \'sice\' field for GX3
vanilla](./CICE_runs/sandbox_gx3_bin/cice6-gx3_sice.png){#fig:gx1_vanilla}

### TX1 (Binary Grid) partial success \[2022-06-06 Mon\]

This run through to completion with the following commands.

``` {.bash org-language="sh"}
./cice.setup -m conda -e macos -c sandbox_tx1 -g tx1
cd sandbox_tx1
./cice.build
./cice.run
```

![This is the \'sice\' field for TX1
vanilla](./CICE_runs/sandbox_tx1_bin/cice6-tx1_sice_BINgrid.png){#fig:tx1_vanilla}

### TX1 (NetCDF Grid -- ACCESS-OM) partial success \[2022-06-06 Mon\]

The grid file can be viewed on the jupyter notebook
[here](../src/py/cice_analysis.ipynb)

``` {.bash org-language="sh"}
./cice.setup -m conda -e macos -c sandbox_tx1_nc -g tx1
cd sandbox_tx1_nc
```

After initial setup I edited
[./CICE_runs/sandbox_tx1_nc/ice_in](./CICE_runs/sandbox_tx1_nc/ice_in)
file changing the following lines

``` text
grid_format  = 'nc'
grid_type    = 'tripole'
grid_file    = '/Users/dpath2o/cice-dirs/input/CICE_data/grid/tx1/grid_tx1.nc'
kmt_file     = 'unknown_kmt_file'
```

This runs through to completion with the following commands but the
results (see [fig:tx1_access-om2](fig:tx1_access-om2), below)

``` {.bash org-language="sh"}
./cice.build
./cice.run
```

![This is the \'sice\' field for ACCESS-OM2
grid](./CICE_runs/sandbox_tx1_nc/cice6-tx1_sice_NCgrid.png){#fig:tx1_access-om2}

### TXp25 (NetCDF Grid) from ACCESS-OM2 forced by ERA5 and BRAN

Here is the cice setup

``` {.bash org-language="sh"}
./cice.setup -m conda -e macos -c sandbox_tx1_nc -g tx1
cd sandbox_tx1_nc
```

After initial setup I edited
[./CICE_runs/sandbox_tx1_nc/ice_in](./CICE_runs/sandbox_tx1_nc/ice_in)
file changing the following lines

``` text
grid_format  = 'nc'
grid_type    = 'tripole'
grid_file    = '/Users/dpath2o/cice-dirs/input/CICE_data/grid/tx1/grid_tx1.nc'
kmt_file     = 'unknown_kmt_file'
```

This runs through to completion with the following commands but the
results (see [fig:tx1_access-om2](fig:tx1_access-om2), below)

``` {.bash org-language="sh"}
./cice.build
./cice.run
```

```{=org}
#+CAPTION: This is the 'cice' field for ACCESS-OM2 grid
```
```{=org}
#+NAME:   fig:tx1_access-om2
```
1.  First considered effort

    \[2022-06-07 Tue 09:00\] See post on [CICE Consortium
    Forum](https://bb.cgd.ucar.edu/cesm/threads/cice6-sandboxing.7412/post-45021)
    for details

2.  CICE setup alternate grid

    Create a new [configuration optional script
    file](file:///Users/dpath2o/src/CICE/configuration/scripts/options/set_nml.afim_0p25)
    that has the name-list parameters that you want start with,
    something like this:

    > dt = 3600.0 runtype = \'initial\' year~init~ = 2010 use~leapyears~
    > = .false. use~restarttime~ = .false. ice~ic~ = \'none\'
    > grid~format~ = \'nc\' grid~type~ = \'tripole\' grid~file~ =
    > \'/Volumes/ioa03/cice-dirs/input/AFIM/grid/0p25/grid.nc\'
    > kmt~file~ =
    > \'/Volumes/ioa03/cice-dirs/input/AFIM/grid/0p25/kmt.nc\'
    > bathymetry~file~ =
    > \'/Volumes/ioa03/cice-dirs/input/AFIM/grid/0p25/topog.nc\'
    > maskhalo~dyn~ = .true. maskhalo~remap~ = .true. maskhalo~bound~ =
    > .true. fyear~init~ = 2010 atm~dataformat~ = \'nc\' atm~datatype~ =
    > \'JRA55~tx1~\' atm~datadir~ =
    > \'*Volumes/ioa03/cice-dirs/input/AFIM/forcing/0p25*\'
    > precip~units~ = \'mm~permonth~\' ocn~datadir~ =
    > \'/Volumes/ioa03/cice-dirs/input/AFIM/forcing/0p25/daily\'
    > bgc~datadir~ = \'\' distribution~wghtfile~ = \'\'

    Then run CICE setup script with this new unique namelist file

    ``` shell
    ./cice.setup -m conda -e macos -c sandbox_gx0p25 -s sandbox_gx0p25
    ```

## AFIM on *gadi*

[NCI support desk
engaged](https://track.nci.org.au/servicedesk/customer/portal/5/HELP-185784)

### CICE files edited before calling cice.setup

1.  ./CICE/configuration/scripts/cice.batch.csh

    These lines were added to the end of the file

    > else if (\${ICE~MACHINE~} =\~ gadi1\*) then cat \>\> \${jobfile}
    > \<\< EOFB #PBS -P jk72 #PBS -l walltime=6:00:00 #PBS -q \${queue}
    > #PBS -l mem=190GB #PBS -l ncpus=48 #PBS -l
    > storage=gdata/cj50+gdata/jk72+gdata/ik11+gdata/hh5+gdata/rt52+gdata/gb6
    > EOFB
    >
    > else if (\${ICE~MACHINE~} =\~ gadi\*) then cat \>\> \${jobfile}
    > \<\< EOFB
    >
    > EOFB

2.  ./CICE/configuration/scripts/cice.launch.csh

    These lines were added to the end of the file

    > \#`=====`{.verbatim} else if (\${ICE~MACHINE~} =\~ gadi\*) then if
    > (\${ICE~COMMDIR~} =\~ serial\*) then cat \>\> \${jobfile} \<\<
    > EOFR ./cice \>&! \\\$ICE~RUNLOGFILE~ EOFR else cat \>\>
    > \${jobfile} \<\< EOFR mpirun -np \${ntasks} ./cice \>&!
    > \\\$ICE~RUNLOGFILE~ EOFR endif

3.  ./CICE/configuration/scripts/options/set~nml~.0p1

    File created (may be edited at any point)

    > dt = 3600.0 runtype = \'initial\' year~init~ = 2005 use~leapyears~
    > = .false. use~restarttime~ = .false. ice~ic~ = \'none\'
    > grid~format~ = \'nc\' grid~type~ = \'tripole\' grid~file~ =
    > \'/g/data/jk72/da1339/cice-dirs/input/AFIM/grid/0p1/grid.nc\'
    > kmt~file~ =
    > \'/g/data/jk72/da1339/cice-dirs/input/AFIM/grid/0p1/kmt.nc\'
    > bathymetry~file~ =
    > \'/g/data/jk72/da1339/cice-dirs/input/AFIM/grid/0p1/topog.nc\'
    > maskhalo~dyn~ = .true. maskhalo~remap~ = .true. maskhalo~bound~ =
    > .true. fyear~init~ = 2005 atm~dataformat~ = \'nc\' atm~datatype~ =
    > \'JRA55~tx1~\' atm~datadir~ =
    > \'*g/data/jk72/da1339/cice-dirs/input/AFIM/forcing/0p1*\'
    > precip~units~ = \'mks\' ocn~datadir~ =
    > \'/g/data/jk72/da1339/cice-dirs/input/AFIM/forcing/0p1/daily\'
    > bgc~datadir~ = \'\' distribution~wghtfile~ = \'\'

4.  ./CICE/configuration/scripts/machines/env.gadi1~intel~

    File created (might be edited at any time)

    > #!/bin/csh
    >
    > set inp = \"undefined\" if (\$#argv == 1) then set inp = \$1 endif
    >
    > module purge module load intel-compiler/2021.7.0 module load
    > openmpi/4.1.3 #module load netcdf/4.9.0p module load netcdf/4.9.0
    >
    > setenv ICE~MACHINEMACHNAME~ gadi setenv ICE~MACHINEMACHINFO~
    > \"gadi at NCI\" setenv ICE~MACHINEENVNAME~ intel setenv
    > ICE~MACHINEENVINFO~ \"ifort 2021.7.0, openmpi/4.1.3 netcdf4.9.0p\"
    > setenv ICE~MACHINEMAKE~ gmake setenv ICE~MACHINEWKDIR~
    > \$HOME/cice-dirs/runs setenv ICE~MACHINEINPUTDATA~
    > \$HOME/cice-dirs/input setenv ICE~MACHINEBASELINE~
    > \$HOME/cice-dirs/baseline setenv ICE~MACHINESUBMIT~ \"qsub \"
    > setenv ICE~MACHINEACCT~ unused setenv ICE~MACHINEQUEUE~ \"normal\"
    > setenv ICE~MACHINETPNODE~ 48 setenv ICE~MACHINEBLDTHRDS~ 1 setenv
    > ICE~MACHINEQSTAT~ \"qstat \"

5.  ./CICE/configuration/scripts/machines/Macros.gadi1~intel~

    File created (might be edited at any time)

    > \#`============================================================================`{.verbatim}
    >
    > #==============================================================================
    >
    > CPP := fpp CPPDEFS := -DFORTRANUNDERSCORE \${ICE~CPPDEFS~} CFLAGS
    > := -c -O2 -fp-model precise
    >
    > FIXEDFLAGS := -132 FREEFLAGS := -FR FFLAGS := -fp-model precise
    > -convert big~endian~ -assume byterecl -ftz -traceback -g
    > FFLAGS~NOOPT~:= -O0
    >
    > ifeq (\$(ICE~BLDDEBUG~), true) FFLAGS += -O0 -g -check uninit
    > -check bounds -check pointers -fpe0 -check noarg~tempcreated~ else
    > FFLAGS += -O2 endif
    >
    > SCC := mpicc SFC := mpif90 MPICC := mpicc MPIFC := mpif90
    >
    > ifeq (\$(ICE~COMMDIR~), mpi) FC := \$(MPIFC) CC := \$(MPICC) else
    > FC := \$(SFC) CC := \$(SCC) endif LD:= \$(FC)
    >
    > INCLDIR := \$(INCLDIR)
    >
    > SLIBS := -lnetcdf -lnetcdff
    >
    > ifeq (\$(ICE~THREADED~), true) LDFLAGS += -qopenmp CFLAGS +=
    > -qopenmp FFLAGS += -qopenmp endif

### `./cice.setup`

After editing the files in the previous section one can issue the
following command on *gadi*

> ./cice.setup -m gadi1 -c afim~test~ -p 48x1

### `./cice.build`

After running `./cice.setup` and prior to running `./cice.build` the
following should occur. Edit `./cice.build` , `./cice.run` , and
`./cice.settings` deleting `-f` in the first line -- i.e. remove from
`#!/bin/csh -f` to be `#!/bin/csh` in each of those files.

Then simply run `./cice.build`

[NOTE]{.underline}: any changes to `./cice.settings` need to be made
prior to running `./cice_build`

### `./cice.submit` or `./cice.run`

### tests with `afim64`

### separating out ACOM2 and BRAN

1.  \[2023-01-26 Thu\]

# compare era5 and JRA55do air temperatures

# be looking at the tendency terms in the output

# monthy ocean output interpolated to daily -- this to move away from JRA55do climatology \"imprinted\" on the ocean

# re-connect with Andrea about ground iceberg code for implementing in CICE6 (could be in github repo)

# re-grid JRA55do atmosphere to 1/10-degree
