
# Table of Contents

1.  [CICE Consortium](#org28fbdde)
    1.  [CICE Github Repository](#org5179e48)
        1.  [Resources Index](#org6d24216)
        2.  [CICE Consortium Test Data](#orgf8b290e)
    2.  [GX3 Description](#orge8e3b2d)
        1.  [GX3 Grid](#org2b3ded1)
        2.  [GX3 Bathymetry](#orgbf2fbdf)
        3.  [GX3 Land Mask](#org7393b2e)
        4.  [GX3 Forcing](#org08df58e)
        5.  [GX3 Initial Conditions](#orgd5aa7e0)
    3.  [GX1 Description](#org05ab192)
        1.  [GX1 Grid](#orgd404fe0)
        2.  [GX1 Bathymetry](#org31c07a2)
        3.  [GX1 Land Mask](#org2268e1c)
        4.  [GX1 Forcing](#org2a6fb4c)
        5.  [GX1 Initial Conditions](#orgfd2e842)
    4.  [TX1 Description](#org79cf8e6)
        1.  [TX1 Grid](#org6421b79)
        2.  [TX1 Bathymetry](#orgc70422f)
        3.  [TX1 Land Mask](#org1ae6ae3)
        4.  [TX1 Forcing](#org0eecc9b)
        5.  [TX1 Initial Conditions](#org9249e4b)
    5.  [Initial Conditions Discussion](#org2e0c1f6)
2.  [Climatological Datasets](#orge2dbf35)
    1.  [COREv2](#orgef674bc)
    2.  [C-GLORS](#orgc8939a8)
    3.  [ERA5](#org150c1e8)
        1.  [10u](#org25f314e)
        2.  [10v](#org703b36c)
        3.  [2d](#org6661811)
        4.  [2t](#org493eaea)
        5.  [metss](#org54397e7)
        6.  [mntss](#org5926abd)
        7.  [mror](#orgd5b196f)
        8.  [msdrswrf](#orgd6a9b35)
        9.  [msdrswrfcs](#orgdda916b)
        10. [msdwlwrf](#org4a5a706)
        11. [msdwlwrfcs](#orgae47d0c)
        12. [msdwswrf](#org939fd6f)
        13. [msdwswrfcs](#org7099012)
        14. [msl](#org36ae245)
        15. [msnlwrf](#org69dda3d)
        16. [msnlwrfcs](#org9e1095b)
        17. [msnswrf](#orgbbbc293)
        18. [msnswrfcs](#org679bb35)
        19. [msr](#org3984d1b)
        20. [msror](#org2f58934)
        21. [mtdwswrf](#orgb3f1f7b)
        22. [mtnlwrf](#orgffb5609)
        23. [mtnlwrfcs](#org8091c1d)
        24. [mtnswrf](#org62bfd96)
        25. [mtnswrfcs](#org61ea099)
        26. [mtpr](#org703fad6)
        27. [sp](#org99b2dfd)
        28. [z](#org64aa233)
    4.  [JRA55do](#org60ab42d)
        1.  [huss](#org9042284)
        2.  [prra](#org41569e1)
        3.  [prsn](#org58abefe)
        4.  [psl](#orgc49e63f)
        5.  [rlds](#orgc84a418)
        6.  [rsds](#orgf231b7e)
        7.  [tas](#org4792b1c)
        8.  [ts](#orgd4d04f7)
        9.  [uas](#org8f639e9)
        10. [vas](#orgbdb0273)
    5.  [MERRAv2](#org6bedea4)
        1.  [Surface Flux Diagnostics ('flx')](#org0abc25e)
        2.  [Ocean Surface Diagnostics ('ocn')](#org03ca400)
        3.  [Radiation Diagnostics ('rad')](#org5063a13)
        4.  [Single-level Diagnostics ('slv')](#org1d519b7)
    6.  [SOSE](#orga65f4b0)
        1.  [Bottom pressure](#org9145a77)
        2.  [Fresh water from atmosphere and land](#org9e5ccb7)
        3.  [Mixed layer depth](#org2a85212)
        4.  [Sea surface height anomaly](#org4281540)
        5.  [Sea ice area](#org5164d2b)
        6.  [Sea ice covered part of SIqnet](#org40cb11a)
        7.  [Ocean surface freshwater flux](#org89f3cbe)
        8.  [Sea ice effective ice thickness](#org2126178)
        9.  [Sea ice surface temperature over sea ice](#org3fd9cb6)
        10. [Sea ice zonal transport effective thickness](#orgec146f6)
        11. [Sea ice zonal velocity](#orgf8a0379)
        12. [Sea ice meridional transport effective thickness](#orgedaedf7)
        13. [Sea ice meridional velocity](#orgc0d5441)
        14. [Sea ice effective snow thickness](#orgb39e948)
        15. [Sea ice snow precipitation over sea ice](#org98336b2)
        16. [Zonal wind stress at surface](#org6a6e72b)
        17. [Meridional wind stress at surface](#orge50e9a6)
        18. [Total salt flux](#orgd114522)
        19. [Total heat flux](#orgcad33e4)
        20. [Atmospheric surface pressure](#org68eff42)
    7.  [BRAN](#orgea1d349)
    8.  [CAWCR](#orgea7c781)
3.  [Pre-processing Forcing and Initial Conditions](#org057ba5d)
    1.  [Mimicing CICE's current options for forcing](#orgcfbc076)
        1.  [Finding which dataset to mimic](#org1b51396)
        2.  [Creating Python module afim.py](#org8a1d213)
    2.  [One-to-One Forcing Discusion](#org08accaa)
        1.  [Table of Stand-alone CICE6 Forcing Fields](#orga0f92be)
4.  [External Drives](#org373a9d3)
    1.  [ioa01 (5000GB capacity &#x2013; 4940GB used)](#org6c2b476)
        1.  [Datasets](#org4552813)
    2.  [ioa02 (5000GB capacity &#x2013; 4170GB used)](#org0a77467)
        1.  [Datasets](#org8c6f352)
    3.  [ioa03 (5000GB capacity &#x2013; varying used)](#orgc15bedc)
    4.  [Tracking Data Downloads](#orgfd9826b)
5.  [Documenting Runs](#org29b629f)
    1.  [Initial Trials](#org8fe0643)
        1.  [CICE Installation](#orgfe74489)
        2.  [Grid, Forcing and Initial Conditions](#org49c392a)
        3.  [Test Run \`\`Vanilla''](#org389d754)
    2.  [Sandbox](#org8150af8)
        1.  [GX1 fail <span class="timestamp-wrapper"><span class="timestamp">[2022-05-23 Mon 06:00]</span></span>](#orgd05f616)
        2.  [GX1 success <span class="timestamp-wrapper"><span class="timestamp">[2022-05-23 Mon 14:00]</span></span>](#org1e539ea)
        3.  [GX3 success <span class="timestamp-wrapper"><span class="timestamp">[2022-05-23 Mon 14:30]</span></span>](#orgc03d9b5)
        4.  [TX1 (Binary Grid) partial success <span class="timestamp-wrapper"><span class="timestamp">[2022-06-06 Mon]</span></span>](#org40b9209)
        5.  [TX1 (NetCDF Grid &#x2013; ACCESS-OM) partial success <span class="timestamp-wrapper"><span class="timestamp">[2022-06-06 Mon]</span></span>](#org21a9a22)
        6.  [TXp25 (NetCDF Grid) from ACCESS-OM2 forced by ERA5 and BRAN](#org5f83a3b)
    3.  [AFIM stand-alone](#org5eb5303)
6.  [Python Notebook and Scripts](#org6fc8564)
    1.  [AFIM CICE Analysis](#org1f7b356)



<a id="org28fbdde"></a>

# CICE Consortium


<a id="org5179e48"></a>

## CICE Github Repository

<https://github.com/CICE-Consortium>


<a id="org6d24216"></a>

### Resources Index

<https://github.com/CICE-Consortium/About-Us/wiki/Resource-Index>


<a id="orgf8b290e"></a>

### CICE Consortium Test Data

The *default* forcing data can be obtained here:
<https://github.com/CICE-Consortium/CICE/wiki/CICE-Input-Data>

[Standalone Forcing](https://cice-consortium-cice.readthedocs.io/en/master/developer_guide/dg_forcing.html) is described in detail and will be helpful when
pulling in more forcing data than the *default*. The following two lists
show the fields that are required for forcing. The documentation
states that all the variables have reasonable static defaults that
will be used in *default* mode. To advance the forcing, the subroutines
`get_forcing_atmo` and `get_forcing_ocn` are called each timestep from the
step loop. That subroutine computes the forcing year (`fyear`), calls
the appropriate forcing data method, and then calls `prepare_forcing`
which converts the input data fields to model forcing fields.

Data was downloaded and unpacked into the following local directory,
<file:///Volumes/ioa03/cice-dirs/input/CICE_data>.


<a id="orge8e3b2d"></a>

## GX3 Description

Global Three-Degree Cartesian


<a id="org2b3ded1"></a>

### GX3 Grid

Header from the [GX3 grid file](file:///Volumes/ioa03/cice-dirs/input/CICE_data/grid/gx3/grid_gx3.nc):

> netcdf grid<sub>gx3</sub> {
> dimensions:
> 	y = 116 ;
> 	x = 100 ;
> variables:
> 	double ulat(y, x) ;
> 		ulat:units = "rad" ;
> 	double ulon(y, x) ;
> 		ulon:units = "rad" ;
> 	double htn(y, x) ;
> 		htn:units = "cm" ;
> 	double hte(y, x) ;
> 		hte:units = "cm" ;
> 	double hus(y, x) ;
> 		hus:units = "cm" ;
> 	double huw(y, x) ;
> 		huw:units = "cm" ;
> 	double angle(y, x) ;
> 		angle:units = "rad" ;
> 
> // global attributes:
> 		:history = "Tue Apr 17 16:11:59 2007: ncatted -astandard<sub>name,,d</sub>,, global<sub>gx3.grid.nc\n</sub>",
> 			"Tue Apr 17 16:11:38 2007: ncatted -aunits,huw,c,c,cm global<sub>gx3.grid.nc\n</sub>",
> 			"Tue Apr 17 16:11:34 2007: ncatted -aunits,hus,c,c,cm global<sub>gx3.grid.nc\n</sub>",
> 			"Tue Apr 17 16:11:30 2007: ncatted -aunits,hte,c,c,cm global<sub>gx3.grid.nc\n</sub>",
> 			"Tue Apr 17 16:11:27 2007: ncatted -aunits,htn,c,c,cm global<sub>gx3.grid.nc\n</sub>",
> 			"Tue Apr 17 16:09:35 2007: ncatted -aunits,ulat,o,c,rad global<sub>gx3.grid.nc\n</sub>",
> 			"Tue Apr 17 16:09:32 2007: ncatted -aunits,ulon,o,c,rad global<sub>gx3.grid.nc\n</sub>",
> 			"Tue Apr 17 16:09:19 2007: ncatted -aunits,angle,c,c,rad global<sub>gx3.grid.nc\n</sub>",
> 			"Tue Apr 17 16:09:04 2007: ncatted -aunits,ulon,c,c,rad global<sub>gx3.grid.nc\n</sub>",
> 			"Tue Apr 17 16:08:59 2007: ncatted -aunits,ulat,c,c,rad global<sub>gx3.grid.nc\n</sub>",
> 			"Tue Apr 17 16:06:06 2007: ncatted -aunits,ulon,c,c,degree<sub>east</sub> global<sub>gx3.grid.nc\n</sub>",
> 			"Tue Apr 17 16:05:53 2007: ncatted -aunits,ulat,c,c,degree<sub>north</sub> global<sub>gx3.grid.nc\n</sub>",
> 			"Tue Apr 17 16:05:15 2007: ncrename -aunits,standard<sub>name</sub> global<sub>gx3.grid.nc\n</sub>",
> 			"Tue Apr 17 16:04:40 2007: ncatted -aunits,ulon,c,c,longitude global<sub>gx3.grid.nc\n</sub>",
> 			"Tue Apr 17 16:04:31 2007: ncatted -aunits,ulat,c,c,latitude global<sub>gx3.grid.nc\n</sub>",
> 			"Tue Apr 17 16:02:44 2007: ncatted -aConventions,global,c,c,CF-1.0 global<sub>gx3.grid.nc\n</sub>",
> 			"Tue Apr 17 15:55:18 2007: ncap -O -sulat=ulat\*0.0174533;ulon=ulon\*0.0174533;htn=htn\*100.0;hte=hte\*100.0;hus=hus\*100.0;huw=huw\*100.0;angle=angle\*0.0174533 global<sub>gx3.grid.nc</sub> temp.nc\n",
> 			"Tue Apr 17 15:34:39 2007: ncrename -vunspecified,ulat -vunspecified<sub>2,ulon</sub> -vunspecified<sub>3,htn</sub> -vunspecified<sub>4,hte</sub> -vunspecified<sub>5,hus</sub> -vunspecified<sub>6,huw</sub> -vunspecified<sub>7,angle</sub> temp.nc\n",
> 			"Tue Apr 17 15:33:13 2007: ncks -x -a -vlongitude,latitude,unspecified<sub>1,t</sub> global<sub>gx3.grid.nc</sub> temp.nc\n",
> 			"Tue Apr 17 15:32:26 2007: ncrename -dlongitude,x -dlatitude,y global<sub>gx3.grid.nc\n</sub>",
> 			"Tue Apr 17 15:31:50 2007: ncwa -O -aunspecified<sub>1</sub> global<sub>gx3.grid.nc</sub> global<sub>gx3.grid.nc\n</sub>",
> 			"Tue Apr 17 15:31:34 2007: ncwa -at global<sub>gx3.grid.nc</sub> global<sub>gx3.grid.nc\n</sub>",
> 			"Tue Apr 17 15:31:16 2007: ncatted -a,,d,, global<sub>gx3.grid.nc\n</sub>",
> 			"Tue Apr 17 15:30:50 BST 2007 - XCONV V1.91 Development" ;
> 		:Conventions = "CF-1.0" ;
> }


<a id="orgbf2fbdf"></a>

### GX3 Bathymetry

Header from the [GX3 bathymetry file](file:///Volumes/ioa03/cice-dirs/input/CICE_data/grid/gx3/global_gx3.bathy.nc):

> netcdf global<sub>gx3.bathy</sub> {
> dimensions:
> 	nj = 116 ;
> 	ni = 100 ;
> variables:
> 	float Bathymetry(nj, ni) ;
> 		Bathymetry:long<sub>name</sub> = "ocean bathymetry for grounding scheme" ;
> 		Bathymetry:units = "m" ;
> 		Bathymetry:coordinates = "TLON TLAT" ;
> 		Bathymetry:missing<sub>value</sub> = 1.e+30f ;
> 		Bathymetry:<sub>FillValue</sub> = 1.e+30f ;
> 	float TLAT(nj, ni) ;
> 		TLAT:long<sub>name</sub> = "T grid center latitude" ;
> 		TLAT:units = "degrees<sub>north</sub>" ;
> 		TLAT:missing<sub>value</sub> = 1.e+30f ;
> 		TLAT:<sub>FillValue</sub> = 1.e+30f ;
> 	float TLON(nj, ni) ;
> 		TLON:long<sub>name</sub> = "T grid center longitude" ;
> 		TLON:units = "degrees<sub>east</sub>" ;
> 		TLON:missing<sub>value</sub> = 1.e+30f ;
> 		TLON:<sub>FillValue</sub> = 1.e+30f ;
> 
> // global attributes:
> 		:title = "ocean bathymetry for grounding scheme" ;
> 		:source = "created from bathy<sub>ORCA025</sub><sub>LIM.std</sub>" ;
> 		:history = "Mon Oct 29 15:57:35 2018: ncatted -a conventions,global,d,, gx3<sub>bathy.nc</sub> gx3<sub>bathyTP12.nc\n</sub>",
> 			"Mon Oct 29 15:56:16 2018: ncatted -a source,global,o,c,created from bathy<sub>ORCA025</sub><sub>LIM.std</sub> gx3<sub>bathyTP10.nc</sub> gx3<sub>bathyTP11.nc\n</sub>",
> 			"Mon Oct 29 15:54:19 2018: ncatted -a contents,global,d,, gx3<sub>bathyTP9.nc</sub> gx3<sub>bathyTP10.nc\n</sub>",
> 			"Mon Oct 29 15:53:44 2018: ncatted -a comment,global,d,, gx3<sub>bathyTP8.nc</sub> gx3<sub>bathyTP9.nc\n</sub>",
> 			"Mon Oct 29 15:53:36 2018: ncatted -a comment2,global,d,, gx3<sub>bathyTP7.nc</sub> gx3<sub>bathyTP8.nc\n</sub>",
> 			"Mon Oct 29 15:53:10 2018: ncatted -a comment3,global,d,, gx3<sub>bathyTP6.nc</sub> gx3<sub>bathyTP7.nc\n</sub>",
> 			"Mon Oct 29 15:27:25 2018: ncatted -a title,global,o,c,ocean bathymetry for grounding scheme gx3<sub>bathyTP5.nc</sub> gx3<sub>bathyTP6.nc\n</sub>",
> 			"Mon Oct 29 15:11:32 2018: ncks -x -v time,time<sub>bounds</sub> gx3<sub>bathyTP4.nc</sub> gx3<sub>bathyTP5.nc\n</sub>",
> 			"Mon Oct 29 15:10:09 2018: ncatted -a units,Bathymetry,o,c,m gx3<sub>bathyTP3.nc</sub> gx3<sub>bathyTP4.nc\n</sub>",
> 			"Mon Oct 29 15:06:47 2018: ncatted -a long<sub>name,Bathymetry,o,c,ocean</sub> bathymetry for grounding scheme gx3<sub>bathyTP2.nc</sub> gx3<sub>bathyTP3.nc\n</sub>",
> 			"Mon Oct 29 14:57:31 2018: ncks -x -v time<sub>bounds</sub> gx3<sub>bathyTP.nc</sub> gx3<sub>bathyTP2.nc\n</sub>",
> 			"Mon Oct 29 14:57:06 2018: ncks -x -v time gx3<sub>bathyTP.nc</sub> gx3<sub>bathyTP2.nc\n</sub>",
> 			"Mon Oct 29 14:55:17 2018: ncks -x -v ULAT,ULON gx3<sub>bathyTP.nc</sub> gx3<sub>bathyTP2.nc\n</sub>",
> 			"Thu Oct 25 18:26:01 2018: ncrename -v tarea,Bathymetry tempgx3.nc\n",
> 			"Thu Oct 25 17:26:26 2018: ncks -x -v hi tempgx3.nc tempgx3out.nc\n",
> 			"This dataset was created on 2018-10-25 at 17:24:50.3" ;
> 		:io<sub>flavor</sub> = "io<sub>netcdf</sub>" ;
> 		:NCO = "4.4.2" ;
> }


<a id="org7393b2e"></a>

### GX3 Land Mask

Header from the [GX3 land mask file](file:///Volumes/ioa03/cice-dirs/input/CICE_data/grid/gx3/kmt_gx3.nc):

> netcdf kmt<sub>gx3</sub> {
> dimensions:
> 	y = 116 ;
> 	x = 100 ;
> variables:
> 	float kmt(y, x) ;
> 		kmt:units = "1" ;
> 
> // global attributes:
> 		:history = "Tue Apr 17 16:12:58 2007: ncatted -astandard<sub>name,,d</sub>,, global<sub>gx3.kmt.nc\n</sub>",
> 			"Tue Apr 17 16:01:24 2007: ncatted -astandard<sub>name,,c,c,model</sub><sub>level</sub><sub>number</sub> global<sub>gx3.kmt.nc\n</sub>",
> 			"Tue Apr 17 16:00:39 2007: ncatted -aunits,,c,c,1 global<sub>gx3.kmt.nc\n</sub>",
> 			"Tue Apr 17 15:58:38 2007: ncatted -aConventions,global,c,c,CF-1.0 global<sub>gx3.kmt.nc\n</sub>",
> 			"Tue Apr 17 15:27:57 2007: ncrename -dlongitude,x -dlatitude,y temp.nc\n",
> 			"Tue Apr 17 15:27:39 2007: ncwa -aunspecified<sub>1</sub> temp.nc temp.nc\n",
> 			"Tue Apr 17 15:27:32 2007: ncwa -at temp.nc temp.nc\n",
> 			"Tue Apr 17 15:26:58 2007: ncrename -vunspecified,kmt temp.nc\n",
> 			"Tue Apr 17 15:26:21 2007: ncks -C -vunspecified global<sub>gx3.kmt.nc</sub> temp.nc\n",
> 			"Tue Apr 17 15:25:37 2007: ncatted -a,,d,, global<sub>gx3.kmt.nc\n</sub>",
> 			"Tue Apr 17 15:25:13 BST 2007 - XCONV V1.91 Development" ;
> 		:Conventions = "CF-1.0" ;
> }


<a id="org08df58e"></a>

### GX3 Forcing

Three datasets are used for forcing: [1.2.4.1](#org12fd10a), [1.2.4.2](#orgbab5221), and [1.2.4.3](#org5f87a81)

1.  GX3 JRA55

    Separated into five yearly NetCDF files that have three-hourly data &#x2013; i.e.
    [8XDAILY](file:///Volumes/ioa03/cice-dirs/input/CICE_data/forcing/gx3/JRA55/8XDAILY).
    
    Dictories contents:
    
    > -rw-r-xr&#x2013;  1 dpath2o  staff   905M 22 Feb  2020 JRA55<sub>gx3</sub><sub>03hr</sub><sub>forcing</sub><sub>2005.nc</sub>
    > -rw-r-xr&#x2013;  1 dpath2o  staff   905M 22 Feb  2020 JRA55<sub>gx3</sub><sub>03hr</sub><sub>forcing</sub><sub>2006.nc</sub>
    > -rw-r-xr&#x2013;  1 dpath2o  staff   905M 22 Feb  2020 JRA55<sub>gx3</sub><sub>03hr</sub><sub>forcing</sub><sub>2007.nc</sub>
    > -rw-r-xr&#x2013;  1 dpath2o  staff   907M 22 Feb  2020 JRA55<sub>gx3</sub><sub>03hr</sub><sub>forcing</sub><sub>2008.nc</sub>
    > -rw-r-xr&#x2013;  1 dpath2o  staff   905M 22 Feb  2020 JRA55<sub>gx3</sub><sub>03hr</sub><sub>forcing</sub><sub>2009.nc</sub>
    
    Here is the header information of the first file
    
    > netcdf JRA55<sub>gx3</sub><sub>03hr</sub><sub>forcing</sub><sub>2005</sub> {
    > dimensions:
    > 	time = UNLIMITED ; // (2920 currently)
    > 	nj = 116 ;
    > 	ni = 100 ;
    > variables:
    > 	float time(time) ;
    > 		time:long<sub>name</sub> = "model time" ;
    > 		time:units = "days since 1900-12-31 00:00:00" ;
    > 		time:calendar = "standard" ;
    > 	float LON(nj, ni) ;
    > 		LON:long<sub>name</sub> = "Longitude" ;
    > 		LON:units = "degrees<sub>east</sub>" ;
    > 	float LAT(nj, ni) ;
    > 		LAT:long<sub>name</sub> = "Latitude" ;
    > 		LAT:units = "degrees<sub>north</sub>" ;
    > 	float glbrad(time, nj, ni) ;
    > 		glbrad:long<sub>name</sub> = "downward surface shortwave" ;
    > 		glbrad:units = "W/m<sup>2</sup>" ;
    > 		glbrad:coordinates = "LON LAT" ;
    > 		glbrad:<sub>FillValue</sub> = 1.e+30f ;
    > 	float dlwsfc(time, nj, ni) ;
    > 		dlwsfc:long<sub>name</sub> = "downward surface longwave" ;
    > 		dlwsfc:units = "W/m<sup>2</sup>" ;
    > 		dlwsfc:coordinates = "LON LAT" ;
    > 		dlwsfc:<sub>FillValue</sub> = 1.e+30f ;
    > 	float wndewd(time, nj, ni) ;
    > 		wndewd:long<sub>name</sub> = "x-ward winds" ;
    > 		wndewd:units = "m/s" ;
    > 		wndewd:coordinates = "LON LAT" ;
    > 		wndewd:<sub>FillValue</sub> = 1.e+30f ;
    > 	float wndnwd(time, nj, ni) ;
    > 		wndnwd:long<sub>name</sub> = "y-ward winds" ;
    > 		wndnwd:units = "m/s" ;
    > 		wndnwd:coordinates = "LON LAT" ;
    > 		wndnwd:<sub>FillValue</sub> = 1.e+30f ;
    > 	float airtmp(time, nj, ni) ;
    > 		airtmp:long<sub>name</sub> = "Air temperature" ;
    > 		airtmp:units = "Kelvin" ;
    > 		airtmp:coordinates = "LON LAT" ;
    > 		airtmp:<sub>FillValue</sub> = 1.e+30f ;
    > 	float spchmd(time, nj, ni) ;
    > 		spchmd:long<sub>name</sub> = "Specific Humidity" ;
    > 		spchmd:units = "kg/kg" ;
    > 		spchmd:coordinates = "LON LAT" ;
    > 		spchmd:<sub>FillValue</sub> = 1.e+30f ;
    > 	float ttlpcp(time, nj, ni) ;
    > 		ttlpcp:long<sub>name</sub> = "Precipitation" ;
    > 		ttlpcp:units = "kg m-2 s-1" ;
    > 		ttlpcp:coordinates = "LON LAT" ;
    > 		ttlpcp:<sub>FillValue</sub> = 1.e+30f ;
    > }

2.  GX3 NCAR

    Separate into [4XDAILY](file:///Volumes/ioa03/cice-dirs/input/CICE_data/forcing/gx3/NCAR_bulk/4XDAILY) (i.e. six-hourly) and [monthly](file:///Volumes/ioa03/cice-dirs/input/CICE_data/forcing/gx3/NCAR_bulk/MONTHLY) files, which are all binary,
    and hence NIL header information is listed here. The contents of the files can
    be dedecued from the file names:
    
    <table border="2" cellspacing="0" cellpadding="6" rules="groups" frame="hsides">
    
    
    <colgroup>
    <col  class="org-left" />
    
    <col  class="org-left" />
    </colgroup>
    <thead>
    <tr>
    <th scope="col" class="org-left">filename</th>
    <th scope="col" class="org-left">my guess</th>
    </tr>
    </thead>
    
    <tbody>
    <tr>
    <td class="org-left">cldf</td>
    <td class="org-left">cloud</td>
    </tr>
    
    
    <tr>
    <td class="org-left">prec</td>
    <td class="org-left">precipitation</td>
    </tr>
    
    
    <tr>
    <td class="org-left">swdn</td>
    <td class="org-left">net downward shortwave</td>
    </tr>
    
    
    <tr>
    <td class="org-left">dn10</td>
    <td class="org-left">dew-point temperature @ 10 m</td>
    </tr>
    
    
    <tr>
    <td class="org-left">q<sub>10</sub></td>
    <td class="org-left">specific humidity @ 10 m</td>
    </tr>
    
    
    <tr>
    <td class="org-left">t<sub>10</sub></td>
    <td class="org-left">air temperature @ 10 m</td>
    </tr>
    
    
    <tr>
    <td class="org-left">u<sub>10</sub></td>
    <td class="org-left">zonal wind speed @ 10 m</td>
    </tr>
    
    
    <tr>
    <td class="org-left">v<sub>10</sub></td>
    <td class="org-left">meridional wind speed @ 10 m</td>
    </tr>
    </tbody>
    </table>
    
    The directory contents are provided as:
    
    1.  GX3 NCAR 4XDAILY
    
        \#+BEGIN<sub>QUOTE</sub>:
        -rw-r&#x2013;r&#x2013;  1 dpath2o  staff   1.1M 12 Dec  2003 cldf.1997.dat
        -rw-r&#x2013;r&#x2013;  1 dpath2o  staff   1.1M 12 Dec  2003 prec.1997.dat
        -rw-r&#x2013;r&#x2013;  1 dpath2o  staff   1.1M 12 Dec  2003 swdn.1997.dat
        \#+END<sub>QUOTE</sub>
    
    2.  GX3 NCAR MONTHLY
    
        > -rw-r&#x2013;r&#x2013;  1 dpath2o  staff   129M 12 Dec  2003 dn10.1997.dat
        > -rw-r&#x2013;r&#x2013;  1 dpath2o  staff   129M 12 Dec  2003 q<sub>10.1997.dat</sub>
        > -rw-r&#x2013;r&#x2013;  1 dpath2o  staff   129M 12 Dec  2003 t<sub>10.1997.dat</sub>
        > -rw-r&#x2013;r&#x2013;  1 dpath2o  staff   129M 12 Dec  2003 u<sub>10.1997.dat</sub>
        > -rw-r&#x2013;r&#x2013;  1 dpath2o  staff   129M 12 Dec  2003 v<sub>10.1997.dat</sub>

3.  GX3 WW3

    [One wave spectral file](file:///Volumes/ioa03/cice-dirs/input/CICE_data/forcing/gx3/WW3/ww3.20100101_efreq_remapgx3.nc) is given with a total of four time steps and 25
    frequencies. The header is provided as:
    
    > 	time = UNLIMITED ; // (4 currently)
    > 	ni = 100 ;
    > 	nj = 116 ;
    > 	f = 25 ;
    > variables:
    > 	double time(time) ;
    > 		time:standard<sub>name</sub> = "time" ;
    > 		time:long<sub>name</sub> = "julian day (UT)" ;
    > 		time:units = "days since 1990-01-01 00:00:00" ;
    > 		time:calendar = "standard" ;
    > 		time:axis = "T" ;
    > 	float TLON(nj, ni) ;
    > 		TLON:standard<sub>name</sub> = "longitude" ;
    > 		TLON:long<sub>name</sub> = "longitude" ;
    > 		TLON:units = "degrees<sub>east</sub>" ;
    > 		TLON:<sub>CoordinateAxisType</sub> = "Lon" ;
    > 	float TLAT(nj, ni) ;
    > 		TLAT:standard<sub>name</sub> = "latitude" ;
    > 		TLAT:long<sub>name</sub> = "latitude" ;
    > 		TLAT:units = "degrees<sub>north</sub>" ;
    > 		TLAT:<sub>CoordinateAxisType</sub> = "Lat" ;
    > 	float f(f) ;
    > 		f:long<sub>name</sub> = "wave<sub>frequency</sub>" ;
    > 		f:units = "s-1" ;
    > 		f:axis = "Hz" ;
    > 		f:standard<sub>name</sub> = "wave<sub>frequency</sub>" ;
    > 	float efreq(time, f, nj, ni) ;
    > 		efreq:standard<sub>name</sub> = "power<sub>spectral</sub><sub>density</sub><sub>of</sub><sub>surface</sub><sub>elevation</sub>" ;
    > 		efreq:long<sub>name</sub> = "wave<sub>elevation</sub><sub>spectrum</sub>" ;
    > 		efreq:units = "m2 s" ;
    > 		efreq:coordinates = "TLAT TLON" ;
    > 		efreq:<sub>FillValue</sub> = 9.96921e+36f ;
    > 		efreq:missing<sub>value</sub> = 9.96921e+36f ;
    > 		efreq:globwave<sub>name</sub> = "power<sub>spectral</sub><sub>density</sub><sub>of</sub><sub>surface</sub><sub>elevation</sub>" ;
    > 
    > // global attributes:
    > 		:CDI = "Climate Data Interface version 1.9.7.1 (<http://mpimet.mpg.de/cdi>)" ;
    > 		:history = "Wed Sep 04 15:05:27 2019: cdo chname,ef,efreq ww3.20100101<sub>ef</sub><sub>remapgx3.nc</sub> ww3.20100101<sub>efreq</sub><sub>remapgx3.nc\nWed</sub> Sep 04 13:37:05 2019: cdo remapnn,../gx3<sub>tgrid</sub><sub>.nc</sub> ww3.20100101<sub>ef</sub><sub>.nc</sub> ww3.20100101<sub>ef</sub><sub>remapgx3.nc\nWed</sub> Sep  4 13:36:54 2019: ncks -x -v MAPSTA ww3.20100101<sub>ef.nc</sub> ww3.20100101<sub>ef</sub><sub>.nc\nWed</sub> Sep  4 13:28:07 2019: ncatted -a coordinates,ef,c,c,longitude latitude ww3.20100101<sub>ef.nc</sub>" ;
    > 		:Conventions = "CF-1.6" ;
    > 		:WAVEWATCH<sub>III</sub><sub>version</sub><sub>number</sub> = "5.16" ;
    > 		:WAVEWATCH<sub>III</sub><sub>switches</sub> = "F90 NOGRB NOPA LRB4 NC4 PR3 UQ FLX0 LN1 ST4 BT1 DB1 MLIM NL1 TR0 BS0 REF0 IS0 IC4 XX0 WNT1 WNX1 CRT1 CRX1 O0 O1 O2 O3 O4 O5 O6 O7 O11 SHRD" ;
    > 		:product<sub>name</sub> = "ww3.20100101<sub>ef.nc</sub>" ;
    > 		:area = "POP 1 degree grid (gx1v6b)" ;
    > 		:latitude<sub>resolution</sub> = "n/a" ;
    > 		:longitude<sub>resolution</sub> = "n/a" ;
    > 		:southernmost<sub>latitude</sub> = "-79.22052" ;
    > 		:northernmost<sub>latitude</sub> = "89.70641" ;
    > 		:westernmost<sub>longitude</sub> = "1.4731102E-02" ;
    > 		:easternmost<sub>longitude</sub> = "359.9960" ;
    > 		:minimum<sub>altitude</sub> = "-12000 m" ;
    > 		:maximum<sub>altitude</sub> = "9000 m" ;
    > 		:altitude<sub>resolution</sub> = "n/a" ;
    > 		:start<sub>date</sub> = "2010-01-01 00:00:00" ;
    > 		:stop<sub>date</sub> = "2010-01-01 18:00:00" ;
    > 		:NCO = "netCDF Operators version 4.7.9 (Homepage = <http://nco.sf.net>, Code = <http://github.com/nco/nco>)" ;
    > 		:CDO = "Climate Data Operators version 1.9.7.1 (<http://mpimet.mpg.de/cdo>)" ;
    > }


<a id="orgd5aa7e0"></a>

### GX3 Initial Conditions

There are a total of 12 initial condition files provided. With each file
containing NIL time step across the This is interesting to
me as I thought only one would be required. Nonetheless, the header of
the first file is given below.

> netcdf iced<sub>gx3</sub><sub>v6.2005</sub>-01-01 {
> dimensions:
> 	ni = 100 ;
> 	nj = 116 ;
> 	ncat = 5 ;
> variables:
> 	double uvel(nj, ni) ;
> 	double vvel(nj, ni) ;
> 	double scale<sub>factor</sub>(nj, ni) ;
> 	double swvdr(nj, ni) ;
> 	double swvdf(nj, ni) ;
> 	double swidr(nj, ni) ;
> 	double swidf(nj, ni) ;
> 	double strocnxT(nj, ni) ;
> 	double strocnyT(nj, ni) ;
> 	double stressp<sub>1</sub>(nj, ni) ;
> 	double stressp<sub>2</sub>(nj, ni) ;
> 	double stressp<sub>3</sub>(nj, ni) ;
> 	double stressp<sub>4</sub>(nj, ni) ;
> 	double stressm<sub>1</sub>(nj, ni) ;
> 	double stressm<sub>2</sub>(nj, ni) ;
> 	double stressm<sub>3</sub>(nj, ni) ;
> 	double stressm<sub>4</sub>(nj, ni) ;
> 	double stress12<sub>1</sub>(nj, ni) ;
> 	double stress12<sub>2</sub>(nj, ni) ;
> 	double stress12<sub>3</sub>(nj, ni) ;
> 	double stress12<sub>4</sub>(nj, ni) ;
> 	double iceumask(nj, ni) ;
> 	double sst(nj, ni) ;
> 	double frzmlt(nj, ni) ;
> 	double frz<sub>onset</sub>(nj, ni) ;
> 	double fsnow(nj, ni) ;
> 	double aicen(ncat, nj, ni) ;
> 	double vicen(ncat, nj, ni) ;
> 	double vsnon(ncat, nj, ni) ;
> 	double Tsfcn(ncat, nj, ni) ;
> 	double iage(ncat, nj, ni) ;
> 	double FY(ncat, nj, ni) ;
> 	double alvl(ncat, nj, ni) ;
> 	double vlvl(ncat, nj, ni) ;
> 	double apnd(ncat, nj, ni) ;
> 	double hpnd(ncat, nj, ni) ;
> 	double ipnd(ncat, nj, ni) ;
> 	double dhs(ncat, nj, ni) ;
> 	double ffrac(ncat, nj, ni) ;
> 	double fbrn(ncat, nj, ni) ;
> 	double first<sub>ice</sub>(ncat, nj, ni) ;
> 	double sice001(ncat, nj, ni) ;
> 	double qice001(ncat, nj, ni) ;
> 	double sice002(ncat, nj, ni) ;
> 	double qice002(ncat, nj, ni) ;
> 	double sice003(ncat, nj, ni) ;
> 	double qice003(ncat, nj, ni) ;
> 	double sice004(ncat, nj, ni) ;
> 	double qice004(ncat, nj, ni) ;
> 	double sice005(ncat, nj, ni) ;
> 	double qice005(ncat, nj, ni) ;
> 	double sice006(ncat, nj, ni) ;
> 	double qice006(ncat, nj, ni) ;
> 	double sice007(ncat, nj, ni) ;
> 	double qice007(ncat, nj, ni) ;
> 	double qsno001(ncat, nj, ni) ;
> 
> // global attributes:
> 		:istep1 = 87672 ;
> 		:time = 315619200. ;
> 		:time<sub>forc</sub> = 315619200. ;
> 		:nyr = 11 ;
> 		:month = 1 ;
> 		:mday = 1 ;
> 		:sec = 0 ;
> 		:created = "2021-03-19" ;
> 		:history = "Wed Apr  7 12:29:43 2021: ncatted &#x2013;attribute created,global,c,c,2021-03-19 iced<sub>gx3</sub><sub>v6.2005</sub>-01-01.nc" ;
> 		:NCO = "netCDF Operators version 4.9.5 (Homepage = <http://nco.sf.net>, Code = <http://github.com/nco/nco>)" ;
> }


<a id="org05ab192"></a>

## GX1 Description

Global One-Degree Cartesian


<a id="orgd404fe0"></a>

### GX1 Grid

For some reason the [GX1 grid file](file:///Volumes/ioa03/cice-dirs/input/CICE_data/grid/gx1/grid_gx1.bin) is in binary format and the header information
is not easily accessible.


<a id="org31c07a2"></a>

### GX1 Bathymetry

Header from the [GX1 bathymetry file](file:///Volumes/ioa03/cice-dirs/input/CICE_data/grid/gx1/grid_gx1.bin):

> netcdf global<sub>gx1.bathy</sub> {
> dimensions:
> 	nj = 384 ;
> 	ni = 320 ;
> variables:
> 	float Bathymetry(nj, ni) ;
> 		Bathymetry:long<sub>name</sub> = "ocean bathymetry for grounding scheme" ;
> 		Bathymetry:units = "m" ;
> 		Bathymetry:coordinates = "TLON TLAT" ;
> 		Bathymetry:missing<sub>value</sub> = 1.e+30f ;
> 		Bathymetry:<sub>FillValue</sub> = 1.e+30f ;
> 	float TLAT(nj, ni) ;
> 		TLAT:long<sub>name</sub> = "T grid center latitude" ;
> 		TLAT:units = "degrees<sub>north</sub>" ;
> 		TLAT:missing<sub>value</sub> = 1.e+30f ;
> 		TLAT:<sub>FillValue</sub> = 1.e+30f ;
> 	float TLON(nj, ni) ;
> 		TLON:long<sub>name</sub> = "T grid center longitude" ;
> 		TLON:units = "degrees<sub>east</sub>" ;
> 		TLON:missing<sub>value</sub> = 1.e+30f ;
> 		TLON:<sub>FillValue</sub> = 1.e+30f ;
> 
> // global attributes:
> 		:title = "ocean bathymetry for grounding scheme" ;
> 		:source = "created from bathy<sub>ORCA025</sub><sub>LIM.std</sub>" ;
> 		:history = "Mon Oct 29 16:01:01 2018: ncatted -a source,global,o,c,created from bathy<sub>ORCA025</sub><sub>LIM.std</sub> gx1<sub>bathyTP10.nc</sub> gx1<sub>bathyTP11.nc\n</sub>",
> 			"Mon Oct 29 16:00:02 2018: ncatted -a contents,global,d,, gx1<sub>bathyTP9.nc</sub> gx1<sub>bathyTP10.nc\n</sub>",
> 			"Mon Oct 29 15:59:40 2018: ncatted -a comment,global,d,, gx1<sub>bathyTP8.nc</sub> gx1<sub>bathyTP9.nc\n</sub>",
> 			"Mon Oct 29 15:59:33 2018: ncatted -a comment2,global,d,, gx1<sub>bathyTP7.nc</sub> gx1<sub>bathyTP8.nc\n</sub>",
> 			"Mon Oct 29 15:59:26 2018: ncatted -a comment3,global,d,, gx1<sub>bathyTP6.nc</sub> gx1<sub>bathyTP7.nc\n</sub>",
> 			"Mon Oct 29 15:59:13 2018: ncatted -a conventions,global,d,, gx1<sub>bathyTP5.nc</sub> gx1<sub>bathyTP6.nc\n</sub>",
> 			"Mon Oct 29 15:27:15 2018: ncatted -a title,global,o,c,ocean bathymetry for grounding scheme gx1<sub>bathyTP4.nc</sub> gx1<sub>bathyTP5.nc\n</sub>",
> 			"Mon Oct 29 15:13:20 2018: ncks -x -v time,time<sub>bounds</sub> gx1<sub>bathyTP3.nc</sub> gx1<sub>bathyTP4.nc\n</sub>",
> 			"Mon Oct 29 15:13:04 2018: ncatted -a units,Bathymetry,o,c,m gx1<sub>bathyTP2.nc</sub> gx1<sub>bathyTP3.nc\n</sub>",
> 			"Mon Oct 29 15:12:46 2018: ncatted -a long<sub>name,Bathymetry,o,c,ocean</sub> bathymetry for grounding scheme gx1<sub>bathy.nc</sub> gx1<sub>bathyTP2.nc\n</sub>",
> 			"Fri Oct 26 13:28:47 2018: ncrename -v tarea,Bathymetry tempgx1.nc\n",
> 			"Thu Oct 25 17:17:42 2018: ncks -x -v ANGLE,ANGLET,NCAT,ULAT,ULON,VGRDa,aice,blkmask,hi,taubx,tauby,time<sub>bounds,tmask,uarea,uvel,vvel</sub> tempgx1.nc tempgx1out.nc\n",
> 			"This dataset was created on 2018-10-22 at 18:12:11.5" ;
> 		:io<sub>flavor</sub> = "io<sub>netcdf</sub>" ;
> 		:NCO = "4.4.2" ;
> }


<a id="org2268e1c"></a>

### GX1 Land Mask

For some reason the [GX1 land mask file](file:///Volumes/ioa03/cice-dirs/input/CICE_data/grid/gx1/kmt_gx1.bin) is in binary format and the header
information is not easily accessible.


<a id="org2a6fb4c"></a>

### GX1 Forcing

Four datasets are used for forcing: [1.3.4.1](#org8da65c5), [1.3.4.2](#orgec47000), [1.3.4.3](#org0bd7e04), and [1.3.4.4](#org01f43ad). I find it interesting to note that WaveWatch III data is not provided.

1.  GX1 CESM

    One file provides monthly summarised oceanographic data in the following format.
    Note that the time dimension is only 12 increments in length.
    
    > netcdf ocean<sub>forcing</sub><sub>clim</sub><sub>2D</sub><sub>gx1</sub> {
    > dimensions:
    > 	time = 12 ;
    > 	nj = 384 ;
    > 	ni = 320 ;
    > variables:
    > 	double area(nj, ni) ;
    > 		area:long<sub>name</sub> = "area of grid cell in radians squared" ;
    > 		area:units = "area" ;
    > 	int mask(nj, ni) ;
    > 		mask:units = "unitless" ;
    > 		mask:long<sub>name</sub> = "domain maskr" ;
    > 		mask:<sub>FillValue</sub> = -2147483647 ;
    > 		mask:missing<sub>value</sub> = -2147483647 ;
    > 	double xc(nj, ni) ;
    > 		xc:missing<sub>value</sub> = 9.96920996838687e+36 ;
    > 		xc:<sub>FillValue</sub> = 9.96920996838687e+36 ;
    > 		xc:units = "degrees east" ;
    > 		xc:long<sub>name</sub> = "longitude of grid cell center" ;
    > 	double yc(nj, ni) ;
    > 		yc:missing<sub>value</sub> = 9.96920996838687e+36 ;
    > 		yc:<sub>FillValue</sub> = 9.96920996838687e+36 ;
    > 		yc:units = "degrees north" ;
    > 		yc:long<sub>name</sub> = "latitude of grid cell center" ;
    > 	float time(time) ;
    > 		time:calendar = "noleap" ;
    > 		time:long<sub>name</sub> = "observation time" ;
    > 		time:units = "days since 0001-01-01 00:00:00" ;
    > 	float S(time, nj, ni) ;
    > 		S:units = "ppt" ;
    > 		S:long<sub>name</sub> = "salinity" ;
    > 		S:<sub>FillValue</sub> = 9.96921e+36f ;
    > 	float T(time, nj, ni) ;
    > 		T:units = "degC" ;
    > 		T:long<sub>name</sub> = "temperature" ;
    > 		T:<sub>FillValue</sub> = 9.96921e+36f ;
    > 	float U(time, nj, ni) ;
    > 		U:units = "m/s" ;
    > 		U:long<sub>name</sub> = "u ocean current" ;
    > 		U:<sub>FillValue</sub> = 9.96921e+36f ;
    > 	float V(time, nj, ni) ;
    > 		V:units = "m/s" ;
    > 		V:long<sub>name</sub> = "v ocean current" ;
    > 		V:<sub>FillValue</sub> = 9.96921e+36f ;
    > 	float dhdx(time, nj, ni) ;
    > 		dhdx:units = "m/m" ;
    > 		dhdx:long<sub>name</sub> = "ocean surface slope: zonal" ;
    > 		dhdx:<sub>FillValue</sub> = 9.96921e+36f ;
    > 	float dhdy(time, nj, ni) ;
    > 		dhdy:units = "m/m" ;
    > 		dhdy:long<sub>name</sub> = "ocean surface slope: meridional" ;
    > 		dhdy:<sub>FillValue</sub> = 9.96921e+36f ;
    > 	float hblt(time, nj, ni) ;
    > 		hblt:units = "m" ;
    > 		hblt:long<sub>name</sub> = "boundary layer depth" ;
    > 		hblt:<sub>FillValue</sub> = 9.96921e+36f ;
    > 	float qdp(time, nj, ni) ;
    > 		qdp:units = "W/m<sup>2</sup>" ;
    > 		qdp:long<sub>name</sub> = "ocean heat flux convergence" ;
    > 		qdp:<sub>FillValue</sub> = 9.96921e+36f ;
    > 
    > // global attributes:
    > 		:creation<sub>date</sub> = "Tue Feb 17 09:45:48 MST 2015" ;
    > 		:comment = "This data is on the displaced pole grid gx1v5" ;
    > 		:calendar = "standard" ;
    > 		:author = "D. Bailey" ;
    > 		:note3 = "qdp is computed from depth summed ocean column" ;
    > 		:note2 = "all fields interpolated to T-grid" ;
    > 		:note1 = "fields computed from years 402 to 1510 monthly means from pop" ;
    > 		:description = "Input data for DOCN7 mixed layer model from b.e11.B1850C5CN.f09<sub>g16.005</sub>" ;
    > 		:source = "pop<sub>frc.ncl</sub>" ;
    > 		:conventions = "CCSM data model domain description" ;
    > 		:title = "Monthly averaged ocean forcing from POP output" ;
    > 		:history = "Wed Mar 24 11:42:29 2021: ncatted &#x2013;glb<sub>att</sub><sub>add</sub> description2=From the CESM-LE control run as in Kay et al. 2015 pop<sub>frc.gx1.nc</sub>" ;
    > 		:description2 = "From the CESM-LE control run as in Kay et al. 2015" ;
    > 		:NCO = "netCDF Operators version 4.9.5 (Homepage = <http://nco.sf.net>, Code = <http://github.com/nco/nco>)" ;
    >     }

2.  GX1 COREII

    These files are given as binary with two sub-directories [4XDAILY](file:///Volumes/ioa03/cice-dirs/input/CICE_data/forcing/gx1/COREII/4XDAILY) and [MONTHLY](file:///Volumes/ioa03/cice-dirs/input/CICE_data/forcing/gx1/COREII/MONTHLY).
    The file names indicate the representative variable contained within.
    
    Here is the directory listing of the [4XDAILY](file:///Volumes/ioa03/cice-dirs/input/CICE_data/forcing/gx1/COREII/4XDAILY/):
    
    > -rwxr-xr-x  1 dpath2o  staff   1.3G 30 Aug  2017 q<sub>10.2005.dat</sub>
    > -rwxr-xr-x  1 dpath2o  staff   1.3G 30 Aug  2017 q<sub>10.2006.dat</sub>
    > -rwxr-xr-x  1 dpath2o  staff   1.3G 30 Aug  2017 q<sub>10.2007.dat</sub>
    > -rwxr-xr-x  1 dpath2o  staff   1.3G 30 Aug  2017 q<sub>10.2008.dat</sub>
    > -rwxr-xr-x  1 dpath2o  staff   1.3G 30 Aug  2017 q<sub>10.2009.dat</sub>
    > -rwxr-xr-x  1 dpath2o  staff   1.3G 30 Aug  2017 t<sub>10.2005.dat</sub>
    > -rwxr-xr-x  1 dpath2o  staff   1.3G 30 Aug  2017 t<sub>10.2006.dat</sub>
    > -rwxr-xr-x  1 dpath2o  staff   1.3G 30 Aug  2017 t<sub>10.2007.dat</sub>
    > -rwxr-xr-x  1 dpath2o  staff   1.3G 30 Aug  2017 t<sub>10.2008.dat</sub>
    > -rwxr-xr-x  1 dpath2o  staff   1.3G 30 Aug  2017 t<sub>10.2009.dat</sub>
    > -rwxr-xr-x  1 dpath2o  staff   1.3G 30 Aug  2017 u<sub>10.2005.dat</sub>
    > -rwxr-xr-x  1 dpath2o  staff   1.3G 30 Aug  2017 u<sub>10.2006.dat</sub>
    > -rwxr-xr-x  1 dpath2o  staff   1.3G 30 Aug  2017 u<sub>10.2007.dat</sub>
    > -rwxr-xr-x  1 dpath2o  staff   1.3G 30 Aug  2017 u<sub>10.2008.dat</sub>
    > -rwxr-xr-x  1 dpath2o  staff   1.3G 31 Aug  2017 u<sub>10.2009.dat</sub>
    > -rwxr-xr-x  1 dpath2o  staff   1.3G 30 Aug  2017 v<sub>10.2005.dat</sub>
    > -rwxr-xr-x  1 dpath2o  staff   1.3G 30 Aug  2017 v<sub>10.2006.dat</sub>
    > -rwxr-xr-x  1 dpath2o  staff   1.3G 30 Aug  2017 v<sub>10.2007.dat</sub>
    > -rwxr-xr-x  1 dpath2o  staff   1.3G 30 Aug  2017 v<sub>10.2008.dat</sub>
    > -rwxr-xr-x  1 dpath2o  staff   1.3G 31 Aug  2017 v<sub>10.2009.dat</sub>
    
    Here is the listing of the [MONTHLY](file:///Volumes/ioa03/cice-dirs/input/CICE_data/forcing/gx1/COREII/MONTHLY):
    
    > -rw-r&#x2013;r&#x2013;  1 dpath2o  staff    11M 12 Jul  2012 cldf.omip.dat
    > -rw-r&#x2013;r&#x2013;  1 dpath2o  staff    11M 12 Jul  2012 prec.nmyr.dat

3.  GX1 JRA55

    The forcing files are nearly identical to [1.2.4.1](#org12fd10a) with the exception of the
    change in grid. This is most notable in comparing the file sizes &#x2013; [1.2.4.1](#org12fd10a)
    files are 0.9GB and these GX1 JRA55 files are 9.5GB. Here is the directory
    listing of the [8XDAILY](file:///Volumes/ioa03/cice-dirs/input/CICE_data/forcing/gx1/JRA55/8XDAILY/):
    
    > -rw-r&#x2013;r&#x2013;@ 1 dpath2o  staff   9.4G 22 Aug  2019 JRA55<sub>03hr</sub><sub>forcing</sub><sub>2005.nc</sub>
    > -rw-r&#x2013;r&#x2013;  1 dpath2o  staff   9.4G 22 Aug  2019 JRA55<sub>03hr</sub><sub>forcing</sub><sub>2006.nc</sub>
    > -rw-r&#x2013;r&#x2013;  1 dpath2o  staff   9.4G 22 Aug  2019 JRA55<sub>03hr</sub><sub>forcing</sub><sub>2007.nc</sub>
    > -rw-r&#x2013;r&#x2013;  1 dpath2o  staff   9.4G 22 Aug  2019 JRA55<sub>03hr</sub><sub>forcing</sub><sub>2008.nc</sub>
    > -rw-r&#x2013;r&#x2013;  1 dpath2o  staff   9.4G 22 Aug  2019 JRA55<sub>03hr</sub><sub>forcing</sub><sub>2009.nc</sub>
    
    Here is the header information of the first file:
    
    > netcdf JRA55<sub>03hr</sub><sub>forcing</sub><sub>2005</sub> {
    > dimensions:
    > 	nj = 384 ;
    > 	ni = 320 ;
    > 	time = UNLIMITED ; // (2920 currently)
    > variables:
    > 	float LAT(nj, ni) ;
    > 		LAT:long<sub>name</sub> = "Latitude" ;
    > 		LAT:units = "degrees<sub>north</sub>" ;
    > 	float LON(nj, ni) ;
    > 		LON:long<sub>name</sub> = "Longitude" ;
    > 		LON:units = "degrees<sub>east</sub>" ;
    > 	float airtmp(time, nj, ni) ;
    > 		airtmp:long<sub>name</sub> = "Air temperature" ;
    > 		airtmp:units = "Kelvin" ;
    > 		airtmp:coordinates = "LON LAT" ;
    > 		airtmp:<sub>FillValue</sub> = 1.e+30f ;
    > 	float dlwsfc(time, nj, ni) ;
    > 		dlwsfc:long<sub>name</sub> = "downward surface longwave" ;
    > 		dlwsfc:units = "W/m<sup>2</sup>" ;
    > 		dlwsfc:coordinates = "LON LAT" ;
    > 		dlwsfc:<sub>FillValue</sub> = 1.e+30f ;
    > 	float glbrad(time, nj, ni) ;
    > 		glbrad:long<sub>name</sub> = "downward surface shortwave" ;
    > 		glbrad:units = "W/m<sup>2</sup>" ;
    > 		glbrad:coordinates = "LON LAT" ;
    > 		glbrad:<sub>FillValue</sub> = 1.e+30f ;
    > 	float spchmd(time, nj, ni) ;
    > 		spchmd:long<sub>name</sub> = "Specific Humidity" ;
    > 		spchmd:units = "kg/kg" ;
    > 		spchmd:coordinates = "LON LAT" ;
    > 		spchmd:<sub>FillValue</sub> = 1.e+30f ;
    > 	float time(time) ;
    > 		time:long<sub>name</sub> = "model time" ;
    > 		time:units = "days since 1900-12-31 00:00:00" ;
    > 		time:calendar = "standard" ;
    > 	float ttlpcp(time, nj, ni) ;
    > 		ttlpcp:long<sub>name</sub> = "Total precipiation" ;
    > 		ttlpcp:units = "kg m-2 s-1" ;
    > 		ttlpcp:coordinates = "LON LAT" ;
    > 		ttlpcp:<sub>FillValue</sub> = 1.e+30f ;
    > 		ttlpcp:comment = "derived from bucket dumps" ;
    > 	float wndewd(time, nj, ni) ;
    > 		wndewd:long<sub>name</sub> = "x-ward winds" ;
    > 		wndewd:units = "m/s" ;
    > 		wndewd:coordinates = "LON LAT" ;
    > 		wndewd:<sub>FillValue</sub> = 1.e+30f ;
    > 	float wndnwd(time, nj, ni) ;
    > 		wndnwd:long<sub>name</sub> = "y-ward winds" ;
    > 		wndnwd:units = "m/s" ;
    > 		wndnwd:coordinates = "LON LAT" ;
    > 		wndnwd:<sub>FillValue</sub> = 1.e+30f ;
    > 
    > // global attributes:
    > 		:history = "Mon Aug  5 15:31:32 2019: ncks -v time,LON,LAT,glbrad,dlwsfc,wndewd,wndnwd,airtmp,spchmd,ttlpcp JRA55<sub>03hr</sub><sub>forcing</sub><sub>20050101.00.netrad.nc</sub> JRA55<sub>03hr</sub><sub>forcing</sub><sub>20050101.00.nc</sub>" ;
    > 		:NCO = "netCDF Operators version 4.7.5 (Homepage = <http://nco.sf.net>, Code = <http://github.com/nco/nco>)" ;
    > }

4.  GX1 WOA

    Monthly statistics for silicate and nitrate are provided across two files. The
    header information for these two files is provided.
    
    [Nitrate](file:///Volumes/ioa03/cice-dirs/input/CICE_data/forcing/gx1/WOA/MONTHLY/nitrate_climatologyWOA_gx1v6f_20150107.nc) file:
    
    > netcdf nitrate<sub>climatologyWOA</sub><sub>gx1v6f</sub><sub>20150107</sub> {
    > dimensions:
    > 	time = 12 ;
    > 	nj = 384 ;
    > 	ni = 320 ;
    > variables:
    > 	float time(time) ;
    > 	int nj(nj) ;
    > 	int ni(ni) ;
    > 	float TLAT(nj, ni) ;
    > 		TLAT:<sub>FillValue</sub> = 1.e+30f ;
    > 	float TLON(nj, ni) ;
    > 		TLON:<sub>FillValue</sub> = 1.e+30f ;
    > 	double nitrate(time, nj, ni) ;
    > 		nitrate:<sub>FillValue</sub> = 1.00000001504747e+30 ;
    > }
    
    [Silicate](file:///Volumes/ioa03/cice-dirs/input/CICE_data/forcing/gx1/WOA/MONTHLY/silicate_climatologyWOA_gx1v6f_20150107.nc) file:
    
    > netcdf silicate<sub>climatologyWOA</sub><sub>gx1v6f</sub><sub>20150107</sub> {
    > dimensions:
    > 	time = 12 ;
    > 	nj = 384 ;
    > 	ni = 320 ;
    > variables:
    > 	float time(time) ;
    > 		time:<sub>FillValue</sub> = 9.96921e+36f ;
    > 	int nj(nj) ;
    > 	int ni(ni) ;
    > 	float TLAT(nj, ni) ;
    > 		TLAT:<sub>FillValue</sub> = 1.e+30f ;
    > 	float TLON(nj, ni) ;
    > 		TLON:<sub>FillValue</sub> = 1.e+30f ;
    > 	float silicate(time, nj, ni) ;
    > 		silicate:<sub>FillValue</sub> = 1.e+30f ;
    > }


<a id="orgfd2e842"></a>

### GX1 Initial Conditions

These files are nearly identical to the [1.2.5](#orgd5aa7e0) **other than the
underlying grid**. So there are a total of 12 initial condition files provided.
With each file containing NIL time step across the This is interesting to me as
I thought only one would be required. Nonetheless, the header of the first file
is given below.

> netcdf iced<sub>gx1</sub><sub>v6.2005</sub>-01-01 {
> dimensions:
> 	ni = 320 ;
> 	nj = 384 ;
> 	ncat = 5 ;
> variables:
> 	double uvel(nj, ni) ;
> 	double vvel(nj, ni) ;
> 	double scale<sub>factor</sub>(nj, ni) ;
> 	double swvdr(nj, ni) ;
> 	double swvdf(nj, ni) ;
> 	double swidr(nj, ni) ;
> 	double swidf(nj, ni) ;
> 	double strocnxT(nj, ni) ;
> 	double strocnyT(nj, ni) ;
> 	double stressp<sub>1</sub>(nj, ni) ;
> 	double stressp<sub>2</sub>(nj, ni) ;
> 	double stressp<sub>3</sub>(nj, ni) ;
> 	double stressp<sub>4</sub>(nj, ni) ;
> 	double stressm<sub>1</sub>(nj, ni) ;
> 	double stressm<sub>2</sub>(nj, ni) ;
> 	double stressm<sub>3</sub>(nj, ni) ;
> 	double stressm<sub>4</sub>(nj, ni) ;
> 	double stress12<sub>1</sub>(nj, ni) ;
> 	double stress12<sub>2</sub>(nj, ni) ;
> 	double stress12<sub>3</sub>(nj, ni) ;
> 	double stress12<sub>4</sub>(nj, ni) ;
> 	double iceumask(nj, ni) ;
> 	double sst(nj, ni) ;
> 	double frzmlt(nj, ni) ;
> 	double frz<sub>onset</sub>(nj, ni) ;
> 	double fsnow(nj, ni) ;
> 	double aicen(ncat, nj, ni) ;
> 	double vicen(ncat, nj, ni) ;
> 	double vsnon(ncat, nj, ni) ;
> 	double Tsfcn(ncat, nj, ni) ;
> 	double iage(ncat, nj, ni) ;
> 	double FY(ncat, nj, ni) ;
> 	double alvl(ncat, nj, ni) ;
> 	double vlvl(ncat, nj, ni) ;
> 	double apnd(ncat, nj, ni) ;
> 	double hpnd(ncat, nj, ni) ;
> 	double ipnd(ncat, nj, ni) ;
> 	double dhs(ncat, nj, ni) ;
> 	double ffrac(ncat, nj, ni) ;
> 	double fbrn(ncat, nj, ni) ;
> 	double first<sub>ice</sub>(ncat, nj, ni) ;
> 	double sice001(ncat, nj, ni) ;
> 	double qice001(ncat, nj, ni) ;
> 	double sice002(ncat, nj, ni) ;
> 	double qice002(ncat, nj, ni) ;
> 	double sice003(ncat, nj, ni) ;
> 	double qice003(ncat, nj, ni) ;
> 	double sice004(ncat, nj, ni) ;
> 	double qice004(ncat, nj, ni) ;
> 	double sice005(ncat, nj, ni) ;
> 	double qice005(ncat, nj, ni) ;
> 	double sice006(ncat, nj, ni) ;
> 	double qice006(ncat, nj, ni) ;
> 	double sice007(ncat, nj, ni) ;
> 	double qice007(ncat, nj, ni) ;
> 	double qsno001(ncat, nj, ni) ;
> 
> // global attributes:
> 		:istep1 = 87672 ;
> 		:time = 315619200. ;
> 		:time<sub>forc</sub> = 315619200. ;
> 		:nyr = 11 ;
> 		:month = 1 ;
> 		:mday = 1 ;
> 		:sec = 0 ;
> 		:created = "2021-03-30" ;
> 		:history = "Fri Apr  2 08:48:44 2021: ncatted &#x2013;attribute created,global,c,c,2021-03-30 iced<sub>gx1</sub><sub>v6.2005</sub>-01-01.nc" ;
> 		:NCO = "netCDF Operators version 4.9.5 (Homepage = <http://nco.sf.net>, Code = <http://github.com/nco/nco>)" ;
> }


<a id="org79cf8e6"></a>

## TX1 Description

Global One-Degree Tri-polar


<a id="org6421b79"></a>

### TX1 Grid

Here this [grid file](file:///Volumes/ioa03/cice-dirs/input/CICE_data/grid/tx1/grid_tx1.nc) is given in both NetCDF and binary format, hence here is the
header contents of this file. This header contents indicates are fair amount
more information than the [1.2.1](#org2b3ded1) file.

> netcdf grid<sub>tx1</sub> {
> dimensions:
> 	nx = 360 ;
> 	ny = 300 ;
> 	nc = 4 ;
> variables:
> 	double ulat(ny, nx) ;
> 		ulat:units = "radians" ;
> 		ulat:title = "Latitude of U points" ;
> 	double ulon(ny, nx) ;
> 		ulon:units = "radians" ;
> 		ulon:title = "Longitude of U points" ;
> 	double tlat(ny, nx) ;
> 		tlat:units = "radians" ;
> 		tlat:title = "Latitude of T points" ;
> 	double tlon(ny, nx) ;
> 		tlon:units = "radians" ;
> 		tlon:title = "Longitude of T points" ;
> 	double htn(ny, nx) ;
> 		htn:units = "cm" ;
> 		htn:title = "Width of T cells on N side" ;
> 	double hte(ny, nx) ;
> 		hte:units = "cm" ;
> 		hte:title = "Width of T cells on E side" ;
> 	double hun(ny, nx) ;
> 		hun:units = "cm" ;
> 		hun:title = "Width of U cells on N side" ;
> 	double hue(ny, nx) ;
> 		hue:units = "cm" ;
> 		hue:title = "Width of U cells on E side" ;
> 	double angle(ny, nx) ;
> 		angle:units = "radians" ;
> 		angle:title = "Rotation Angle of U cells" ;
> 	double angleT(ny, nx) ;
> 		angleT:units = "radians" ;
> 		angleT:title = "Rotation Angle of T cells" ;
> 	double tarea(ny, nx) ;
> 		tarea:units = "m<sup>2</sup>" ;
> 		tarea:title = "area in m<sup>2</sup> for T cells" ;
> 	double uarea(ny, nx) ;
> 		uarea:units = "m<sup>2</sup>" ;
> 		uarea:title = "area in m<sup>2</sup> for U cells" ;
> 	double lont<sub>bonds</sub>(nc, ny, nx) ;
> 		lont<sub>bonds</sub>:units = "radians" ;
> 		lont<sub>bonds</sub>:title = "long. of T cell vertices" ;
> 	double latt<sub>bonds</sub>(nc, ny, nx) ;
> 		latt<sub>bonds</sub>:units = "radians" ;
> 		latt<sub>bonds</sub>:title = "lat. of T cell vertices" ;
> 	double lonu<sub>bonds</sub>(nc, ny, nx) ;
> 		lonu<sub>bonds</sub>:units = "radians" ;
> 		lonu<sub>bonds</sub>:title = "long. of U cell vertices" ;
> 	double latu<sub>bonds</sub>(nc, ny, nx) ;
> 		latu<sub>bonds</sub>:units = "radians" ;
> 		latu<sub>bonds</sub>:title = "lat. of U cell vertices" ;
> }


<a id="orgc70422f"></a>

### TX1 Bathymetry

The [bathymetry file](file:///Volumes/ioa03/cice-dirs/input/CICE_data/grid/tx1/tx1_bathy.nc) for the tri-polar grid is pretty straight forward:

> netcdf tx1<sub>bathy</sub> {
> dimensions:
> 	dim<sub>j</sub> = 240 ;
> 	dim<sub>i</sub> = 360 ;
> variables:
> 	double Bathymetry(dim<sub>j</sub>, dim<sub>i</sub>) ;
> 		Bathymetry:long<sub>name</sub> = "ocean bathymetry for grounding scheme" ;
> 		Bathymetry:units = "m" ;
> 		Bathymetry:<sub>FillValue</sub> = 1.e+30 ;
> 	double lat(dim<sub>j</sub>, dim<sub>i</sub>) ;
> 		lat:long<sub>name</sub> = "Latitude" ;
> 		lat:units = "degrees<sub>north</sub>" ;
> 		lat:<sub>FillValue</sub> = 1.e+30 ;
> 	double lon(dim<sub>j</sub>, dim<sub>i</sub>) ;
> 		lon:long<sub>name</sub> = "Longitude" ;
> 		lon:units = "degrees<sub>east</sub>" ;
> 		lon:<sub>FillValue</sub> = 1.e+30 ;
> 
> // global attributes:
> 		:history = "Wed Feb 10 21:25:30 2021: ncatted -a units,Bathymetry,o,c,m temptx1.nc tx1<sub>bathy.nc\n</sub>",
> 			"Wed Feb 10 21:24:26 2021: ncatted -a long<sub>name,Bathymetry,o,c,ocean</sub> bathymetry for grounding scheme temptx1.nc temptx1.nc\n",
> 			"Wed Feb 10 21:23:52 2021: ncks -x -v dx,dy,mask temptx1.nc temptx1.nc\n",
> 			"Wed Feb 10 21:22:19 2021: ncrename -v ang,Bathymetry temptx1.nc" ;
> 		:NCO = "4.4.2" ;
> }


<a id="org1ae6ae3"></a>

### TX1 Land Mask

The [land mask file](file:///Volumes/ioa03/cice-dirs/input/CICE_data/grid/tx1/kmt_tx1.bin) is binary and cannot be easily shown here.


<a id="org0eecc9b"></a>

### TX1 Forcing

The only dataset utilised for the tri-polar forcing is JRA55. It is equivalent
in it's structure and content to the other two JRA55 ([1.2.4.1](#org12fd10a) and [1.3.4.3](#org0bd7e04)),
but just re-gridded onto the [1.4.1](#org6421b79).


<a id="org9249e4b"></a>

### TX1 Initial Conditions

There is no version 6 files provided [1](#org28fbdde)


<a id="org2e0c1f6"></a>

## Initial Conditions Discussion

All three (GX3, GX1 and TX1) cases of [1.1.2](#orgf8b290e) have the same
structure of the initial condition files. Each case provides twelve files named
for each month. The contents of those files show 56 gridded parameters over two
and three dimensions (one categorical dimension and two spatial dimensions).
There is not a time dimension. Hence time-step is implied as monthly with 12
files named for each month &#x2013; it is not known wether this is a snap snot or an
average over the month. The table below is made in order to consolidate the
initial condition parameters should at some point it be determined that initial
conditions should be provided. The short name and description (the first two
columns) are obtained from correlating the [1.1.2](#orgf8b290e) information
and the information found in the documentation, online [here](https://cice-consortium-cice.readthedocs.io/en/main/cice_index.html).

<table border="2" cellspacing="0" cellpadding="6" rules="groups" frame="hsides">


<colgroup>
<col  class="org-left" />

<col  class="org-left" />

<col  class="org-left" />

<col  class="org-left" />
</colgroup>
<thead>
<tr>
<th scope="col" class="org-left">CICE6 short</th>
<th scope="col" class="org-left">CICE6 description</th>
<th scope="col" class="org-left">Dims</th>
<th scope="col" class="org-left">Units</th>
</tr>
</thead>

<tbody>
<tr>
<td class="org-left">uvel</td>
<td class="org-left">sea ice velocity, eastward</td>
<td class="org-left">ni, nj</td>
<td class="org-left">m/s</td>
</tr>


<tr>
<td class="org-left">vvel</td>
<td class="org-left">sea ice velocity, northward</td>
<td class="org-left">ni, nj</td>
<td class="org-left">m/s</td>
</tr>


<tr>
<td class="org-left">scale<sub>factor</sub></td>
<td class="org-left">timestep offsets in shortwave radiation components</td>
<td class="org-left">ni, nj</td>
<td class="org-left">-</td>
</tr>


<tr>
<td class="org-left">swvdr</td>
<td class="org-left">incoming shortwave radiation, visible direct</td>
<td class="org-left">ni, nj</td>
<td class="org-left">W/m<sup>2</sup></td>
</tr>


<tr>
<td class="org-left">swvdf</td>
<td class="org-left">incoming shortwave radiation, visible diffuse</td>
<td class="org-left">ni, nj</td>
<td class="org-left">W/m<sup>2</sup></td>
</tr>


<tr>
<td class="org-left">swidr</td>
<td class="org-left">incoming shortwave radiation, infrared direct</td>
<td class="org-left">ni, nj</td>
<td class="org-left">W/m<sup>2</sup></td>
</tr>


<tr>
<td class="org-left">swidf</td>
<td class="org-left">incoming shortwave radiation, infrared diffuse</td>
<td class="org-left">ni, nj</td>
<td class="org-left">W/m<sup>2</sup></td>
</tr>


<tr>
<td class="org-left">strocnxT</td>
<td class="org-left">ice-ocean stress</td>
<td class="org-left">ni, nj</td>
<td class="org-left">N/m</td>
</tr>


<tr>
<td class="org-left">strocnyT</td>
<td class="org-left">ice-ocean stress</td>
<td class="org-left">ni, nj</td>
<td class="org-left">N/m</td>
</tr>


<tr>
<td class="org-left">stressp_[1-4]</td>
<td class="org-left">internal strain tensor</td>
<td class="org-left">ni, nj</td>
<td class="org-left">N/m</td>
</tr>


<tr>
<td class="org-left">stressm_[1-4]</td>
<td class="org-left">internal strain tensor</td>
<td class="org-left">ni, nj</td>
<td class="org-left">N/m</td>
</tr>


<tr>
<td class="org-left">stress12_[1-4]</td>
<td class="org-left">internal strain tensor</td>
<td class="org-left">ni, nj</td>
<td class="org-left">N/m</td>
</tr>


<tr>
<td class="org-left">iceumask</td>
<td class="org-left">the mask where there is ice in a gridcell</td>
<td class="org-left">ni, nj</td>
<td class="org-left">-</td>
</tr>


<tr>
<td class="org-left">sst</td>
<td class="org-left">sea surface temperature</td>
<td class="org-left">ni, nj</td>
<td class="org-left">C</td>
</tr>


<tr>
<td class="org-left">frzmlt</td>
<td class="org-left">ice ocean heat exhange</td>
<td class="org-left">ni, nj</td>
<td class="org-left">W/m<sup>2</sup></td>
</tr>


<tr>
<td class="org-left">frz<sub>onset</sub></td>
<td class="org-left">Julian day in which the onset of freezing occurs</td>
<td class="org-left">ni, nj</td>
<td class="org-left">-</td>
</tr>


<tr>
<td class="org-left">ulat</td>
<td class="org-left">u-grid latitudes</td>
<td class="org-left">ni, nj</td>
<td class="org-left">-</td>
</tr>


<tr>
<td class="org-left">ulon</td>
<td class="org-left">u-grid longitudes</td>
<td class="org-left">ni, nj</td>
<td class="org-left">-</td>
</tr>


<tr>
<td class="org-left">tlat</td>
<td class="org-left">t-grid latitudes</td>
<td class="org-left">ni, nj</td>
<td class="org-left">-</td>
</tr>


<tr>
<td class="org-left">tlon</td>
<td class="org-left">t-grid latitudes</td>
<td class="org-left">ni, nj</td>
<td class="org-left">-</td>
</tr>


<tr>
<td class="org-left">aicen</td>
<td class="org-left">ice fraction per category</td>
<td class="org-left">ncat, ni, nj</td>
<td class="org-left">-</td>
</tr>


<tr>
<td class="org-left">vicen</td>
<td class="org-left">ice volume per unit area per category</td>
<td class="org-left">ncat, ni, nj</td>
<td class="org-left">-</td>
</tr>


<tr>
<td class="org-left">vsnon</td>
<td class="org-left">snow volume per unit area per category</td>
<td class="org-left">ncat, ni, nj</td>
<td class="org-left">-</td>
</tr>


<tr>
<td class="org-left">Tsfcn</td>
<td class="org-left">surface snow/ice temperature per category</td>
<td class="org-left">ncat, ni, nj</td>
<td class="org-left">C</td>
</tr>


<tr>
<td class="org-left">iage</td>
<td class="org-left">ice age</td>
<td class="org-left">ncat, ni, nj</td>
<td class="org-left">-</td>
</tr>


<tr>
<td class="org-left">FY</td>
<td class="org-left">first year ice area</td>
<td class="org-left">ncat, ni, nj</td>
<td class="org-left">-</td>
</tr>


<tr>
<td class="org-left">alvl</td>
<td class="org-left">fraction of level ice area</td>
<td class="org-left">ncat, ni, nj</td>
<td class="org-left">-</td>
</tr>


<tr>
<td class="org-left">vlvl</td>
<td class="org-left">volume of the level ice area</td>
<td class="org-left">ncat, ni, nj</td>
<td class="org-left">-</td>
</tr>


<tr>
<td class="org-left">apnd</td>
<td class="org-left">fraction of ponds per cell</td>
<td class="org-left">ncat, ni, nj</td>
<td class="org-left">-</td>
</tr>


<tr>
<td class="org-left">hpnd</td>
<td class="org-left">depth of ponds per cell</td>
<td class="org-left">ncat, ni, nj</td>
<td class="org-left">-</td>
</tr>


<tr>
<td class="org-left">ipnd</td>
<td class="org-left">other pond characteristics</td>
<td class="org-left">ncat, ni, nj</td>
<td class="org-left">-</td>
</tr>


<tr>
<td class="org-left">dhs</td>
<td class="org-left">-</td>
<td class="org-left">ncat, ni, nj</td>
<td class="org-left">-</td>
</tr>


<tr>
<td class="org-left">ffrac</td>
<td class="org-left">-</td>
<td class="org-left">ncat, ni, nj</td>
<td class="org-left">-</td>
</tr>


<tr>
<td class="org-left">fbrs</td>
<td class="org-left">-</td>
<td class="org-left">ncat, ni, nj</td>
<td class="org-left">-</td>
</tr>


<tr>
<td class="org-left">sice00[1-5]</td>
<td class="org-left">-</td>
<td class="org-left">ncat, ni, nj</td>
<td class="org-left">-</td>
</tr>


<tr>
<td class="org-left">qice00[1-5]</td>
<td class="org-left">-</td>
<td class="org-left">ncat, ni, nj</td>
<td class="org-left">-</td>
</tr>


<tr>
<td class="org-left">qsno00[1-5]</td>
<td class="org-left">-</td>
<td class="org-left">ncat, ni, nj</td>
<td class="org-left">-</td>
</tr>
</tbody>
</table>


<a id="orge2dbf35"></a>

# Climatological Datasets

The following external datasets are relevant to stand-alone [CICE6](#org28fbdde). Below is
given information that is helpful to reference for this purpose.


<a id="orgef674bc"></a>

## COREv2

[COREv2](https://data1.gfdl.noaa.gov/nomads/forms/core/COREv2/CIAF_v2.html) is a high quality, but now aging, global one degree and one-hour
atmospheric climate dataset with from years 1948-2009. I have downloaded the
following files: 

> -rw-r&#x2013;r&#x2013;  1 dpath2o  staff   102M 24 Oct  2012 ncar<sub>precip.1948</sub>-2009.23OCT2012.nc
> -rw-r&#x2013;r&#x2013;  1 dpath2o  staff   3.0G 24 Oct  2012 ncar<sub>rad.1948</sub>-2009.23OCT2012.nc
> -rw-r&#x2013;r&#x2013;  1 dpath2o  staff   6.1G 24 Oct  2012 q<sub>10.1948</sub>-2009.23OCT2012.nc
> -rw-r&#x2013;r&#x2013;  1 dpath2o  staff   509K 18 Jun  2009 runoff.15JUNE2009.nc
> -rw-r&#x2013;r&#x2013;  1 dpath2o  staff   386M 26 Jul  2016 runoff.daitren.iaf.20120419.nc
> -rw-r&#x2013;r&#x2013;  1 dpath2o  staff   6.1G 24 Oct  2012 slp.1948-2009.23OCT2012.nc
> -rw-r&#x2013;r&#x2013;  1 dpath2o  staff   6.1G 24 Oct  2012 u<sub>10.1948</sub>-2009.23OCT2012.nc
> -rw-r&#x2013;r&#x2013;  1 dpath2o  staff   6.1G 24 Oct  2012 v<sub>10.1948</sub>-2009.23OCT2012.nc

Header of u10 file:

> netcdf u<sub>10.1948</sub>-2009.23OCT2012 {
> dimensions:
> 	LAT = 94 ;
> 	LON = 192 ;
> 	TIME = UNLIMITED ; // (90520 currently)
> 	bnds = 2 ;
> variables:
> 	double LAT(LAT) ;
> 		LAT:units = "degrees<sub>north</sub>" ;
> 		LAT:point<sub>spacing</sub> = "uneven" ;
> 		LAT:axis = "Y" ;
> 	double LON(LON) ;
> 		LON:units = "degrees<sub>east</sub>" ;
> 		LON:modulo = 360. ;
> 		LON:point<sub>spacing</sub> = "even" ;
> 		LON:axis = "X" ;
> 	double TIME(TIME) ;
> 		TIME:units = "days since 1948-01-01 00:00:00" ;
> 		TIME:axis = "T" ;
> 		TIME:bounds = "TIME<sub>bnds</sub>" ;
> 		TIME:time<sub>origin</sub> = "1-JAN-1948" ;
> 		TIME:calendar = "NOLEAP" ;
> 	double TIME<sub>bnds</sub>(TIME, bnds) ;
> 	float U<sub>10</sub><sub>MOD</sub>(TIME, LAT, LON) ;
> 		U<sub>10</sub><sub>MOD</sub>:missing<sub>value</sub> = -1.e+34f ;
> 		U<sub>10</sub><sub>MOD</sub>:<sub>FillValue</sub> = -1.e+34f ;
> 		U<sub>10</sub><sub>MOD</sub>:long<sub>name</sub> = "10m U Wind" ;
> 		U<sub>10</sub><sub>MOD</sub>:units = "m/s" ;
> 
> // global attributes:
> 		:history = "Thu Sep 27 07:10:36 2012: ncatted -O -a bounds,LAT,d,c,LAT<sub>bnds</sub> u<sub>10</sub><sub>mod.1948</sub>-2009.nomads.nc\n",
> 			"Thu Sep 27 07:10:08 2012: ncks -x -v LAT<sub>bnds</sub> u<sub>10</sub><sub>mod.1948</sub>-2009.new.nc u<sub>10</sub><sub>mod.1948</sub>-2009.nomads.nc\n",
> 			"FERRET V6.725   18-Sep-12" ;
> 		:Conventions = "CF-1.0" ;
> 		:NCO = "4.0.3" ;
> }


<a id="orgc8939a8"></a>

## C-GLORS

The CMCC Global Ocean Physical Reanalysis System ([C-GLORS](http://c-glors.cmcc.it/index/index.html)) is used  to simulate
the state of the ocean in the last decades. It consists of a variational data
assimilation system (OceanVar), capable of assimilating all in-situ observations
along with altimetry data, and a forecast step performed by the ocean model NEMO
coupled with the LIM2 sea-ice model. 

I have downloaded the MLD files for 2010-2020 and here is the header
information:

> netcdf CGv7<sub>2010</sub><sub>MLD</sub> {
> dimensions:
> 	y = 1050 ;
> 	x = 1442 ;
> 	time<sub>counter</sub> = UNLIMITED ; // (365 currently)
> 	axis<sub>nbounds</sub> = 2 ;
> variables:
> 	float nav<sub>lat</sub>(y, x) ;
> 		nav<sub>lat</sub>:standard<sub>name</sub> = "latitude" ;
> 		nav<sub>lat</sub>:long<sub>name</sub> = "Latitude" ;
> 		nav<sub>lat</sub>:units = "degrees<sub>north</sub>" ;
> 		nav<sub>lat</sub>:nav<sub>model</sub> = "grid<sub>T</sub>" ;
> 	float nav<sub>lon</sub>(y, x) ;
> 		nav<sub>lon</sub>:standard<sub>name</sub> = "longitude" ;
> 		nav<sub>lon</sub>:long<sub>name</sub> = "Longitude" ;
> 		nav<sub>lon</sub>:units = "degrees<sub>east</sub>" ;
> 		nav<sub>lon</sub>:nav<sub>model</sub> = "grid<sub>T</sub>" ;
> 	float somxl010(time<sub>counter</sub>, y, x) ;
> 		somxl010:standard<sub>name</sub> = "ocean<sub>mixed</sub><sub>layer</sub><sub>thickness</sub><sub>defined</sub><sub>by</sub><sub>sigma</sub><sub>theta</sub>" ;
> 		somxl010:long<sub>name</sub> = "Mixed Layer Depth (dsigma = 0.01 wrt 10m)" ;
> 		somxl010:units = "m" ;
> 		somxl010:online<sub>operation</sub> = "average" ;
> 		somxl010:interval<sub>operation</sub> = "1080 s" ;
> 		somxl010:interval<sub>write</sub> = "1 d" ;
> 		somxl010:cell<sub>methods</sub> = "time: mean (interval: 1080 s)" ;
> 		somxl010:<sub>FillValue</sub> = 1.e+20f ;
> 		somxl010:missing<sub>value</sub> = 1.e+20f ;
> 		somxl010:coordinates = "time<sub>centered</sub> nav<sub>lon</sub> nav<sub>lat</sub>" ;
> 	double time<sub>centered</sub>(time<sub>counter</sub>) ;
> 		time<sub>centered</sub>:standard<sub>name</sub> = "time" ;
> 		time<sub>centered</sub>:long<sub>name</sub> = "Time axis" ;
> 		time<sub>centered</sub>:calendar = "gregorian" ;
> 		time<sub>centered</sub>:units = "seconds since 1950-01-01 00:00:00" ;
> 		time<sub>centered</sub>:time<sub>origin</sub> = "1950-01-01 00:00:00" ;
> 		time<sub>centered</sub>:bounds = "time<sub>centered</sub><sub>bounds</sub>" ;
> 	double time<sub>centered</sub><sub>bounds</sub>(time<sub>counter</sub>, axis<sub>nbounds</sub>) ;
> 	double time<sub>counter</sub>(time<sub>counter</sub>) ;
> 		time<sub>counter</sub>:axis = "T" ;
> 		time<sub>counter</sub>:standard<sub>name</sub> = "time" ;
> 		time<sub>counter</sub>:long<sub>name</sub> = "Time axis" ;
> 		time<sub>counter</sub>:calendar = "gregorian" ;
> 		time<sub>counter</sub>:units = "seconds since 1950-01-01 00:00:00" ;
> 		time<sub>counter</sub>:time<sub>origin</sub> = "1950-01-01 00:00:00" ;
> 		time<sub>counter</sub>:bounds = "time<sub>counter</sub><sub>bounds</sub>" ;
> 	double time<sub>counter</sub><sub>bounds</sub>(time<sub>counter</sub>, axis<sub>nbounds</sub>) ;
> 
> // global attributes:
> 		:name = "NEMO<sub>1d</sub><sub>20091227</sub><sub>20100106</sub>" ;
> 		:description = "ocean T grid variables" ;
> 		:title = "ocean T grid variables" ;
> 		:Conventions = "CF-1.5" ;
> 		:production = "An IPSL model" ;
> 		:timeStamp = "2016-Aug-14 13:15:38 CEST" ;
> }


<a id="org150c1e8"></a>

## ERA5

$1/4$ degree spatial and 1 hour temporal resolutions. A thorough description of
all the fields in ERA5 is given [here](https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview) 

[COSIMA GitHub Working Group/Issue](https://github.com/COSIMA/access-om2/issues/242). For a comprehensive listing on gadi of mapped
ERA5 fields [see Andrew Kiss's](https://github.com/COSIMA/access-om2/issues/242#issuecomment-913908910) helpful table &#x2013; that thread also contains other
helpful for information re. ERA5 setup for CICE6. 

-   [Latitudes and Longitudes (grid) in text (ascii) format](https://rda.ucar.edu/datasets/ds633.0/docs/ERA5.025deg_lats_lons.txt)
-   [Pressure Levels (levels) in text (ascii) format](https://rda.ucar.edu/datasets/ds633.0/docs/ERA5.025deg_std_pres_levels.txt)
-   [ERA5 Atmospheric Surface Analysis](https://rda.ucar.edu/datasets/ds633.0/docs/ds633.0.e5.oper.an.sfc.grib1.table.web.txt)
-   [ERA5 Atmospheric Pressure Level Analysis](https://rda.ucar.edu/datasets/ds633.0/docs/ds633.0.e5.oper.an.pl.grib1.table.web.txt)
-   [ERA5 Accummulated Fields](https://rda.ucar.edu/datasets/ds633.0/docs/ds633.0.e5.oper.fc.sfc.accumu.grib1.table.web.txt)
-   [ERA5 Instantaneous Fields](https://rda.ucar.edu/datasets/ds633.0/docs/ds633.0.e5.oper.fc.sfc.instan.grib1.table.web.txt)

The following ERA5 fields have been downloaded from gadi to my [local repository](file:///Volumes/ioa02/reanalysis/ERA5) and
organised in same year sub-directory structure for the years 2010-2019.


<a id="org25f314e"></a>

### [10u](file:///Volumes/ioa02/reanalysis/ERA5/10u)

Ten-metre zonal surface winds

> netcdf \\10u<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131 {
> dimensions:
> 	longitude = 1440 ;
> 	latitude = 721 ;
> 	time = 744 ;
> variables:
> 	float longitude(longitude) ;
> 		longitude:units = "degrees<sub>east</sub>" ;
> 		longitude:long<sub>name</sub> = "longitude" ;
> 	float latitude(latitude) ;
> 		latitude:units = "degrees<sub>north</sub>" ;
> 		latitude:long<sub>name</sub> = "latitude" ;
> 	int time(time) ;
> 		time:units = "hours since 1900-01-01 00:00:00.0" ;
> 		time:long<sub>name</sub> = "time" ;
> 		time:calendar = "gregorian" ;
> 	short u10(time, latitude, longitude) ;
> 		u10:scale<sub>factor</sub> = 0.00113594702823864 ;
> 		u10:add<sub>offset</sub> = 7.62408107433744 ;
> 		u10:<sub>FillValue</sub> = -32767s ;
> 		u10:missing<sub>value</sub> = -32767s ;
> 		u10:units = "m s\*\*-1" ;
> 		u10:long<sub>name</sub> = "10 metre U wind component" ;
> 
> // global attributes:
> 		:Conventions = "CF-1.6" ;
> 		:history = "2020-09-28 12:42:00 UTC+1000 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: mv *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/10u/2010*.10u<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp *g/data/rt52/era5/single-levels/reanalysis/10u/2010/10u<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc\n",
> 			"2020-09-28 12:40:58 UTC+1000 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: nccopy -k4 -d5 -s /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/10u/2010/10u<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc /g/data/rt52/admin/incoming/era5/single-levels/reanalysis/10u/2010*.10u<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp\n",
> 			"2020-09-28 12:36:04 UTC+1000 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: curl &#x2013;connect-timeout 20 &#x2013;show-error &#x2013;silent &#x2013;max-time 36000 -o /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/10u/2010/10u<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc <http://136.156.133.25/cache-compute-0008/cache/data0/adaptor.mars.internal-1601260358.5388222-19803-23-f152e830-9d4e-468b-a146-f2132486120c.nc\n>",
> 			"2020-09-28 02:33:37 GMT by grib<sub>to</sub><sub>netcdf</sub>-2.16.0: *opt/ecmwf/eccodes/bin/grib<sub>to</sub><sub>netcdf</sub> -S param -o /cache/data0/adaptor.mars.internal-1601260358.5388222-19803-23-f152e830-9d4e-468b-a146-f2132486120c.nc /cache/tmp/f152e830-9d4e-468b-a146-f2132486120c-adaptor.mars.internal-1601260358.5393667-19803-6-tmp.grib" ;
> 		:license = "Licence to use Copernicus Products: <https://apps.ecmwf.int/datasets/licences/copernicus>*" ;
> 		:summary = "ERA5 is the fifth generation ECMWF atmospheric reanalysis of the global climate. This file is part of the ERA5 replica hosted at NCI Australia. For more information please see <http://dx.doi.org/10.25914/5f48874388857>" ;
> 		:title = "ERA5 single-levels reanalysis 10m<sub>u</sub><sub>component</sub><sub>of</sub><sub>wind</sub> 20100101-20100131" ;
> }


<a id="org703b36c"></a>

### [10v](file:///Volumes/ioa02/reanalysis/ERA5/10v)

Ten-metre meridional surface winds 

> netcdf \\10v<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131 {
> dimensions:
> 	longitude = 1440 ;
> 	latitude = 721 ;
> 	time = 744 ;
> variables:
> 	float longitude(longitude) ;
> 		longitude:units = "degrees<sub>east</sub>" ;
> 		longitude:long<sub>name</sub> = "longitude" ;
> 	float latitude(latitude) ;
> 		latitude:units = "degrees<sub>north</sub>" ;
> 		latitude:long<sub>name</sub> = "latitude" ;
> 	int time(time) ;
> 		time:units = "hours since 1900-01-01 00:00:00.0" ;
> 		time:long<sub>name</sub> = "time" ;
> 		time:calendar = "gregorian" ;
> 	short v10(time, latitude, longitude) ;
> 		v10:scale<sub>factor</sub> = 0.00113620990606891 ;
> 		v10:add<sub>offset</sub> = -6.70077104684757 ;
> 		v10:<sub>FillValue</sub> = -32767s ;
> 		v10:missing<sub>value</sub> = -32767s ;
> 		v10:units = "m s\*\*-1" ;
> 		v10:long<sub>name</sub> = "10 metre V wind component" ;
> 
> // global attributes:
> 		:Conventions = "CF-1.6" ;
> 		:history = "2020-09-28 12:40:15 UTC+1000 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: mv *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/10v/2010*.10v<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp *g/data/rt52/era5/single-levels/reanalysis/10v/2010/10v<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc\n",
> 			"2020-09-28 12:39:12 UTC+1000 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: nccopy -k4 -d5 -s /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/10v/2010/10v<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc /g/data/rt52/admin/incoming/era5/single-levels/reanalysis/10v/2010*.10v<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp\n",
> 			"2020-09-28 12:35:24 UTC+1000 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: curl &#x2013;connect-timeout 20 &#x2013;show-error &#x2013;silent &#x2013;max-time 36000 -o /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/10v/2010/10v<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc <http://136.156.132.153/cache-compute-0002/cache/data1/adaptor.mars.internal-1601260363.3011663-8591-19-7ca36757-cc7c-4890-9fc2-43150cac7864.nc\n>",
> 			"2020-09-28 02:33:37 GMT by grib<sub>to</sub><sub>netcdf</sub>-2.16.0: *opt/ecmwf/eccodes/bin/grib<sub>to</sub><sub>netcdf</sub> -S param -o /cache/data1/adaptor.mars.internal-1601260363.3011663-8591-19-7ca36757-cc7c-4890-9fc2-43150cac7864.nc /cache/tmp/7ca36757-cc7c-4890-9fc2-43150cac7864-adaptor.mars.internal-1601260363.3017175-8591-7-tmp.grib" ;
> 		:license = "Licence to use Copernicus Products: <https://apps.ecmwf.int/datasets/licences/copernicus>*" ;
> 		:summary = "ERA5 is the fifth generation ECMWF atmospheric reanalysis of the global climate. This file is part of the ERA5 replica hosted at NCI Australia. For more information please see <http://dx.doi.org/10.25914/5f48874388857>" ;
> 		:title = "ERA5 single-levels reanalysis 10m<sub>v</sub><sub>component</sub><sub>of</sub><sub>wind</sub> 20100101-20100131" ;
> }


<a id="org6661811"></a>

### [2d](file:///Volumes/ioa02/reanalysis/ERA5/2d)

Two-metre dew-point temperature

> netcdf \\2d<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131 {
> dimensions:
> 	longitude = 1440 ;
> 	latitude = 721 ;
> 	time = 744 ;
> variables:
> 	float longitude(longitude) ;
> 		longitude:units = "degrees<sub>east</sub>" ;
> 		longitude:long<sub>name</sub> = "longitude" ;
> 	float latitude(latitude) ;
> 		latitude:units = "degrees<sub>north</sub>" ;
> 		latitude:long<sub>name</sub> = "latitude" ;
> 	int time(time) ;
> 		time:units = "hours since 1900-01-01 00:00:00.0" ;
> 		time:long<sub>name</sub> = "time" ;
> 		time:calendar = "gregorian" ;
> 	short d2m(time, latitude, longitude) ;
> 		d2m:scale<sub>factor</sub> = 0.00140529551830638 ;
> 		d2m:add<sub>offset</sub> = 257.339492054389 ;
> 		d2m:<sub>FillValue</sub> = -32767s ;
> 		d2m:missing<sub>value</sub> = -32767s ;
> 		d2m:units = "K" ;
> 		d2m:long<sub>name</sub> = "2 metre dewpoint temperature" ;
> 
> // global attributes:
> 		:Conventions = "CF-1.6" ;
> 		:history = "2020-09-28 12:43:43 UTC+1000 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: mv *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/2d/2010*.2d<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp *g/data/rt52/era5/single-levels/reanalysis/2d/2010/2d<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc\n",
> 			"2020-09-28 12:42:44 UTC+1000 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: nccopy -k4 -d5 -s /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/2d/2010/2d<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc /g/data/rt52/admin/incoming/era5/single-levels/reanalysis/2d/2010*.2d<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp\n",
> 			"2020-09-28 12:37:40 UTC+1000 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: curl &#x2013;connect-timeout 20 &#x2013;show-error &#x2013;silent &#x2013;max-time 36000 -o /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/2d/2010/2d<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc <http://136.156.132.201/cache-compute-0004/cache/data5/adaptor.mars.internal-1601260378.023689-11063-21-0f2ccb24-24bb-4785-bfe1-620949bb0f34.nc\n>",
> 			"2020-09-28 02:33:43 GMT by grib<sub>to</sub><sub>netcdf</sub>-2.16.0: *opt/ecmwf/eccodes/bin/grib<sub>to</sub><sub>netcdf</sub> -S param -o /cache/data5/adaptor.mars.internal-1601260378.023689-11063-21-0f2ccb24-24bb-4785-bfe1-620949bb0f34.nc /cache/tmp/0f2ccb24-24bb-4785-bfe1-620949bb0f34-adaptor.mars.internal-1601260378.0242786-11063-3-tmp.grib" ;
> 		:license = "Licence to use Copernicus Products: <https://apps.ecmwf.int/datasets/licences/copernicus>*" ;
> 		:summary = "ERA5 is the fifth generation ECMWF atmospheric reanalysis of the global climate. This file is part of the ERA5 replica hosted at NCI Australia. For more information please see <http://dx.doi.org/10.25914/5f48874388857>" ;
> 		:title = "ERA5 single-levels reanalysis 2m<sub>dewpoint</sub><sub>temperature</sub> 20100101-20100131" ;
> }


<a id="org493eaea"></a>

### [2t](file:///Volumes/ioa02/reanalysis/ERA5/2t)

Two-metre air temperature

> netcdf \\2t<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131 {
> dimensions:
> 	longitude = 1440 ;
> 	latitude = 721 ;
> 	time = 744 ;
> variables:
> 	float longitude(longitude) ;
> 		longitude:units = "degrees<sub>east</sub>" ;
> 		longitude:long<sub>name</sub> = "longitude" ;
> 	float latitude(latitude) ;
> 		latitude:units = "degrees<sub>north</sub>" ;
> 		latitude:long<sub>name</sub> = "latitude" ;
> 	int time(time) ;
> 		time:units = "hours since 1900-01-01 00:00:00.0" ;
> 		time:long<sub>name</sub> = "time" ;
> 		time:calendar = "gregorian" ;
> 	short t2m(time, latitude, longitude) ;
> 		t2m:scale<sub>factor</sub> = 0.00163849236925778 ;
> 		t2m:add<sub>offset</sub> = 268.255307157624 ;
> 		t2m:<sub>FillValue</sub> = -32767s ;
> 		t2m:missing<sub>value</sub> = -32767s ;
> 		t2m:units = "K" ;
> 		t2m:long<sub>name</sub> = "2 metre temperature" ;
> 
> // global attributes:
> 		:Conventions = "CF-1.6" ;
> 		:history = "2020-09-28 12:42:00 UTC+1000 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: mv *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/2t/2010*.2t<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp *g/data/rt52/era5/single-levels/reanalysis/2t/2010/2t<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc\n",
> 			"2020-09-28 12:41:06 UTC+1000 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: nccopy -k4 -d5 -s /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/2t/2010/2t<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc /g/data/rt52/admin/incoming/era5/single-levels/reanalysis/2t/2010*.2t<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp\n",
> 			"2020-09-28 12:37:30 UTC+1000 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: curl &#x2013;connect-timeout 20 &#x2013;show-error &#x2013;silent &#x2013;max-time 36000 -o /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/2t/2010/2t<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc <http://136.156.133.36/cache-compute-0010/cache/data6/adaptor.mars.internal-1601260376.2580519-20638-15-e6b74d84-3516-47c1-900e-48a2e567b012.nc\n>",
> 			"2020-09-28 02:33:41 GMT by grib<sub>to</sub><sub>netcdf</sub>-2.16.0: *opt/ecmwf/eccodes/bin/grib<sub>to</sub><sub>netcdf</sub> -S param -o /cache/data6/adaptor.mars.internal-1601260376.2580519-20638-15-e6b74d84-3516-47c1-900e-48a2e567b012.nc /cache/tmp/e6b74d84-3516-47c1-900e-48a2e567b012-adaptor.mars.internal-1601260376.258729-20638-3-tmp.grib" ;
> 		:license = "Licence to use Copernicus Products: <https://apps.ecmwf.int/datasets/licences/copernicus>*" ;
> 		:summary = "ERA5 is the fifth generation ECMWF atmospheric reanalysis of the global climate. This file is part of the ERA5 replica hosted at NCI Australia. For more information please see <http://dx.doi.org/10.25914/5f48874388857>" ;
> 		:title = "ERA5 single-levels reanalysis 2m<sub>temperature</sub> 20100101-20100131" ;
> }


<a id="org54397e7"></a>

### [metss](file:///Volumes/ioa02/reanalysis/ERA5/metss)

Mean eastward turbulent surface stress

> netcdf metss<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131 {
> dimensions:
> 	longitude = 1440 ;
> 	latitude = 721 ;
> 	time = 744 ;
> variables:
> 	float longitude(longitude) ;
> 		longitude:units = "degrees<sub>east</sub>" ;
> 		longitude:long<sub>name</sub> = "longitude" ;
> 	float latitude(latitude) ;
> 		latitude:units = "degrees<sub>north</sub>" ;
> 		latitude:long<sub>name</sub> = "latitude" ;
> 	int time(time) ;
> 		time:units = "hours since 1900-01-01 00:00:00.0" ;
> 		time:long<sub>name</sub> = "time" ;
> 		time:calendar = "gregorian" ;
> 	short metss(time, latitude, longitude) ;
> 		metss:scale<sub>factor</sub> = 0.000366911687951802 ;
> 		metss:add<sub>offset</sub> = -0.80675398578172 ;
> 		metss:<sub>FillValue</sub> = -32767s ;
> 		metss:missing<sub>value</sub> = -32767s ;
> 		metss:units = "N m\*\*-2" ;
> 		metss:long<sub>name</sub> = "Mean eastward turbulent surface stress" ;
> 
> // global attributes:
> 		:Conventions = "CF-1.6" ;
> 		:history = "2020-10-06 09:10:44 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: mv *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/metss/2010*.metss<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp *g/data/rt52/era5/single-levels/reanalysis/metss/2010/metss<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc\n",
> 			"2020-10-06 09:09:59 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: nccopy -k4 -d5 -s /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/metss/2010/metss<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc /g/data/rt52/admin/incoming/era5/single-levels/reanalysis/metss/2010*.metss<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp\n",
> 			"2020-10-06 09:04:24 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: curl &#x2013;connect-timeout 20 &#x2013;show-error &#x2013;silent &#x2013;max-time 36000 -o /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/metss/2010/metss<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc <http://136.156.132.110/cache-compute-0001/cache/data2/adaptor.mars.internal-1601935170.1456351-18982-33-a6b9dca0-74b0-4ce5-8017-3660c3fbe2d5.nc\n>",
> 			"2020-10-05 22:00:31 GMT by grib<sub>to</sub><sub>netcdf</sub>-2.16.0: *opt/ecmwf/eccodes/bin/grib<sub>to</sub><sub>netcdf</sub> -S param -o /cache/data2/adaptor.mars.internal-1601935170.1456351-18982-33-a6b9dca0-74b0-4ce5-8017-3660c3fbe2d5.nc /cache/tmp/a6b9dca0-74b0-4ce5-8017-3660c3fbe2d5-adaptor.mars.internal-1601935170.1462164-18982-8-tmp.grib" ;
> 		:license = "Licence to use Copernicus Products: <https://apps.ecmwf.int/datasets/licences/copernicus>*" ;
> 		:summary = "ERA5 is the fifth generation ECMWF atmospheric reanalysis of the global climate. This file is part of the ERA5 replica hosted at NCI Australia. For more information please see <http://dx.doi.org/10.25914/5f48874388857>" ;
> 		:title = "ERA5 single-levels reanalysis mean<sub>eastward</sub><sub>turbulent</sub><sub>surface</sub><sub>stress</sub> 20100101-20100131" ;
> }


<a id="org5926abd"></a>

### [mntss](file:///Volumes/ioa02/reanalysis/ERA5/mntss)

Mean northward turbulent surface stress

> netcdf mntss<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131 {
> dimensions:
> 	longitude = 1440 ;
> 	latitude = 721 ;
> 	time = 744 ;
> variables:
> 	float longitude(longitude) ;
> 		longitude:units = "degrees<sub>east</sub>" ;
> 		longitude:long<sub>name</sub> = "longitude" ;
> 	float latitude(latitude) ;
> 		latitude:units = "degrees<sub>north</sub>" ;
> 		latitude:long<sub>name</sub> = "latitude" ;
> 	int time(time) ;
> 		time:units = "hours since 1900-01-01 00:00:00.0" ;
> 		time:long<sub>name</sub> = "time" ;
> 		time:calendar = "gregorian" ;
> 	short mntss(time, latitude, longitude) ;
> 		mntss:scale<sub>factor</sub> = 0.000288605498999811 ;
> 		mntss:add<sub>offset</sub> = 0.0627023578645138 ;
> 		mntss:<sub>FillValue</sub> = -32767s ;
> 		mntss:missing<sub>value</sub> = -32767s ;
> 		mntss:units = "N m\*\*-2" ;
> 		mntss:long<sub>name</sub> = "Mean northward turbulent surface stress" ;
> 
> // global attributes:
> 		:Conventions = "CF-1.6" ;
> 		:history = "2020-10-06 09:07:45 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: mv *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/mntss/2010*.mntss<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp *g/data/rt52/era5/single-levels/reanalysis/mntss/2010/mntss<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc\n",
> 			"2020-10-06 09:06:59 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: nccopy -k4 -d5 -s /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/mntss/2010/mntss<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc /g/data/rt52/admin/incoming/era5/single-levels/reanalysis/mntss/2010*.mntss<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp\n",
> 			"2020-10-06 09:03:25 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: curl &#x2013;connect-timeout 20 &#x2013;show-error &#x2013;silent &#x2013;max-time 36000 -o /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/mntss/2010/mntss<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc <http://136.156.133.37/cache-compute-0011/cache/data9/adaptor.mars.internal-1601935174.829866-29967-28-db022499-88a8-4c1e-9b82-746968321737.nc\n>",
> 			"2020-10-05 22:00:30 GMT by grib<sub>to</sub><sub>netcdf</sub>-2.16.0: *opt/ecmwf/eccodes/bin/grib<sub>to</sub><sub>netcdf</sub> -S param -o /cache/data9/adaptor.mars.internal-1601935174.829866-29967-28-db022499-88a8-4c1e-9b82-746968321737.nc /cache/tmp/db022499-88a8-4c1e-9b82-746968321737-adaptor.mars.internal-1601935174.8304174-29967-7-tmp.grib" ;
> 		:license = "Licence to use Copernicus Products: <https://apps.ecmwf.int/datasets/licences/copernicus>*" ;
> 		:summary = "ERA5 is the fifth generation ECMWF atmospheric reanalysis of the global climate. This file is part of the ERA5 replica hosted at NCI Australia. For more information please see <http://dx.doi.org/10.25914/5f48874388857>" ;
> 		:title = "ERA5 single-levels reanalysis mean<sub>northward</sub><sub>turbulent</sub><sub>surface</sub><sub>stress</sub> 20100101-20100131" ;
> }


<a id="orgd5b196f"></a>

### [mror](file:///Volumes/ioa02/reanalysis/ERA5/mror)

Mean runoff rate

> netcdf mror<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131 {
> dimensions:
> 	longitude = 1440 ;
> 	latitude = 721 ;
> 	time = 744 ;
> variables:
> 	float longitude(longitude) ;
> 		longitude:units = "degrees<sub>east</sub>" ;
> 		longitude:long<sub>name</sub> = "longitude" ;
> 	float latitude(latitude) ;
> 		latitude:units = "degrees<sub>north</sub>" ;
> 		latitude:long<sub>name</sub> = "latitude" ;
> 	int time(time) ;
> 		time:units = "hours since 1900-01-01 00:00:00.0" ;
> 		time:long<sub>name</sub> = "time" ;
> 		time:calendar = "gregorian" ;
> 	short mror(time, latitude, longitude) ;
> 		mror:scale<sub>factor</sub> = 3.5196145727257e-07 ;
> 		mror:add<sub>offset</sub> = 0.011532369108993 ;
> 		mror:<sub>FillValue</sub> = -32767s ;
> 		mror:missing<sub>value</sub> = -32767s ;
> 		mror:units = "kg m\*\*-2 s\*\*-1" ;
> 		mror:long<sub>name</sub> = "Mean runoff rate" ;
> 
> // global attributes:
> 		:Conventions = "CF-1.6" ;
> 		:history = "2020-10-06 09:09:27 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: mv *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/mror/2010*.mror<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp *g/data/rt52/era5/single-levels/reanalysis/mror/2010/mror<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc\n",
> 			"2020-10-06 09:09:11 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: nccopy -k4 -d5 -s /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/mror/2010/mror<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc /g/data/rt52/admin/incoming/era5/single-levels/reanalysis/mror/2010*.mror<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp\n",
> 			"2020-10-06 09:04:54 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: curl &#x2013;connect-timeout 20 &#x2013;show-error &#x2013;silent &#x2013;max-time 36000 -o /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/mror/2010/mror<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc <http://136.156.132.198/cache-compute-0003/cache/data2/adaptor.mars.internal-1601935221.2261062-30150-20-349424a8-a7ab-4eae-8d46-b300f17ce397.nc\n>",
> 			"2020-10-05 22:01:06 GMT by grib<sub>to</sub><sub>netcdf</sub>-2.16.0: *opt/ecmwf/eccodes/bin/grib<sub>to</sub><sub>netcdf</sub> -S param -o /cache/data2/adaptor.mars.internal-1601935221.2261062-30150-20-349424a8-a7ab-4eae-8d46-b300f17ce397.nc /cache/tmp/349424a8-a7ab-4eae-8d46-b300f17ce397-adaptor.mars.internal-1601935221.2266717-30150-7-tmp.grib" ;
> 		:license = "Licence to use Copernicus Products: <https://apps.ecmwf.int/datasets/licences/copernicus>*" ;
> 		:summary = "ERA5 is the fifth generation ECMWF atmospheric reanalysis of the global climate. This file is part of the ERA5 replica hosted at NCI Australia. For more information please see <http://dx.doi.org/10.25914/5f48874388857>" ;
> 		:title = "ERA5 single-levels reanalysis mean<sub>runoff</sub><sub>rate</sub> 20100101-20100131" ;
> }


<a id="orgd6a9b35"></a>

### [msdrswrf](file:///Volumes/ioa02/reanalysis/ERA5/msdrswrf)

Mean surface direct short-wave radiation flux

> netcdf msdrswrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131 {
> dimensions:
> 	longitude = 1440 ;
> 	latitude = 721 ;
> 	time = 744 ;
> variables:
> 	float longitude(longitude) ;
> 		longitude:units = "degrees<sub>east</sub>" ;
> 		longitude:long<sub>name</sub> = "longitude" ;
> 	float latitude(latitude) ;
> 		latitude:units = "degrees<sub>north</sub>" ;
> 		latitude:long<sub>name</sub> = "latitude" ;
> 	int time(time) ;
> 		time:units = "hours since 1900-01-01 00:00:00.0" ;
> 		time:long<sub>name</sub> = "time" ;
> 		time:calendar = "gregorian" ;
> 	short msdrswrf(time, latitude, longitude) ;
> 		msdrswrf:scale<sub>factor</sub> = 0.0188225779378328 ;
> 		msdrswrf:add<sub>offset</sub> = 616.740588711031 ;
> 		msdrswrf:<sub>FillValue</sub> = -32767s ;
> 		msdrswrf:missing<sub>value</sub> = -32767s ;
> 		msdrswrf:units = "W m\*\*-2" ;
> 		msdrswrf:long<sub>name</sub> = "Mean surface direct short-wave radiation flux" ;
> 
> // global attributes:
> 		:Conventions = "CF-1.6" ;
> 		:history = "2020-10-06 09:14:19 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: mv *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msdrswrf/2010*.msdrswrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp *g/data/rt52/era5/single-levels/reanalysis/msdrswrf/2010/msdrswrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc\n",
> 			"2020-10-06 09:13:44 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: nccopy -k4 -d5 -s /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msdrswrf/2010/msdrswrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc /g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msdrswrf/2010*.msdrswrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp\n",
> 			"2020-10-06 09:08:14 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: curl &#x2013;connect-timeout 20 &#x2013;show-error &#x2013;silent &#x2013;max-time 36000 -o /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msdrswrf/2010/msdrswrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc <http://136.156.133.36/cache-compute-0010/cache/data1/adaptor.mars.internal-1601935408.957139-15597-14-b57f2409-d32e-448d-bf0e-ee6d5e9e4503.nc\n>",
> 			"2020-10-05 22:04:10 GMT by grib<sub>to</sub><sub>netcdf</sub>-2.16.0: *opt/ecmwf/eccodes/bin/grib<sub>to</sub><sub>netcdf</sub> -S param -o /cache/data1/adaptor.mars.internal-1601935408.957139-15597-14-b57f2409-d32e-448d-bf0e-ee6d5e9e4503.nc /cache/tmp/b57f2409-d32e-448d-bf0e-ee6d5e9e4503-adaptor.mars.internal-1601935408.9578357-15597-6-tmp.grib" ;
> 		:license = "Licence to use Copernicus Products: <https://apps.ecmwf.int/datasets/licences/copernicus>*" ;
> 		:summary = "ERA5 is the fifth generation ECMWF atmospheric reanalysis of the global climate. This file is part of the ERA5 replica hosted at NCI Australia. For more information please see <http://dx.doi.org/10.25914/5f48874388857>" ;
> 		:title = "ERA5 single-levels reanalysis mean<sub>surface</sub><sub>direct</sub><sub>short</sub><sub>wave</sub><sub>radiation</sub><sub>flux</sub> 20100101-20100131" ;
> }


<a id="orgdda916b"></a>

### [msdrswrfcs](file:///Volumes/ioa02/reanalysis/ERA5/msdrswrfcs)

Mean surface direct short-wave radiation flux, clear sky

> netcdf msdrswrfcs<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131 {
> dimensions:
> 	longitude = 1440 ;
> 	latitude = 721 ;
> 	time = 744 ;
> variables:
> 	float longitude(longitude) ;
> 		longitude:units = "degrees<sub>east</sub>" ;
> 		longitude:long<sub>name</sub> = "longitude" ;
> 	float latitude(latitude) ;
> 		latitude:units = "degrees<sub>north</sub>" ;
> 		latitude:long<sub>name</sub> = "latitude" ;
> 	int time(time) ;
> 		time:units = "hours since 1900-01-01 00:00:00.0" ;
> 		time:long<sub>name</sub> = "time" ;
> 		time:calendar = "gregorian" ;
> 	short msdrswrfcs(time, latitude, longitude) ;
> 		msdrswrfcs:scale<sub>factor</sub> = 0.0182069529855187 ;
> 		msdrswrfcs:add<sub>offset</sub> = 596.569021523507 ;
> 		msdrswrfcs:<sub>FillValue</sub> = -32767s ;
> 		msdrswrfcs:missing<sub>value</sub> = -32767s ;
> 		msdrswrfcs:units = "W m\*\*-2" ;
> 		msdrswrfcs:long<sub>name</sub> = "Mean surface direct short-wave radiation flux, clear sky" ;
> 
> // global attributes:
> 		:Conventions = "CF-1.6" ;
> 		:history = "2020-10-06 09:12:53 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: mv *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msdrswrfcs/2010*.msdrswrfcs<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp *g/data/rt52/era5/single-levels/reanalysis/msdrswrfcs/2010/msdrswrfcs<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc\n",
> 			"2020-10-06 09:12:21 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: nccopy -k4 -d5 -s /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msdrswrfcs/2010/msdrswrfcs<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc /g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msdrswrfcs/2010*.msdrswrfcs<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp\n",
> 			"2020-10-06 09:08:21 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: curl &#x2013;connect-timeout 20 &#x2013;show-error &#x2013;silent &#x2013;max-time 36000 -o /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msdrswrfcs/2010/msdrswrfcs<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc <http://136.156.132.153/cache-compute-0002/cache/data3/adaptor.mars.internal-1601935426.7869933-12103-12-b44cb567-fe87-4047-bde0-d85c747fde81.nc\n>",
> 			"2020-10-05 22:04:43 GMT by grib<sub>to</sub><sub>netcdf</sub>-2.16.0: *opt/ecmwf/eccodes/bin/grib<sub>to</sub><sub>netcdf</sub> -S param -o /cache/data3/adaptor.mars.internal-1601935426.7869933-12103-12-b44cb567-fe87-4047-bde0-d85c747fde81.nc /cache/tmp/b44cb567-fe87-4047-bde0-d85c747fde81-adaptor.mars.internal-1601935426.787668-12103-5-tmp.grib" ;
> 		:license = "Licence to use Copernicus Products: <https://apps.ecmwf.int/datasets/licences/copernicus>*" ;
> 		:summary = "ERA5 is the fifth generation ECMWF atmospheric reanalysis of the global climate. This file is part of the ERA5 replica hosted at NCI Australia. For more information please see <http://dx.doi.org/10.25914/5f48874388857>" ;
> 		:title = "ERA5 single-levels reanalysis mean<sub>surface</sub><sub>direct</sub><sub>short</sub><sub>wave</sub><sub>radiation</sub><sub>flux</sub><sub>clear</sub><sub>sky</sub> 20100101-20100131" ;
> }


<a id="org4a5a706"></a>

### [msdwlwrf](file:///Volumes/ioa02/reanalysis/ERA5/msdwlwrf)

Mean surface downward long-wave radiation flux

> netcdf msdwlwrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131 {
> dimensions:
> 	longitude = 1440 ;
> 	latitude = 721 ;
> 	time = 744 ;
> variables:
> 	float longitude(longitude) ;
> 		longitude:units = "degrees<sub>east</sub>" ;
> 		longitude:long<sub>name</sub> = "longitude" ;
> 	float latitude(latitude) ;
> 		latitude:units = "degrees<sub>north</sub>" ;
> 		latitude:long<sub>name</sub> = "latitude" ;
> 	int time(time) ;
> 		time:units = "hours since 1900-01-01 00:00:00.0" ;
> 		time:long<sub>name</sub> = "time" ;
> 		time:calendar = "gregorian" ;
> 	short msdwlwrf(time, latitude, longitude) ;
> 		msdwlwrf:scale<sub>factor</sub> = 0.00672298478721479 ;
> 		msdwlwrf:add<sub>offset</sub> = 289.4882766912 ;
> 		msdwlwrf:<sub>FillValue</sub> = -32767s ;
> 		msdwlwrf:missing<sub>value</sub> = -32767s ;
> 		msdwlwrf:units = "W m\*\*-2" ;
> 		msdwlwrf:long<sub>name</sub> = "Mean surface downward long-wave radiation flux" ;
> 
> // global attributes:
> 		:Conventions = "CF-1.6" ;
> 		:history = "2020-10-06 09:09:32 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: mv *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msdwlwrf/2010*.msdwlwrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp *g/data/rt52/era5/single-levels/reanalysis/msdwlwrf/2010/msdwlwrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc\n",
> 			"2020-10-06 09:08:27 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: nccopy -k4 -d5 -s /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msdwlwrf/2010/msdwlwrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc /g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msdwlwrf/2010*.msdwlwrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp\n",
> 			"2020-10-06 09:01:48 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: curl &#x2013;connect-timeout 20 &#x2013;show-error &#x2013;silent &#x2013;max-time 36000 -o /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msdwlwrf/2010/msdwlwrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc <http://136.156.133.46/cache-compute-0015/cache/data3/adaptor.mars.internal-1601935022.9836493-31252-13-a31306dc-69ab-4a06-bc9e-ad88c60a42ca.nc\n>",
> 			"2020-10-05 21:57:45 GMT by grib<sub>to</sub><sub>netcdf</sub>-2.16.0: *opt/ecmwf/eccodes/bin/grib<sub>to</sub><sub>netcdf</sub> -S param -o /cache/data3/adaptor.mars.internal-1601935022.9836493-31252-13-a31306dc-69ab-4a06-bc9e-ad88c60a42ca.nc /cache/tmp/a31306dc-69ab-4a06-bc9e-ad88c60a42ca-adaptor.mars.internal-1601935022.9842257-31252-5-tmp.grib" ;
> 		:license = "Licence to use Copernicus Products: <https://apps.ecmwf.int/datasets/licences/copernicus>*" ;
> 		:summary = "ERA5 is the fifth generation ECMWF atmospheric reanalysis of the global climate. This file is part of the ERA5 replica hosted at NCI Australia. For more information please see <http://dx.doi.org/10.25914/5f48874388857>" ;
> 		:title = "ERA5 single-levels reanalysis mean<sub>surface</sub><sub>downward</sub><sub>long</sub><sub>wave</sub><sub>radiation</sub><sub>flux</sub> 20100101-20100131" ;
> }


<a id="orgae47d0c"></a>

### [msdwlwrfcs](file:///Volumes/ioa02/reanalysis/ERA5/msdwlwrfcs)

Mean surface downward long-wave radiation flux, clear sky

> netcdf msdwlwrfcs<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131 {
> dimensions:
> 	longitude = 1440 ;
> 	latitude = 721 ;
> 	time = 744 ;
> variables:
> 	float longitude(longitude) ;
> 		longitude:units = "degrees<sub>east</sub>" ;
> 		longitude:long<sub>name</sub> = "longitude" ;
> 	float latitude(latitude) ;
> 		latitude:units = "degrees<sub>north</sub>" ;
> 		latitude:long<sub>name</sub> = "latitude" ;
> 	int time(time) ;
> 		time:units = "hours since 1900-01-01 00:00:00.0" ;
> 		time:long<sub>name</sub> = "time" ;
> 		time:calendar = "gregorian" ;
> 	short msdwlwrfcs(time, latitude, longitude) ;
> 		msdwlwrfcs:scale<sub>factor</sub> = 0.00657778122565039 ;
> 		msdwlwrfcs:add<sub>offset</sub> = 287.159476612317 ;
> 		msdwlwrfcs:<sub>FillValue</sub> = -32767s ;
> 		msdwlwrfcs:missing<sub>value</sub> = -32767s ;
> 		msdwlwrfcs:units = "W m\*\*-2" ;
> 		msdwlwrfcs:long<sub>name</sub> = "Mean surface downward long-wave radiation flux, clear sky" ;
> 
> // global attributes:
> 		:Conventions = "CF-1.6" ;
> 		:history = "2020-10-06 09:16:35 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: mv *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msdwlwrfcs/2010*.msdwlwrfcs<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp *g/data/rt52/era5/single-levels/reanalysis/msdwlwrfcs/2010/msdwlwrfcs<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc\n",
> 			"2020-10-06 09:15:39 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: nccopy -k4 -d5 -s /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msdwlwrfcs/2010/msdwlwrfcs<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc /g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msdwlwrfcs/2010*.msdwlwrfcs<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp\n",
> 			"2020-10-06 09:09:11 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: curl &#x2013;connect-timeout 20 &#x2013;show-error &#x2013;silent &#x2013;max-time 36000 -o /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msdwlwrfcs/2010/msdwlwrfcs<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc <http://136.156.132.105/cache-compute-0000/cache/data4/adaptor.mars.internal-1601935467.8587792-11179-31-6f536992-db19-476c-8205-d5ebf4a7e4c2.nc\n>",
> 			"2020-10-05 22:05:11 GMT by grib<sub>to</sub><sub>netcdf</sub>-2.16.0: *opt/ecmwf/eccodes/bin/grib<sub>to</sub><sub>netcdf</sub> -S param -o /cache/data4/adaptor.mars.internal-1601935467.8587792-11179-31-6f536992-db19-476c-8205-d5ebf4a7e4c2.nc /cache/tmp/6f536992-db19-476c-8205-d5ebf4a7e4c2-adaptor.mars.internal-1601935467.859394-11179-9-tmp.grib" ;
> 		:license = "Licence to use Copernicus Products: <https://apps.ecmwf.int/datasets/licences/copernicus>*" ;
> 		:summary = "ERA5 is the fifth generation ECMWF atmospheric reanalysis of the global climate. This file is part of the ERA5 replica hosted at NCI Australia. For more information please see <http://dx.doi.org/10.25914/5f48874388857>" ;
> 		:title = "ERA5 single-levels reanalysis mean<sub>surface</sub><sub>downward</sub><sub>long</sub><sub>wave</sub><sub>radiation</sub><sub>flux</sub><sub>clear</sub><sub>sky</sub> 20100101-20100131" ;
> }


<a id="org939fd6f"></a>

### [msdwswrf](file:///Volumes/ioa02/reanalysis/ERA5/msdwswrf)

Mean surface downward short-wave radiation flux

> netcdf msdwswrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131 {
> dimensions:
> 	longitude = 1440 ;
> 	latitude = 721 ;
> 	time = 744 ;
> variables:
> 	float longitude(longitude) ;
> 		longitude:units = "degrees<sub>east</sub>" ;
> 		longitude:long<sub>name</sub> = "longitude" ;
> 	float latitude(latitude) ;
> 		latitude:units = "degrees<sub>north</sub>" ;
> 		latitude:long<sub>name</sub> = "latitude" ;
> 	int time(time) ;
> 		time:units = "hours since 1900-01-01 00:00:00.0" ;
> 		time:long<sub>name</sub> = "time" ;
> 		time:calendar = "gregorian" ;
> 	short msdwswrf(time, latitude, longitude) ;
> 		msdwswrf:scale<sub>factor</sub> = 0.0201353707292509 ;
> 		msdwswrf:add<sub>offset</sub> = 659.755557314635 ;
> 		msdwswrf:<sub>FillValue</sub> = -32767s ;
> 		msdwswrf:missing<sub>value</sub> = -32767s ;
> 		msdwswrf:units = "W m\*\*-2" ;
> 		msdwswrf:long<sub>name</sub> = "Mean surface downward short-wave radiation flux" ;
> 
> // global attributes:
> 		:Conventions = "CF-1.6" ;
> 		:history = "2020-10-06 09:05:31 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: mv *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msdwswrf/2010*.msdwswrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp *g/data/rt52/era5/single-levels/reanalysis/msdwswrf/2010/msdwswrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc\n",
> 			"2020-10-06 09:04:54 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: nccopy -k4 -d5 -s /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msdwswrf/2010/msdwswrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc /g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msdwswrf/2010*.msdwswrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp\n",
> 			"2020-10-06 09:01:28 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: curl &#x2013;connect-timeout 20 &#x2013;show-error &#x2013;silent &#x2013;max-time 36000 -o /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msdwswrf/2010/msdwswrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc <http://136.156.132.105/cache-compute-0000/cache/data1/adaptor.mars.internal-1601935012.7223494-5877-21-fad3c536-1dd2-422a-9b5d-47a7db1e9840.nc\n>",
> 			"2020-10-05 21:57:36 GMT by grib<sub>to</sub><sub>netcdf</sub>-2.16.0: *opt/ecmwf/eccodes/bin/grib<sub>to</sub><sub>netcdf</sub> -S param -o /cache/data1/adaptor.mars.internal-1601935012.7223494-5877-21-fad3c536-1dd2-422a-9b5d-47a7db1e9840.nc /cache/tmp/fad3c536-1dd2-422a-9b5d-47a7db1e9840-adaptor.mars.internal-1601935012.7229025-5877-7-tmp.grib" ;
> 		:license = "Licence to use Copernicus Products: <https://apps.ecmwf.int/datasets/licences/copernicus>*" ;
> 		:summary = "ERA5 is the fifth generation ECMWF atmospheric reanalysis of the global climate. This file is part of the ERA5 replica hosted at NCI Australia. For more information please see <http://dx.doi.org/10.25914/5f48874388857>" ;
> 		:title = "ERA5 single-levels reanalysis mean<sub>surface</sub><sub>downward</sub><sub>short</sub><sub>wave</sub><sub>radiation</sub><sub>flux</sub> 20100101-20100131" ;
> }


<a id="org7099012"></a>

### [msdwswrfcs](file:///Volumes/ioa02/reanalysis/ERA5/msdwswrfcs)

Mean surface downward short-wave radiation flux, clear sky

> netcdf msdwswrfcs<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131 {
> dimensions:
> 	longitude = 1440 ;
> 	latitude = 721 ;
> 	time = 744 ;
> variables:
> 	float longitude(longitude) ;
> 		longitude:units = "degrees<sub>east</sub>" ;
> 		longitude:long<sub>name</sub> = "longitude" ;
> 	float latitude(latitude) ;
> 		latitude:units = "degrees<sub>north</sub>" ;
> 		latitude:long<sub>name</sub> = "latitude" ;
> 	int time(time) ;
> 		time:units = "hours since 1900-01-01 00:00:00.0" ;
> 		time:long<sub>name</sub> = "time" ;
> 		time:calendar = "gregorian" ;
> 	short msdwswrfcs(time, latitude, longitude) ;
> 		msdwswrfcs:scale<sub>factor</sub> = 0.0194467863519143 ;
> 		msdwswrfcs:add<sub>offset</sub> = 637.193401606824 ;
> 		msdwswrfcs:<sub>FillValue</sub> = -32767s ;
> 		msdwswrfcs:missing<sub>value</sub> = -32767s ;
> 		msdwswrfcs:units = "W m\*\*-2" ;
> 		msdwswrfcs:long<sub>name</sub> = "Mean surface downward short-wave radiation flux, clear sky" ;
> 
> // global attributes:
> 		:Conventions = "CF-1.6" ;
> 		:history = "2020-10-06 09:12:58 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: mv *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msdwswrfcs/2010*.msdwswrfcs<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp *g/data/rt52/era5/single-levels/reanalysis/msdwswrfcs/2010/msdwswrfcs<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc\n",
> 			"2020-10-06 09:12:25 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: nccopy -k4 -d5 -s /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msdwswrfcs/2010/msdwswrfcs<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc /g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msdwswrfcs/2010*.msdwswrfcs<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp\n",
> 			"2020-10-06 09:08:27 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: curl &#x2013;connect-timeout 20 &#x2013;show-error &#x2013;silent &#x2013;max-time 36000 -o /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msdwswrfcs/2010/msdwswrfcs<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc <http://136.156.132.235/cache-compute-0006/cache/data8/adaptor.mars.internal-1601935429.6491008-27899-16-e96963ed-62ea-4560-b279-6b1ed587aa92.nc\n>",
> 			"2020-10-05 22:04:51 GMT by grib<sub>to</sub><sub>netcdf</sub>-2.16.0: *opt/ecmwf/eccodes/bin/grib<sub>to</sub><sub>netcdf</sub> -S param -o /cache/data8/adaptor.mars.internal-1601935429.6491008-27899-16-e96963ed-62ea-4560-b279-6b1ed587aa92.nc /cache/tmp/e96963ed-62ea-4560-b279-6b1ed587aa92-adaptor.mars.internal-1601935429.6498413-27899-8-tmp.grib" ;
> 		:license = "Licence to use Copernicus Products: <https://apps.ecmwf.int/datasets/licences/copernicus>*" ;
> 		:summary = "ERA5 is the fifth generation ECMWF atmospheric reanalysis of the global climate. This file is part of the ERA5 replica hosted at NCI Australia. For more information please see <http://dx.doi.org/10.25914/5f48874388857>" ;
> 		:title = "ERA5 single-levels reanalysis mean<sub>surface</sub><sub>downward</sub><sub>short</sub><sub>wave</sub><sub>radiation</sub><sub>flux</sub><sub>clear</sub><sub>sky</sub> 20100101-20100131" ;
> }


<a id="org36ae245"></a>

### [msl](file:///Volumes/ioa02/reanalysis/ERA5/msl)

Mean sea level pressure

> netcdf msl<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131 {
> dimensions:
> 	longitude = 1440 ;
> 	latitude = 721 ;
> 	time = 744 ;
> variables:
> 	float longitude(longitude) ;
> 		longitude:units = "degrees<sub>east</sub>" ;
> 		longitude:long<sub>name</sub> = "longitude" ;
> 	float latitude(latitude) ;
> 		latitude:units = "degrees<sub>north</sub>" ;
> 		latitude:long<sub>name</sub> = "latitude" ;
> 	int time(time) ;
> 		time:units = "hours since 1900-01-01 00:00:00.0" ;
> 		time:long<sub>name</sub> = "time" ;
> 		time:calendar = "gregorian" ;
> 	short msl(time, latitude, longitude) ;
> 		msl:scale<sub>factor</sub> = 0.196378160621366 ;
> 		msl:add<sub>offset</sub> = 100651.52681092 ;
> 		msl:<sub>FillValue</sub> = -32767s ;
> 		msl:missing<sub>value</sub> = -32767s ;
> 		msl:units = "Pa" ;
> 		msl:long<sub>name</sub> = "Mean sea level pressure" ;
> 		msl:standard<sub>name</sub> = "air<sub>pressure</sub><sub>at</sub><sub>mean</sub><sub>sea</sub><sub>level</sub>" ;
> 
> // global attributes:
> 		:Conventions = "CF-1.6" ;
> 		:history = "2020-09-28 12:39:13 UTC+1000 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: mv *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msl/2010*.msl<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp *g/data/rt52/era5/single-levels/reanalysis/msl/2010/msl<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc\n",
> 			"2020-09-28 12:38:23 UTC+1000 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: nccopy -k4 -d5 -s /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msl/2010/msl<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc /g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msl/2010*.msl<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp\n",
> 			"2020-09-28 12:35:03 UTC+1000 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: curl &#x2013;connect-timeout 20 &#x2013;show-error &#x2013;silent &#x2013;max-time 36000 -o /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msl/2010/msl<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc <http://136.156.133.39/cache-compute-0012/cache/data5/adaptor.mars.internal-1601260296.1220453-31017-5-1e4f2d61-ed30-44f1-bdfb-012d72c51c74.nc\n>",
> 			"2020-09-28 02:32:49 GMT by grib<sub>to</sub><sub>netcdf</sub>-2.16.0: *opt/ecmwf/eccodes/bin/grib<sub>to</sub><sub>netcdf</sub> -S param -o /cache/data5/adaptor.mars.internal-1601260296.1220453-31017-5-1e4f2d61-ed30-44f1-bdfb-012d72c51c74.nc /cache/tmp/1e4f2d61-ed30-44f1-bdfb-012d72c51c74-adaptor.mars.internal-1601260296.12292-31017-3-tmp.grib" ;
> 		:license = "Licence to use Copernicus Products: <https://apps.ecmwf.int/datasets/licences/copernicus>*" ;
> 		:summary = "ERA5 is the fifth generation ECMWF atmospheric reanalysis of the global climate. This file is part of the ERA5 replica hosted at NCI Australia. For more information please see <http://dx.doi.org/10.25914/5f48874388857>" ;
> 		:title = "ERA5 single-levels reanalysis mean<sub>sea</sub><sub>level</sub><sub>pressure</sub> 20100101-20100131" ;
> }


<a id="org69dda3d"></a>

### [msnlwrf](file:///Volumes/ioa02/reanalysis/ERA5/msnlwrf)

Mean surface net long-wave radiation flux

> netcdf msnlwrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131 {
> dimensions:
> 	longitude = 1440 ;
> 	latitude = 721 ;
> 	time = 744 ;
> variables:
> 	float longitude(longitude) ;
> 		longitude:units = "degrees<sub>east</sub>" ;
> 		longitude:long<sub>name</sub> = "longitude" ;
> 	float latitude(latitude) ;
> 		latitude:units = "degrees<sub>north</sub>" ;
> 		latitude:long<sub>name</sub> = "latitude" ;
> 	int time(time) ;
> 		time:units = "hours since 1900-01-01 00:00:00.0" ;
> 		time:long<sub>name</sub> = "time" ;
> 		time:calendar = "gregorian" ;
> 	short msnlwrf(time, latitude, longitude) ;
> 		msnlwrf:scale<sub>factor</sub> = 0.00765515516156936 ;
> 		msnlwrf:add<sub>offset</sub> = -185.368207460393 ;
> 		msnlwrf:<sub>FillValue</sub> = -32767s ;
> 		msnlwrf:missing<sub>value</sub> = -32767s ;
> 		msnlwrf:units = "W m\*\*-2" ;
> 		msnlwrf:long<sub>name</sub> = "Mean surface net long-wave radiation flux" ;
> 
> // global attributes:
> 		:Conventions = "CF-1.6" ;
> 		:history = "2020-10-06 09:09:18 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: mv *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msnlwrf/2010*.msnlwrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp *g/data/rt52/era5/single-levels/reanalysis/msnlwrf/2010/msnlwrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc\n",
> 			"2020-10-06 09:08:14 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: nccopy -k4 -d5 -s /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msnlwrf/2010/msnlwrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc /g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msnlwrf/2010*.msnlwrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp\n",
> 			"2020-10-06 09:02:18 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: curl &#x2013;connect-timeout 20 &#x2013;show-error &#x2013;silent &#x2013;max-time 36000 -o /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msnlwrf/2010/msnlwrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc <http://136.156.132.153/cache-compute-0002/cache/data6/adaptor.mars.internal-1601935105.2208984-18608-5-8a196cc4-a7d0-47dd-b326-528e0e8ffb31.nc\n>",
> 			"2020-10-05 21:59:08 GMT by grib<sub>to</sub><sub>netcdf</sub>-2.16.0: *opt/ecmwf/eccodes/bin/grib<sub>to</sub><sub>netcdf</sub> -S param -o /cache/data6/adaptor.mars.internal-1601935105.2208984-18608-5-8a196cc4-a7d0-47dd-b326-528e0e8ffb31.nc /cache/tmp/8a196cc4-a7d0-47dd-b326-528e0e8ffb31-adaptor.mars.internal-1601935105.2214694-18608-1-tmp.grib" ;
> 		:license = "Licence to use Copernicus Products: <https://apps.ecmwf.int/datasets/licences/copernicus>*" ;
> 		:summary = "ERA5 is the fifth generation ECMWF atmospheric reanalysis of the global climate. This file is part of the ERA5 replica hosted at NCI Australia. For more information please see <http://dx.doi.org/10.25914/5f48874388857>" ;
> 		:title = "ERA5 single-levels reanalysis mean<sub>surface</sub><sub>net</sub><sub>long</sub><sub>wave</sub><sub>radiation</sub><sub>flux</sub> 20100101-20100131" ;
> }


<a id="org9e1095b"></a>

### [msnlwrfcs](file:///Volumes/ioa02/reanalysis/ERA5/msnlwrfcs)

Mean surface net long-wave radiation flux, clear sky

> netcdf msnlwrfcs<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131 {
> dimensions:
> 	longitude = 1440 ;
> 	latitude = 721 ;
> 	time = 744 ;
> variables:
> 	float longitude(longitude) ;
> 		longitude:units = "degrees<sub>east</sub>" ;
> 		longitude:long<sub>name</sub> = "longitude" ;
> 	float latitude(latitude) ;
> 		latitude:units = "degrees<sub>north</sub>" ;
> 		latitude:long<sub>name</sub> = "latitude" ;
> 	int time(time) ;
> 		time:units = "hours since 1900-01-01 00:00:00.0" ;
> 		time:long<sub>name</sub> = "time" ;
> 		time:calendar = "gregorian" ;
> 	short msnlwrfcs(time, latitude, longitude) ;
> 		msnlwrfcs:scale<sub>factor</sub> = 0.00686127063368837 ;
> 		msnlwrfcs:add<sub>offset</sub> = -207.683118135317 ;
> 		msnlwrfcs:<sub>FillValue</sub> = -32767s ;
> 		msnlwrfcs:missing<sub>value</sub> = -32767s ;
> 		msnlwrfcs:units = "W m\*\*-2" ;
> 		msnlwrfcs:long<sub>name</sub> = "Mean surface net long-wave radiation flux, clear sky" ;
> 
> // global attributes:
> 		:Conventions = "CF-1.6" ;
> 		:history = "2020-10-06 09:11:31 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: mv *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msnlwrfcs/2010*.msnlwrfcs<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp *g/data/rt52/era5/single-levels/reanalysis/msnlwrfcs/2010/msnlwrfcs<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc\n",
> 			"2020-10-06 09:10:36 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: nccopy -k4 -d5 -s /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msnlwrfcs/2010/msnlwrfcs<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc /g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msnlwrfcs/2010*.msnlwrfcs<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp\n",
> 			"2020-10-06 09:06:58 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: curl &#x2013;connect-timeout 20 &#x2013;show-error &#x2013;silent &#x2013;max-time 36000 -o /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msnlwrfcs/2010/msnlwrfcs<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc <http://136.156.133.25/cache-compute-0008/cache/data7/adaptor.mars.internal-1601935291.563169-20885-28-16218919-c599-4310-b83c-d3e5b2d7ef7d.nc\n>",
> 			"2020-10-05 22:02:20 GMT by grib<sub>to</sub><sub>netcdf</sub>-2.16.0: *opt/ecmwf/eccodes/bin/grib<sub>to</sub><sub>netcdf</sub> -S param -o /cache/data7/adaptor.mars.internal-1601935291.563169-20885-28-16218919-c599-4310-b83c-d3e5b2d7ef7d.nc /cache/tmp/16218919-c599-4310-b83c-d3e5b2d7ef7d-adaptor.mars.internal-1601935291.5818343-20885-10-tmp.grib" ;
> 		:license = "Licence to use Copernicus Products: <https://apps.ecmwf.int/datasets/licences/copernicus>*" ;
> 		:summary = "ERA5 is the fifth generation ECMWF atmospheric reanalysis of the global climate. This file is part of the ERA5 replica hosted at NCI Australia. For more information please see <http://dx.doi.org/10.25914/5f48874388857>" ;
> 		:title = "ERA5 single-levels reanalysis mean<sub>surface</sub><sub>net</sub><sub>long</sub><sub>wave</sub><sub>radiation</sub><sub>flux</sub><sub>clear</sub><sub>sky</sub> 20100101-20100131" ;
> }


<a id="orgbbbc293"></a>

### [msnswrf](file:///Volumes/ioa02/reanalysis/ERA5/msnswrf)

Mean surface net short-wave radiation flux

> netcdf msnswrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131 {
> dimensions:
> 	longitude = 1440 ;
> 	latitude = 721 ;
> 	time = 744 ;
> variables:
> 	float longitude(longitude) ;
> 		longitude:units = "degrees<sub>east</sub>" ;
> 		longitude:long<sub>name</sub> = "longitude" ;
> 	float latitude(latitude) ;
> 		latitude:units = "degrees<sub>north</sub>" ;
> 		latitude:long<sub>name</sub> = "latitude" ;
> 	int time(time) ;
> 		time:units = "hours since 1900-01-01 00:00:00.0" ;
> 		time:long<sub>name</sub> = "time" ;
> 		time:calendar = "gregorian" ;
> 	short msnswrf(time, latitude, longitude) ;
> 		msnswrf:scale<sub>factor</sub> = 0.0171468954572505 ;
> 		msnswrf:add<sub>offset</sub> = 561.835176552271 ;
> 		msnswrf:<sub>FillValue</sub> = -32767s ;
> 		msnswrf:missing<sub>value</sub> = -32767s ;
> 		msnswrf:units = "W m\*\*-2" ;
> 		msnswrf:long<sub>name</sub> = "Mean surface net short-wave radiation flux" ;
> 
> // global attributes:
> 		:Conventions = "CF-1.6" ;
> 		:history = "2020-10-06 09:07:34 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: mv *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msnswrf/2010*.msnswrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp *g/data/rt52/era5/single-levels/reanalysis/msnswrf/2010/msnswrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc\n",
> 			"2020-10-06 09:06:58 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: nccopy -k4 -d5 -s /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msnswrf/2010/msnswrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc /g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msnswrf/2010*.msnswrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp\n",
> 			"2020-10-06 09:02:13 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: curl &#x2013;connect-timeout 20 &#x2013;show-error &#x2013;silent &#x2013;max-time 36000 -o /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msnswrf/2010/msnswrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc <http://136.156.133.39/cache-compute-0012/cache/data8/adaptor.mars.internal-1601935069.7217205-6425-1-6f222a72-e8c4-4c60-acb7-d3b89ffc8e8f.nc\n>",
> 			"2020-10-05 21:58:59 GMT by grib<sub>to</sub><sub>netcdf</sub>-2.16.0: *opt/ecmwf/eccodes/bin/grib<sub>to</sub><sub>netcdf</sub> -S param -o /cache/data8/adaptor.mars.internal-1601935069.7217205-6425-1-6f222a72-e8c4-4c60-acb7-d3b89ffc8e8f.nc /cache/tmp/6f222a72-e8c4-4c60-acb7-d3b89ffc8e8f-adaptor.mars.internal-1601935069.723566-6425-1-tmp.grib" ;
> 		:license = "Licence to use Copernicus Products: <https://apps.ecmwf.int/datasets/licences/copernicus>*" ;
> 		:summary = "ERA5 is the fifth generation ECMWF atmospheric reanalysis of the global climate. This file is part of the ERA5 replica hosted at NCI Australia. For more information please see <http://dx.doi.org/10.25914/5f48874388857>" ;
> 		:title = "ERA5 single-levels reanalysis mean<sub>surface</sub><sub>net</sub><sub>short</sub><sub>wave</sub><sub>radiation</sub><sub>flux</sub> 20100101-20100131" ;
> }


<a id="org679bb35"></a>

### [msnswrfcs](file:///Volumes/ioa02/reanalysis/ERA5/msnswrfcs)

Mean surface net short-wave radiation flux, clear sky

> netcdf msnswrfcs<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131 {
> dimensions:
> 	longitude = 1440 ;
> 	latitude = 721 ;
> 	time = 744 ;
> variables:
> 	float longitude(longitude) ;
> 		longitude:units = "degrees<sub>east</sub>" ;
> 		longitude:long<sub>name</sub> = "longitude" ;
> 	float latitude(latitude) ;
> 		latitude:units = "degrees<sub>north</sub>" ;
> 		latitude:long<sub>name</sub> = "latitude" ;
> 	int time(time) ;
> 		time:units = "hours since 1900-01-01 00:00:00.0" ;
> 		time:long<sub>name</sub> = "time" ;
> 		time:calendar = "gregorian" ;
> 	short msnswrfcs(time, latitude, longitude) ;
> 		msnswrfcs:scale<sub>factor</sub> = 0.0169985923122702 ;
> 		msnswrfcs:add<sub>offset</sub> = 556.975875703844 ;
> 		msnswrfcs:<sub>FillValue</sub> = -32767s ;
> 		msnswrfcs:missing<sub>value</sub> = -32767s ;
> 		msnswrfcs:units = "W m\*\*-2" ;
> 		msnswrfcs:long<sub>name</sub> = "Mean surface net short-wave radiation flux, clear sky" ;
> 
> // global attributes:
> 		:Conventions = "CF-1.6" ;
> 		:history = "2020-10-06 09:12:36 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: mv *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msnswrfcs/2010*.msnswrfcs<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp *g/data/rt52/era5/single-levels/reanalysis/msnswrfcs/2010/msnswrfcs<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc\n",
> 			"2020-10-06 09:12:01 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: nccopy -k4 -d5 -s /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msnswrfcs/2010/msnswrfcs<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc /g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msnswrfcs/2010*.msnswrfcs<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp\n",
> 			"2020-10-06 09:06:10 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: curl &#x2013;connect-timeout 20 &#x2013;show-error &#x2013;silent &#x2013;max-time 36000 -o /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msnswrfcs/2010/msnswrfcs<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc <http://136.156.133.39/cache-compute-0012/cache/data6/adaptor.mars.internal-1601935253.1640558-4989-18-bababa4c-ec9e-4ba6-98cd-9ad90e6be27b.nc\n>",
> 			"2020-10-05 22:01:38 GMT by grib<sub>to</sub><sub>netcdf</sub>-2.16.0: *opt/ecmwf/eccodes/bin/grib<sub>to</sub><sub>netcdf</sub> -S param -o /cache/data6/adaptor.mars.internal-1601935253.1640558-4989-18-bababa4c-ec9e-4ba6-98cd-9ad90e6be27b.nc /cache/tmp/bababa4c-ec9e-4ba6-98cd-9ad90e6be27b-adaptor.mars.internal-1601935253.164695-4989-6-tmp.grib" ;
> 		:license = "Licence to use Copernicus Products: <https://apps.ecmwf.int/datasets/licences/copernicus>*" ;
> 		:summary = "ERA5 is the fifth generation ECMWF atmospheric reanalysis of the global climate. This file is part of the ERA5 replica hosted at NCI Australia. For more information please see <http://dx.doi.org/10.25914/5f48874388857>" ;
> 		:title = "ERA5 single-levels reanalysis mean<sub>surface</sub><sub>net</sub><sub>short</sub><sub>wave</sub><sub>radiation</sub><sub>flux</sub><sub>clear</sub><sub>sky</sub> 20100101-20100131" ;
> }


<a id="org3984d1b"></a>

### [msr](file:///Volumes/ioa02/reanalysis/ERA5/msr)

Mean snowfall rate

> netcdf msr<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131 {
> dimensions:
> 	longitude = 1440 ;
> 	latitude = 721 ;
> 	time = 744 ;
> variables:
> 	float longitude(longitude) ;
> 		longitude:units = "degrees<sub>east</sub>" ;
> 		longitude:long<sub>name</sub> = "longitude" ;
> 	float latitude(latitude) ;
> 		latitude:units = "degrees<sub>north</sub>" ;
> 		latitude:long<sub>name</sub> = "latitude" ;
> 	int time(time) ;
> 		time:units = "hours since 1900-01-01 00:00:00.0" ;
> 		time:long<sub>name</sub> = "time" ;
> 		time:calendar = "gregorian" ;
> 	short msr(time, latitude, longitude) ;
> 		msr:scale<sub>factor</sub> = 3.16955222732488e-08 ;
> 		msr:add<sub>offset</sub> = 0.00103853548280527 ;
> 		msr:<sub>FillValue</sub> = -32767s ;
> 		msr:missing<sub>value</sub> = -32767s ;
> 		msr:units = "kg m\*\*-2 s\*\*-1" ;
> 		msr:long<sub>name</sub> = "Mean snowfall rate" ;
> 
> // global attributes:
> 		:Conventions = "CF-1.6" ;
> 		:history = "2020-10-06 09:04:04 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: mv *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msr/2010*.msr<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp *g/data/rt52/era5/single-levels/reanalysis/msr/2010/msr<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc\n",
> 			"2020-10-06 09:03:42 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: nccopy -k4 -d5 -s /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msr/2010/msr<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc /g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msr/2010*.msr<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp\n",
> 			"2020-10-06 09:00:18 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: curl &#x2013;connect-timeout 20 &#x2013;show-error &#x2013;silent &#x2013;max-time 36000 -o /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msr/2010/msr<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc <http://136.156.133.25/cache-compute-0008/cache/data8/adaptor.mars.internal-1601934958.8138766-19540-22-829a84f4-b413-4b71-b15a-99bb8c146d46.nc\n>",
> 			"2020-10-05 21:56:43 GMT by grib<sub>to</sub><sub>netcdf</sub>-2.16.0: *opt/ecmwf/eccodes/bin/grib<sub>to</sub><sub>netcdf</sub> -S param -o /cache/data8/adaptor.mars.internal-1601934958.8138766-19540-22-829a84f4-b413-4b71-b15a-99bb8c146d46.nc /cache/tmp/829a84f4-b413-4b71-b15a-99bb8c146d46-adaptor.mars.internal-1601934958.8144617-19540-8-tmp.grib" ;
> 		:license = "Licence to use Copernicus Products: <https://apps.ecmwf.int/datasets/licences/copernicus>*" ;
> 		:summary = "ERA5 is the fifth generation ECMWF atmospheric reanalysis of the global climate. This file is part of the ERA5 replica hosted at NCI Australia. For more information please see <http://dx.doi.org/10.25914/5f48874388857>" ;
> 		:title = "ERA5 single-levels reanalysis mean<sub>snowfall</sub><sub>rate</sub> 20100101-20100131" ;
> }


<a id="org2f58934"></a>

### [msror](file:///Volumes/ioa02/reanalysis/ERA5/msror)

Mean surface runoff rate

> netcdf msror<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131 {
> dimensions:
> 	longitude = 1440 ;
> 	latitude = 721 ;
> 	time = 744 ;
> variables:
> 	float longitude(longitude) ;
> 		longitude:units = "degrees<sub>east</sub>" ;
> 		longitude:long<sub>name</sub> = "longitude" ;
> 	float latitude(latitude) ;
> 		latitude:units = "degrees<sub>north</sub>" ;
> 		latitude:long<sub>name</sub> = "latitude" ;
> 	int time(time) ;
> 		time:units = "hours since 1900-01-01 00:00:00.0" ;
> 		time:long<sub>name</sub> = "time" ;
> 		time:calendar = "gregorian" ;
> 	short msror(time, latitude, longitude) ;
> 		msror:scale<sub>factor</sub> = 3.49596662796261e-07 ;
> 		msror:add<sub>offset</sub> = 0.0114548842531823 ;
> 		msror:<sub>FillValue</sub> = -32767s ;
> 		msror:missing<sub>value</sub> = -32767s ;
> 		msror:units = "kg m\*\*-2 s\*\*-1" ;
> 		msror:long<sub>name</sub> = "Mean surface runoff rate" ;
> 
> // global attributes:
> 		:Conventions = "CF-1.6" ;
> 		:history = "2020-10-06 09:02:38 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: mv *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msror/2010*.msror<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp *g/data/rt52/era5/single-levels/reanalysis/msror/2010/msror<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc\n",
> 			"2020-10-06 09:02:24 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: nccopy -k4 -d5 -s /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msror/2010/msror<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc /g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msror/2010*.msror<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp\n",
> 			"2020-10-06 08:58:36 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: curl &#x2013;connect-timeout 20 &#x2013;show-error &#x2013;silent &#x2013;max-time 36000 -o /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msror/2010/msror<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc <http://136.156.133.46/cache-compute-0015/cache/data2/adaptor.mars.internal-1601934851.1730347-28836-22-a7717ea4-1d51-4b28-a876-87ce3356d63c.nc\n>",
> 			"2020-10-05 21:55:17 GMT by grib<sub>to</sub><sub>netcdf</sub>-2.16.0: *opt/ecmwf/eccodes/bin/grib<sub>to</sub><sub>netcdf</sub> -S param -o /cache/data2/adaptor.mars.internal-1601934851.1730347-28836-22-a7717ea4-1d51-4b28-a876-87ce3356d63c.nc /cache/tmp/a7717ea4-1d51-4b28-a876-87ce3356d63c-adaptor.mars.internal-1601934851.173735-28836-8-tmp.grib" ;
> 		:license = "Licence to use Copernicus Products: <https://apps.ecmwf.int/datasets/licences/copernicus>*" ;
> 		:summary = "ERA5 is the fifth generation ECMWF atmospheric reanalysis of the global climate. This file is part of the ERA5 replica hosted at NCI Australia. For more information please see <http://dx.doi.org/10.25914/5f48874388857>" ;
> 		:title = "ERA5 single-levels reanalysis mean<sub>surface</sub><sub>runoff</sub><sub>rate</sub> 20100101-20100131" ;
> }


<a id="orgb3f1f7b"></a>

### [mtdwswrf](file:///Volumes/ioa02/reanalysis/ERA5/mtdwswrf)

Mean top downward short-wave radiation flux

> netcdf mtdwswrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131 {
> dimensions:
> 	longitude = 1440 ;
> 	latitude = 721 ;
> 	time = 744 ;
> variables:
> 	float longitude(longitude) ;
> 		longitude:units = "degrees<sub>east</sub>" ;
> 		longitude:long<sub>name</sub> = "longitude" ;
> 	float latitude(latitude) ;
> 		latitude:units = "degrees<sub>north</sub>" ;
> 		latitude:long<sub>name</sub> = "latitude" ;
> 	int time(time) ;
> 		time:units = "hours since 1900-01-01 00:00:00.0" ;
> 		time:long<sub>name</sub> = "time" ;
> 		time:calendar = "gregorian" ;
> 	short mtdwswrf(time, latitude, longitude) ;
> 		mtdwswrf:scale<sub>factor</sub> = 0.021422889994354 ;
> 		mtdwswrf:add<sub>offset</sub> = 701.942413555003 ;
> 		mtdwswrf:<sub>FillValue</sub> = -32767s ;
> 		mtdwswrf:missing<sub>value</sub> = -32767s ;
> 		mtdwswrf:units = "W m\*\*-2" ;
> 		mtdwswrf:long<sub>name</sub> = "Mean top downward short-wave radiation flux" ;
> 
> // global attributes:
> 		:Conventions = "CF-1.6" ;
> 		:history = "2020-10-06 09:11:01 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: mv *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/mtdwswrf/2010*.mtdwswrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp *g/data/rt52/era5/single-levels/reanalysis/mtdwswrf/2010/mtdwswrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc\n",
> 			"2020-10-06 09:10:32 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: nccopy -k4 -d5 -s /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/mtdwswrf/2010/mtdwswrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc /g/data/rt52/admin/incoming/era5/single-levels/reanalysis/mtdwswrf/2010*.mtdwswrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp\n",
> 			"2020-10-06 09:06:38 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: curl &#x2013;connect-timeout 20 &#x2013;show-error &#x2013;silent &#x2013;max-time 36000 -o /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/mtdwswrf/2010/mtdwswrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc <http://136.156.132.236/cache-compute-0007/cache/data7/adaptor.mars.internal-1601935311.62434-12055-32-fb7e93c7-424f-4538-bf13-ba57fce74df4.nc\n>",
> 			"2020-10-05 22:02:35 GMT by grib<sub>to</sub><sub>netcdf</sub>-2.16.0: *opt/ecmwf/eccodes/bin/grib<sub>to</sub><sub>netcdf</sub> -S param -o /cache/data7/adaptor.mars.internal-1601935311.62434-12055-32-fb7e93c7-424f-4538-bf13-ba57fce74df4.nc /cache/tmp/fb7e93c7-424f-4538-bf13-ba57fce74df4-adaptor.mars.internal-1601935311.6249025-12055-14-tmp.grib" ;
> 		:license = "Licence to use Copernicus Products: <https://apps.ecmwf.int/datasets/licences/copernicus>*" ;
> 		:summary = "ERA5 is the fifth generation ECMWF atmospheric reanalysis of the global climate. This file is part of the ERA5 replica hosted at NCI Australia. For more information please see <http://dx.doi.org/10.25914/5f48874388857>" ;
> 		:title = "ERA5 single-levels reanalysis mean<sub>top</sub><sub>downward</sub><sub>short</sub><sub>wave</sub><sub>radiation</sub><sub>flux</sub> 20100101-20100131" ;
> }


<a id="orgffb5609"></a>

### [mtnlwrf](file:///Volumes/ioa02/reanalysis/ERA5/mtnlwrf)

Mean top net long-wave radiation flux

> netcdf mtnlwrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131 {
> dimensions:
> 	longitude = 1440 ;
> 	latitude = 721 ;
> 	time = 744 ;
> variables:
> 	float longitude(longitude) ;
> 		longitude:units = "degrees<sub>east</sub>" ;
> 		longitude:long<sub>name</sub> = "longitude" ;
> 	float latitude(latitude) ;
> 		latitude:units = "degrees<sub>north</sub>" ;
> 		latitude:long<sub>name</sub> = "latitude" ;
> 	int time(time) ;
> 		time:units = "hours since 1900-01-01 00:00:00.0" ;
> 		time:long<sub>name</sub> = "time" ;
> 		time:calendar = "gregorian" ;
> 	short mtnlwrf(time, latitude, longitude) ;
> 		mtnlwrf:scale<sub>factor</sub> = 0.00547009042080898 ;
> 		mtnlwrf:add<sub>offset</sub> = -239.618579771773 ;
> 		mtnlwrf:<sub>FillValue</sub> = -32767s ;
> 		mtnlwrf:missing<sub>value</sub> = -32767s ;
> 		mtnlwrf:units = "W m\*\*-2" ;
> 		mtnlwrf:long<sub>name</sub> = "Mean top net long-wave radiation flux" ;
> 
> // global attributes:
> 		:Conventions = "CF-1.6" ;
> 		:history = "2020-10-06 09:09:23 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: mv *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/mtnlwrf/2010*.mtnlwrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp *g/data/rt52/era5/single-levels/reanalysis/mtnlwrf/2010/mtnlwrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc\n",
> 			"2020-10-06 09:08:21 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: nccopy -k4 -d5 -s /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/mtnlwrf/2010/mtnlwrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc /g/data/rt52/admin/incoming/era5/single-levels/reanalysis/mtnlwrf/2010*.mtnlwrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp\n",
> 			"2020-10-06 09:02:24 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: curl &#x2013;connect-timeout 20 &#x2013;show-error &#x2013;silent &#x2013;max-time 36000 -o /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/mtnlwrf/2010/mtnlwrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc <http://136.156.133.41/cache-compute-0013/cache/data1/adaptor.mars.internal-1601935118.9578016-15806-13-c7cd5c57-4899-4528-8267-f01896333fba.nc\n>",
> 			"2020-10-05 21:59:22 GMT by grib<sub>to</sub><sub>netcdf</sub>-2.16.0: *opt/ecmwf/eccodes/bin/grib<sub>to</sub><sub>netcdf</sub> -S param -o /cache/data1/adaptor.mars.internal-1601935118.9578016-15806-13-c7cd5c57-4899-4528-8267-f01896333fba.nc /cache/tmp/c7cd5c57-4899-4528-8267-f01896333fba-adaptor.mars.internal-1601935118.9584107-15806-4-tmp.grib" ;
> 		:license = "Licence to use Copernicus Products: <https://apps.ecmwf.int/datasets/licences/copernicus>*" ;
> 		:summary = "ERA5 is the fifth generation ECMWF atmospheric reanalysis of the global climate. This file is part of the ERA5 replica hosted at NCI Australia. For more information please see <http://dx.doi.org/10.25914/5f48874388857>" ;
> 		:title = "ERA5 single-levels reanalysis mean<sub>top</sub><sub>net</sub><sub>long</sub><sub>wave</sub><sub>radiation</sub><sub>flux</sub> 20100101-20100131" ;
> }


<a id="org8091c1d"></a>

### [mtnlwrfcs](file:///Volumes/ioa02/reanalysis/ERA5/mtnlwrfcs)

Mean top net long-wave radiation flux, clear sky

> netcdf mtnlwrfcs<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131 {
> dimensions:
> 	longitude = 1440 ;
> 	latitude = 721 ;
> 	time = 744 ;
> variables:
> 	float longitude(longitude) ;
> 		longitude:units = "degrees<sub>east</sub>" ;
> 		longitude:long<sub>name</sub> = "longitude" ;
> 	float latitude(latitude) ;
> 		latitude:units = "degrees<sub>north</sub>" ;
> 		latitude:long<sub>name</sub> = "latitude" ;
> 	int time(time) ;
> 		time:units = "hours since 1900-01-01 00:00:00.0" ;
> 		time:long<sub>name</sub> = "time" ;
> 		time:calendar = "gregorian" ;
> 	short mtnlwrfcs(time, latitude, longitude) ;
> 		mtnlwrfcs:scale<sub>factor</sub> = 0.00447323906109327 ;
> 		mtnlwrfcs:add<sub>offset</sub> = -272.350747361718 ;
> 		mtnlwrfcs:<sub>FillValue</sub> = -32767s ;
> 		mtnlwrfcs:missing<sub>value</sub> = -32767s ;
> 		mtnlwrfcs:units = "W m\*\*-2" ;
> 		mtnlwrfcs:long<sub>name</sub> = "Mean top net long-wave radiation flux, clear sky" ;
> 
> // global attributes:
> 		:Conventions = "CF-1.6" ;
> 		:history = "2020-10-06 09:12:24 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: mv *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/mtnlwrfcs/2010*.mtnlwrfcs<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp *g/data/rt52/era5/single-levels/reanalysis/mtnlwrfcs/2010/mtnlwrfcs<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc\n",
> 			"2020-10-06 09:11:30 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: nccopy -k4 -d5 -s /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/mtnlwrfcs/2010/mtnlwrfcs<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc /g/data/rt52/admin/incoming/era5/single-levels/reanalysis/mtnlwrfcs/2010*.mtnlwrfcs<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp\n",
> 			"2020-10-06 09:06:09 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: curl &#x2013;connect-timeout 20 &#x2013;show-error &#x2013;silent &#x2013;max-time 36000 -o /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/mtnlwrfcs/2010/mtnlwrfcs<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc <http://136.156.132.110/cache-compute-0001/cache/data6/adaptor.mars.internal-1601935245.2007327-31712-13-9fe894bb-6f03-4d13-b576-5f138e35f11b.nc\n>",
> 			"2020-10-05 22:01:29 GMT by grib<sub>to</sub><sub>netcdf</sub>-2.16.0: *opt/ecmwf/eccodes/bin/grib<sub>to</sub><sub>netcdf</sub> -S param -o /cache/data6/adaptor.mars.internal-1601935245.2007327-31712-13-9fe894bb-6f03-4d13-b576-5f138e35f11b.nc /cache/tmp/9fe894bb-6f03-4d13-b576-5f138e35f11b-adaptor.mars.internal-1601935245.2013144-31712-6-tmp.grib" ;
> 		:license = "Licence to use Copernicus Products: <https://apps.ecmwf.int/datasets/licences/copernicus>*" ;
> 		:summary = "ERA5 is the fifth generation ECMWF atmospheric reanalysis of the global climate. This file is part of the ERA5 replica hosted at NCI Australia. For more information please see <http://dx.doi.org/10.25914/5f48874388857>" ;
> 		:title = "ERA5 single-levels reanalysis mean<sub>top</sub><sub>net</sub><sub>long</sub><sub>wave</sub><sub>radiation</sub><sub>flux</sub><sub>clear</sub><sub>sky</sub> 20100101-20100131" ;
> }


<a id="org62bfd96"></a>

### [mtnswrf](file:///Volumes/ioa02/reanalysis/ERA5/mtnswrf)

Mean top net short-wave radiation flux

> netcdf mtnswrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131 {
> dimensions:
> 	longitude = 1440 ;
> 	latitude = 721 ;
> 	time = 744 ;
> variables:
> 	float longitude(longitude) ;
> 		longitude:units = "degrees<sub>east</sub>" ;
> 		longitude:long<sub>name</sub> = "longitude" ;
> 	float latitude(latitude) ;
> 		latitude:units = "degrees<sub>north</sub>" ;
> 		latitude:long<sub>name</sub> = "latitude" ;
> 	int time(time) ;
> 		time:units = "hours since 1900-01-01 00:00:00.0" ;
> 		time:long<sub>name</sub> = "time" ;
> 		time:calendar = "gregorian" ;
> 	short mtnswrf(time, latitude, longitude) ;
> 		mtnswrf:scale<sub>factor</sub> = 0.0207843758106603 ;
> 		mtnswrf:add<sub>offset</sub> = 681.020857812095 ;
> 		mtnswrf:<sub>FillValue</sub> = -32767s ;
> 		mtnswrf:missing<sub>value</sub> = -32767s ;
> 		mtnswrf:units = "W m\*\*-2" ;
> 		mtnswrf:long<sub>name</sub> = "Mean top net short-wave radiation flux" ;
> 
> // global attributes:
> 		:Conventions = "CF-1.6" ;
> 		:history = "2020-10-06 09:08:09 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: mv *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/mtnswrf/2010*.mtnswrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp *g/data/rt52/era5/single-levels/reanalysis/mtnswrf/2010/mtnswrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc\n",
> 			"2020-10-06 09:07:33 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: nccopy -k4 -d5 -s /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/mtnswrf/2010/mtnswrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc /g/data/rt52/admin/incoming/era5/single-levels/reanalysis/mtnswrf/2010*.mtnswrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp\n",
> 			"2020-10-06 09:03:06 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: curl &#x2013;connect-timeout 20 &#x2013;show-error &#x2013;silent &#x2013;max-time 36000 -o /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/mtnswrf/2010/mtnswrf<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc <http://136.156.132.201/cache-compute-0004/cache/data9/adaptor.mars.internal-1601935110.2694445-4501-11-b5809659-6840-445d-a99a-233947a79209.nc\n>",
> 			"2020-10-05 21:59:15 GMT by grib<sub>to</sub><sub>netcdf</sub>-2.16.0: *opt/ecmwf/eccodes/bin/grib<sub>to</sub><sub>netcdf</sub> -S param -o /cache/data9/adaptor.mars.internal-1601935110.2694445-4501-11-b5809659-6840-445d-a99a-233947a79209.nc /cache/tmp/b5809659-6840-445d-a99a-233947a79209-adaptor.mars.internal-1601935110.270101-4501-4-tmp.grib" ;
> 		:license = "Licence to use Copernicus Products: <https://apps.ecmwf.int/datasets/licences/copernicus>*" ;
> 		:summary = "ERA5 is the fifth generation ECMWF atmospheric reanalysis of the global climate. This file is part of the ERA5 replica hosted at NCI Australia. For more information please see <http://dx.doi.org/10.25914/5f48874388857>" ;
> 		:title = "ERA5 single-levels reanalysis mean<sub>top</sub><sub>net</sub><sub>short</sub><sub>wave</sub><sub>radiation</sub><sub>flux</sub> 20100101-20100131" ;
> }


<a id="org61ea099"></a>

### [mtnswrfcs](file:///Volumes/ioa02/reanalysis/ERA5/mtnswrfcs)

Mean top net short-wave radiation flux, clear sky

> netcdf mtnswrfcs<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131 {
> dimensions:
> 	longitude = 1440 ;
> 	latitude = 721 ;
> 	time = 744 ;
> variables:
> 	float longitude(longitude) ;
> 		longitude:units = "degrees<sub>east</sub>" ;
> 		longitude:long<sub>name</sub> = "longitude" ;
> 	float latitude(latitude) ;
> 		latitude:units = "degrees<sub>north</sub>" ;
> 		latitude:long<sub>name</sub> = "latitude" ;
> 	int time(time) ;
> 		time:units = "hours since 1900-01-01 00:00:00.0" ;
> 		time:long<sub>name</sub> = "time" ;
> 		time:calendar = "gregorian" ;
> 	short mtnswrfcs(time, latitude, longitude) ;
> 		mtnswrfcs:scale<sub>factor</sub> = 0.0203156234263653 ;
> 		mtnswrfcs:add<sub>offset</sub> = 665.661717188287 ;
> 		mtnswrfcs:<sub>FillValue</sub> = -32767s ;
> 		mtnswrfcs:missing<sub>value</sub> = -32767s ;
> 		mtnswrfcs:units = "W m\*\*-2" ;
> 		mtnswrfcs:long<sub>name</sub> = "Mean top net short-wave radiation flux, clear sky" ;
> 
> // global attributes:
> 		:Conventions = "CF-1.6" ;
> 		:history = "2020-10-06 09:10:27 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: mv *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/mtnswrfcs/2010*.mtnswrfcs<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp *g/data/rt52/era5/single-levels/reanalysis/mtnswrfcs/2010/mtnswrfcs<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc\n",
> 			"2020-10-06 09:09:55 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: nccopy -k4 -d5 -s /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/mtnswrfcs/2010/mtnswrfcs<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc /g/data/rt52/admin/incoming/era5/single-levels/reanalysis/mtnswrfcs/2010*.mtnswrfcs<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp\n",
> 			"2020-10-06 09:06:03 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: curl &#x2013;connect-timeout 20 &#x2013;show-error &#x2013;silent &#x2013;max-time 36000 -o /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/mtnswrfcs/2010/mtnswrfcs<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc <http://136.156.133.37/cache-compute-0011/cache/data8/adaptor.mars.internal-1601935241.8915317-1334-3-4c3516f9-0313-410d-9c58-6f6f823a7e3f.nc\n>",
> 			"2020-10-05 22:01:23 GMT by grib<sub>to</sub><sub>netcdf</sub>-2.16.0: *opt/ecmwf/eccodes/bin/grib<sub>to</sub><sub>netcdf</sub> -S param -o /cache/data8/adaptor.mars.internal-1601935241.8915317-1334-3-4c3516f9-0313-410d-9c58-6f6f823a7e3f.nc /cache/tmp/4c3516f9-0313-410d-9c58-6f6f823a7e3f-adaptor.mars.internal-1601935241.8921776-1334-1-tmp.grib" ;
> 		:license = "Licence to use Copernicus Products: <https://apps.ecmwf.int/datasets/licences/copernicus>*" ;
> 		:summary = "ERA5 is the fifth generation ECMWF atmospheric reanalysis of the global climate. This file is part of the ERA5 replica hosted at NCI Australia. For more information please see <http://dx.doi.org/10.25914/5f48874388857>" ;
> 		:title = "ERA5 single-levels reanalysis mean<sub>top</sub><sub>net</sub><sub>short</sub><sub>wave</sub><sub>radiation</sub><sub>flux</sub><sub>clear</sub><sub>sky</sub> 20100101-20100131" ;
> }


<a id="org703fad6"></a>

### [mtpr](file:///Volumes/ioa02/reanalysis/ERA5/mtpr)

Mean total precipitation rate

> netcdf mtpr<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131 {
> dimensions:
> 	longitude = 1440 ;
> 	latitude = 721 ;
> 	time = 744 ;
> variables:
> 	float longitude(longitude) ;
> 		longitude:units = "degrees<sub>east</sub>" ;
> 		longitude:long<sub>name</sub> = "longitude" ;
> 	float latitude(latitude) ;
> 		latitude:units = "degrees<sub>north</sub>" ;
> 		latitude:long<sub>name</sub> = "latitude" ;
> 	int time(time) ;
> 		time:units = "hours since 1900-01-01 00:00:00.0" ;
> 		time:long<sub>name</sub> = "time" ;
> 		time:calendar = "gregorian" ;
> 	short mtpr(time, latitude, longitude) ;
> 		mtpr:scale<sub>factor</sub> = 3.71643823606165e-07 ;
> 		mtpr:add<sub>offset</sub> = 0.0121772815242796 ;
> 		mtpr:<sub>FillValue</sub> = -32767s ;
> 		mtpr:missing<sub>value</sub> = -32767s ;
> 		mtpr:units = "kg m\*\*-2 s\*\*-1" ;
> 		mtpr:long<sub>name</sub> = "Mean total precipitation rate" ;
> 
> // global attributes:
> 		:Conventions = "CF-1.6" ;
> 		:history = "2020-10-06 09:14:14 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: mv *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/mtpr/2010*.mtpr<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp *g/data/rt52/era5/single-levels/reanalysis/mtpr/2010/mtpr<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc\n",
> 			"2020-10-06 09:13:40 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: nccopy -k4 -d5 -s /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/mtpr/2010/mtpr<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc /g/data/rt52/admin/incoming/era5/single-levels/reanalysis/mtpr/2010*.mtpr<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp\n",
> 			"2020-10-06 09:07:33 UTC+1100 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: curl &#x2013;connect-timeout 20 &#x2013;show-error &#x2013;silent &#x2013;max-time 36000 -o /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/mtpr/2010/mtpr<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc <http://136.156.133.46/cache-compute-0015/cache/data8/adaptor.mars.internal-1601935342.505529-1432-7-5cc4d452-44f9-435e-8afc-b8d12664c4a7.nc\n>",
> 			"2020-10-05 22:03:07 GMT by grib<sub>to</sub><sub>netcdf</sub>-2.16.0: *opt/ecmwf/eccodes/bin/grib<sub>to</sub><sub>netcdf</sub> -S param -o /cache/data8/adaptor.mars.internal-1601935342.505529-1432-7-5cc4d452-44f9-435e-8afc-b8d12664c4a7.nc /cache/tmp/5cc4d452-44f9-435e-8afc-b8d12664c4a7-adaptor.mars.internal-1601935342.5060759-1432-3-tmp.grib" ;
> 		:license = "Licence to use Copernicus Products: <https://apps.ecmwf.int/datasets/licences/copernicus>*" ;
> 		:summary = "ERA5 is the fifth generation ECMWF atmospheric reanalysis of the global climate. This file is part of the ERA5 replica hosted at NCI Australia. For more information please see <http://dx.doi.org/10.25914/5f48874388857>" ;
> 		:title = "ERA5 single-levels reanalysis mean<sub>total</sub><sub>precipitation</sub><sub>rate</sub> 20100101-20100131" ;
> }


<a id="org99b2dfd"></a>

### [sp](file:///Volumes/ioa02/reanalysis/ERA5/sp)

Surface pressure

> netcdf sp<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131 {
> dimensions:
> 	longitude = 1440 ;
> 	latitude = 721 ;
> 	time = 744 ;
> variables:
> 	float longitude(longitude) ;
> 		longitude:units = "degrees<sub>east</sub>" ;
> 		longitude:long<sub>name</sub> = "longitude" ;
> 	float latitude(latitude) ;
> 		latitude:units = "degrees<sub>north</sub>" ;
> 		latitude:long<sub>name</sub> = "latitude" ;
> 	int time(time) ;
> 		time:units = "hours since 1900-01-01 00:00:00.0" ;
> 		time:long<sub>name</sub> = "time" ;
> 		time:calendar = "gregorian" ;
> 	short sp(time, latitude, longitude) ;
> 		sp:scale<sub>factor</sub> = 0.882383148089512 ;
> 		sp:add<sub>offset</sub> = 77345.396699051 ;
> 		sp:<sub>FillValue</sub> = -32767s ;
> 		sp:missing<sub>value</sub> = -32767s ;
> 		sp:units = "Pa" ;
> 		sp:long<sub>name</sub> = "Surface pressure" ;
> 		sp:standard<sub>name</sub> = "surface<sub>air</sub><sub>pressure</sub>" ;
> 
> // global attributes:
> 		:Conventions = "CF-1.6" ;
> 		:history = "2020-09-28 12:39:13 UTC+1000 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: mv *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/sp/2010*.sp<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp *g/data/rt52/era5/single-levels/reanalysis/sp/2010/sp<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc\n",
> 			"2020-09-28 12:38:26 UTC+1000 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: nccopy -k4 -d5 -s /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/sp/2010/sp<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc /g/data/rt52/admin/incoming/era5/single-levels/reanalysis/sp/2010*.sp<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp\n",
> 			"2020-09-28 12:34:44 UTC+1000 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: curl &#x2013;connect-timeout 20 &#x2013;show-error &#x2013;silent &#x2013;max-time 36000 -o /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/sp/2010/sp<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc <http://136.156.132.198/cache-compute-0003/cache/data4/adaptor.mars.internal-1601260248.6544905-25194-9-c80c2295-8d73-4f03-8c30-c28361fb17ff.nc\n>",
> 			"2020-09-28 02:32:17 GMT by grib<sub>to</sub><sub>netcdf</sub>-2.16.0: *opt/ecmwf/eccodes/bin/grib<sub>to</sub><sub>netcdf</sub> -S param -o /cache/data4/adaptor.mars.internal-1601260248.6544905-25194-9-c80c2295-8d73-4f03-8c30-c28361fb17ff.nc /cache/tmp/c80c2295-8d73-4f03-8c30-c28361fb17ff-adaptor.mars.internal-1601260248.6552908-25194-4-tmp.grib" ;
> 		:license = "Licence to use Copernicus Products: <https://apps.ecmwf.int/datasets/licences/copernicus>*" ;
> 		:summary = "ERA5 is the fifth generation ECMWF atmospheric reanalysis of the global climate. This file is part of the ERA5 replica hosted at NCI Australia. For more information please see <http://dx.doi.org/10.25914/5f48874388857>" ;
> 		:title = "ERA5 single-levels reanalysis surface<sub>pressure</sub> 20100101-20100131" ;
> }


<a id="org64aa233"></a>

### [z](file:///Volumes/ioa02/reanalysis/ERA5/z)

Geopotential

> netcdf z<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131 {
> dimensions:
> 	longitude = 1440 ;
> 	latitude = 721 ;
> 	time = 744 ;
> variables:
> 	float longitude(longitude) ;
> 		longitude:units = "degrees<sub>east</sub>" ;
> 		longitude:long<sub>name</sub> = "longitude" ;
> 	float latitude(latitude) ;
> 		latitude:units = "degrees<sub>north</sub>" ;
> 		latitude:long<sub>name</sub> = "latitude" ;
> 	int time(time) ;
> 		time:units = "hours since 1900-01-01 00:00:00.0" ;
> 		time:long<sub>name</sub> = "time" ;
> 		time:calendar = "gregorian" ;
> 	short z(time, latitude, longitude) ;
> 		z:scale<sub>factor</sub> = 0.897939420788 ;
> 		z:add<sub>offset</sub> = 28035.3894091959 ;
> 		z:<sub>FillValue</sub> = -32767s ;
> 		z:missing<sub>value</sub> = -32767s ;
> 		z:units = "m\*\*2 s\*\*-2" ;
> 		z:long<sub>name</sub> = "Geopotential" ;
> 		z:standard<sub>name</sub> = "geopotential" ;
> 
> // global attributes:
> 		:Conventions = "CF-1.6" ;
> 		:history = "2020-09-24 12:45:20 UTC+1000 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: mv *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/z/2010*.z<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp *g/data/rt52/era5/single-levels/reanalysis/z/2010/z<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc\n",
> 			"2020-09-24 12:45:06 UTC+1000 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: nccopy -k4 -d5 -s /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/z/2010/z<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc /g/data/rt52/admin/incoming/era5/single-levels/reanalysis/z/2010*.z<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.tmp\n",
> 			"2020-09-24 12:41:04 UTC+1000 by era5<sub>replication</sub><sub>tools</sub>-1.2.1: curl &#x2013;connect-timeout 20 &#x2013;show-error &#x2013;silent &#x2013;max-time 36000 -o /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/z/2010/z<sub>era5</sub><sub>oper</sub><sub>sfc</sub><sub>20100101</sub>-20100131.nc <http://136.156.133.41/cache-compute-0013/cache/data7/adaptor.mars.internal-1600915045.6589952-3753-13-24833832-41d7-4def-96fb-c55f6af6e05a.nc\n>",
> 			"2020-09-24 02:38:37 GMT by grib<sub>to</sub><sub>netcdf</sub>-2.16.0: *opt/ecmwf/eccodes/bin/grib<sub>to</sub><sub>netcdf</sub> -S param -o /cache/data7/adaptor.mars.internal-1600915045.6589952-3753-13-24833832-41d7-4def-96fb-c55f6af6e05a.nc /cache/tmp/24833832-41d7-4def-96fb-c55f6af6e05a-adaptor.mars.internal-1600915045.660192-3753-4-tmp.grib" ;
> 		:license = "Licence to use Copernicus Products: <https://apps.ecmwf.int/datasets/licences/copernicus>*" ;
> 		:summary = "ERA5 is the fifth generation ECMWF atmospheric reanalysis of the global climate. This file is part of the ERA5 replica hosted at NCI Australia. For more information please see <http://dx.doi.org/10.25914/5f48874388857>" ;
> 		:title = "ERA5 single-levels reanalysis geopotential 20100101-20100131" ;
> }


<a id="org60ab42d"></a>

## JRA55do

[JRA55do](https://climate.mri-jma.go.jp/pub/ocean/JRA55-do/) is a $1/2$ degree spatial and three-hour temporal resolution atmospheric
climate dataset. Version 1.5 has been downloaded to my [local repository](file:///Volumes/ioa01/reanalysis/JRA55do/) for the
following parameters. 


<a id="org9042284"></a>

### [huss](file:///Volumes/ioa01/reanalysis/JRA55do/huss)

Near-surface specific humidity

> netcdf huss<sub>input4MIPs</sub><sub>atmosphericState</sub><sub>OMIP</sub><sub>MRI</sub>-JRA55-do-1-5-0<sub>gr</sub><sub>201001010000</sub>-201012312100 {
> dimensions:
> 	time = UNLIMITED ; // (2920 currently)
> 	lat = 320 ;
> 	lon = 640 ;
> 	bnds = 2 ;
> variables:
> 	double time(time) ;
> 		time:bounds = "time<sub>bnds</sub>" ;
> 		time:units = "days since 1900-01-01 00:00:00" ;
> 		time:calendar = "gregorian" ;
> 		time:axis = "T" ;
> 		time:long<sub>name</sub> = "time" ;
> 		time:standard<sub>name</sub> = "time" ;
> 	double time<sub>bnds</sub>(time, bnds) ;
> 	double lat(lat) ;
> 		lat:bounds = "lat<sub>bnds</sub>" ;
> 		lat:units = "degrees<sub>north</sub>" ;
> 		lat:axis = "Y" ;
> 		lat:long<sub>name</sub> = "Latitude" ;
> 		lat:standard<sub>name</sub> = "latitude" ;
> 	double lat<sub>bnds</sub>(lat, bnds) ;
> 	double lon(lon) ;
> 		lon:bounds = "lon<sub>bnds</sub>" ;
> 		lon:units = "degrees<sub>east</sub>" ;
> 		lon:axis = "X" ;
> 		lon:long<sub>name</sub> = "Longitude" ;
> 		lon:standard<sub>name</sub> = "longitude" ;
> 	double lon<sub>bnds</sub>(lon, bnds) ;
> 	double height ;
> 		height:units = "m" ;
> 		height:axis = "Z" ;
> 		height:positive = "up" ;
> 		height:long<sub>name</sub> = "height" ;
> 		height:standard<sub>name</sub> = "height" ;
> 	float huss(time, lat, lon) ;
> 		huss:standard<sub>name</sub> = "specific<sub>humidity</sub>" ;
> 		huss:long<sub>name</sub> = "Near-Surface Specific Humidity" ;
> 		huss:comment = "Near-surface (usually, 2 meter) specific humidity" ;
> 		huss:units = "1" ;
> 		huss:original<sub>units</sub> = "1.0" ;
> 		huss:history = "2020-09-15T15:01:10Z altered by CMOR: Converted units from \\'1.0\\' to \\'1\\'. 2020-09-15T15:01:10Z altered by CMOR: Treated scalar dimension: \\'height\\'." ;
> 		huss:cell<sub>methods</sub> = "area: mean time: point" ;
> 		huss:cell<sub>measures</sub> = "area: areacella" ;
> 		huss:coordinates = "height" ;
> 		huss:missing<sub>value</sub> = 1.e+20f ;
> 		huss:<sub>FillValue</sub> = 1.e+20f ;
> 
> // global attributes:
> 		:Conventions = "CF-1.7 CMIP-6.2" ;
> 		:activity<sub>id</sub> = "input4MIPs" ;
> 		:cell<sub>measures</sub> = "area: areacella" ;
> 		:comment = "Based on JRA-55 reanalysis (1958-01 to 2020-07)" ;
> 		:contact = "Hiroyuki Tsujino (htsujino@mri-jma.go.jp)" ;
> 		:creation<sub>date</sub> = "2020-09-15T15:01:18Z" ;
> 		:data<sub>specs</sub><sub>version</sub> = "01.00.32" ;
> 		:dataset<sub>category</sub> = "atmosphericState" ;
> 		:external<sub>variables</sub> = "areacella" ;
> 		:frequency = "3hrPt" ;
> 		:further<sub>info</sub><sub>url</sub> = "<http://climate.mri-jma.go.jp/~htsujino/jra55do.html>" ;
> 		:grid = "data regridded to the normal atmosphere TL319 gaussian grid (320x640 latxlon) from a reduced TL319 gaussian grid" ;
> 		:grid<sub>label</sub> = "gr" ;
> 		:history = "2020-09-15T15:01:18Z; CMOR rewrote data to be consistent with input4MIPs, CMIP6, and CF-1.7 CMIP-6.2 standards" ;
> 		:institution = "Meteorological Research Institute, Tsukuba, Ibaraki 305-0052, Japan" ;
> 		:institution<sub>id</sub> = "MRI" ;
> 		:mip<sub>era</sub> = "CMIP6" ;
> 		:nominal<sub>resolution</sub> = "50 km" ;
> 		:product = "reanalysis" ;
> 		:realm = "atmos" ;
> 		:references = "Tsujino et al., 2018: JRA-55 based surface dataset for driving ocean-sea-ice models (JRA55-do), Ocean Modelling, 130(1), pp 79-139. <https://doi.org/10.1016/j.ocemod.2018.07.002>" ;
> 		:region = "global<sub>ocean</sub>" ;
> 		:release<sub>year</sub> = "2020" ;
> 		:source = "MRI JRA55-do 1.5.0: Atmospheric state generated for OMIP based on the JRA-55 reanalysis" ;
> 		:source<sub>description</sub> = "Atmospheric state and terrestrial runoff datasets produced by MRI for the OMIP experiment of CMIP6" ;
> 		:source<sub>id</sub> = "MRI-JRA55-do-1-5-0" ;
> 		:source<sub>type</sub> = "satellite<sub>blended</sub>" ;
> 		:source<sub>version</sub> = "1.5.0" ;
> 		:table<sub>id</sub> = "input4MIPs<sub>A3hrPt</sub>" ;
> 		:table<sub>info</sub> = "Creation Date:(14 September 2020) MD5:d10deb7a45fa93934ec0a16c9d534e45" ;
> 		:target<sub>mip</sub> = "OMIP" ;
> 		:title = "MRI JRA55-do 1.5.0 dataset prepared for input4MIPs" ;
> 		:tracking<sub>id</sub> = "hdl:21.14100/369cad98-ed94-4040-a47d-72258fd68584" ;
> 		:variable<sub>id</sub> = "huss" ;
> 		:license = "OMIP boundary condition data produced by MRI is licensed under a Creative Commons Attribution-[NonCommercial-]ShareAlike 4.0 International License (<https://creativecommons.org/licenses>). Consult <https://pcmdi.llnl.gov/CMIP6/TermsOfUse> for terms of use governing input4MIPs output, including citation requirements and proper acknowledgment. Further information about this data, including some limitations, can be found via the further<sub>info</sub><sub>url</sub> (recorded as a global attribute in this file). The data producers and data providers make no warranty, either express or implied, including, but not limited to, warranties of merchantability and fitness for a particular purpose. All liabilities arising from the supply of the information (including any liability arising in negligence) are excluded to the fullest extent permitted by law." ;
> 		:cmor<sub>version</sub> = "3.6.0" ;
> }


<a id="org41569e1"></a>

### [prra](file:///Volumes/ioa01/reanalysis/JRA55do/prra)

Rainfall flux

> netcdf prra<sub>input4MIPs</sub><sub>atmosphericState</sub><sub>OMIP</sub><sub>MRI</sub>-JRA55-do-1-5-0<sub>gr</sub><sub>201001010130</sub>-201012312230 {
> dimensions:
> 	time = UNLIMITED ; // (2920 currently)
> 	lat = 320 ;
> 	lon = 640 ;
> 	bnds = 2 ;
> variables:
> 	double time(time) ;
> 		time:bounds = "time<sub>bnds</sub>" ;
> 		time:units = "days since 1900-01-01 00:00:00" ;
> 		time:calendar = "gregorian" ;
> 		time:axis = "T" ;
> 		time:long<sub>name</sub> = "time" ;
> 		time:standard<sub>name</sub> = "time" ;
> 	double time<sub>bnds</sub>(time, bnds) ;
> 	double lat(lat) ;
> 		lat:bounds = "lat<sub>bnds</sub>" ;
> 		lat:units = "degrees<sub>north</sub>" ;
> 		lat:axis = "Y" ;
> 		lat:long<sub>name</sub> = "Latitude" ;
> 		lat:standard<sub>name</sub> = "latitude" ;
> 	double lat<sub>bnds</sub>(lat, bnds) ;
> 	double lon(lon) ;
> 		lon:bounds = "lon<sub>bnds</sub>" ;
> 		lon:units = "degrees<sub>east</sub>" ;
> 		lon:axis = "X" ;
> 		lon:long<sub>name</sub> = "Longitude" ;
> 		lon:standard<sub>name</sub> = "longitude" ;
> 	double lon<sub>bnds</sub>(lon, bnds) ;
> 	float prra(time, lat, lon) ;
> 		prra:standard<sub>name</sub> = "rainfall<sub>flux</sub>" ;
> 		prra:long<sub>name</sub> = "Rainfall Flux" ;
> 		prra:comment = "In accordance with common usage in geophysical disciplines, \\'flux\\' implies per unit area, called \\'flux density\\' in physics" ;
> 		prra:units = "kg m-2 s-1" ;
> 		prra:cell<sub>methods</sub> = "area: time: mean" ;
> 		prra:cell<sub>measures</sub> = "area: areacella" ;
> 		prra:missing<sub>value</sub> = 1.e+20f ;
> 		prra:<sub>FillValue</sub> = 1.e+20f ;
> 
> // global attributes:
> 		:Conventions = "CF-1.7 CMIP-6.2" ;
> 		:activity<sub>id</sub> = "input4MIPs" ;
> 		:cell<sub>measures</sub> = "area: areacella" ;
> 		:comment = "Based on JRA-55 reanalysis (1958-01 to 2020-07)" ;
> 		:contact = "Hiroyuki Tsujino (htsujino@mri-jma.go.jp)" ;
> 		:creation<sub>date</sub> = "2020-09-15T06:03:44Z" ;
> 		:data<sub>specs</sub><sub>version</sub> = "01.00.32" ;
> 		:dataset<sub>category</sub> = "atmosphericState" ;
> 		:external<sub>variables</sub> = "areacella" ;
> 		:frequency = "3hr" ;
> 		:further<sub>info</sub><sub>url</sub> = "<http://climate.mri-jma.go.jp/~htsujino/jra55do.html>" ;
> 		:grid = "data regridded to the normal atmosphere TL319 gaussian grid (320x640 latxlon) from a reduced TL319 gaussian grid" ;
> 		:grid<sub>label</sub> = "gr" ;
> 		:history = "2020-09-15T06:03:44Z; CMOR rewrote data to be consistent with input4MIPs, CMIP6, and CF-1.7 CMIP-6.2 standards" ;
> 		:institution = "Meteorological Research Institute, Tsukuba, Ibaraki 305-0052, Japan" ;
> 		:institution<sub>id</sub> = "MRI" ;
> 		:mip<sub>era</sub> = "CMIP6" ;
> 		:nominal<sub>resolution</sub> = "50 km" ;
> 		:product = "reanalysis" ;
> 		:realm = "atmos" ;
> 		:references = "Tsujino et al., 2018: JRA-55 based surface dataset for driving ocean-sea-ice models (JRA55-do), Ocean Modelling, 130(1), pp 79-139. <https://doi.org/10.1016/j.ocemod.2018.07.002>" ;
> 		:region = "global<sub>ocean</sub>" ;
> 		:release<sub>year</sub> = "2020" ;
> 		:source = "MRI JRA55-do 1.5.0: Atmospheric state generated for OMIP based on the JRA-55 reanalysis" ;
> 		:source<sub>description</sub> = "Atmospheric state and terrestrial runoff datasets produced by MRI for the OMIP experiment of CMIP6" ;
> 		:source<sub>id</sub> = "MRI-JRA55-do-1-5-0" ;
> 		:source<sub>type</sub> = "satellite<sub>blended</sub>" ;
> 		:source<sub>version</sub> = "1.5.0" ;
> 		:table<sub>id</sub> = "input4MIPs<sub>A3hr</sub>" ;
> 		:table<sub>info</sub> = "Creation Date:(14 September 2020) MD5:d10deb7a45fa93934ec0a16c9d534e45" ;
> 		:target<sub>mip</sub> = "OMIP" ;
> 		:title = "MRI JRA55-do 1.5.0 dataset prepared for input4MIPs" ;
> 		:tracking<sub>id</sub> = "hdl:21.14100/0ebc0818-9a24-47ad-a6d5-112506908d87" ;
> 		:variable<sub>id</sub> = "prra" ;
> 		:license = "OMIP boundary condition data produced by MRI is licensed under a Creative Commons Attribution-[NonCommercial-]ShareAlike 4.0 International License (<https://creativecommons.org/licenses>). Consult <https://pcmdi.llnl.gov/CMIP6/TermsOfUse> for terms of use governing input4MIPs output, including citation requirements and proper acknowledgment. Further information about this data, including some limitations, can be found via the further<sub>info</sub><sub>url</sub> (recorded as a global attribute in this file). The data producers and data providers make no warranty, either express or implied, including, but not limited to, warranties of merchantability and fitness for a particular purpose. All liabilities arising from the supply of the information (including any liability arising in negligence) are excluded to the fullest extent permitted by law." ;
> 		:cmor<sub>version</sub> = "3.6.0" ;
> }


<a id="org58abefe"></a>

### [prsn](file:///Volumes/ioa01/reanalysis/JRA55do/prsn/)

Snowfall flux

> netcdf prsn<sub>input4MIPs</sub><sub>atmosphericState</sub><sub>OMIP</sub><sub>MRI</sub>-JRA55-do-1-5-0<sub>gr</sub><sub>201001010130</sub>-201012312230 {
> dimensions:
> 	time = UNLIMITED ; // (2920 currently)
> 	lat = 320 ;
> 	lon = 640 ;
> 	bnds = 2 ;
> variables:
> 	double time(time) ;
> 		time:bounds = "time<sub>bnds</sub>" ;
> 		time:units = "days since 1900-01-01 00:00:00" ;
> 		time:calendar = "gregorian" ;
> 		time:axis = "T" ;
> 		time:long<sub>name</sub> = "time" ;
> 		time:standard<sub>name</sub> = "time" ;
> 	double time<sub>bnds</sub>(time, bnds) ;
> 	double lat(lat) ;
> 		lat:bounds = "lat<sub>bnds</sub>" ;
> 		lat:units = "degrees<sub>north</sub>" ;
> 		lat:axis = "Y" ;
> 		lat:long<sub>name</sub> = "Latitude" ;
> 		lat:standard<sub>name</sub> = "latitude" ;
> 	double lat<sub>bnds</sub>(lat, bnds) ;
> 	double lon(lon) ;
> 		lon:bounds = "lon<sub>bnds</sub>" ;
> 		lon:units = "degrees<sub>east</sub>" ;
> 		lon:axis = "X" ;
> 		lon:long<sub>name</sub> = "Longitude" ;
> 		lon:standard<sub>name</sub> = "longitude" ;
> 	double lon<sub>bnds</sub>(lon, bnds) ;
> 	float prsn(time, lat, lon) ;
> 		prsn:standard<sub>name</sub> = "snowfall<sub>flux</sub>" ;
> 		prsn:long<sub>name</sub> = "Snowfall Flux" ;
> 		prsn:comment = "At surface; includes precipitation of all forms of water in the solid phase" ;
> 		prsn:units = "kg m-2 s-1" ;
> 		prsn:cell<sub>methods</sub> = "area: time: mean" ;
> 		prsn:cell<sub>measures</sub> = "area: areacella" ;
> 		prsn:missing<sub>value</sub> = 1.e+20f ;
> 		prsn:<sub>FillValue</sub> = 1.e+20f ;
> 
> // global attributes:
> 		:Conventions = "CF-1.7 CMIP-6.2" ;
> 		:activity<sub>id</sub> = "input4MIPs" ;
> 		:cell<sub>measures</sub> = "area: areacella" ;
> 		:comment = "Based on JRA-55 reanalysis (1958-01 to 2020-07)" ;
> 		:contact = "Hiroyuki Tsujino (htsujino@mri-jma.go.jp)" ;
> 		:creation<sub>date</sub> = "2020-09-15T07:53:10Z" ;
> 		:data<sub>specs</sub><sub>version</sub> = "01.00.32" ;
> 		:dataset<sub>category</sub> = "atmosphericState" ;
> 		:external<sub>variables</sub> = "areacella" ;
> 		:frequency = "3hr" ;
> 		:further<sub>info</sub><sub>url</sub> = "<http://climate.mri-jma.go.jp/~htsujino/jra55do.html>" ;
> 		:grid = "data regridded to the normal atmosphere TL319 gaussian grid (320x640 latxlon) from a reduced TL319 gaussian grid" ;
> 		:grid<sub>label</sub> = "gr" ;
> 		:history = "2020-09-15T07:53:10Z; CMOR rewrote data to be consistent with input4MIPs, CMIP6, and CF-1.7 CMIP-6.2 standards" ;
> 		:institution = "Meteorological Research Institute, Tsukuba, Ibaraki 305-0052, Japan" ;
> 		:institution<sub>id</sub> = "MRI" ;
> 		:mip<sub>era</sub> = "CMIP6" ;
> 		:nominal<sub>resolution</sub> = "50 km" ;
> 		:product = "reanalysis" ;
> 		:realm = "atmos" ;
> 		:references = "Tsujino et al., 2018: JRA-55 based surface dataset for driving ocean-sea-ice models (JRA55-do), Ocean Modelling, 130(1), pp 79-139. <https://doi.org/10.1016/j.ocemod.2018.07.002>" ;
> 		:region = "global<sub>ocean</sub>" ;
> 		:release<sub>year</sub> = "2020" ;
> 		:source = "MRI JRA55-do 1.5.0: Atmospheric state generated for OMIP based on the JRA-55 reanalysis" ;
> 		:source<sub>description</sub> = "Atmospheric state and terrestrial runoff datasets produced by MRI for the OMIP experiment of CMIP6" ;
> 		:source<sub>id</sub> = "MRI-JRA55-do-1-5-0" ;
> 		:source<sub>type</sub> = "satellite<sub>blended</sub>" ;
> 		:source<sub>version</sub> = "1.5.0" ;
> 		:table<sub>id</sub> = "input4MIPs<sub>A3hr</sub>" ;
> 		:table<sub>info</sub> = "Creation Date:(14 September 2020) MD5:d10deb7a45fa93934ec0a16c9d534e45" ;
> 		:target<sub>mip</sub> = "OMIP" ;
> 		:title = "MRI JRA55-do 1.5.0 dataset prepared for input4MIPs" ;
> 		:tracking<sub>id</sub> = "hdl:21.14100/11edbd94-09f9-4e92-b307-c69ee1238c4e" ;
> 		:variable<sub>id</sub> = "prsn" ;
> 		:license = "OMIP boundary condition data produced by MRI is licensed under a Creative Commons Attribution-[NonCommercial-]ShareAlike 4.0 International License (<https://creativecommons.org/licenses>). Consult <https://pcmdi.llnl.gov/CMIP6/TermsOfUse> for terms of use governing input4MIPs output, including citation requirements and proper acknowledgment. Further information about this data, including some limitations, can be found via the further<sub>info</sub><sub>url</sub> (recorded as a global attribute in this file). The data producers and data providers make no warranty, either express or implied, including, but not limited to, warranties of merchantability and fitness for a particular purpose. All liabilities arising from the supply of the information (including any liability arising in negligence) are excluded to the fullest extent permitted by law." ;
> 		:cmor<sub>version</sub> = "3.6.0" ;
> }


<a id="orgc49e63f"></a>

### [psl](file:///Volumes/ioa01/reanalysis/JRA55do/psl)

Mean surface pressure

> netcdf psl<sub>input4MIPs</sub><sub>atmosphericState</sub><sub>OMIP</sub><sub>MRI</sub>-JRA55-do-1-5-0<sub>gr</sub><sub>201001010000</sub>-201012312100 {
> dimensions:
> 	time = UNLIMITED ; // (2920 currently)
> 	lat = 320 ;
> 	lon = 640 ;
> 	bnds = 2 ;
> variables:
> 	double time(time) ;
> 		time:bounds = "time<sub>bnds</sub>" ;
> 		time:units = "days since 1900-01-01 00:00:00" ;
> 		time:calendar = "gregorian" ;
> 		time:axis = "T" ;
> 		time:long<sub>name</sub> = "time" ;
> 		time:standard<sub>name</sub> = "time" ;
> 	double time<sub>bnds</sub>(time, bnds) ;
> 	double lat(lat) ;
> 		lat:bounds = "lat<sub>bnds</sub>" ;
> 		lat:units = "degrees<sub>north</sub>" ;
> 		lat:axis = "Y" ;
> 		lat:long<sub>name</sub> = "Latitude" ;
> 		lat:standard<sub>name</sub> = "latitude" ;
> 	double lat<sub>bnds</sub>(lat, bnds) ;
> 	double lon(lon) ;
> 		lon:bounds = "lon<sub>bnds</sub>" ;
> 		lon:units = "degrees<sub>east</sub>" ;
> 		lon:axis = "X" ;
> 		lon:long<sub>name</sub> = "Longitude" ;
> 		lon:standard<sub>name</sub> = "longitude" ;
> 	double lon<sub>bnds</sub>(lon, bnds) ;
> 	float psl(time, lat, lon) ;
> 		psl:standard<sub>name</sub> = "air<sub>pressure</sub><sub>at</sub><sub>mean</sub><sub>sea</sub><sub>level</sub>" ;
> 		psl:long<sub>name</sub> = "Sea Level Pressure" ;
> 		psl:comment = "Sea Level Pressure" ;
> 		psl:units = "Pa" ;
> 		psl:cell<sub>methods</sub> = "area: mean time: point" ;
> 		psl:cell<sub>measures</sub> = "area: areacella" ;
> 		psl:missing<sub>value</sub> = 1.e+20f ;
> 		psl:<sub>FillValue</sub> = 1.e+20f ;
> 
> // global attributes:
> 		:Conventions = "CF-1.7 CMIP-6.2" ;
> 		:activity<sub>id</sub> = "input4MIPs" ;
> 		:cell<sub>measures</sub> = "area: areacella" ;
> 		:comment = "Based on JRA-55 reanalysis (1958-01 to 2020-07)" ;
> 		:contact = "Hiroyuki Tsujino (htsujino@mri-jma.go.jp)" ;
> 		:creation<sub>date</sub> = "2020-09-15T17:20:36Z" ;
> 		:data<sub>specs</sub><sub>version</sub> = "01.00.32" ;
> 		:dataset<sub>category</sub> = "atmosphericState" ;
> 		:external<sub>variables</sub> = "areacella" ;
> 		:frequency = "3hrPt" ;
> 		:further<sub>info</sub><sub>url</sub> = "<http://climate.mri-jma.go.jp/~htsujino/jra55do.html>" ;
> 		:grid = "data regridded to the normal atmosphere TL319 gaussian grid (320x640 latxlon) from a reduced TL319 gaussian grid" ;
> 		:grid<sub>label</sub> = "gr" ;
> 		:history = "2020-09-15T17:20:36Z; CMOR rewrote data to be consistent with input4MIPs, CMIP6, and CF-1.7 CMIP-6.2 standards" ;
> 		:institution = "Meteorological Research Institute, Tsukuba, Ibaraki 305-0052, Japan" ;
> 		:institution<sub>id</sub> = "MRI" ;
> 		:mip<sub>era</sub> = "CMIP6" ;
> 		:nominal<sub>resolution</sub> = "50 km" ;
> 		:product = "reanalysis" ;
> 		:realm = "atmos" ;
> 		:references = "Tsujino et al., 2018: JRA-55 based surface dataset for driving ocean-sea-ice models (JRA55-do), Ocean Modelling, 130(1), pp 79-139. <https://doi.org/10.1016/j.ocemod.2018.07.002>" ;
> 		:region = "global<sub>ocean</sub>" ;
> 		:release<sub>year</sub> = "2020" ;
> 		:source = "MRI JRA55-do 1.5.0: Atmospheric state generated for OMIP based on the JRA-55 reanalysis" ;
> 		:source<sub>description</sub> = "Atmospheric state and terrestrial runoff datasets produced by MRI for the OMIP experiment of CMIP6" ;
> 		:source<sub>id</sub> = "MRI-JRA55-do-1-5-0" ;
> 		:source<sub>type</sub> = "satellite<sub>blended</sub>" ;
> 		:source<sub>version</sub> = "1.5.0" ;
> 		:table<sub>id</sub> = "input4MIPs<sub>A3hrPt</sub>" ;
> 		:table<sub>info</sub> = "Creation Date:(14 September 2020) MD5:d10deb7a45fa93934ec0a16c9d534e45" ;
> 		:target<sub>mip</sub> = "OMIP" ;
> 		:title = "MRI JRA55-do 1.5.0 dataset prepared for input4MIPs" ;
> 		:tracking<sub>id</sub> = "hdl:21.14100/4e5646ac-45d5-40cd-a999-1554bb4cd75a" ;
> 		:variable<sub>id</sub> = "psl" ;
> 		:license = "OMIP boundary condition data produced by MRI is licensed under a Creative Commons Attribution-[NonCommercial-]ShareAlike 4.0 International License (<https://creativecommons.org/licenses>). Consult <https://pcmdi.llnl.gov/CMIP6/TermsOfUse> for terms of use governing input4MIPs output, including citation requirements and proper acknowledgment. Further information about this data, including some limitations, can be found via the further<sub>info</sub><sub>url</sub> (recorded as a global attribute in this file). The data producers and data providers make no warranty, either express or implied, including, but not limited to, warranties of merchantability and fitness for a particular purpose. All liabilities arising from the supply of the information (including any liability arising in negligence) are excluded to the fullest extent permitted by law." ;
> 		:cmor<sub>version</sub> = "3.6.0" ;
> }


<a id="orgc84a418"></a>

### [rlds](file:///Volumes/ioa01/reanalysis/JRA55do/rlds)

Surface downwelling longwave flux

> netcdf rlds<sub>input4MIPs</sub><sub>atmosphericState</sub><sub>OMIP</sub><sub>MRI</sub>-JRA55-do-1-5-0<sub>gr</sub><sub>201001010130</sub>-201012312230 {
> dimensions:
> 	time = UNLIMITED ; // (2920 currently)
> 	lat = 320 ;
> 	lon = 640 ;
> 	bnds = 2 ;
> variables:
> 	double time(time) ;
> 		time:bounds = "time<sub>bnds</sub>" ;
> 		time:units = "days since 1900-01-01 00:00:00" ;
> 		time:calendar = "gregorian" ;
> 		time:axis = "T" ;
> 		time:long<sub>name</sub> = "time" ;
> 		time:standard<sub>name</sub> = "time" ;
> 	double time<sub>bnds</sub>(time, bnds) ;
> 	double lat(lat) ;
> 		lat:bounds = "lat<sub>bnds</sub>" ;
> 		lat:units = "degrees<sub>north</sub>" ;
> 		lat:axis = "Y" ;
> 		lat:long<sub>name</sub> = "Latitude" ;
> 		lat:standard<sub>name</sub> = "latitude" ;
> 	double lat<sub>bnds</sub>(lat, bnds) ;
> 	double lon(lon) ;
> 		lon:bounds = "lon<sub>bnds</sub>" ;
> 		lon:units = "degrees<sub>east</sub>" ;
> 		lon:axis = "X" ;
> 		lon:long<sub>name</sub> = "Longitude" ;
> 		lon:standard<sub>name</sub> = "longitude" ;
> 	double lon<sub>bnds</sub>(lon, bnds) ;
> 	float rlds(time, lat, lon) ;
> 		rlds:standard<sub>name</sub> = "surface<sub>downwelling</sub><sub>longwave</sub><sub>flux</sub><sub>in</sub><sub>air</sub>" ;
> 		rlds:long<sub>name</sub> = "Surface Downwelling Longwave Radiation" ;
> 		rlds:comment = "The surface called \\'surface\\' means the lower boundary of the atmosphere. \\'longwave\\' means longwave radiation. Downwelling radiation is radiation from above. It does not mean \\'net downward\\'. When thought of as being incident on a surface, a radiative flux is sometimes called \\'irradiance\\'. In addition, it is identical with the quantity measured by a cosine-collector light-meter and sometimes called \\'vector irradiance\\'. In accordance with common usage in geophysical disciplines, \\'flux\\' implies per unit area, called \\'flux density\\' in physics." ;
> 		rlds:units = "W m-2" ;
> 		rlds:cell<sub>methods</sub> = "area: time: mean" ;
> 		rlds:cell<sub>measures</sub> = "area: areacella" ;
> 		rlds:missing<sub>value</sub> = 1.e+20f ;
> 		rlds:<sub>FillValue</sub> = 1.e+20f ;
> 
> // global attributes:
> 		:Conventions = "CF-1.7 CMIP-6.2" ;
> 		:activity<sub>id</sub> = "input4MIPs" ;
> 		:cell<sub>measures</sub> = "area: areacella" ;
> 		:comment = "Based on JRA-55 reanalysis (1958-01 to 2020-07)" ;
> 		:contact = "Hiroyuki Tsujino (htsujino@mri-jma.go.jp)" ;
> 		:creation<sub>date</sub> = "2020-09-15T10:15:08Z" ;
> 		:data<sub>specs</sub><sub>version</sub> = "01.00.32" ;
> 		:dataset<sub>category</sub> = "atmosphericState" ;
> 		:external<sub>variables</sub> = "areacella" ;
> 		:frequency = "3hr" ;
> 		:further<sub>info</sub><sub>url</sub> = "<http://climate.mri-jma.go.jp/~htsujino/jra55do.html>" ;
> 		:grid = "data regridded to the normal atmosphere TL319 gaussian grid (320x640 latxlon) from a reduced TL319 gaussian grid" ;
> 		:grid<sub>label</sub> = "gr" ;
> 		:history = "2020-09-15T10:15:08Z; CMOR rewrote data to be consistent with input4MIPs, CMIP6, and CF-1.7 CMIP-6.2 standards" ;
> 		:institution = "Meteorological Research Institute, Tsukuba, Ibaraki 305-0052, Japan" ;
> 		:institution<sub>id</sub> = "MRI" ;
> 		:mip<sub>era</sub> = "CMIP6" ;
> 		:nominal<sub>resolution</sub> = "50 km" ;
> 		:product = "reanalysis" ;
> 		:realm = "atmos" ;
> 		:references = "Tsujino et al., 2018: JRA-55 based surface dataset for driving ocean-sea-ice models (JRA55-do), Ocean Modelling, 130(1), pp 79-139. <https://doi.org/10.1016/j.ocemod.2018.07.002>" ;
> 		:region = "global<sub>ocean</sub>" ;
> 		:release<sub>year</sub> = "2020" ;
> 		:source = "MRI JRA55-do 1.5.0: Atmospheric state generated for OMIP based on the JRA-55 reanalysis" ;
> 		:source<sub>description</sub> = "Atmospheric state and terrestrial runoff datasets produced by MRI for the OMIP experiment of CMIP6" ;
> 		:source<sub>id</sub> = "MRI-JRA55-do-1-5-0" ;
> 		:source<sub>type</sub> = "satellite<sub>blended</sub>" ;
> 		:source<sub>version</sub> = "1.5.0" ;
> 		:table<sub>id</sub> = "input4MIPs<sub>A3hr</sub>" ;
> 		:table<sub>info</sub> = "Creation Date:(14 September 2020) MD5:d10deb7a45fa93934ec0a16c9d534e45" ;
> 		:target<sub>mip</sub> = "OMIP" ;
> 		:title = "MRI JRA55-do 1.5.0 dataset prepared for input4MIPs" ;
> 		:tracking<sub>id</sub> = "hdl:21.14100/9a1230a5-bb91-4cb5-9c99-179ba031866e" ;
> 		:variable<sub>id</sub> = "rlds" ;
> 		:license = "OMIP boundary condition data produced by MRI is licensed under a Creative Commons Attribution-[NonCommercial-]ShareAlike 4.0 International License (<https://creativecommons.org/licenses>). Consult <https://pcmdi.llnl.gov/CMIP6/TermsOfUse> for terms of use governing input4MIPs output, including citation requirements and proper acknowledgment. Further information about this data, including some limitations, can be found via the further<sub>info</sub><sub>url</sub> (recorded as a global attribute in this file). The data producers and data providers make no warranty, either express or implied, including, but not limited to, warranties of merchantability and fitness for a particular purpose. All liabilities arising from the supply of the information (including any liability arising in negligence) are excluded to the fullest extent permitted by law." ;
> 		:cmor<sub>version</sub> = "3.6.0" ;
> }


<a id="orgf231b7e"></a>

### [rsds](file:///Volumes/ioa01/reanalysis/JRA55do/rsds)

Surface downwelling shortwave flux

> netcdf rsds<sub>input4MIPs</sub><sub>atmosphericState</sub><sub>OMIP</sub><sub>MRI</sub>-JRA55-do-1-5-0<sub>gr</sub><sub>201001010130</sub>-201012312230 {
> dimensions:
> 	time = UNLIMITED ; // (2920 currently)
> 	lat = 320 ;
> 	lon = 640 ;
> 	bnds = 2 ;
> variables:
> 	double time(time) ;
> 		time:bounds = "time<sub>bnds</sub>" ;
> 		time:units = "days since 1900-01-01 00:00:00" ;
> 		time:calendar = "gregorian" ;
> 		time:axis = "T" ;
> 		time:long<sub>name</sub> = "time" ;
> 		time:standard<sub>name</sub> = "time" ;
> 	double time<sub>bnds</sub>(time, bnds) ;
> 	double lat(lat) ;
> 		lat:bounds = "lat<sub>bnds</sub>" ;
> 		lat:units = "degrees<sub>north</sub>" ;
> 		lat:axis = "Y" ;
> 		lat:long<sub>name</sub> = "Latitude" ;
> 		lat:standard<sub>name</sub> = "latitude" ;
> 	double lat<sub>bnds</sub>(lat, bnds) ;
> 	double lon(lon) ;
> 		lon:bounds = "lon<sub>bnds</sub>" ;
> 		lon:units = "degrees<sub>east</sub>" ;
> 		lon:axis = "X" ;
> 		lon:long<sub>name</sub> = "Longitude" ;
> 		lon:standard<sub>name</sub> = "longitude" ;
> 	double lon<sub>bnds</sub>(lon, bnds) ;
> 	float rsds(time, lat, lon) ;
> 		rsds:standard<sub>name</sub> = "surface<sub>downwelling</sub><sub>shortwave</sub><sub>flux</sub><sub>in</sub><sub>air</sub>" ;
> 		rsds:long<sub>name</sub> = "Surface Downwelling Shortwave Radiation" ;
> 		rsds:comment = "Surface solar irradiance for UV calculations." ;
> 		rsds:units = "W m-2" ;
> 		rsds:cell<sub>methods</sub> = "area: time: mean" ;
> 		rsds:cell<sub>measures</sub> = "area: areacella" ;
> 		rsds:missing<sub>value</sub> = 1.e+20f ;
> 		rsds:<sub>FillValue</sub> = 1.e+20f ;
> 
> // global attributes:
> 		:Conventions = "CF-1.7 CMIP-6.2" ;
> 		:activity<sub>id</sub> = "input4MIPs" ;
> 		:cell<sub>measures</sub> = "area: areacella" ;
> 		:comment = "Based on JRA-55 reanalysis (1958-01 to 2020-07)" ;
> 		:contact = "Hiroyuki Tsujino (htsujino@mri-jma.go.jp)" ;
> 		:creation<sub>date</sub> = "2020-09-15T12:28:57Z" ;
> 		:data<sub>specs</sub><sub>version</sub> = "01.00.32" ;
> 		:dataset<sub>category</sub> = "atmosphericState" ;
> 		:external<sub>variables</sub> = "areacella" ;
> 		:frequency = "3hr" ;
> 		:further<sub>info</sub><sub>url</sub> = "<http://climate.mri-jma.go.jp/~htsujino/jra55do.html>" ;
> 		:grid = "data regridded to the normal atmosphere TL319 gaussian grid (320x640 latxlon) from a reduced TL319 gaussian grid" ;
> 		:grid<sub>label</sub> = "gr" ;
> 		:history = "2020-09-15T12:28:57Z; CMOR rewrote data to be consistent with input4MIPs, CMIP6, and CF-1.7 CMIP-6.2 standards" ;
> 		:institution = "Meteorological Research Institute, Tsukuba, Ibaraki 305-0052, Japan" ;
> 		:institution<sub>id</sub> = "MRI" ;
> 		:mip<sub>era</sub> = "CMIP6" ;
> 		:nominal<sub>resolution</sub> = "50 km" ;
> 		:product = "reanalysis" ;
> 		:realm = "atmos" ;
> 		:references = "Tsujino et al., 2018: JRA-55 based surface dataset for driving ocean-sea-ice models (JRA55-do), Ocean Modelling, 130(1), pp 79-139. <https://doi.org/10.1016/j.ocemod.2018.07.002>" ;
> 		:region = "global<sub>ocean</sub>" ;
> 		:release<sub>year</sub> = "2020" ;
> 		:source = "MRI JRA55-do 1.5.0: Atmospheric state generated for OMIP based on the JRA-55 reanalysis" ;
> 		:source<sub>description</sub> = "Atmospheric state and terrestrial runoff datasets produced by MRI for the OMIP experiment of CMIP6" ;
> 		:source<sub>id</sub> = "MRI-JRA55-do-1-5-0" ;
> 		:source<sub>type</sub> = "satellite<sub>blended</sub>" ;
> 		:source<sub>version</sub> = "1.5.0" ;
> 		:table<sub>id</sub> = "input4MIPs<sub>A3hr</sub>" ;
> 		:table<sub>info</sub> = "Creation Date:(14 September 2020) MD5:d10deb7a45fa93934ec0a16c9d534e45" ;
> 		:target<sub>mip</sub> = "OMIP" ;
> 		:title = "MRI JRA55-do 1.5.0 dataset prepared for input4MIPs" ;
> 		:tracking<sub>id</sub> = "hdl:21.14100/f376d4e4-6396-44b4-ab8b-c7220a1afcb3" ;
> 		:variable<sub>id</sub> = "rsds" ;
> 		:license = "OMIP boundary condition data produced by MRI is licensed under a Creative Commons Attribution-[NonCommercial-]ShareAlike 4.0 International License (<https://creativecommons.org/licenses>). Consult <https://pcmdi.llnl.gov/CMIP6/TermsOfUse> for terms of use governing input4MIPs output, including citation requirements and proper acknowledgment. Further information about this data, including some limitations, can be found via the further<sub>info</sub><sub>url</sub> (recorded as a global attribute in this file). The data producers and data providers make no warranty, either express or implied, including, but not limited to, warranties of merchantability and fitness for a particular purpose. All liabilities arising from the supply of the information (including any liability arising in negligence) are excluded to the fullest extent permitted by law." ;
> 		:cmor<sub>version</sub> = "3.6.0" ;
> }


<a id="org4792b1c"></a>

### [tas](file:///Volumes/ioa01/reanalysis/JRA55do/tas)

Near-surface air temperature

> netcdf tas<sub>input4MIPs</sub><sub>atmosphericState</sub><sub>OMIP</sub><sub>MRI</sub>-JRA55-do-1-5-0<sub>gr</sub><sub>201001010000</sub>-201012312100 {
> dimensions:
> 	time = UNLIMITED ; // (2920 currently)
> 	lat = 320 ;
> 	lon = 640 ;
> 	bnds = 2 ;
> variables:
> 	double time(time) ;
> 		time:bounds = "time<sub>bnds</sub>" ;
> 		time:units = "days since 1900-01-01 00:00:00" ;
> 		time:calendar = "gregorian" ;
> 		time:axis = "T" ;
> 		time:long<sub>name</sub> = "time" ;
> 		time:standard<sub>name</sub> = "time" ;
> 	double time<sub>bnds</sub>(time, bnds) ;
> 	double lat(lat) ;
> 		lat:bounds = "lat<sub>bnds</sub>" ;
> 		lat:units = "degrees<sub>north</sub>" ;
> 		lat:axis = "Y" ;
> 		lat:long<sub>name</sub> = "Latitude" ;
> 		lat:standard<sub>name</sub> = "latitude" ;
> 	double lat<sub>bnds</sub>(lat, bnds) ;
> 	double lon(lon) ;
> 		lon:bounds = "lon<sub>bnds</sub>" ;
> 		lon:units = "degrees<sub>east</sub>" ;
> 		lon:axis = "X" ;
> 		lon:long<sub>name</sub> = "Longitude" ;
> 		lon:standard<sub>name</sub> = "longitude" ;
> 	double lon<sub>bnds</sub>(lon, bnds) ;
> 	double height ;
> 		height:units = "m" ;
> 		height:axis = "Z" ;
> 		height:positive = "up" ;
> 		height:long<sub>name</sub> = "height" ;
> 		height:standard<sub>name</sub> = "height" ;
> 	float tas(time, lat, lon) ;
> 		tas:standard<sub>name</sub> = "air<sub>temperature</sub>" ;
> 		tas:long<sub>name</sub> = "Near-Surface Air Temperature" ;
> 		tas:comment = "near-surface (usually, 2 meter) air temperature" ;
> 		tas:units = "K" ;
> 		tas:cell<sub>methods</sub> = "area: mean time: point" ;
> 		tas:cell<sub>measures</sub> = "area: areacella" ;
> 		tas:history = "2020-09-15T19:35:02Z altered by CMOR: Treated scalar dimension: \\'height\\'." ;
> 		tas:coordinates = "height" ;
> 		tas:missing<sub>value</sub> = 1.e+20f ;
> 		tas:<sub>FillValue</sub> = 1.e+20f ;
> 
> // global attributes:
> 		:Conventions = "CF-1.7 CMIP-6.2" ;
> 		:activity<sub>id</sub> = "input4MIPs" ;
> 		:cell<sub>measures</sub> = "area: areacella" ;
> 		:comment = "Based on JRA-55 reanalysis (1958-01 to 2020-07)" ;
> 		:contact = "Hiroyuki Tsujino (htsujino@mri-jma.go.jp)" ;
> 		:creation<sub>date</sub> = "2020-09-15T19:35:07Z" ;
> 		:data<sub>specs</sub><sub>version</sub> = "01.00.32" ;
> 		:dataset<sub>category</sub> = "atmosphericState" ;
> 		:external<sub>variables</sub> = "areacella" ;
> 		:frequency = "3hrPt" ;
> 		:further<sub>info</sub><sub>url</sub> = "<http://climate.mri-jma.go.jp/~htsujino/jra55do.html>" ;
> 		:grid = "data regridded to the normal atmosphere TL319 gaussian grid (320x640 latxlon) from a reduced TL319 gaussian grid" ;
> 		:grid<sub>label</sub> = "gr" ;
> 		:history = "2020-09-15T19:35:07Z; CMOR rewrote data to be consistent with input4MIPs, CMIP6, and CF-1.7 CMIP-6.2 standards" ;
> 		:institution = "Meteorological Research Institute, Tsukuba, Ibaraki 305-0052, Japan" ;
> 		:institution<sub>id</sub> = "MRI" ;
> 		:mip<sub>era</sub> = "CMIP6" ;
> 		:nominal<sub>resolution</sub> = "50 km" ;
> 		:product = "reanalysis" ;
> 		:realm = "atmos" ;
> 		:references = "Tsujino et al., 2018: JRA-55 based surface dataset for driving ocean-sea-ice models (JRA55-do), Ocean Modelling, 130(1), pp 79-139. <https://doi.org/10.1016/j.ocemod.2018.07.002>" ;
> 		:region = "global<sub>ocean</sub>" ;
> 		:release<sub>year</sub> = "2020" ;
> 		:source = "MRI JRA55-do 1.5.0: Atmospheric state generated for OMIP based on the JRA-55 reanalysis" ;
> 		:source<sub>description</sub> = "Atmospheric state and terrestrial runoff datasets produced by MRI for the OMIP experiment of CMIP6" ;
> 		:source<sub>id</sub> = "MRI-JRA55-do-1-5-0" ;
> 		:source<sub>type</sub> = "satellite<sub>blended</sub>" ;
> 		:source<sub>version</sub> = "1.5.0" ;
> 		:table<sub>id</sub> = "input4MIPs<sub>A3hrPt</sub>" ;
> 		:table<sub>info</sub> = "Creation Date:(14 September 2020) MD5:d10deb7a45fa93934ec0a16c9d534e45" ;
> 		:target<sub>mip</sub> = "OMIP" ;
> 		:title = "MRI JRA55-do 1.5.0 dataset prepared for input4MIPs" ;
> 		:tracking<sub>id</sub> = "hdl:21.14100/a0582304-8b18-4297-b548-f33f0b40eccb" ;
> 		:variable<sub>id</sub> = "tas" ;
> 		:license = "OMIP boundary condition data produced by MRI is licensed under a Creative Commons Attribution-[NonCommercial-]ShareAlike 4.0 International License (<https://creativecommons.org/licenses>). Consult <https://pcmdi.llnl.gov/CMIP6/TermsOfUse> for terms of use governing input4MIPs output, including citation requirements and proper acknowledgment. Further information about this data, including some limitations, can be found via the further<sub>info</sub><sub>url</sub> (recorded as a global attribute in this file). The data producers and data providers make no warranty, either express or implied, including, but not limited to, warranties of merchantability and fitness for a particular purpose. All liabilities arising from the supply of the information (including any liability arising in negligence) are excluded to the fullest extent permitted by law." ;
> 		:cmor<sub>version</sub> = "3.6.0" ;
> }


<a id="orgd4d04f7"></a>

### [ts](file:///Volumes/ioa01/reanalysis/JRA55do/ts)

Surface temperature

> netcdf ts<sub>input4MIPs</sub><sub>atmosphericState</sub><sub>OMIP</sub><sub>MRI</sub>-JRA55-do-1-5-0<sub>gr</sub><sub>201001010000</sub>-201012312100 {
> dimensions:
> 	time = UNLIMITED ; // (2920 currently)
> 	lat = 320 ;
> 	lon = 640 ;
> 	bnds = 2 ;
> variables:
> 	double time(time) ;
> 		time:bounds = "time<sub>bnds</sub>" ;
> 		time:units = "days since 1900-01-01 00:00:00" ;
> 		time:calendar = "gregorian" ;
> 		time:axis = "T" ;
> 		time:long<sub>name</sub> = "time" ;
> 		time:standard<sub>name</sub> = "time" ;
> 	double time<sub>bnds</sub>(time, bnds) ;
> 	double lat(lat) ;
> 		lat:bounds = "lat<sub>bnds</sub>" ;
> 		lat:units = "degrees<sub>north</sub>" ;
> 		lat:axis = "Y" ;
> 		lat:long<sub>name</sub> = "Latitude" ;
> 		lat:standard<sub>name</sub> = "latitude" ;
> 	double lat<sub>bnds</sub>(lat, bnds) ;
> 	double lon(lon) ;
> 		lon:bounds = "lon<sub>bnds</sub>" ;
> 		lon:units = "degrees<sub>east</sub>" ;
> 		lon:axis = "X" ;
> 		lon:long<sub>name</sub> = "Longitude" ;
> 		lon:standard<sub>name</sub> = "longitude" ;
> 	double lon<sub>bnds</sub>(lon, bnds) ;
> 	float ts(time, lat, lon) ;
> 		ts:standard<sub>name</sub> = "surface<sub>temperature</sub>" ;
> 		ts:long<sub>name</sub> = "Surface Temperature" ;
> 		ts:comment = "Temperature of the lower boundary of the atmosphere" ;
> 		ts:units = "K" ;
> 		ts:cell<sub>methods</sub> = "area: mean time: point" ;
> 		ts:cell<sub>measures</sub> = "area: areacella" ;
> 		ts:missing<sub>value</sub> = 1.e+20f ;
> 		ts:<sub>FillValue</sub> = 1.e+20f ;
> 
> // global attributes:
> 		:Conventions = "CF-1.7 CMIP-6.2" ;
> 		:activity<sub>id</sub> = "input4MIPs" ;
> 		:cell<sub>measures</sub> = "area: areacella" ;
> 		:comment = "Based on JRA-55 reanalysis (1958-01 to 2020-07)" ;
> 		:contact = "Hiroyuki Tsujino (htsujino@mri-jma.go.jp)" ;
> 		:creation<sub>date</sub> = "2020-09-15T21:48:05Z" ;
> 		:data<sub>specs</sub><sub>version</sub> = "01.00.32" ;
> 		:dataset<sub>category</sub> = "atmosphericState" ;
> 		:external<sub>variables</sub> = "areacella" ;
> 		:frequency = "3hrPt" ;
> 		:further<sub>info</sub><sub>url</sub> = "<http://climate.mri-jma.go.jp/~htsujino/jra55do.html>" ;
> 		:grid = "data regridded to the normal atmosphere TL319 gaussian grid (320x640 latxlon) from a reduced TL319 gaussian grid" ;
> 		:grid<sub>label</sub> = "gr" ;
> 		:history = "2020-09-15T21:48:05Z; CMOR rewrote data to be consistent with input4MIPs, CMIP6, and CF-1.7 CMIP-6.2 standards" ;
> 		:institution = "Meteorological Research Institute, Tsukuba, Ibaraki 305-0052, Japan" ;
> 		:institution<sub>id</sub> = "MRI" ;
> 		:mip<sub>era</sub> = "CMIP6" ;
> 		:nominal<sub>resolution</sub> = "50 km" ;
> 		:product = "reanalysis" ;
> 		:realm = "atmos" ;
> 		:references = "Tsujino et al., 2018: JRA-55 based surface dataset for driving ocean-sea-ice models (JRA55-do), Ocean Modelling, 130(1), pp 79-139. <https://doi.org/10.1016/j.ocemod.2018.07.002>" ;
> 		:region = "global<sub>ocean</sub>" ;
> 		:release<sub>year</sub> = "2020" ;
> 		:source = "MRI JRA55-do 1.5.0: Atmospheric state generated for OMIP based on the JRA-55 reanalysis" ;
> 		:source<sub>description</sub> = "Atmospheric state and terrestrial runoff datasets produced by MRI for the OMIP experiment of CMIP6" ;
> 		:source<sub>id</sub> = "MRI-JRA55-do-1-5-0" ;
> 		:source<sub>type</sub> = "satellite<sub>blended</sub>" ;
> 		:source<sub>version</sub> = "1.5.0" ;
> 		:table<sub>id</sub> = "input4MIPs<sub>A3hrPt</sub>" ;
> 		:table<sub>info</sub> = "Creation Date:(14 September 2020) MD5:d10deb7a45fa93934ec0a16c9d534e45" ;
> 		:target<sub>mip</sub> = "OMIP" ;
> 		:title = "MRI JRA55-do 1.5.0 dataset prepared for input4MIPs" ;
> 		:tracking<sub>id</sub> = "hdl:21.14100/d725c7b9-ce5a-4fe7-8f65-3003acf15323" ;
> 		:variable<sub>id</sub> = "ts" ;
> 		:license = "OMIP boundary condition data produced by MRI is licensed under a Creative Commons Attribution-[NonCommercial-]ShareAlike 4.0 International License (<https://creativecommons.org/licenses>). Consult <https://pcmdi.llnl.gov/CMIP6/TermsOfUse> for terms of use governing input4MIPs output, including citation requirements and proper acknowledgment. Further information about this data, including some limitations, can be found via the further<sub>info</sub><sub>url</sub> (recorded as a global attribute in this file). The data producers and data providers make no warranty, either express or implied, including, but not limited to, warranties of merchantability and fitness for a particular purpose. All liabilities arising from the supply of the information (including any liability arising in negligence) are excluded to the fullest extent permitted by law." ;
> 		:cmor<sub>version</sub> = "3.6.0" ;
> }


<a id="org8f639e9"></a>

### [uas](file:///Volumes/ioa01/reanalysis/JRA55do/uas)

Eastward near-surface wind

> netcdf uas<sub>input4MIPs</sub><sub>atmosphericState</sub><sub>OMIP</sub><sub>MRI</sub>-JRA55-do-1-5-0<sub>gr</sub><sub>201001010000</sub>-201012312100 {
> dimensions:
> 	time = UNLIMITED ; // (2920 currently)
> 	lat = 320 ;
> 	lon = 640 ;
> 	bnds = 2 ;
> variables:
> 	double time(time) ;
> 		time:bounds = "time<sub>bnds</sub>" ;
> 		time:units = "days since 1900-01-01 00:00:00" ;
> 		time:calendar = "gregorian" ;
> 		time:axis = "T" ;
> 		time:long<sub>name</sub> = "time" ;
> 		time:standard<sub>name</sub> = "time" ;
> 	double time<sub>bnds</sub>(time, bnds) ;
> 	double lat(lat) ;
> 		lat:bounds = "lat<sub>bnds</sub>" ;
> 		lat:units = "degrees<sub>north</sub>" ;
> 		lat:axis = "Y" ;
> 		lat:long<sub>name</sub> = "Latitude" ;
> 		lat:standard<sub>name</sub> = "latitude" ;
> 	double lat<sub>bnds</sub>(lat, bnds) ;
> 	double lon(lon) ;
> 		lon:bounds = "lon<sub>bnds</sub>" ;
> 		lon:units = "degrees<sub>east</sub>" ;
> 		lon:axis = "X" ;
> 		lon:long<sub>name</sub> = "Longitude" ;
> 		lon:standard<sub>name</sub> = "longitude" ;
> 	double lon<sub>bnds</sub>(lon, bnds) ;
> 	double height ;
> 		height:units = "m" ;
> 		height:axis = "Z" ;
> 		height:positive = "up" ;
> 		height:long<sub>name</sub> = "height" ;
> 		height:standard<sub>name</sub> = "height" ;
> 	float uas(time, lat, lon) ;
> 		uas:standard<sub>name</sub> = "eastward<sub>wind</sub>" ;
> 		uas:long<sub>name</sub> = "Eastward Near-Surface Wind" ;
> 		uas:comment = "Eastward component of the near-surface wind" ;
> 		uas:units = "m s-1" ;
> 		uas:cell<sub>methods</sub> = "area: mean time: point" ;
> 		uas:cell<sub>measures</sub> = "area: areacella" ;
> 		uas:history = "2020-09-16T01:30:40Z altered by CMOR: Treated scalar dimension: \\'height\\'." ;
> 		uas:coordinates = "height" ;
> 		uas:missing<sub>value</sub> = 1.e+20f ;
> 		uas:<sub>FillValue</sub> = 1.e+20f ;
> 
> // global attributes:
> 		:Conventions = "CF-1.7 CMIP-6.2" ;
> 		:activity<sub>id</sub> = "input4MIPs" ;
> 		:cell<sub>measures</sub> = "area: areacella" ;
> 		:comment = "Based on JRA-55 reanalysis (1958-01 to 2020-07)" ;
> 		:contact = "Hiroyuki Tsujino (htsujino@mri-jma.go.jp)" ;
> 		:creation<sub>date</sub> = "2020-09-16T01:30:45Z" ;
> 		:data<sub>specs</sub><sub>version</sub> = "01.00.32" ;
> 		:dataset<sub>category</sub> = "atmosphericState" ;
> 		:external<sub>variables</sub> = "areacella" ;
> 		:frequency = "3hrPt" ;
> 		:further<sub>info</sub><sub>url</sub> = "<http://climate.mri-jma.go.jp/~htsujino/jra55do.html>" ;
> 		:grid = "data regridded to the normal atmosphere TL319 gaussian grid (320x640 latxlon) from a reduced TL319 gaussian grid" ;
> 		:grid<sub>label</sub> = "gr" ;
> 		:history = "2020-09-16T01:30:45Z; CMOR rewrote data to be consistent with input4MIPs, CMIP6, and CF-1.7 CMIP-6.2 standards" ;
> 		:institution = "Meteorological Research Institute, Tsukuba, Ibaraki 305-0052, Japan" ;
> 		:institution<sub>id</sub> = "MRI" ;
> 		:mip<sub>era</sub> = "CMIP6" ;
> 		:nominal<sub>resolution</sub> = "50 km" ;
> 		:product = "reanalysis" ;
> 		:realm = "atmos" ;
> 		:references = "Tsujino et al., 2018: JRA-55 based surface dataset for driving ocean-sea-ice models (JRA55-do), Ocean Modelling, 130(1), pp 79-139. <https://doi.org/10.1016/j.ocemod.2018.07.002>" ;
> 		:region = "global<sub>ocean</sub>" ;
> 		:release<sub>year</sub> = "2020" ;
> 		:source = "MRI JRA55-do 1.5.0: Atmospheric state generated for OMIP based on the JRA-55 reanalysis" ;
> 		:source<sub>description</sub> = "Atmospheric state and terrestrial runoff datasets produced by MRI for the OMIP experiment of CMIP6" ;
> 		:source<sub>id</sub> = "MRI-JRA55-do-1-5-0" ;
> 		:source<sub>type</sub> = "satellite<sub>blended</sub>" ;
> 		:source<sub>version</sub> = "1.5.0" ;
> 		:table<sub>id</sub> = "input4MIPs<sub>A3hrPt</sub>" ;
> 		:table<sub>info</sub> = "Creation Date:(14 September 2020) MD5:d10deb7a45fa93934ec0a16c9d534e45" ;
> 		:target<sub>mip</sub> = "OMIP" ;
> 		:title = "MRI JRA55-do 1.5.0 dataset prepared for input4MIPs" ;
> 		:tracking<sub>id</sub> = "hdl:21.14100/9e2256e3-29b7-4e8a-906a-47a98a256308" ;
> 		:variable<sub>id</sub> = "uas" ;
> 		:license = "OMIP boundary condition data produced by MRI is licensed under a Creative Commons Attribution-[NonCommercial-]ShareAlike 4.0 International License (<https://creativecommons.org/licenses>). Consult <https://pcmdi.llnl.gov/CMIP6/TermsOfUse> for terms of use governing input4MIPs output, including citation requirements and proper acknowledgment. Further information about this data, including some limitations, can be found via the further<sub>info</sub><sub>url</sub> (recorded as a global attribute in this file). The data producers and data providers make no warranty, either express or implied, including, but not limited to, warranties of merchantability and fitness for a particular purpose. All liabilities arising from the supply of the information (including any liability arising in negligence) are excluded to the fullest extent permitted by law." ;
> 		:cmor<sub>version</sub> = "3.6.0" ;
> }


<a id="orgbdb0273"></a>

### [vas](file:///Volumes/ioa01/reanalysis/JRA55do/vas)

Northward near-surface wind

> netcdf vas<sub>input4MIPs</sub><sub>atmosphericState</sub><sub>OMIP</sub><sub>MRI</sub>-JRA55-do-1-5-0<sub>gr</sub><sub>201001010000</sub>-201012312100 {
> dimensions:
> 	time = UNLIMITED ; // (2920 currently)
> 	lat = 320 ;
> 	lon = 640 ;
> 	bnds = 2 ;
> variables:
> 	double time(time) ;
> 		time:bounds = "time<sub>bnds</sub>" ;
> 		time:units = "days since 1900-01-01 00:00:00" ;
> 		time:calendar = "gregorian" ;
> 		time:axis = "T" ;
> 		time:long<sub>name</sub> = "time" ;
> 		time:standard<sub>name</sub> = "time" ;
> 	double time<sub>bnds</sub>(time, bnds) ;
> 	double lat(lat) ;
> 		lat:bounds = "lat<sub>bnds</sub>" ;
> 		lat:units = "degrees<sub>north</sub>" ;
> 		lat:axis = "Y" ;
> 		lat:long<sub>name</sub> = "Latitude" ;
> 		lat:standard<sub>name</sub> = "latitude" ;
> 	double lat<sub>bnds</sub>(lat, bnds) ;
> 	double lon(lon) ;
> 		lon:bounds = "lon<sub>bnds</sub>" ;
> 		lon:units = "degrees<sub>east</sub>" ;
> 		lon:axis = "X" ;
> 		lon:long<sub>name</sub> = "Longitude" ;
> 		lon:standard<sub>name</sub> = "longitude" ;
> 	double lon<sub>bnds</sub>(lon, bnds) ;
> 	double height ;
> 		height:units = "m" ;
> 		height:axis = "Z" ;
> 		height:positive = "up" ;
> 		height:long<sub>name</sub> = "height" ;
> 		height:standard<sub>name</sub> = "height" ;
> 	float vas(time, lat, lon) ;
> 		vas:standard<sub>name</sub> = "northward<sub>wind</sub>" ;
> 		vas:long<sub>name</sub> = "Northward Near-Surface Wind" ;
> 		vas:comment = "Northward component of the near surface wind" ;
> 		vas:units = "m s-1" ;
> 		vas:cell<sub>methods</sub> = "area: mean time: point" ;
> 		vas:cell<sub>measures</sub> = "area: areacella" ;
> 		vas:history = "2020-09-16T04:10:06Z altered by CMOR: Treated scalar dimension: \\'height\\'." ;
> 		vas:coordinates = "height" ;
> 		vas:missing<sub>value</sub> = 1.e+20f ;
> 		vas:<sub>FillValue</sub> = 1.e+20f ;
> 
> // global attributes:
> 		:Conventions = "CF-1.7 CMIP-6.2" ;
> 		:activity<sub>id</sub> = "input4MIPs" ;
> 		:cell<sub>measures</sub> = "area: areacella" ;
> 		:comment = "Based on JRA-55 reanalysis (1958-01 to 2020-07)" ;
> 		:contact = "Hiroyuki Tsujino (htsujino@mri-jma.go.jp)" ;
> 		:creation<sub>date</sub> = "2020-09-16T04:10:11Z" ;
> 		:data<sub>specs</sub><sub>version</sub> = "01.00.32" ;
> 		:dataset<sub>category</sub> = "atmosphericState" ;
> 		:external<sub>variables</sub> = "areacella" ;
> 		:frequency = "3hrPt" ;
> 		:further<sub>info</sub><sub>url</sub> = "<http://climate.mri-jma.go.jp/~htsujino/jra55do.html>" ;
> 		:grid = "data regridded to the normal atmosphere TL319 gaussian grid (320x640 latxlon) from a reduced TL319 gaussian grid" ;
> 		:grid<sub>label</sub> = "gr" ;
> 		:history = "2020-09-16T04:10:11Z; CMOR rewrote data to be consistent with input4MIPs, CMIP6, and CF-1.7 CMIP-6.2 standards" ;
> 		:institution = "Meteorological Research Institute, Tsukuba, Ibaraki 305-0052, Japan" ;
> 		:institution<sub>id</sub> = "MRI" ;
> 		:mip<sub>era</sub> = "CMIP6" ;
> 		:nominal<sub>resolution</sub> = "50 km" ;
> 		:product = "reanalysis" ;
> 		:realm = "atmos" ;
> 		:references = "Tsujino et al., 2018: JRA-55 based surface dataset for driving ocean-sea-ice models (JRA55-do), Ocean Modelling, 130(1), pp 79-139. <https://doi.org/10.1016/j.ocemod.2018.07.002>" ;
> 		:region = "global<sub>ocean</sub>" ;
> 		:release<sub>year</sub> = "2020" ;
> 		:source = "MRI JRA55-do 1.5.0: Atmospheric state generated for OMIP based on the JRA-55 reanalysis" ;
> 		:source<sub>description</sub> = "Atmospheric state and terrestrial runoff datasets produced by MRI for the OMIP experiment of CMIP6" ;
> 		:source<sub>id</sub> = "MRI-JRA55-do-1-5-0" ;
> 		:source<sub>type</sub> = "satellite<sub>blended</sub>" ;
> 		:source<sub>version</sub> = "1.5.0" ;
> 		:table<sub>id</sub> = "input4MIPs<sub>A3hrPt</sub>" ;
> 		:table<sub>info</sub> = "Creation Date:(14 September 2020) MD5:d10deb7a45fa93934ec0a16c9d534e45" ;
> 		:target<sub>mip</sub> = "OMIP" ;
> 		:title = "MRI JRA55-do 1.5.0 dataset prepared for input4MIPs" ;
> 		:tracking<sub>id</sub> = "hdl:21.14100/b80b6346-812b-4171-9be7-75718695dd0e" ;
> 		:variable<sub>id</sub> = "vas" ;
> 		:license = "OMIP boundary condition data produced by MRI is licensed under a Creative Commons Attribution-[NonCommercial-]ShareAlike 4.0 International License (<https://creativecommons.org/licenses>). Consult <https://pcmdi.llnl.gov/CMIP6/TermsOfUse> for terms of use governing input4MIPs output, including citation requirements and proper acknowledgment. Further information about this data, including some limitations, can be found via the further<sub>info</sub><sub>url</sub> (recorded as a global attribute in this file). The data producers and data providers make no warranty, either express or implied, including, but not limited to, warranties of merchantability and fitness for a particular purpose. All liabilities arising from the supply of the information (including any liability arising in negligence) are excluded to the fullest extent permitted by law." ;
> 		:cmor<sub>version</sub> = "3.6.0" ;
> }


<a id="org6bedea4"></a>

## MERRAv2

[MERRAv2](https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/) is a $1/2$ degree spatial and 1 hour temporal resolution atmospheric
climate dataset maintained by NASA. gadi has a repository described [here](http://climate-cms.wikis.unsw.edu.au/MERRA2) and an
[excel spreadsheet](file:///Users/dpath2o/PHD/references/data/MERRA2_variables.xlsx) has also been provided to further describe this repository.


<a id="org0abc25e"></a>

### Surface Flux Diagnostics ('flx')

I have downloaded years 2010-2019 files from gadi to my [local repository](file:///Volumes/ioa01/reanalysis/MERRAv2/flx). Here
is the header of the first file:

> netcdf MERRA2<sub>300.tavg1</sub><sub>2d</sub><sub>flx</sub><sub>Nx.20100101</sub> {
> dimensions:
> 	lon = 576 ;
> 	lat = 361 ;
> 	time = UNLIMITED ; // (24 currently)
> variables:
> 	double lon(lon) ;
> 		lon:long<sub>name</sub> = "longitude" ;
> 		lon:units = "degrees<sub>east</sub>" ;
> 		lon:vmax = 1.e+15f ;
> 		lon:vmin = -1.e+15f ;
> 		lon:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	double lat(lat) ;
> 		lat:long<sub>name</sub> = "latitude" ;
> 		lat:units = "degrees<sub>north</sub>" ;
> 		lat:vmax = 1.e+15f ;
> 		lat:vmin = -1.e+15f ;
> 		lat:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	int time(time) ;
> 		time:long<sub>name</sub> = "time" ;
> 		time:units = "minutes since 2010-01-01 00:30:00" ;
> 		time:time<sub>increment</sub> = 10000 ;
> 		time:begin<sub>date</sub> = 20100101 ;
> 		time:begin<sub>time</sub> = 3000 ;
> 		time:vmax = 1.e+15f ;
> 		time:vmin = -1.e+15f ;
> 		time:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float BSTAR(time, lat, lon) ;
> 		BSTAR:long<sub>name</sub> = "surface<sub>bouyancy</sub><sub>scale</sub>" ;
> 		BSTAR:units = "m s-2" ;
> 		BSTAR:<sub>FillValue</sub> = 1.e+15f ;
> 		BSTAR:missing<sub>value</sub> = 1.e+15f ;
> 		BSTAR:fmissing<sub>value</sub> = 1.e+15f ;
> 		BSTAR:scale<sub>factor</sub> = 1.f ;
> 		BSTAR:add<sub>offset</sub> = 0.f ;
> 		BSTAR:standard<sub>name</sub> = "surface<sub>bouyancy</sub><sub>scale</sub>" ;
> 		BSTAR:vmax = 1.e+15f ;
> 		BSTAR:vmin = -1.e+15f ;
> 		BSTAR:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float CDH(time, lat, lon) ;
> 		CDH:long<sub>name</sub> = "surface<sub>exchange</sub><sub>coefficient</sub><sub>for</sub><sub>heat</sub>" ;
> 		CDH:units = "kg m-2 s-1" ;
> 		CDH:<sub>FillValue</sub> = 1.e+15f ;
> 		CDH:missing<sub>value</sub> = 1.e+15f ;
> 		CDH:fmissing<sub>value</sub> = 1.e+15f ;
> 		CDH:scale<sub>factor</sub> = 1.f ;
> 		CDH:add<sub>offset</sub> = 0.f ;
> 		CDH:standard<sub>name</sub> = "surface<sub>exchange</sub><sub>coefficient</sub><sub>for</sub><sub>heat</sub>" ;
> 		CDH:vmax = 1.e+15f ;
> 		CDH:vmin = -1.e+15f ;
> 		CDH:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float CDM(time, lat, lon) ;
> 		CDM:long<sub>name</sub> = "surface<sub>exchange</sub><sub>coefficient</sub><sub>for</sub><sub>momentum</sub>" ;
> 		CDM:units = "kg m-2 s-1" ;
> 		CDM:<sub>FillValue</sub> = 1.e+15f ;
> 		CDM:missing<sub>value</sub> = 1.e+15f ;
> 		CDM:fmissing<sub>value</sub> = 1.e+15f ;
> 		CDM:scale<sub>factor</sub> = 1.f ;
> 		CDM:add<sub>offset</sub> = 0.f ;
> 		CDM:standard<sub>name</sub> = "surface<sub>exchange</sub><sub>coefficient</sub><sub>for</sub><sub>momentum</sub>" ;
> 		CDM:vmax = 1.e+15f ;
> 		CDM:vmin = -1.e+15f ;
> 		CDM:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float CDQ(time, lat, lon) ;
> 		CDQ:long<sub>name</sub> = "surface<sub>exchange</sub><sub>coefficient</sub><sub>for</sub><sub>moisture</sub>" ;
> 		CDQ:units = "kg m-2 s-1" ;
> 		CDQ:<sub>FillValue</sub> = 1.e+15f ;
> 		CDQ:missing<sub>value</sub> = 1.e+15f ;
> 		CDQ:fmissing<sub>value</sub> = 1.e+15f ;
> 		CDQ:scale<sub>factor</sub> = 1.f ;
> 		CDQ:add<sub>offset</sub> = 0.f ;
> 		CDQ:standard<sub>name</sub> = "surface<sub>exchange</sub><sub>coefficient</sub><sub>for</sub><sub>moisture</sub>" ;
> 		CDQ:vmax = 1.e+15f ;
> 		CDQ:vmin = -1.e+15f ;
> 		CDQ:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float CN(time, lat, lon) ;
> 		CN:long<sub>name</sub> = "surface<sub>neutral</sub><sub>drag</sub><sub>coefficient</sub>" ;
> 		CN:units = "1" ;
> 		CN:<sub>FillValue</sub> = 1.e+15f ;
> 		CN:missing<sub>value</sub> = 1.e+15f ;
> 		CN:fmissing<sub>value</sub> = 1.e+15f ;
> 		CN:scale<sub>factor</sub> = 1.f ;
> 		CN:add<sub>offset</sub> = 0.f ;
> 		CN:standard<sub>name</sub> = "surface<sub>neutral</sub><sub>drag</sub><sub>coefficient</sub>" ;
> 		CN:vmax = 1.e+15f ;
> 		CN:vmin = -1.e+15f ;
> 		CN:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float DISPH(time, lat, lon) ;
> 		DISPH:long<sub>name</sub> = "zero<sub>plane</sub><sub>displacement</sub><sub>height</sub>" ;
> 		DISPH:units = "m" ;
> 		DISPH:<sub>FillValue</sub> = 1.e+15f ;
> 		DISPH:missing<sub>value</sub> = 1.e+15f ;
> 		DISPH:fmissing<sub>value</sub> = 1.e+15f ;
> 		DISPH:scale<sub>factor</sub> = 1.f ;
> 		DISPH:add<sub>offset</sub> = 0.f ;
> 		DISPH:standard<sub>name</sub> = "zero<sub>plane</sub><sub>displacement</sub><sub>height</sub>" ;
> 		DISPH:vmax = 1.e+15f ;
> 		DISPH:vmin = -1.e+15f ;
> 		DISPH:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float EFLUX(time, lat, lon) ;
> 		EFLUX:long<sub>name</sub> = "total<sub>latent</sub><sub>energy</sub><sub>flux</sub>" ;
> 		EFLUX:units = "W m-2" ;
> 		EFLUX:<sub>FillValue</sub> = 1.e+15f ;
> 		EFLUX:missing<sub>value</sub> = 1.e+15f ;
> 		EFLUX:fmissing<sub>value</sub> = 1.e+15f ;
> 		EFLUX:scale<sub>factor</sub> = 1.f ;
> 		EFLUX:add<sub>offset</sub> = 0.f ;
> 		EFLUX:standard<sub>name</sub> = "total<sub>latent</sub><sub>energy</sub><sub>flux</sub>" ;
> 		EFLUX:vmax = 1.e+15f ;
> 		EFLUX:vmin = -1.e+15f ;
> 		EFLUX:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float EVAP(time, lat, lon) ;
> 		EVAP:long<sub>name</sub> = "evaporation<sub>from</sub><sub>turbulence</sub>" ;
> 		EVAP:units = "kg m-2 s-1" ;
> 		EVAP:<sub>FillValue</sub> = 1.e+15f ;
> 		EVAP:missing<sub>value</sub> = 1.e+15f ;
> 		EVAP:fmissing<sub>value</sub> = 1.e+15f ;
> 		EVAP:scale<sub>factor</sub> = 1.f ;
> 		EVAP:add<sub>offset</sub> = 0.f ;
> 		EVAP:standard<sub>name</sub> = "evaporation<sub>from</sub><sub>turbulence</sub>" ;
> 		EVAP:vmax = 1.e+15f ;
> 		EVAP:vmin = -1.e+15f ;
> 		EVAP:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float FRCAN(time, lat, lon) ;
> 		FRCAN:long<sub>name</sub> = "areal<sub>fraction</sub><sub>of</sub><sub>anvil</sub><sub>showers</sub>" ;
> 		FRCAN:units = "1" ;
> 		FRCAN:<sub>FillValue</sub> = 1.e+15f ;
> 		FRCAN:missing<sub>value</sub> = 1.e+15f ;
> 		FRCAN:fmissing<sub>value</sub> = 1.e+15f ;
> 		FRCAN:scale<sub>factor</sub> = 1.f ;
> 		FRCAN:add<sub>offset</sub> = 0.f ;
> 		FRCAN:standard<sub>name</sub> = "areal<sub>fraction</sub><sub>of</sub><sub>anvil</sub><sub>showers</sub>" ;
> 		FRCAN:vmax = 1.e+15f ;
> 		FRCAN:vmin = -1.e+15f ;
> 		FRCAN:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float FRCCN(time, lat, lon) ;
> 		FRCCN:long<sub>name</sub> = "areal<sub>fraction</sub><sub>of</sub><sub>convective</sub><sub>showers</sub>" ;
> 		FRCCN:units = "1" ;
> 		FRCCN:<sub>FillValue</sub> = 1.e+15f ;
> 		FRCCN:missing<sub>value</sub> = 1.e+15f ;
> 		FRCCN:fmissing<sub>value</sub> = 1.e+15f ;
> 		FRCCN:scale<sub>factor</sub> = 1.f ;
> 		FRCCN:add<sub>offset</sub> = 0.f ;
> 		FRCCN:standard<sub>name</sub> = "areal<sub>fraction</sub><sub>of</sub><sub>convective</sub><sub>showers</sub>" ;
> 		FRCCN:vmax = 1.e+15f ;
> 		FRCCN:vmin = -1.e+15f ;
> 		FRCCN:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float FRCLS(time, lat, lon) ;
> 		FRCLS:long<sub>name</sub> = "areal<sub>fraction</sub><sub>of</sub><sub>nonanvil</sub><sub>large</sub><sub>scale</sub><sub>showers</sub>" ;
> 		FRCLS:units = "1" ;
> 		FRCLS:<sub>FillValue</sub> = 1.e+15f ;
> 		FRCLS:missing<sub>value</sub> = 1.e+15f ;
> 		FRCLS:fmissing<sub>value</sub> = 1.e+15f ;
> 		FRCLS:scale<sub>factor</sub> = 1.f ;
> 		FRCLS:add<sub>offset</sub> = 0.f ;
> 		FRCLS:standard<sub>name</sub> = "areal<sub>fraction</sub><sub>of</sub><sub>nonanvil</sub><sub>large</sub><sub>scale</sub><sub>showers</sub>" ;
> 		FRCLS:vmax = 1.e+15f ;
> 		FRCLS:vmin = -1.e+15f ;
> 		FRCLS:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float FRSEAICE(time, lat, lon) ;
> 		FRSEAICE:long<sub>name</sub> = "ice<sub>covered</sub><sub>fraction</sub><sub>of</sub><sub>tile</sub>" ;
> 		FRSEAICE:units = "1" ;
> 		FRSEAICE:<sub>FillValue</sub> = 1.e+15f ;
> 		FRSEAICE:missing<sub>value</sub> = 1.e+15f ;
> 		FRSEAICE:fmissing<sub>value</sub> = 1.e+15f ;
> 		FRSEAICE:scale<sub>factor</sub> = 1.f ;
> 		FRSEAICE:add<sub>offset</sub> = 0.f ;
> 		FRSEAICE:standard<sub>name</sub> = "ice<sub>covered</sub><sub>fraction</sub><sub>of</sub><sub>tile</sub>" ;
> 		FRSEAICE:vmax = 1.e+15f ;
> 		FRSEAICE:vmin = -1.e+15f ;
> 		FRSEAICE:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float GHTSKIN(time, lat, lon) ;
> 		GHTSKIN:long<sub>name</sub> = "Ground<sub>heating</sub><sub>for</sub><sub>skin</sub><sub>temp</sub>" ;
> 		GHTSKIN:units = "W m-2" ;
> 		GHTSKIN:<sub>FillValue</sub> = 1.e+15f ;
> 		GHTSKIN:missing<sub>value</sub> = 1.e+15f ;
> 		GHTSKIN:fmissing<sub>value</sub> = 1.e+15f ;
> 		GHTSKIN:scale<sub>factor</sub> = 1.f ;
> 		GHTSKIN:add<sub>offset</sub> = 0.f ;
> 		GHTSKIN:standard<sub>name</sub> = "Ground<sub>heating</sub><sub>for</sub><sub>skin</sub><sub>temp</sub>" ;
> 		GHTSKIN:vmax = 1.e+15f ;
> 		GHTSKIN:vmin = -1.e+15f ;
> 		GHTSKIN:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float HFLUX(time, lat, lon) ;
> 		HFLUX:long<sub>name</sub> = "sensible<sub>heat</sub><sub>flux</sub><sub>from</sub><sub>turbulence</sub>" ;
> 		HFLUX:units = "W m-2" ;
> 		HFLUX:<sub>FillValue</sub> = 1.e+15f ;
> 		HFLUX:missing<sub>value</sub> = 1.e+15f ;
> 		HFLUX:fmissing<sub>value</sub> = 1.e+15f ;
> 		HFLUX:scale<sub>factor</sub> = 1.f ;
> 		HFLUX:add<sub>offset</sub> = 0.f ;
> 		HFLUX:standard<sub>name</sub> = "sensible<sub>heat</sub><sub>flux</sub><sub>from</sub><sub>turbulence</sub>" ;
> 		HFLUX:vmax = 1.e+15f ;
> 		HFLUX:vmin = -1.e+15f ;
> 		HFLUX:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float HLML(time, lat, lon) ;
> 		HLML:long<sub>name</sub> = "surface<sub>layer</sub><sub>height</sub>" ;
> 		HLML:units = "m" ;
> 		HLML:<sub>FillValue</sub> = 1.e+15f ;
> 		HLML:missing<sub>value</sub> = 1.e+15f ;
> 		HLML:fmissing<sub>value</sub> = 1.e+15f ;
> 		HLML:scale<sub>factor</sub> = 1.f ;
> 		HLML:add<sub>offset</sub> = 0.f ;
> 		HLML:standard<sub>name</sub> = "surface<sub>layer</sub><sub>height</sub>" ;
> 		HLML:vmax = 1.e+15f ;
> 		HLML:vmin = -1.e+15f ;
> 		HLML:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float NIRDF(time, lat, lon) ;
> 		NIRDF:long<sub>name</sub> = "surface<sub>downwelling</sub><sub>nearinfrared</sub><sub>diffuse</sub><sub>flux</sub>" ;
> 		NIRDF:units = "W m-2" ;
> 		NIRDF:<sub>FillValue</sub> = 1.e+15f ;
> 		NIRDF:missing<sub>value</sub> = 1.e+15f ;
> 		NIRDF:fmissing<sub>value</sub> = 1.e+15f ;
> 		NIRDF:scale<sub>factor</sub> = 1.f ;
> 		NIRDF:add<sub>offset</sub> = 0.f ;
> 		NIRDF:standard<sub>name</sub> = "surface<sub>downwelling</sub><sub>nearinfrared</sub><sub>diffuse</sub><sub>flux</sub>" ;
> 		NIRDF:vmax = 1.e+15f ;
> 		NIRDF:vmin = -1.e+15f ;
> 		NIRDF:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float NIRDR(time, lat, lon) ;
> 		NIRDR:long<sub>name</sub> = "surface<sub>downwelling</sub><sub>nearinfrared</sub><sub>beam</sub><sub>flux</sub>" ;
> 		NIRDR:units = "W m-2" ;
> 		NIRDR:<sub>FillValue</sub> = 1.e+15f ;
> 		NIRDR:missing<sub>value</sub> = 1.e+15f ;
> 		NIRDR:fmissing<sub>value</sub> = 1.e+15f ;
> 		NIRDR:scale<sub>factor</sub> = 1.f ;
> 		NIRDR:add<sub>offset</sub> = 0.f ;
> 		NIRDR:standard<sub>name</sub> = "surface<sub>downwelling</sub><sub>nearinfrared</sub><sub>beam</sub><sub>flux</sub>" ;
> 		NIRDR:vmax = 1.e+15f ;
> 		NIRDR:vmin = -1.e+15f ;
> 		NIRDR:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float PBLH(time, lat, lon) ;
> 		PBLH:long<sub>name</sub> = "planetary<sub>boundary</sub><sub>layer</sub><sub>height</sub>" ;
> 		PBLH:units = "m" ;
> 		PBLH:<sub>FillValue</sub> = 1.e+15f ;
> 		PBLH:missing<sub>value</sub> = 1.e+15f ;
> 		PBLH:fmissing<sub>value</sub> = 1.e+15f ;
> 		PBLH:scale<sub>factor</sub> = 1.f ;
> 		PBLH:add<sub>offset</sub> = 0.f ;
> 		PBLH:standard<sub>name</sub> = "planetary<sub>boundary</sub><sub>layer</sub><sub>height</sub>" ;
> 		PBLH:vmax = 1.e+15f ;
> 		PBLH:vmin = -1.e+15f ;
> 		PBLH:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float PGENTOT(time, lat, lon) ;
> 		PGENTOT:long<sub>name</sub> = "Total<sub>column</sub><sub>production</sub><sub>of</sub><sub>precipitation</sub>" ;
> 		PGENTOT:units = "kg m-2 s-1" ;
> 		PGENTOT:<sub>FillValue</sub> = 1.e+15f ;
> 		PGENTOT:missing<sub>value</sub> = 1.e+15f ;
> 		PGENTOT:fmissing<sub>value</sub> = 1.e+15f ;
> 		PGENTOT:scale<sub>factor</sub> = 1.f ;
> 		PGENTOT:add<sub>offset</sub> = 0.f ;
> 		PGENTOT:standard<sub>name</sub> = "Total<sub>column</sub><sub>production</sub><sub>of</sub><sub>precipitation</sub>" ;
> 		PGENTOT:vmax = 1.e+15f ;
> 		PGENTOT:vmin = -1.e+15f ;
> 		PGENTOT:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float PRECANV(time, lat, lon) ;
> 		PRECANV:long<sub>name</sub> = "anvil<sub>precipitation</sub>" ;
> 		PRECANV:units = "kg m-2 s-1" ;
> 		PRECANV:<sub>FillValue</sub> = 1.e+15f ;
> 		PRECANV:missing<sub>value</sub> = 1.e+15f ;
> 		PRECANV:fmissing<sub>value</sub> = 1.e+15f ;
> 		PRECANV:scale<sub>factor</sub> = 1.f ;
> 		PRECANV:add<sub>offset</sub> = 0.f ;
> 		PRECANV:standard<sub>name</sub> = "anvil<sub>precipitation</sub>" ;
> 		PRECANV:vmax = 1.e+15f ;
> 		PRECANV:vmin = -1.e+15f ;
> 		PRECANV:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float PRECCON(time, lat, lon) ;
> 		PRECCON:long<sub>name</sub> = "convective<sub>precipitation</sub>" ;
> 		PRECCON:units = "kg m-2 s-1" ;
> 		PRECCON:<sub>FillValue</sub> = 1.e+15f ;
> 		PRECCON:missing<sub>value</sub> = 1.e+15f ;
> 		PRECCON:fmissing<sub>value</sub> = 1.e+15f ;
> 		PRECCON:scale<sub>factor</sub> = 1.f ;
> 		PRECCON:add<sub>offset</sub> = 0.f ;
> 		PRECCON:standard<sub>name</sub> = "convective<sub>precipitation</sub>" ;
> 		PRECCON:vmax = 1.e+15f ;
> 		PRECCON:vmin = -1.e+15f ;
> 		PRECCON:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float PRECLSC(time, lat, lon) ;
> 		PRECLSC:long<sub>name</sub> = "nonanvil<sub>large</sub><sub>scale</sub><sub>precipitation</sub>" ;
> 		PRECLSC:units = "kg m-2 s-1" ;
> 		PRECLSC:<sub>FillValue</sub> = 1.e+15f ;
> 		PRECLSC:missing<sub>value</sub> = 1.e+15f ;
> 		PRECLSC:fmissing<sub>value</sub> = 1.e+15f ;
> 		PRECLSC:scale<sub>factor</sub> = 1.f ;
> 		PRECLSC:add<sub>offset</sub> = 0.f ;
> 		PRECLSC:standard<sub>name</sub> = "nonanvil<sub>large</sub><sub>scale</sub><sub>precipitation</sub>" ;
> 		PRECLSC:vmax = 1.e+15f ;
> 		PRECLSC:vmin = -1.e+15f ;
> 		PRECLSC:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float PRECSNO(time, lat, lon) ;
> 		PRECSNO:long<sub>name</sub> = "snowfall" ;
> 		PRECSNO:units = "kg m-2 s-1" ;
> 		PRECSNO:<sub>FillValue</sub> = 1.e+15f ;
> 		PRECSNO:missing<sub>value</sub> = 1.e+15f ;
> 		PRECSNO:fmissing<sub>value</sub> = 1.e+15f ;
> 		PRECSNO:scale<sub>factor</sub> = 1.f ;
> 		PRECSNO:add<sub>offset</sub> = 0.f ;
> 		PRECSNO:standard<sub>name</sub> = "snowfall" ;
> 		PRECSNO:vmax = 1.e+15f ;
> 		PRECSNO:vmin = -1.e+15f ;
> 		PRECSNO:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float PRECTOT(time, lat, lon) ;
> 		PRECTOT:long<sub>name</sub> = "total<sub>precipitation</sub>" ;
> 		PRECTOT:units = "kg m-2 s-1" ;
> 		PRECTOT:<sub>FillValue</sub> = 1.e+15f ;
> 		PRECTOT:missing<sub>value</sub> = 1.e+15f ;
> 		PRECTOT:fmissing<sub>value</sub> = 1.e+15f ;
> 		PRECTOT:scale<sub>factor</sub> = 1.f ;
> 		PRECTOT:add<sub>offset</sub> = 0.f ;
> 		PRECTOT:standard<sub>name</sub> = "total<sub>precipitation</sub>" ;
> 		PRECTOT:vmax = 1.e+15f ;
> 		PRECTOT:vmin = -1.e+15f ;
> 		PRECTOT:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float PRECTOTCORR(time, lat, lon) ;
> 		PRECTOTCORR:long<sub>name</sub> = "total<sub>precipitation</sub>" ;
> 		PRECTOTCORR:units = "kg m-2 s-1" ;
> 		PRECTOTCORR:<sub>FillValue</sub> = 1.e+15f ;
> 		PRECTOTCORR:missing<sub>value</sub> = 1.e+15f ;
> 		PRECTOTCORR:fmissing<sub>value</sub> = 1.e+15f ;
> 		PRECTOTCORR:scale<sub>factor</sub> = 1.f ;
> 		PRECTOTCORR:add<sub>offset</sub> = 0.f ;
> 		PRECTOTCORR:standard<sub>name</sub> = "total<sub>precipitation</sub>" ;
> 		PRECTOTCORR:vmax = 1.e+15f ;
> 		PRECTOTCORR:vmin = -1.e+15f ;
> 		PRECTOTCORR:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float PREVTOT(time, lat, lon) ;
> 		PREVTOT:long<sub>name</sub> = "Total<sub>column</sub><sub>re</sub>-evap/subl<sub>of</sub><sub>precipitation</sub>" ;
> 		PREVTOT:units = "kg m-2 s-1" ;
> 		PREVTOT:<sub>FillValue</sub> = 1.e+15f ;
> 		PREVTOT:missing<sub>value</sub> = 1.e+15f ;
> 		PREVTOT:fmissing<sub>value</sub> = 1.e+15f ;
> 		PREVTOT:scale<sub>factor</sub> = 1.f ;
> 		PREVTOT:add<sub>offset</sub> = 0.f ;
> 		PREVTOT:standard<sub>name</sub> = "Total<sub>column</sub><sub>re</sub>-evap/subl<sub>of</sub><sub>precipitation</sub>" ;
> 		PREVTOT:vmax = 1.e+15f ;
> 		PREVTOT:vmin = -1.e+15f ;
> 		PREVTOT:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float QLML(time, lat, lon) ;
> 		QLML:long<sub>name</sub> = "surface<sub>specific</sub><sub>humidity</sub>" ;
> 		QLML:units = "1" ;
> 		QLML:<sub>FillValue</sub> = 1.e+15f ;
> 		QLML:missing<sub>value</sub> = 1.e+15f ;
> 		QLML:fmissing<sub>value</sub> = 1.e+15f ;
> 		QLML:scale<sub>factor</sub> = 1.f ;
> 		QLML:add<sub>offset</sub> = 0.f ;
> 		QLML:standard<sub>name</sub> = "surface<sub>specific</sub><sub>humidity</sub>" ;
> 		QLML:vmax = 1.e+15f ;
> 		QLML:vmin = -1.e+15f ;
> 		QLML:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float QSH(time, lat, lon) ;
> 		QSH:long<sub>name</sub> = "effective<sub>surface</sub><sub>specific</sub><sub>humidity</sub>" ;
> 		QSH:units = "kg kg-1" ;
> 		QSH:<sub>FillValue</sub> = 1.e+15f ;
> 		QSH:missing<sub>value</sub> = 1.e+15f ;
> 		QSH:fmissing<sub>value</sub> = 1.e+15f ;
> 		QSH:scale<sub>factor</sub> = 1.f ;
> 		QSH:add<sub>offset</sub> = 0.f ;
> 		QSH:standard<sub>name</sub> = "effective<sub>surface</sub><sub>specific</sub><sub>humidity</sub>" ;
> 		QSH:vmax = 1.e+15f ;
> 		QSH:vmin = -1.e+15f ;
> 		QSH:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float QSTAR(time, lat, lon) ;
> 		QSTAR:long<sub>name</sub> = "surface<sub>moisture</sub><sub>scale</sub>" ;
> 		QSTAR:units = "kg kg-1" ;
> 		QSTAR:<sub>FillValue</sub> = 1.e+15f ;
> 		QSTAR:missing<sub>value</sub> = 1.e+15f ;
> 		QSTAR:fmissing<sub>value</sub> = 1.e+15f ;
> 		QSTAR:scale<sub>factor</sub> = 1.f ;
> 		QSTAR:add<sub>offset</sub> = 0.f ;
> 		QSTAR:standard<sub>name</sub> = "surface<sub>moisture</sub><sub>scale</sub>" ;
> 		QSTAR:vmax = 1.e+15f ;
> 		QSTAR:vmin = -1.e+15f ;
> 		QSTAR:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float RHOA(time, lat, lon) ;
> 		RHOA:long<sub>name</sub> = "air<sub>density</sub><sub>at</sub><sub>surface</sub>" ;
> 		RHOA:units = "kg m-3" ;
> 		RHOA:<sub>FillValue</sub> = 1.e+15f ;
> 		RHOA:missing<sub>value</sub> = 1.e+15f ;
> 		RHOA:fmissing<sub>value</sub> = 1.e+15f ;
> 		RHOA:scale<sub>factor</sub> = 1.f ;
> 		RHOA:add<sub>offset</sub> = 0.f ;
> 		RHOA:standard<sub>name</sub> = "air<sub>density</sub><sub>at</sub><sub>surface</sub>" ;
> 		RHOA:vmax = 1.e+15f ;
> 		RHOA:vmin = -1.e+15f ;
> 		RHOA:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float RISFC(time, lat, lon) ;
> 		RISFC:long<sub>name</sub> = "surface<sub>bulk</sub><sub>richardson</sub><sub>number</sub>" ;
> 		RISFC:units = "1" ;
> 		RISFC:<sub>FillValue</sub> = 1.e+15f ;
> 		RISFC:missing<sub>value</sub> = 1.e+15f ;
> 		RISFC:fmissing<sub>value</sub> = 1.e+15f ;
> 		RISFC:scale<sub>factor</sub> = 1.f ;
> 		RISFC:add<sub>offset</sub> = 0.f ;
> 		RISFC:standard<sub>name</sub> = "surface<sub>bulk</sub><sub>richardson</sub><sub>number</sub>" ;
> 		RISFC:vmax = 1.e+15f ;
> 		RISFC:vmin = -1.e+15f ;
> 		RISFC:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float SPEED(time, lat, lon) ;
> 		SPEED:long<sub>name</sub> = "surface<sub>wind</sub><sub>speed</sub>" ;
> 		SPEED:units = "m s-1" ;
> 		SPEED:<sub>FillValue</sub> = 1.e+15f ;
> 		SPEED:missing<sub>value</sub> = 1.e+15f ;
> 		SPEED:fmissing<sub>value</sub> = 1.e+15f ;
> 		SPEED:scale<sub>factor</sub> = 1.f ;
> 		SPEED:add<sub>offset</sub> = 0.f ;
> 		SPEED:standard<sub>name</sub> = "surface<sub>wind</sub><sub>speed</sub>" ;
> 		SPEED:vmax = 1.e+15f ;
> 		SPEED:vmin = -1.e+15f ;
> 		SPEED:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float SPEEDMAX(time, lat, lon) ;
> 		SPEEDMAX:long<sub>name</sub> = "surface<sub>wind</sub><sub>speed</sub>" ;
> 		SPEEDMAX:units = "m s-1" ;
> 		SPEEDMAX:<sub>FillValue</sub> = 1.e+15f ;
> 		SPEEDMAX:missing<sub>value</sub> = 1.e+15f ;
> 		SPEEDMAX:fmissing<sub>value</sub> = 1.e+15f ;
> 		SPEEDMAX:scale<sub>factor</sub> = 1.f ;
> 		SPEEDMAX:add<sub>offset</sub> = 0.f ;
> 		SPEEDMAX:standard<sub>name</sub> = "surface<sub>wind</sub><sub>speed</sub>" ;
> 		SPEEDMAX:vmax = 1.e+15f ;
> 		SPEEDMAX:vmin = -1.e+15f ;
> 		SPEEDMAX:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float TAUGWX(time, lat, lon) ;
> 		TAUGWX:long<sub>name</sub> = "surface<sub>eastward</sub><sub>gravity</sub><sub>wave</sub><sub>stress</sub>" ;
> 		TAUGWX:units = "N m-2" ;
> 		TAUGWX:<sub>FillValue</sub> = 1.e+15f ;
> 		TAUGWX:missing<sub>value</sub> = 1.e+15f ;
> 		TAUGWX:fmissing<sub>value</sub> = 1.e+15f ;
> 		TAUGWX:scale<sub>factor</sub> = 1.f ;
> 		TAUGWX:add<sub>offset</sub> = 0.f ;
> 		TAUGWX:standard<sub>name</sub> = "surface<sub>eastward</sub><sub>gravity</sub><sub>wave</sub><sub>stress</sub>" ;
> 		TAUGWX:vmax = 1.e+15f ;
> 		TAUGWX:vmin = -1.e+15f ;
> 		TAUGWX:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float TAUGWY(time, lat, lon) ;
> 		TAUGWY:long<sub>name</sub> = "surface<sub>northward</sub><sub>gravity</sub><sub>wave</sub><sub>stress</sub>" ;
> 		TAUGWY:units = "N m-2" ;
> 		TAUGWY:<sub>FillValue</sub> = 1.e+15f ;
> 		TAUGWY:missing<sub>value</sub> = 1.e+15f ;
> 		TAUGWY:fmissing<sub>value</sub> = 1.e+15f ;
> 		TAUGWY:scale<sub>factor</sub> = 1.f ;
> 		TAUGWY:add<sub>offset</sub> = 0.f ;
> 		TAUGWY:standard<sub>name</sub> = "surface<sub>northward</sub><sub>gravity</sub><sub>wave</sub><sub>stress</sub>" ;
> 		TAUGWY:vmax = 1.e+15f ;
> 		TAUGWY:vmin = -1.e+15f ;
> 		TAUGWY:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float TAUX(time, lat, lon) ;
> 		TAUX:long<sub>name</sub> = "eastward<sub>surface</sub><sub>stress</sub>" ;
> 		TAUX:units = "N m-2" ;
> 		TAUX:<sub>FillValue</sub> = 1.e+15f ;
> 		TAUX:missing<sub>value</sub> = 1.e+15f ;
> 		TAUX:fmissing<sub>value</sub> = 1.e+15f ;
> 		TAUX:scale<sub>factor</sub> = 1.f ;
> 		TAUX:add<sub>offset</sub> = 0.f ;
> 		TAUX:standard<sub>name</sub> = "eastward<sub>surface</sub><sub>stress</sub>" ;
> 		TAUX:vmax = 1.e+15f ;
> 		TAUX:vmin = -1.e+15f ;
> 		TAUX:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float TAUY(time, lat, lon) ;
> 		TAUY:long<sub>name</sub> = "northward<sub>surface</sub><sub>stress</sub>" ;
> 		TAUY:units = "N m-2" ;
> 		TAUY:<sub>FillValue</sub> = 1.e+15f ;
> 		TAUY:missing<sub>value</sub> = 1.e+15f ;
> 		TAUY:fmissing<sub>value</sub> = 1.e+15f ;
> 		TAUY:scale<sub>factor</sub> = 1.f ;
> 		TAUY:add<sub>offset</sub> = 0.f ;
> 		TAUY:standard<sub>name</sub> = "northward<sub>surface</sub><sub>stress</sub>" ;
> 		TAUY:vmax = 1.e+15f ;
> 		TAUY:vmin = -1.e+15f ;
> 		TAUY:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float TCZPBL(time, lat, lon) ;
> 		TCZPBL:long<sub>name</sub> = "transcom<sub>planetary</sub><sub>boundary</sub><sub>layer</sub><sub>height</sub>" ;
> 		TCZPBL:units = "m" ;
> 		TCZPBL:<sub>FillValue</sub> = 1.e+15f ;
> 		TCZPBL:missing<sub>value</sub> = 1.e+15f ;
> 		TCZPBL:fmissing<sub>value</sub> = 1.e+15f ;
> 		TCZPBL:scale<sub>factor</sub> = 1.f ;
> 		TCZPBL:add<sub>offset</sub> = 0.f ;
> 		TCZPBL:standard<sub>name</sub> = "transcom<sub>planetary</sub><sub>boundary</sub><sub>layer</sub><sub>height</sub>" ;
> 		TCZPBL:vmax = 1.e+15f ;
> 		TCZPBL:vmin = -1.e+15f ;
> 		TCZPBL:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float TLML(time, lat, lon) ;
> 		TLML:long<sub>name</sub> = "surface<sub>air</sub><sub>temperature</sub>" ;
> 		TLML:units = "K" ;
> 		TLML:<sub>FillValue</sub> = 1.e+15f ;
> 		TLML:missing<sub>value</sub> = 1.e+15f ;
> 		TLML:fmissing<sub>value</sub> = 1.e+15f ;
> 		TLML:scale<sub>factor</sub> = 1.f ;
> 		TLML:add<sub>offset</sub> = 0.f ;
> 		TLML:standard<sub>name</sub> = "surface<sub>air</sub><sub>temperature</sub>" ;
> 		TLML:vmax = 1.e+15f ;
> 		TLML:vmin = -1.e+15f ;
> 		TLML:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float TSH(time, lat, lon) ;
> 		TSH:long<sub>name</sub> = "effective<sub>surface</sub><sub>skin</sub><sub>temperature</sub>" ;
> 		TSH:units = "K" ;
> 		TSH:<sub>FillValue</sub> = 1.e+15f ;
> 		TSH:missing<sub>value</sub> = 1.e+15f ;
> 		TSH:fmissing<sub>value</sub> = 1.e+15f ;
> 		TSH:scale<sub>factor</sub> = 1.f ;
> 		TSH:add<sub>offset</sub> = 0.f ;
> 		TSH:standard<sub>name</sub> = "effective<sub>surface</sub><sub>skin</sub><sub>temperature</sub>" ;
> 		TSH:vmax = 1.e+15f ;
> 		TSH:vmin = -1.e+15f ;
> 		TSH:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float TSTAR(time, lat, lon) ;
> 		TSTAR:long<sub>name</sub> = "surface<sub>temperature</sub><sub>scale</sub>" ;
> 		TSTAR:units = "K" ;
> 		TSTAR:<sub>FillValue</sub> = 1.e+15f ;
> 		TSTAR:missing<sub>value</sub> = 1.e+15f ;
> 		TSTAR:fmissing<sub>value</sub> = 1.e+15f ;
> 		TSTAR:scale<sub>factor</sub> = 1.f ;
> 		TSTAR:add<sub>offset</sub> = 0.f ;
> 		TSTAR:standard<sub>name</sub> = "surface<sub>temperature</sub><sub>scale</sub>" ;
> 		TSTAR:vmax = 1.e+15f ;
> 		TSTAR:vmin = -1.e+15f ;
> 		TSTAR:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float ULML(time, lat, lon) ;
> 		ULML:long<sub>name</sub> = "surface<sub>eastward</sub><sub>wind</sub>" ;
> 		ULML:units = "m s-1" ;
> 		ULML:<sub>FillValue</sub> = 1.e+15f ;
> 		ULML:missing<sub>value</sub> = 1.e+15f ;
> 		ULML:fmissing<sub>value</sub> = 1.e+15f ;
> 		ULML:scale<sub>factor</sub> = 1.f ;
> 		ULML:add<sub>offset</sub> = 0.f ;
> 		ULML:standard<sub>name</sub> = "surface<sub>eastward</sub><sub>wind</sub>" ;
> 		ULML:vmax = 1.e+15f ;
> 		ULML:vmin = -1.e+15f ;
> 		ULML:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float USTAR(time, lat, lon) ;
> 		USTAR:long<sub>name</sub> = "surface<sub>velocity</sub><sub>scale</sub>" ;
> 		USTAR:units = "m s-1" ;
> 		USTAR:<sub>FillValue</sub> = 1.e+15f ;
> 		USTAR:missing<sub>value</sub> = 1.e+15f ;
> 		USTAR:fmissing<sub>value</sub> = 1.e+15f ;
> 		USTAR:scale<sub>factor</sub> = 1.f ;
> 		USTAR:add<sub>offset</sub> = 0.f ;
> 		USTAR:standard<sub>name</sub> = "surface<sub>velocity</sub><sub>scale</sub>" ;
> 		USTAR:vmax = 1.e+15f ;
> 		USTAR:vmin = -1.e+15f ;
> 		USTAR:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float VLML(time, lat, lon) ;
> 		VLML:long<sub>name</sub> = "surface<sub>northward</sub><sub>wind</sub>" ;
> 		VLML:units = "m s-1" ;
> 		VLML:<sub>FillValue</sub> = 1.e+15f ;
> 		VLML:missing<sub>value</sub> = 1.e+15f ;
> 		VLML:fmissing<sub>value</sub> = 1.e+15f ;
> 		VLML:scale<sub>factor</sub> = 1.f ;
> 		VLML:add<sub>offset</sub> = 0.f ;
> 		VLML:standard<sub>name</sub> = "surface<sub>northward</sub><sub>wind</sub>" ;
> 		VLML:vmax = 1.e+15f ;
> 		VLML:vmin = -1.e+15f ;
> 		VLML:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float Z0H(time, lat, lon) ;
> 		Z0H:long<sub>name</sub> = "surface<sub>roughness</sub><sub>for</sub><sub>heat</sub>" ;
> 		Z0H:units = "m" ;
> 		Z0H:<sub>FillValue</sub> = 1.e+15f ;
> 		Z0H:missing<sub>value</sub> = 1.e+15f ;
> 		Z0H:fmissing<sub>value</sub> = 1.e+15f ;
> 		Z0H:scale<sub>factor</sub> = 1.f ;
> 		Z0H:add<sub>offset</sub> = 0.f ;
> 		Z0H:standard<sub>name</sub> = "surface<sub>roughness</sub><sub>for</sub><sub>heat</sub>" ;
> 		Z0H:vmax = 1.e+15f ;
> 		Z0H:vmin = -1.e+15f ;
> 		Z0H:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float Z0M(time, lat, lon) ;
> 		Z0M:long<sub>name</sub> = "surface<sub>roughness</sub>" ;
> 		Z0M:units = "m" ;
> 		Z0M:<sub>FillValue</sub> = 1.e+15f ;
> 		Z0M:missing<sub>value</sub> = 1.e+15f ;
> 		Z0M:fmissing<sub>value</sub> = 1.e+15f ;
> 		Z0M:scale<sub>factor</sub> = 1.f ;
> 		Z0M:add<sub>offset</sub> = 0.f ;
> 		Z0M:standard<sub>name</sub> = "surface<sub>roughness</sub>" ;
> 		Z0M:vmax = 1.e+15f ;
> 		Z0M:vmin = -1.e+15f ;
> 		Z0M:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 
> // global attributes:
> 		:History = "Original file generated: Mon Mar 23 03:29:05 2015 GMT" ;
> 		:Comment = "GMAO filename: d5124<sub>m2</sub><sub>jan00.tavg1</sub><sub>2d</sub><sub>flx</sub><sub>Nx.20100101.nc4</sub>" ;
> 		:Filename = "MERRA2<sub>300.tavg1</sub><sub>2d</sub><sub>flx</sub><sub>Nx.20100101.nc4</sub>" ;
> 		:Conventions = "CF-1" ;
> 		:Institution = "NASA Global Modeling and Assimilation Office" ;
> 		:References = "<http://gmao.gsfc.nasa.gov>" ;
> 		:Format = "NetCDF-4/HDF-5" ;
> 		:SpatialCoverage = "global" ;
> 		:VersionID = "5.12.4" ;
> 		:TemporalRange = "1980-01-01 -> 2016-12-31" ;
> 		:identifier<sub>product</sub><sub>doi</sub><sub>authority</sub> = "<http://dx.doi.org/>" ;
> 		:ShortName = "M2T1NXFLX" ;
> 		:GranuleID = "MERRA2<sub>300.tavg1</sub><sub>2d</sub><sub>flx</sub><sub>Nx.20100101.nc4</sub>" ;
> 		:ProductionDateTime = "Original file generated: Mon Mar 23 03:29:05 2015 GMT" ;
> 		:LongName = "MERRA2 tavg1<sub>2d</sub><sub>flx</sub><sub>Nx</sub>: 2d,1-Hourly,Time-Averaged,Single-Level,Assimilation,Surface Flux Diagnostics" ;
> 		:Title = "MERRA2 tavg1<sub>2d</sub><sub>flx</sub><sub>Nx</sub>: 2d,1-Hourly,Time-Averaged,Single-Level,Assimilation,Surface Flux Diagnostics" ;
> 		:SouthernmostLatitude = "-90.0" ;
> 		:NorthernmostLatitude = "90.0" ;
> 		:WesternmostLongitude = "-180.0" ;
> 		:EasternmostLongitude = "179.375" ;
> 		:LatitudeResolution = "0.5" ;
> 		:LongitudeResolution = "0.625" ;
> 		:DataResolution = "0.5 x 0.625" ;
> 		:Source = "CVS tag: GEOSadas-5<sub>12</sub><sub>4</sub>" ;
> 		:Contact = "<http://gmao.gsfc.nasa.gov>" ;
> 		:identifier<sub>product</sub><sub>doi</sub> = "10.5067/7MCPBJ41Y0K6" ;
> 		:RangeBeginningDate = "2010-01-01" ;
> 		:RangeBeginningTime = "00:00:00.000000" ;
> 		:RangeEndingDate = "2010-01-01" ;
> 		:RangeEndingTime = "23:59:59.000000" ;
> }


<a id="org03ca400"></a>

### Ocean Surface Diagnostics ('ocn')

Not presently on any projects associated with gadi. I have inquired with the one
of the project leaders at UTas if there is any scope to 'host' (for the lack of
a better term) this portion of MERRAv2 in gadi project 'ua8'


<a id="org5063a13"></a>

### Radiation Diagnostics ('rad')

I have downloaded years 2010-2019 files from gadi to my [local repository](file:///volumes/ioa01/reanalysis/MERRAv2/rad). Here
is the header of the first file:

> netcdf MERRA2<sub>300.tavg1</sub><sub>2d</sub><sub>rad</sub><sub>Nx.20100101</sub> {
> dimensions:
> 	lon = 576 ;
> 	lat = 361 ;
> 	time = UNLIMITED ; // (24 currently)
> variables:
> 	double lon(lon) ;
> 		lon:long<sub>name</sub> = "longitude" ;
> 		lon:units = "degrees<sub>east</sub>" ;
> 		lon:vmax = 1.e+15f ;
> 		lon:vmin = -1.e+15f ;
> 		lon:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	double lat(lat) ;
> 		lat:long<sub>name</sub> = "latitude" ;
> 		lat:units = "degrees<sub>north</sub>" ;
> 		lat:vmax = 1.e+15f ;
> 		lat:vmin = -1.e+15f ;
> 		lat:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	int time(time) ;
> 		time:long<sub>name</sub> = "time" ;
> 		time:units = "minutes since 2010-01-01 00:30:00" ;
> 		time:time<sub>increment</sub> = 10000 ;
> 		time:begin<sub>date</sub> = 20100101 ;
> 		time:begin<sub>time</sub> = 3000 ;
> 		time:vmax = 1.e+15f ;
> 		time:vmin = -1.e+15f ;
> 		time:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float ALBEDO(time, lat, lon) ;
> 		ALBEDO:long<sub>name</sub> = "surface<sub>albedo</sub>" ;
> 		ALBEDO:units = "1" ;
> 		ALBEDO:<sub>FillValue</sub> = 1.e+15f ;
> 		ALBEDO:missing<sub>value</sub> = 1.e+15f ;
> 		ALBEDO:fmissing<sub>value</sub> = 1.e+15f ;
> 		ALBEDO:scale<sub>factor</sub> = 1.f ;
> 		ALBEDO:add<sub>offset</sub> = 0.f ;
> 		ALBEDO:standard<sub>name</sub> = "surface<sub>albedo</sub>" ;
> 		ALBEDO:vmax = 1.e+15f ;
> 		ALBEDO:vmin = -1.e+15f ;
> 		ALBEDO:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float ALBNIRDF(time, lat, lon) ;
> 		ALBNIRDF:long<sub>name</sub> = "surface<sub>albedo</sub><sub>for</sub><sub>near</sub><sub>infrared</sub><sub>diffuse</sub>" ;
> 		ALBNIRDF:units = "1" ;
> 		ALBNIRDF:<sub>FillValue</sub> = 1.e+15f ;
> 		ALBNIRDF:missing<sub>value</sub> = 1.e+15f ;
> 		ALBNIRDF:fmissing<sub>value</sub> = 1.e+15f ;
> 		ALBNIRDF:scale<sub>factor</sub> = 1.f ;
> 		ALBNIRDF:add<sub>offset</sub> = 0.f ;
> 		ALBNIRDF:standard<sub>name</sub> = "surface<sub>albedo</sub><sub>for</sub><sub>near</sub><sub>infrared</sub><sub>diffuse</sub>" ;
> 		ALBNIRDF:vmax = 1.e+15f ;
> 		ALBNIRDF:vmin = -1.e+15f ;
> 		ALBNIRDF:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float ALBNIRDR(time, lat, lon) ;
> 		ALBNIRDR:long<sub>name</sub> = "surface<sub>albedo</sub><sub>for</sub><sub>near</sub><sub>infrared</sub><sub>beam</sub>" ;
> 		ALBNIRDR:units = "1" ;
> 		ALBNIRDR:<sub>FillValue</sub> = 1.e+15f ;
> 		ALBNIRDR:missing<sub>value</sub> = 1.e+15f ;
> 		ALBNIRDR:fmissing<sub>value</sub> = 1.e+15f ;
> 		ALBNIRDR:scale<sub>factor</sub> = 1.f ;
> 		ALBNIRDR:add<sub>offset</sub> = 0.f ;
> 		ALBNIRDR:standard<sub>name</sub> = "surface<sub>albedo</sub><sub>for</sub><sub>near</sub><sub>infrared</sub><sub>beam</sub>" ;
> 		ALBNIRDR:vmax = 1.e+15f ;
> 		ALBNIRDR:vmin = -1.e+15f ;
> 		ALBNIRDR:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float ALBVISDF(time, lat, lon) ;
> 		ALBVISDF:long<sub>name</sub> = "surface<sub>albedo</sub><sub>for</sub><sub>visible</sub><sub>diffuse</sub>" ;
> 		ALBVISDF:units = "1" ;
> 		ALBVISDF:<sub>FillValue</sub> = 1.e+15f ;
> 		ALBVISDF:missing<sub>value</sub> = 1.e+15f ;
> 		ALBVISDF:fmissing<sub>value</sub> = 1.e+15f ;
> 		ALBVISDF:scale<sub>factor</sub> = 1.f ;
> 		ALBVISDF:add<sub>offset</sub> = 0.f ;
> 		ALBVISDF:standard<sub>name</sub> = "surface<sub>albedo</sub><sub>for</sub><sub>visible</sub><sub>diffuse</sub>" ;
> 		ALBVISDF:vmax = 1.e+15f ;
> 		ALBVISDF:vmin = -1.e+15f ;
> 		ALBVISDF:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float ALBVISDR(time, lat, lon) ;
> 		ALBVISDR:long<sub>name</sub> = "surface<sub>albedo</sub><sub>for</sub><sub>visible</sub><sub>beam</sub>" ;
> 		ALBVISDR:units = "1" ;
> 		ALBVISDR:<sub>FillValue</sub> = 1.e+15f ;
> 		ALBVISDR:missing<sub>value</sub> = 1.e+15f ;
> 		ALBVISDR:fmissing<sub>value</sub> = 1.e+15f ;
> 		ALBVISDR:scale<sub>factor</sub> = 1.f ;
> 		ALBVISDR:add<sub>offset</sub> = 0.f ;
> 		ALBVISDR:standard<sub>name</sub> = "surface<sub>albedo</sub><sub>for</sub><sub>visible</sub><sub>beam</sub>" ;
> 		ALBVISDR:vmax = 1.e+15f ;
> 		ALBVISDR:vmin = -1.e+15f ;
> 		ALBVISDR:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float CLDHGH(time, lat, lon) ;
> 		CLDHGH:long<sub>name</sub> = "cloud<sub>area</sub><sub>fraction</sub><sub>for</sub><sub>high</sub><sub>clouds</sub>" ;
> 		CLDHGH:units = "1" ;
> 		CLDHGH:<sub>FillValue</sub> = 1.e+15f ;
> 		CLDHGH:missing<sub>value</sub> = 1.e+15f ;
> 		CLDHGH:fmissing<sub>value</sub> = 1.e+15f ;
> 		CLDHGH:scale<sub>factor</sub> = 1.f ;
> 		CLDHGH:add<sub>offset</sub> = 0.f ;
> 		CLDHGH:standard<sub>name</sub> = "cloud<sub>area</sub><sub>fraction</sub><sub>for</sub><sub>high</sub><sub>clouds</sub>" ;
> 		CLDHGH:vmax = 1.e+15f ;
> 		CLDHGH:vmin = -1.e+15f ;
> 		CLDHGH:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float CLDLOW(time, lat, lon) ;
> 		CLDLOW:long<sub>name</sub> = "cloud<sub>area</sub><sub>fraction</sub><sub>for</sub><sub>low</sub><sub>clouds</sub>" ;
> 		CLDLOW:units = "1" ;
> 		CLDLOW:<sub>FillValue</sub> = 1.e+15f ;
> 		CLDLOW:missing<sub>value</sub> = 1.e+15f ;
> 		CLDLOW:fmissing<sub>value</sub> = 1.e+15f ;
> 		CLDLOW:scale<sub>factor</sub> = 1.f ;
> 		CLDLOW:add<sub>offset</sub> = 0.f ;
> 		CLDLOW:standard<sub>name</sub> = "cloud<sub>area</sub><sub>fraction</sub><sub>for</sub><sub>low</sub><sub>clouds</sub>" ;
> 		CLDLOW:vmax = 1.e+15f ;
> 		CLDLOW:vmin = -1.e+15f ;
> 		CLDLOW:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float CLDMID(time, lat, lon) ;
> 		CLDMID:long<sub>name</sub> = "cloud<sub>area</sub><sub>fraction</sub><sub>for</sub><sub>middle</sub><sub>clouds</sub>" ;
> 		CLDMID:units = "1" ;
> 		CLDMID:<sub>FillValue</sub> = 1.e+15f ;
> 		CLDMID:missing<sub>value</sub> = 1.e+15f ;
> 		CLDMID:fmissing<sub>value</sub> = 1.e+15f ;
> 		CLDMID:scale<sub>factor</sub> = 1.f ;
> 		CLDMID:add<sub>offset</sub> = 0.f ;
> 		CLDMID:standard<sub>name</sub> = "cloud<sub>area</sub><sub>fraction</sub><sub>for</sub><sub>middle</sub><sub>clouds</sub>" ;
> 		CLDMID:vmax = 1.e+15f ;
> 		CLDMID:vmin = -1.e+15f ;
> 		CLDMID:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float CLDTOT(time, lat, lon) ;
> 		CLDTOT:long<sub>name</sub> = "total<sub>cloud</sub><sub>area</sub><sub>fraction</sub>" ;
> 		CLDTOT:units = "1" ;
> 		CLDTOT:<sub>FillValue</sub> = 1.e+15f ;
> 		CLDTOT:missing<sub>value</sub> = 1.e+15f ;
> 		CLDTOT:fmissing<sub>value</sub> = 1.e+15f ;
> 		CLDTOT:scale<sub>factor</sub> = 1.f ;
> 		CLDTOT:add<sub>offset</sub> = 0.f ;
> 		CLDTOT:standard<sub>name</sub> = "total<sub>cloud</sub><sub>area</sub><sub>fraction</sub>" ;
> 		CLDTOT:vmax = 1.e+15f ;
> 		CLDTOT:vmin = -1.e+15f ;
> 		CLDTOT:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float EMIS(time, lat, lon) ;
> 		EMIS:long<sub>name</sub> = "surface<sub>emissivity</sub>" ;
> 		EMIS:units = "1" ;
> 		EMIS:<sub>FillValue</sub> = 1.e+15f ;
> 		EMIS:missing<sub>value</sub> = 1.e+15f ;
> 		EMIS:fmissing<sub>value</sub> = 1.e+15f ;
> 		EMIS:scale<sub>factor</sub> = 1.f ;
> 		EMIS:add<sub>offset</sub> = 0.f ;
> 		EMIS:standard<sub>name</sub> = "surface<sub>emissivity</sub>" ;
> 		EMIS:vmax = 1.e+15f ;
> 		EMIS:vmin = -1.e+15f ;
> 		EMIS:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float LWGAB(time, lat, lon) ;
> 		LWGAB:long<sub>name</sub> = "surface<sub>absorbed</sub><sub>longwave</sub><sub>radiation</sub>" ;
> 		LWGAB:units = "W m-2" ;
> 		LWGAB:<sub>FillValue</sub> = 1.e+15f ;
> 		LWGAB:missing<sub>value</sub> = 1.e+15f ;
> 		LWGAB:fmissing<sub>value</sub> = 1.e+15f ;
> 		LWGAB:scale<sub>factor</sub> = 1.f ;
> 		LWGAB:add<sub>offset</sub> = 0.f ;
> 		LWGAB:standard<sub>name</sub> = "surface<sub>absorbed</sub><sub>longwave</sub><sub>radiation</sub>" ;
> 		LWGAB:vmax = 1.e+15f ;
> 		LWGAB:vmin = -1.e+15f ;
> 		LWGAB:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float LWGABCLR(time, lat, lon) ;
> 		LWGABCLR:long<sub>name</sub> = "surface<sub>absorbed</sub><sub>longwave</sub><sub>radiation</sub><sub>assuming</sub><sub>clear</sub><sub>sky</sub>" ;
> 		LWGABCLR:units = "W m-2" ;
> 		LWGABCLR:<sub>FillValue</sub> = 1.e+15f ;
> 		LWGABCLR:missing<sub>value</sub> = 1.e+15f ;
> 		LWGABCLR:fmissing<sub>value</sub> = 1.e+15f ;
> 		LWGABCLR:scale<sub>factor</sub> = 1.f ;
> 		LWGABCLR:add<sub>offset</sub> = 0.f ;
> 		LWGABCLR:standard<sub>name</sub> = "surface<sub>absorbed</sub><sub>longwave</sub><sub>radiation</sub><sub>assuming</sub><sub>clear</sub><sub>sky</sub>" ;
> 		LWGABCLR:vmax = 1.e+15f ;
> 		LWGABCLR:vmin = -1.e+15f ;
> 		LWGABCLR:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float LWGABCLRCLN(time, lat, lon) ;
> 		LWGABCLRCLN:long<sub>name</sub> = "surface<sub>absorbed</sub><sub>longwave</sub><sub>radiation</sub><sub>assuming</sub><sub>clear</sub><sub>sky</sub><sub>and</sub><sub>no</sub><sub>aerosol</sub>" ;
> 		LWGABCLRCLN:units = "W m-2" ;
> 		LWGABCLRCLN:<sub>FillValue</sub> = 1.e+15f ;
> 		LWGABCLRCLN:missing<sub>value</sub> = 1.e+15f ;
> 		LWGABCLRCLN:fmissing<sub>value</sub> = 1.e+15f ;
> 		LWGABCLRCLN:scale<sub>factor</sub> = 1.f ;
> 		LWGABCLRCLN:add<sub>offset</sub> = 0.f ;
> 		LWGABCLRCLN:standard<sub>name</sub> = "surface<sub>absorbed</sub><sub>longwave</sub><sub>radiation</sub><sub>assuming</sub><sub>clear</sub><sub>sky</sub><sub>and</sub><sub>no</sub><sub>aerosol</sub>" ;
> 		LWGABCLRCLN:vmax = 1.e+15f ;
> 		LWGABCLRCLN:vmin = -1.e+15f ;
> 		LWGABCLRCLN:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float LWGEM(time, lat, lon) ;
> 		LWGEM:long<sub>name</sub> = "longwave<sub>flux</sub><sub>emitted</sub><sub>from</sub><sub>surface</sub>" ;
> 		LWGEM:units = "W m-2" ;
> 		LWGEM:<sub>FillValue</sub> = 1.e+15f ;
> 		LWGEM:missing<sub>value</sub> = 1.e+15f ;
> 		LWGEM:fmissing<sub>value</sub> = 1.e+15f ;
> 		LWGEM:scale<sub>factor</sub> = 1.f ;
> 		LWGEM:add<sub>offset</sub> = 0.f ;
> 		LWGEM:standard<sub>name</sub> = "longwave<sub>flux</sub><sub>emitted</sub><sub>from</sub><sub>surface</sub>" ;
> 		LWGEM:vmax = 1.e+15f ;
> 		LWGEM:vmin = -1.e+15f ;
> 		LWGEM:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float LWGNT(time, lat, lon) ;
> 		LWGNT:long<sub>name</sub> = "surface<sub>net</sub><sub>downward</sub><sub>longwave</sub><sub>flux</sub>" ;
> 		LWGNT:units = "W m-2" ;
> 		LWGNT:<sub>FillValue</sub> = 1.e+15f ;
> 		LWGNT:missing<sub>value</sub> = 1.e+15f ;
> 		LWGNT:fmissing<sub>value</sub> = 1.e+15f ;
> 		LWGNT:scale<sub>factor</sub> = 1.f ;
> 		LWGNT:add<sub>offset</sub> = 0.f ;
> 		LWGNT:standard<sub>name</sub> = "surface<sub>net</sub><sub>downward</sub><sub>longwave</sub><sub>flux</sub>" ;
> 		LWGNT:vmax = 1.e+15f ;
> 		LWGNT:vmin = -1.e+15f ;
> 		LWGNT:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float LWGNTCLR(time, lat, lon) ;
> 		LWGNTCLR:long<sub>name</sub> = "surface<sub>net</sub><sub>downward</sub><sub>longwave</sub><sub>flux</sub><sub>assuming</sub><sub>clear</sub><sub>sky</sub>" ;
> 		LWGNTCLR:units = "W m-2" ;
> 		LWGNTCLR:<sub>FillValue</sub> = 1.e+15f ;
> 		LWGNTCLR:missing<sub>value</sub> = 1.e+15f ;
> 		LWGNTCLR:fmissing<sub>value</sub> = 1.e+15f ;
> 		LWGNTCLR:scale<sub>factor</sub> = 1.f ;
> 		LWGNTCLR:add<sub>offset</sub> = 0.f ;
> 		LWGNTCLR:standard<sub>name</sub> = "surface<sub>net</sub><sub>downward</sub><sub>longwave</sub><sub>flux</sub><sub>assuming</sub><sub>clear</sub><sub>sky</sub>" ;
> 		LWGNTCLR:vmax = 1.e+15f ;
> 		LWGNTCLR:vmin = -1.e+15f ;
> 		LWGNTCLR:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float LWGNTCLRCLN(time, lat, lon) ;
> 		LWGNTCLRCLN:long<sub>name</sub> = "surface<sub>net</sub><sub>downward</sub><sub>longwave</sub><sub>flux</sub><sub>assuming</sub><sub>clear</sub><sub>sky</sub><sub>and</sub><sub>no</sub><sub>aerosol</sub>" ;
> 		LWGNTCLRCLN:units = "W m-2" ;
> 		LWGNTCLRCLN:<sub>FillValue</sub> = 1.e+15f ;
> 		LWGNTCLRCLN:missing<sub>value</sub> = 1.e+15f ;
> 		LWGNTCLRCLN:fmissing<sub>value</sub> = 1.e+15f ;
> 		LWGNTCLRCLN:scale<sub>factor</sub> = 1.f ;
> 		LWGNTCLRCLN:add<sub>offset</sub> = 0.f ;
> 		LWGNTCLRCLN:standard<sub>name</sub> = "surface<sub>net</sub><sub>downward</sub><sub>longwave</sub><sub>flux</sub><sub>assuming</sub><sub>clear</sub><sub>sky</sub><sub>and</sub><sub>no</sub><sub>aerosol</sub>" ;
> 		LWGNTCLRCLN:vmax = 1.e+15f ;
> 		LWGNTCLRCLN:vmin = -1.e+15f ;
> 		LWGNTCLRCLN:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float LWTUP(time, lat, lon) ;
> 		LWTUP:long<sub>name</sub> = "upwelling<sub>longwave</sub><sub>flux</sub><sub>at</sub><sub>toa</sub>" ;
> 		LWTUP:units = "W m-2" ;
> 		LWTUP:<sub>FillValue</sub> = 1.e+15f ;
> 		LWTUP:missing<sub>value</sub> = 1.e+15f ;
> 		LWTUP:fmissing<sub>value</sub> = 1.e+15f ;
> 		LWTUP:scale<sub>factor</sub> = 1.f ;
> 		LWTUP:add<sub>offset</sub> = 0.f ;
> 		LWTUP:standard<sub>name</sub> = "upwelling<sub>longwave</sub><sub>flux</sub><sub>at</sub><sub>toa</sub>" ;
> 		LWTUP:vmax = 1.e+15f ;
> 		LWTUP:vmin = -1.e+15f ;
> 		LWTUP:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float LWTUPCLR(time, lat, lon) ;
> 		LWTUPCLR:long<sub>name</sub> = "upwelling<sub>longwave</sub><sub>flux</sub><sub>at</sub><sub>toa</sub><sub>assuming</sub><sub>clear</sub><sub>sky</sub>" ;
> 		LWTUPCLR:units = "W m-2" ;
> 		LWTUPCLR:<sub>FillValue</sub> = 1.e+15f ;
> 		LWTUPCLR:missing<sub>value</sub> = 1.e+15f ;
> 		LWTUPCLR:fmissing<sub>value</sub> = 1.e+15f ;
> 		LWTUPCLR:scale<sub>factor</sub> = 1.f ;
> 		LWTUPCLR:add<sub>offset</sub> = 0.f ;
> 		LWTUPCLR:standard<sub>name</sub> = "upwelling<sub>longwave</sub><sub>flux</sub><sub>at</sub><sub>toa</sub><sub>assuming</sub><sub>clear</sub><sub>sky</sub>" ;
> 		LWTUPCLR:vmax = 1.e+15f ;
> 		LWTUPCLR:vmin = -1.e+15f ;
> 		LWTUPCLR:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float LWTUPCLRCLN(time, lat, lon) ;
> 		LWTUPCLRCLN:long<sub>name</sub> = "upwelling<sub>longwave</sub><sub>flux</sub><sub>at</sub><sub>toa</sub><sub>assuming</sub><sub>clear</sub><sub>sky</sub><sub>and</sub><sub>no</sub><sub>aerosol</sub>" ;
> 		LWTUPCLRCLN:units = "W m-2" ;
> 		LWTUPCLRCLN:<sub>FillValue</sub> = 1.e+15f ;
> 		LWTUPCLRCLN:missing<sub>value</sub> = 1.e+15f ;
> 		LWTUPCLRCLN:fmissing<sub>value</sub> = 1.e+15f ;
> 		LWTUPCLRCLN:scale<sub>factor</sub> = 1.f ;
> 		LWTUPCLRCLN:add<sub>offset</sub> = 0.f ;
> 		LWTUPCLRCLN:standard<sub>name</sub> = "upwelling<sub>longwave</sub><sub>flux</sub><sub>at</sub><sub>toa</sub><sub>assuming</sub><sub>clear</sub><sub>sky</sub><sub>and</sub><sub>no</sub><sub>aerosol</sub>" ;
> 		LWTUPCLRCLN:vmax = 1.e+15f ;
> 		LWTUPCLRCLN:vmin = -1.e+15f ;
> 		LWTUPCLRCLN:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float SWGDN(time, lat, lon) ;
> 		SWGDN:long<sub>name</sub> = "surface<sub>incoming</sub><sub>shortwave</sub><sub>flux</sub>" ;
> 		SWGDN:units = "W m-2" ;
> 		SWGDN:<sub>FillValue</sub> = 1.e+15f ;
> 		SWGDN:missing<sub>value</sub> = 1.e+15f ;
> 		SWGDN:fmissing<sub>value</sub> = 1.e+15f ;
> 		SWGDN:scale<sub>factor</sub> = 1.f ;
> 		SWGDN:add<sub>offset</sub> = 0.f ;
> 		SWGDN:standard<sub>name</sub> = "surface<sub>incoming</sub><sub>shortwave</sub><sub>flux</sub>" ;
> 		SWGDN:vmax = 1.e+15f ;
> 		SWGDN:vmin = -1.e+15f ;
> 		SWGDN:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float SWGDNCLR(time, lat, lon) ;
> 		SWGDNCLR:long<sub>name</sub> = "surface<sub>incoming</sub><sub>shortwave</sub><sub>flux</sub><sub>assuming</sub><sub>clear</sub><sub>sky</sub>" ;
> 		SWGDNCLR:units = "W m-2" ;
> 		SWGDNCLR:<sub>FillValue</sub> = 1.e+15f ;
> 		SWGDNCLR:missing<sub>value</sub> = 1.e+15f ;
> 		SWGDNCLR:fmissing<sub>value</sub> = 1.e+15f ;
> 		SWGDNCLR:scale<sub>factor</sub> = 1.f ;
> 		SWGDNCLR:add<sub>offset</sub> = 0.f ;
> 		SWGDNCLR:standard<sub>name</sub> = "surface<sub>incoming</sub><sub>shortwave</sub><sub>flux</sub><sub>assuming</sub><sub>clear</sub><sub>sky</sub>" ;
> 		SWGDNCLR:vmax = 1.e+15f ;
> 		SWGDNCLR:vmin = -1.e+15f ;
> 		SWGDNCLR:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float SWGNT(time, lat, lon) ;
> 		SWGNT:long<sub>name</sub> = "surface<sub>net</sub><sub>downward</sub><sub>shortwave</sub><sub>flux</sub>" ;
> 		SWGNT:units = "W m-2" ;
> 		SWGNT:<sub>FillValue</sub> = 1.e+15f ;
> 		SWGNT:missing<sub>value</sub> = 1.e+15f ;
> 		SWGNT:fmissing<sub>value</sub> = 1.e+15f ;
> 		SWGNT:scale<sub>factor</sub> = 1.f ;
> 		SWGNT:add<sub>offset</sub> = 0.f ;
> 		SWGNT:standard<sub>name</sub> = "surface<sub>net</sub><sub>downward</sub><sub>shortwave</sub><sub>flux</sub>" ;
> 		SWGNT:vmax = 1.e+15f ;
> 		SWGNT:vmin = -1.e+15f ;
> 		SWGNT:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float SWGNTCLN(time, lat, lon) ;
> 		SWGNTCLN:long<sub>name</sub> = "surface<sub>net</sub><sub>downward</sub><sub>shortwave</sub><sub>flux</sub><sub>assuming</sub><sub>no</sub><sub>aerosol</sub>" ;
> 		SWGNTCLN:units = "W m-2" ;
> 		SWGNTCLN:<sub>FillValue</sub> = 1.e+15f ;
> 		SWGNTCLN:missing<sub>value</sub> = 1.e+15f ;
> 		SWGNTCLN:fmissing<sub>value</sub> = 1.e+15f ;
> 		SWGNTCLN:scale<sub>factor</sub> = 1.f ;
> 		SWGNTCLN:add<sub>offset</sub> = 0.f ;
> 		SWGNTCLN:standard<sub>name</sub> = "surface<sub>net</sub><sub>downward</sub><sub>shortwave</sub><sub>flux</sub><sub>assuming</sub><sub>no</sub><sub>aerosol</sub>" ;
> 		SWGNTCLN:vmax = 1.e+15f ;
> 		SWGNTCLN:vmin = -1.e+15f ;
> 		SWGNTCLN:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float SWGNTCLR(time, lat, lon) ;
> 		SWGNTCLR:long<sub>name</sub> = "surface<sub>net</sub><sub>downward</sub><sub>shortwave</sub><sub>flux</sub><sub>assuming</sub><sub>clear</sub><sub>sky</sub>" ;
> 		SWGNTCLR:units = "W m-2" ;
> 		SWGNTCLR:<sub>FillValue</sub> = 1.e+15f ;
> 		SWGNTCLR:missing<sub>value</sub> = 1.e+15f ;
> 		SWGNTCLR:fmissing<sub>value</sub> = 1.e+15f ;
> 		SWGNTCLR:scale<sub>factor</sub> = 1.f ;
> 		SWGNTCLR:add<sub>offset</sub> = 0.f ;
> 		SWGNTCLR:standard<sub>name</sub> = "surface<sub>net</sub><sub>downward</sub><sub>shortwave</sub><sub>flux</sub><sub>assuming</sub><sub>clear</sub><sub>sky</sub>" ;
> 		SWGNTCLR:vmax = 1.e+15f ;
> 		SWGNTCLR:vmin = -1.e+15f ;
> 		SWGNTCLR:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float SWGNTCLRCLN(time, lat, lon) ;
> 		SWGNTCLRCLN:long<sub>name</sub> = "surface<sub>net</sub><sub>downward</sub><sub>shortwave</sub><sub>flux</sub><sub>assuming</sub><sub>clear</sub><sub>sky</sub><sub>and</sub><sub>no</sub><sub>aerosol</sub>" ;
> 		SWGNTCLRCLN:units = "W m-2" ;
> 		SWGNTCLRCLN:<sub>FillValue</sub> = 1.e+15f ;
> 		SWGNTCLRCLN:missing<sub>value</sub> = 1.e+15f ;
> 		SWGNTCLRCLN:fmissing<sub>value</sub> = 1.e+15f ;
> 		SWGNTCLRCLN:scale<sub>factor</sub> = 1.f ;
> 		SWGNTCLRCLN:add<sub>offset</sub> = 0.f ;
> 		SWGNTCLRCLN:standard<sub>name</sub> = "surface<sub>net</sub><sub>downward</sub><sub>shortwave</sub><sub>flux</sub><sub>assuming</sub><sub>clear</sub><sub>sky</sub><sub>and</sub><sub>no</sub><sub>aerosol</sub>" ;
> 		SWGNTCLRCLN:vmax = 1.e+15f ;
> 		SWGNTCLRCLN:vmin = -1.e+15f ;
> 		SWGNTCLRCLN:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float SWTDN(time, lat, lon) ;
> 		SWTDN:long<sub>name</sub> = "toa<sub>incoming</sub><sub>shortwave</sub><sub>flux</sub>" ;
> 		SWTDN:units = "W m-2" ;
> 		SWTDN:<sub>FillValue</sub> = 1.e+15f ;
> 		SWTDN:missing<sub>value</sub> = 1.e+15f ;
> 		SWTDN:fmissing<sub>value</sub> = 1.e+15f ;
> 		SWTDN:scale<sub>factor</sub> = 1.f ;
> 		SWTDN:add<sub>offset</sub> = 0.f ;
> 		SWTDN:standard<sub>name</sub> = "toa<sub>incoming</sub><sub>shortwave</sub><sub>flux</sub>" ;
> 		SWTDN:vmax = 1.e+15f ;
> 		SWTDN:vmin = -1.e+15f ;
> 		SWTDN:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float SWTNT(time, lat, lon) ;
> 		SWTNT:long<sub>name</sub> = "toa<sub>net</sub><sub>downward</sub><sub>shortwave</sub><sub>flux</sub>" ;
> 		SWTNT:units = "W m-2" ;
> 		SWTNT:<sub>FillValue</sub> = 1.e+15f ;
> 		SWTNT:missing<sub>value</sub> = 1.e+15f ;
> 		SWTNT:fmissing<sub>value</sub> = 1.e+15f ;
> 		SWTNT:scale<sub>factor</sub> = 1.f ;
> 		SWTNT:add<sub>offset</sub> = 0.f ;
> 		SWTNT:standard<sub>name</sub> = "toa<sub>net</sub><sub>downward</sub><sub>shortwave</sub><sub>flux</sub>" ;
> 		SWTNT:vmax = 1.e+15f ;
> 		SWTNT:vmin = -1.e+15f ;
> 		SWTNT:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float SWTNTCLN(time, lat, lon) ;
> 		SWTNTCLN:long<sub>name</sub> = "toa<sub>net</sub><sub>downward</sub><sub>shortwave</sub><sub>flux</sub><sub>assuming</sub><sub>no</sub><sub>aerosol</sub>" ;
> 		SWTNTCLN:units = "W m-2" ;
> 		SWTNTCLN:<sub>FillValue</sub> = 1.e+15f ;
> 		SWTNTCLN:missing<sub>value</sub> = 1.e+15f ;
> 		SWTNTCLN:fmissing<sub>value</sub> = 1.e+15f ;
> 		SWTNTCLN:scale<sub>factor</sub> = 1.f ;
> 		SWTNTCLN:add<sub>offset</sub> = 0.f ;
> 		SWTNTCLN:standard<sub>name</sub> = "toa<sub>net</sub><sub>downward</sub><sub>shortwave</sub><sub>flux</sub><sub>assuming</sub><sub>no</sub><sub>aerosol</sub>" ;
> 		SWTNTCLN:vmax = 1.e+15f ;
> 		SWTNTCLN:vmin = -1.e+15f ;
> 		SWTNTCLN:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float SWTNTCLR(time, lat, lon) ;
> 		SWTNTCLR:long<sub>name</sub> = "toa<sub>net</sub><sub>downward</sub><sub>shortwave</sub><sub>flux</sub><sub>assuming</sub><sub>clear</sub><sub>sky</sub>" ;
> 		SWTNTCLR:units = "W m-2" ;
> 		SWTNTCLR:<sub>FillValue</sub> = 1.e+15f ;
> 		SWTNTCLR:missing<sub>value</sub> = 1.e+15f ;
> 		SWTNTCLR:fmissing<sub>value</sub> = 1.e+15f ;
> 		SWTNTCLR:scale<sub>factor</sub> = 1.f ;
> 		SWTNTCLR:add<sub>offset</sub> = 0.f ;
> 		SWTNTCLR:standard<sub>name</sub> = "toa<sub>net</sub><sub>downward</sub><sub>shortwave</sub><sub>flux</sub><sub>assuming</sub><sub>clear</sub><sub>sky</sub>" ;
> 		SWTNTCLR:vmax = 1.e+15f ;
> 		SWTNTCLR:vmin = -1.e+15f ;
> 		SWTNTCLR:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float SWTNTCLRCLN(time, lat, lon) ;
> 		SWTNTCLRCLN:long<sub>name</sub> = "toa<sub>net</sub><sub>downward</sub><sub>shortwave</sub><sub>flux</sub><sub>assuming</sub><sub>clear</sub><sub>sky</sub><sub>and</sub><sub>no</sub><sub>aerosol</sub>" ;
> 		SWTNTCLRCLN:units = "W m-2" ;
> 		SWTNTCLRCLN:<sub>FillValue</sub> = 1.e+15f ;
> 		SWTNTCLRCLN:missing<sub>value</sub> = 1.e+15f ;
> 		SWTNTCLRCLN:fmissing<sub>value</sub> = 1.e+15f ;
> 		SWTNTCLRCLN:scale<sub>factor</sub> = 1.f ;
> 		SWTNTCLRCLN:add<sub>offset</sub> = 0.f ;
> 		SWTNTCLRCLN:standard<sub>name</sub> = "toa<sub>net</sub><sub>downward</sub><sub>shortwave</sub><sub>flux</sub><sub>assuming</sub><sub>clear</sub><sub>sky</sub><sub>and</sub><sub>no</sub><sub>aerosol</sub>" ;
> 		SWTNTCLRCLN:vmax = 1.e+15f ;
> 		SWTNTCLRCLN:vmin = -1.e+15f ;
> 		SWTNTCLRCLN:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float TAUHGH(time, lat, lon) ;
> 		TAUHGH:long<sub>name</sub> = "in<sub>cloud</sub><sub>optical</sub><sub>thickness</sub><sub>of</sub><sub>high</sub><sub>clouds</sub>(EXPORT)" ;
> 		TAUHGH:units = "1" ;
> 		TAUHGH:<sub>FillValue</sub> = 1.e+15f ;
> 		TAUHGH:missing<sub>value</sub> = 1.e+15f ;
> 		TAUHGH:fmissing<sub>value</sub> = 1.e+15f ;
> 		TAUHGH:scale<sub>factor</sub> = 1.f ;
> 		TAUHGH:add<sub>offset</sub> = 0.f ;
> 		TAUHGH:standard<sub>name</sub> = "in<sub>cloud</sub><sub>optical</sub><sub>thickness</sub><sub>of</sub><sub>high</sub><sub>clouds</sub>(EXPORT)" ;
> 		TAUHGH:vmax = 1.e+15f ;
> 		TAUHGH:vmin = -1.e+15f ;
> 		TAUHGH:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float TAULOW(time, lat, lon) ;
> 		TAULOW:long<sub>name</sub> = "in<sub>cloud</sub><sub>optical</sub><sub>thickness</sub><sub>of</sub><sub>low</sub><sub>clouds</sub>" ;
> 		TAULOW:units = "1" ;
> 		TAULOW:<sub>FillValue</sub> = 1.e+15f ;
> 		TAULOW:missing<sub>value</sub> = 1.e+15f ;
> 		TAULOW:fmissing<sub>value</sub> = 1.e+15f ;
> 		TAULOW:scale<sub>factor</sub> = 1.f ;
> 		TAULOW:add<sub>offset</sub> = 0.f ;
> 		TAULOW:standard<sub>name</sub> = "in<sub>cloud</sub><sub>optical</sub><sub>thickness</sub><sub>of</sub><sub>low</sub><sub>clouds</sub>" ;
> 		TAULOW:vmax = 1.e+15f ;
> 		TAULOW:vmin = -1.e+15f ;
> 		TAULOW:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float TAUMID(time, lat, lon) ;
> 		TAUMID:long<sub>name</sub> = "in<sub>cloud</sub><sub>optical</sub><sub>thickness</sub><sub>of</sub><sub>middle</sub><sub>clouds</sub>" ;
> 		TAUMID:units = "1" ;
> 		TAUMID:<sub>FillValue</sub> = 1.e+15f ;
> 		TAUMID:missing<sub>value</sub> = 1.e+15f ;
> 		TAUMID:fmissing<sub>value</sub> = 1.e+15f ;
> 		TAUMID:scale<sub>factor</sub> = 1.f ;
> 		TAUMID:add<sub>offset</sub> = 0.f ;
> 		TAUMID:standard<sub>name</sub> = "in<sub>cloud</sub><sub>optical</sub><sub>thickness</sub><sub>of</sub><sub>middle</sub><sub>clouds</sub>" ;
> 		TAUMID:vmax = 1.e+15f ;
> 		TAUMID:vmin = -1.e+15f ;
> 		TAUMID:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float TAUTOT(time, lat, lon) ;
> 		TAUTOT:long<sub>name</sub> = "in<sub>cloud</sub><sub>optical</sub><sub>thickness</sub><sub>of</sub><sub>all</sub><sub>clouds</sub>" ;
> 		TAUTOT:units = "1" ;
> 		TAUTOT:<sub>FillValue</sub> = 1.e+15f ;
> 		TAUTOT:missing<sub>value</sub> = 1.e+15f ;
> 		TAUTOT:fmissing<sub>value</sub> = 1.e+15f ;
> 		TAUTOT:scale<sub>factor</sub> = 1.f ;
> 		TAUTOT:add<sub>offset</sub> = 0.f ;
> 		TAUTOT:standard<sub>name</sub> = "in<sub>cloud</sub><sub>optical</sub><sub>thickness</sub><sub>of</sub><sub>all</sub><sub>clouds</sub>" ;
> 		TAUTOT:vmax = 1.e+15f ;
> 		TAUTOT:vmin = -1.e+15f ;
> 		TAUTOT:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 	float TS(time, lat, lon) ;
> 		TS:long<sub>name</sub> = "surface<sub>skin</sub><sub>temperature</sub>" ;
> 		TS:units = "K" ;
> 		TS:<sub>FillValue</sub> = 1.e+15f ;
> 		TS:missing<sub>value</sub> = 1.e+15f ;
> 		TS:fmissing<sub>value</sub> = 1.e+15f ;
> 		TS:scale<sub>factor</sub> = 1.f ;
> 		TS:add<sub>offset</sub> = 0.f ;
> 		TS:standard<sub>name</sub> = "surface<sub>skin</sub><sub>temperature</sub>" ;
> 		TS:vmax = 1.e+15f ;
> 		TS:vmin = -1.e+15f ;
> 		TS:valid<sub>range</sub> = -1.e+15f, 1.e+15f ;
> 
> // global attributes:
> 		:History = "Original file generated: Mon Mar 23 03:29:51 2015 GMT" ;
> 		:Comment = "GMAO filename: d5124<sub>m2</sub><sub>jan00.tavg1</sub><sub>2d</sub><sub>rad</sub><sub>Nx.20100101.nc4</sub>" ;
> 		:Filename = "MERRA2<sub>300.tavg1</sub><sub>2d</sub><sub>rad</sub><sub>Nx.20100101.nc4</sub>" ;
> 		:Conventions = "CF-1" ;
> 		:Institution = "NASA Global Modeling and Assimilation Office" ;
> 		:References = "<http://gmao.gsfc.nasa.gov>" ;
> 		:Format = "NetCDF-4/HDF-5" ;
> 		:SpatialCoverage = "global" ;
> 		:VersionID = "5.12.4" ;
> 		:TemporalRange = "1980-01-01 -> 2016-12-31" ;
> 		:identifier<sub>product</sub><sub>doi</sub><sub>authority</sub> = "<http://dx.doi.org/>" ;
> 		:ShortName = "M2T1NXRAD" ;
> 		:GranuleID = "MERRA2<sub>300.tavg1</sub><sub>2d</sub><sub>rad</sub><sub>Nx.20100101.nc4</sub>" ;
> 		:ProductionDateTime = "Original file generated: Mon Mar 23 03:29:51 2015 GMT" ;
> 		:LongName = "MERRA2 tavg1<sub>2d</sub><sub>rad</sub><sub>Nx</sub>: 2d,1-Hourly,Time-Averaged,Single-Level,Assimilation,Radiation Diagnostics" ;
> 		:Title = "MERRA2 tavg1<sub>2d</sub><sub>rad</sub><sub>Nx</sub>: 2d,1-Hourly,Time-Averaged,Single-Level,Assimilation,Radiation Diagnostics" ;
> 		:SouthernmostLatitude = "-90.0" ;
> 		:NorthernmostLatitude = "90.0" ;
> 		:WesternmostLongitude = "-180.0" ;
> 		:EasternmostLongitude = "179.375" ;
> 		:LatitudeResolution = "0.5" ;
> 		:LongitudeResolution = "0.625" ;
> 		:DataResolution = "0.5 x 0.625" ;
> 		:Source = "CVS tag: GEOSadas-5<sub>12</sub><sub>4</sub>" ;
> 		:Contact = "<http://gmao.gsfc.nasa.gov>" ;
> 		:identifier<sub>product</sub><sub>doi</sub> = "10.5067/Q9QMY5PBNV1T" ;
> 		:RangeBeginningDate = "2010-01-01" ;
> 		:RangeBeginningTime = "00:00:00.000000" ;
> 		:RangeEndingDate = "2010-01-01" ;
> 		:RangeEndingTime = "23:59:59.000000" ;
> }


<a id="org1d519b7"></a>

### Single-level Diagnostics ('slv')


<a id="orga65f4b0"></a>

## SOSE

The following files were downloaded from SOSE [Iteration 122](http://sose.ucsd.edu/BSOSE6_iter122_solution.html) (2013-2017).


<a id="org9145a77"></a>

### Bottom pressure

> netcdf bsose<sub>i122</sub><sub>2013to2017</sub><sub>1day</sub><sub>BottomPres</sub> {
> dimensions:
> 	time = 1826 ;
> 	XC = 2160 ;
> 	YC = 588 ;
> variables:
> 	int64 iter(time) ;
> 		iter:long<sub>name</sub> = "model timestep number" ;
> 		iter:standard<sub>name</sub> = "timestep" ;
> 	int64 time(time) ;
> 		time:long<sub>name</sub> = "Time" ;
> 		time:standard<sub>name</sub> = "time" ;
> 		time:axis = "T" ;
> 		time:units = "seconds since 2012-12-01" ;
> 		time:calendar = "proleptic<sub>gregorian</sub>" ;
> 	float XC(XC) ;
> 		XC:<sub>FillValue</sub> = NaNf ;
> 		XC:coordinate = "YC XC" ;
> 		XC:units = "degrees<sub>east</sub>" ;
> 		XC:standard<sub>name</sub> = "longitude" ;
> 		XC:long<sub>name</sub> = "longitude" ;
> 		XC:axis = "X" ;
> 	float YC(YC) ;
> 		YC:<sub>FillValue</sub> = NaNf ;
> 		YC:coordinate = "YC XC" ;
> 		YC:units = "degrees<sub>north</sub>" ;
> 		YC:standard<sub>name</sub> = "latitude" ;
> 		YC:long<sub>name</sub> = "latitude" ;
> 		YC:axis = "Y" ;
> 	float Depth(YC, XC) ;
> 		Depth:<sub>FillValue</sub> = NaNf ;
> 		Depth:coordinate = "XC YC" ;
> 		Depth:units = "m" ;
> 		Depth:standard<sub>name</sub> = "ocean<sub>depth</sub>" ;
> 		Depth:long<sub>name</sub> = "ocean depth" ;
> 	float rA(YC, XC) ;
> 		rA:<sub>FillValue</sub> = NaNf ;
> 		rA:coordinate = "YC XC" ;
> 		rA:units = "m2" ;
> 		rA:standard<sub>name</sub> = "cell<sub>area</sub>" ;
> 		rA:long<sub>name</sub> = "cell area" ;
> 	float PHIBOT(time, YC, XC) ;
> 		PHIBOT:<sub>FillValue</sub> = NaNf ;
> 		PHIBOT:units = "m<sup>2</sup>/s<sup>2</sup>" ;
> 		PHIBOT:long<sub>name</sub> = "Bottom Pressure Pot.(p/rho) Anomaly" ;
> 		PHIBOT:standard<sub>name</sub> = "PHIBOT" ;
> 		PHIBOT:coordinates = "Depth rA iter" ;
> }


<a id="org9e5ccb7"></a>

### Fresh water from atmosphere and land

> netcdf bsose<sub>i122</sub><sub>2013to2017</sub><sub>1day</sub><sub>FWfrmAtmLnd</sub> {
> dimensions:
> 	time = 1826 ;
> 	XC = 2160 ;
> 	YC = 588 ;
> variables:
> 	int64 iter(time) ;
> 		iter:long<sub>name</sub> = "model timestep number" ;
> 		iter:standard<sub>name</sub> = "timestep" ;
> 	int64 time(time) ;
> 		time:long<sub>name</sub> = "Time" ;
> 		time:standard<sub>name</sub> = "time" ;
> 		time:axis = "T" ;
> 		time:units = "seconds since 2012-12-01" ;
> 		time:calendar = "proleptic<sub>gregorian</sub>" ;
> 	float XC(XC) ;
> 		XC:<sub>FillValue</sub> = NaNf ;
> 		XC:coordinate = "YC XC" ;
> 		XC:units = "degrees<sub>east</sub>" ;
> 		XC:standard<sub>name</sub> = "longitude" ;
> 		XC:long<sub>name</sub> = "longitude" ;
> 		XC:axis = "X" ;
> 	float YC(YC) ;
> 		YC:<sub>FillValue</sub> = NaNf ;
> 		YC:coordinate = "YC XC" ;
> 		YC:units = "degrees<sub>north</sub>" ;
> 		YC:standard<sub>name</sub> = "latitude" ;
> 		YC:long<sub>name</sub> = "latitude" ;
> 		YC:axis = "Y" ;
> 	float Depth(YC, XC) ;
> 		Depth:<sub>FillValue</sub> = NaNf ;
> 		Depth:coordinate = "XC YC" ;
> 		Depth:units = "m" ;
> 		Depth:standard<sub>name</sub> = "ocean<sub>depth</sub>" ;
> 		Depth:long<sub>name</sub> = "ocean depth" ;
> 	float rA(YC, XC) ;
> 		rA:<sub>FillValue</sub> = NaNf ;
> 		rA:coordinate = "YC XC" ;
> 		rA:units = "m2" ;
> 		rA:standard<sub>name</sub> = "cell<sub>area</sub>" ;
> 		rA:long<sub>name</sub> = "cell area" ;
> 	float SIatmFW(time, YC, XC) ;
> 		SIatmFW:<sub>FillValue</sub> = NaNf ;
> 		SIatmFW:units = "kg/m<sup>2</sup>/s" ;
> 		SIatmFW:long<sub>name</sub> = "Net freshwater flux from atmosphere & land (+=down)" ;
> 		SIatmFW:standard<sub>name</sub> = "SIatmFW" ;
> 		SIatmFW:coordinates = "Depth rA iter" ;
> }


<a id="org2a85212"></a>

### Mixed layer depth

> netcdf bsose<sub>i122</sub><sub>2013to2017</sub><sub>1day</sub><sub>MLD</sub> {
> dimensions:
> 	time = 1826 ;
> 	XC = 2160 ;
> 	YC = 588 ;
> variables:
> 	int64 iter(time) ;
> 		iter:long<sub>name</sub> = "model timestep number" ;
> 		iter:standard<sub>name</sub> = "timestep" ;
> 	int64 time(time) ;
> 		time:long<sub>name</sub> = "Time" ;
> 		time:standard<sub>name</sub> = "time" ;
> 		time:axis = "T" ;
> 		time:units = "seconds since 2012-12-01" ;
> 		time:calendar = "proleptic<sub>gregorian</sub>" ;
> 	float XC(XC) ;
> 		XC:<sub>FillValue</sub> = NaNf ;
> 		XC:coordinate = "YC XC" ;
> 		XC:units = "degrees<sub>east</sub>" ;
> 		XC:standard<sub>name</sub> = "longitude" ;
> 		XC:long<sub>name</sub> = "longitude" ;
> 		XC:axis = "X" ;
> 	float YC(YC) ;
> 		YC:<sub>FillValue</sub> = NaNf ;
> 		YC:coordinate = "YC XC" ;
> 		YC:units = "degrees<sub>north</sub>" ;
> 		YC:standard<sub>name</sub> = "latitude" ;
> 		YC:long<sub>name</sub> = "latitude" ;
> 		YC:axis = "Y" ;
> 	float Depth(YC, XC) ;
> 		Depth:<sub>FillValue</sub> = NaNf ;
> 		Depth:coordinate = "XC YC" ;
> 		Depth:units = "m" ;
> 		Depth:standard<sub>name</sub> = "ocean<sub>depth</sub>" ;
> 		Depth:long<sub>name</sub> = "ocean depth" ;
> 	float rA(YC, XC) ;
> 		rA:<sub>FillValue</sub> = NaNf ;
> 		rA:coordinate = "YC XC" ;
> 		rA:units = "m2" ;
> 		rA:standard<sub>name</sub> = "cell<sub>area</sub>" ;
> 		rA:long<sub>name</sub> = "cell area" ;
> 	float BLGMLD(time, YC, XC) ;
> 		BLGMLD:<sub>FillValue</sub> = NaNf ;
> 		BLGMLD:units = "m" ;
> 		BLGMLD:long<sub>name</sub> = "Diagnosed mixed layer depth" ;
> 		BLGMLD:standard<sub>name</sub> = "BLGMLD" ;
> 		BLGMLD:coordinates = "Depth rA iter" ;
> }


<a id="org4281540"></a>

### Sea surface height anomaly

> netcdf bsose<sub>i122</sub><sub>2013to2017</sub><sub>1day</sub><sub>SSH</sub> {
> dimensions:
> 	time = 1826 ;
> 	XC = 2160 ;
> 	YC = 588 ;
> variables:
> 	int64 iter(time) ;
> 		iter:long<sub>name</sub> = "model timestep number" ;
> 		iter:standard<sub>name</sub> = "timestep" ;
> 	int64 time(time) ;
> 		time:long<sub>name</sub> = "Time" ;
> 		time:standard<sub>name</sub> = "time" ;
> 		time:axis = "T" ;
> 		time:units = "seconds since 2012-12-01" ;
> 		time:calendar = "proleptic<sub>gregorian</sub>" ;
> 	float XC(XC) ;
> 		XC:<sub>FillValue</sub> = NaNf ;
> 		XC:coordinate = "YC XC" ;
> 		XC:units = "degrees<sub>east</sub>" ;
> 		XC:standard<sub>name</sub> = "longitude" ;
> 		XC:long<sub>name</sub> = "longitude" ;
> 		XC:axis = "X" ;
> 	float YC(YC) ;
> 		YC:<sub>FillValue</sub> = NaNf ;
> 		YC:coordinate = "YC XC" ;
> 		YC:units = "degrees<sub>north</sub>" ;
> 		YC:standard<sub>name</sub> = "latitude" ;
> 		YC:long<sub>name</sub> = "latitude" ;
> 		YC:axis = "Y" ;
> 	float Depth(YC, XC) ;
> 		Depth:<sub>FillValue</sub> = NaNf ;
> 		Depth:coordinate = "XC YC" ;
> 		Depth:units = "m" ;
> 		Depth:standard<sub>name</sub> = "ocean<sub>depth</sub>" ;
> 		Depth:long<sub>name</sub> = "ocean depth" ;
> 	float rA(YC, XC) ;
> 		rA:<sub>FillValue</sub> = NaNf ;
> 		rA:coordinate = "YC XC" ;
> 		rA:units = "m2" ;
> 		rA:standard<sub>name</sub> = "cell<sub>area</sub>" ;
> 		rA:long<sub>name</sub> = "cell area" ;
> 	float ETAN(time, YC, XC) ;
> 		ETAN:<sub>FillValue</sub> = NaNf ;
> 		ETAN:units = "m" ;
> 		ETAN:long<sub>name</sub> = "Surface Height Anomaly" ;
> 		ETAN:standard<sub>name</sub> = "ETAN" ;
> 		ETAN:coordinates = "Depth rA iter" ;
> }


<a id="org5164d2b"></a>

### Sea ice area

> netcdf bsose<sub>i122</sub><sub>2013to2017</sub><sub>1day</sub><sub>SeaIceArea</sub> {
> dimensions:
> 	time = 1826 ;
> 	XC = 2160 ;
> 	YC = 588 ;
> variables:
> 	int64 iter(time) ;
> 		iter:long<sub>name</sub> = "model timestep number" ;
> 		iter:standard<sub>name</sub> = "timestep" ;
> 	int64 time(time) ;
> 		time:long<sub>name</sub> = "Time" ;
> 		time:standard<sub>name</sub> = "time" ;
> 		time:axis = "T" ;
> 		time:units = "seconds since 2012-12-01" ;
> 		time:calendar = "proleptic<sub>gregorian</sub>" ;
> 	float XC(XC) ;
> 		XC:<sub>FillValue</sub> = NaNf ;
> 		XC:coordinate = "YC XC" ;
> 		XC:units = "degrees<sub>east</sub>" ;
> 		XC:standard<sub>name</sub> = "longitude" ;
> 		XC:long<sub>name</sub> = "longitude" ;
> 		XC:axis = "X" ;
> 	float YC(YC) ;
> 		YC:<sub>FillValue</sub> = NaNf ;
> 		YC:coordinate = "YC XC" ;
> 		YC:units = "degrees<sub>north</sub>" ;
> 		YC:standard<sub>name</sub> = "latitude" ;
> 		YC:long<sub>name</sub> = "latitude" ;
> 		YC:axis = "Y" ;
> 	float Depth(YC, XC) ;
> 		Depth:<sub>FillValue</sub> = NaNf ;
> 		Depth:coordinate = "XC YC" ;
> 		Depth:units = "m" ;
> 		Depth:standard<sub>name</sub> = "ocean<sub>depth</sub>" ;
> 		Depth:long<sub>name</sub> = "ocean depth" ;
> 	float rA(YC, XC) ;
> 		rA:<sub>FillValue</sub> = NaNf ;
> 		rA:coordinate = "YC XC" ;
> 		rA:units = "m2" ;
> 		rA:standard<sub>name</sub> = "cell<sub>area</sub>" ;
> 		rA:long<sub>name</sub> = "cell area" ;
> 	float SIarea(time, YC, XC) ;
> 		SIarea:<sub>FillValue</sub> = NaNf ;
> 		SIarea:units = "m<sup>2</sup>/m<sup>2</sup>" ;
> 		SIarea:long<sub>name</sub> = "SEAICE fractional ice-covered area [0 to 1]" ;
> 		SIarea:standard<sub>name</sub> = "SIarea" ;
> 		SIarea:coordinates = "Depth rA iter" ;
> }


<a id="org40cb11a"></a>

### Sea ice covered part of SIqnet

> netcdf bsose<sub>i122</sub><sub>2013to2017</sub><sub>1day</sub><sub>SeaIceCvrdQnet</sub> {
> dimensions:
> 	time = 1826 ;
> 	XC = 2160 ;
> 	YC = 588 ;
> variables:
> 	int64 iter(time) ;
> 		iter:long<sub>name</sub> = "model timestep number" ;
> 		iter:standard<sub>name</sub> = "timestep" ;
> 	int64 time(time) ;
> 		time:long<sub>name</sub> = "Time" ;
> 		time:standard<sub>name</sub> = "time" ;
> 		time:axis = "T" ;
> 		time:units = "seconds since 2012-12-01" ;
> 		time:calendar = "proleptic<sub>gregorian</sub>" ;
> 	float XC(XC) ;
> 		XC:<sub>FillValue</sub> = NaNf ;
> 		XC:coordinate = "YC XC" ;
> 		XC:units = "degrees<sub>east</sub>" ;
> 		XC:standard<sub>name</sub> = "longitude" ;
> 		XC:long<sub>name</sub> = "longitude" ;
> 		XC:axis = "X" ;
> 	float YC(YC) ;
> 		YC:<sub>FillValue</sub> = NaNf ;
> 		YC:coordinate = "YC XC" ;
> 		YC:units = "degrees<sub>north</sub>" ;
> 		YC:standard<sub>name</sub> = "latitude" ;
> 		YC:long<sub>name</sub> = "latitude" ;
> 		YC:axis = "Y" ;
> 	float Depth(YC, XC) ;
> 		Depth:<sub>FillValue</sub> = NaNf ;
> 		Depth:coordinate = "XC YC" ;
> 		Depth:units = "m" ;
> 		Depth:standard<sub>name</sub> = "ocean<sub>depth</sub>" ;
> 		Depth:long<sub>name</sub> = "ocean depth" ;
> 	float rA(YC, XC) ;
> 		rA:<sub>FillValue</sub> = NaNf ;
> 		rA:coordinate = "YC XC" ;
> 		rA:units = "m2" ;
> 		rA:standard<sub>name</sub> = "cell<sub>area</sub>" ;
> 		rA:long<sub>name</sub> = "cell area" ;
> 	float SIqneti(time, YC, XC) ;
> 		SIqneti:<sub>FillValue</sub> = NaNf ;
> 		SIqneti:units = "W/m<sup>2</sup>" ;
> 		SIqneti:long<sub>name</sub> = "Ice Covered Part of SIqnet, turb+rad, >0 decr theta" ;
> 		SIqneti:standard<sub>name</sub> = "SIqneti" ;
> 		SIqneti:coordinates = "Depth rA iter" ;
> }


<a id="org89f3cbe"></a>

### Ocean surface freshwater flux

> netcdf bsose<sub>i122</sub><sub>2013to2017</sub><sub>1day</sub><sub>SeaIceEmPmR</sub> {
> dimensions:
> 	time = 1826 ;
> 	XC = 2160 ;
> 	YC = 588 ;
> variables:
> 	int64 iter(time) ;
> 		iter:long<sub>name</sub> = "model timestep number" ;
> 		iter:standard<sub>name</sub> = "timestep" ;
> 	int64 time(time) ;
> 		time:long<sub>name</sub> = "Time" ;
> 		time:standard<sub>name</sub> = "time" ;
> 		time:axis = "T" ;
> 		time:units = "seconds since 2012-12-01" ;
> 		time:calendar = "proleptic<sub>gregorian</sub>" ;
> 	float XC(XC) ;
> 		XC:<sub>FillValue</sub> = NaNf ;
> 		XC:coordinate = "YC XC" ;
> 		XC:units = "degrees<sub>east</sub>" ;
> 		XC:standard<sub>name</sub> = "longitude" ;
> 		XC:long<sub>name</sub> = "longitude" ;
> 		XC:axis = "X" ;
> 	float YC(YC) ;
> 		YC:<sub>FillValue</sub> = NaNf ;
> 		YC:coordinate = "YC XC" ;
> 		YC:units = "degrees<sub>north</sub>" ;
> 		YC:standard<sub>name</sub> = "latitude" ;
> 		YC:long<sub>name</sub> = "latitude" ;
> 		YC:axis = "Y" ;
> 	float Depth(YC, XC) ;
> 		Depth:<sub>FillValue</sub> = NaNf ;
> 		Depth:coordinate = "XC YC" ;
> 		Depth:units = "m" ;
> 		Depth:standard<sub>name</sub> = "ocean<sub>depth</sub>" ;
> 		Depth:long<sub>name</sub> = "ocean depth" ;
> 	float rA(YC, XC) ;
> 		rA:<sub>FillValue</sub> = NaNf ;
> 		rA:coordinate = "YC XC" ;
> 		rA:units = "m2" ;
> 		rA:standard<sub>name</sub> = "cell<sub>area</sub>" ;
> 		rA:long<sub>name</sub> = "cell area" ;
> 	float SIempmr(time, YC, XC) ;
> 		SIempmr:<sub>FillValue</sub> = NaNf ;
> 		SIempmr:units = "kg/m<sup>2</sup>/s" ;
> 		SIempmr:long<sub>name</sub> = "Ocean surface freshwater flux, > 0 increases salt" ;
> 		SIempmr:standard<sub>name</sub> = "SIempmr" ;
> 		SIempmr:coordinates = "Depth rA iter" ;
> }


<a id="org2126178"></a>

### Sea ice effective ice thickness

> netcdf bsose<sub>i122</sub><sub>2013to2017</sub><sub>1day</sub><sub>SeaIceHeff</sub> {
> dimensions:
> 	time = 1826 ;
> 	XC = 2160 ;
> 	YC = 588 ;
> variables:
> 	int64 iter(time) ;
> 		iter:long<sub>name</sub> = "model timestep number" ;
> 		iter:standard<sub>name</sub> = "timestep" ;
> 	int64 time(time) ;
> 		time:long<sub>name</sub> = "Time" ;
> 		time:standard<sub>name</sub> = "time" ;
> 		time:axis = "T" ;
> 		time:units = "seconds since 2012-12-01" ;
> 		time:calendar = "proleptic<sub>gregorian</sub>" ;
> 	float XC(XC) ;
> 		XC:<sub>FillValue</sub> = NaNf ;
> 		XC:coordinate = "YC XC" ;
> 		XC:units = "degrees<sub>east</sub>" ;
> 		XC:standard<sub>name</sub> = "longitude" ;
> 		XC:long<sub>name</sub> = "longitude" ;
> 		XC:axis = "X" ;
> 	float YC(YC) ;
> 		YC:<sub>FillValue</sub> = NaNf ;
> 		YC:coordinate = "YC XC" ;
> 		YC:units = "degrees<sub>north</sub>" ;
> 		YC:standard<sub>name</sub> = "latitude" ;
> 		YC:long<sub>name</sub> = "latitude" ;
> 		YC:axis = "Y" ;
> 	float Depth(YC, XC) ;
> 		Depth:<sub>FillValue</sub> = NaNf ;
> 		Depth:coordinate = "XC YC" ;
> 		Depth:units = "m" ;
> 		Depth:standard<sub>name</sub> = "ocean<sub>depth</sub>" ;
> 		Depth:long<sub>name</sub> = "ocean depth" ;
> 	float rA(YC, XC) ;
> 		rA:<sub>FillValue</sub> = NaNf ;
> 		rA:coordinate = "YC XC" ;
> 		rA:units = "m2" ;
> 		rA:standard<sub>name</sub> = "cell<sub>area</sub>" ;
> 		rA:long<sub>name</sub> = "cell area" ;
> 	float SIheff(time, YC, XC) ;
> 		SIheff:<sub>FillValue</sub> = NaNf ;
> 		SIheff:units = "m" ;
> 		SIheff:long<sub>name</sub> = "SEAICE effective ice thickness" ;
> 		SIheff:standard<sub>name</sub> = "SIheff" ;
> 		SIheff:coordinates = "Depth rA iter" ;
> }


<a id="org3fd9cb6"></a>

### Sea ice surface temperature over sea ice

> netcdf bsose<sub>i122</sub><sub>2013to2017</sub><sub>1day</sub><sub>SeaIceSurfT</sub> {
> dimensions:
> 	time = 1826 ;
> 	XC = 2160 ;
> 	YC = 588 ;
> variables:
> 	int64 iter(time) ;
> 		iter:long<sub>name</sub> = "model timestep number" ;
> 		iter:standard<sub>name</sub> = "timestep" ;
> 	int64 time(time) ;
> 		time:long<sub>name</sub> = "Time" ;
> 		time:standard<sub>name</sub> = "time" ;
> 		time:axis = "T" ;
> 		time:units = "seconds since 2012-12-01" ;
> 		time:calendar = "proleptic<sub>gregorian</sub>" ;
> 	float XC(XC) ;
> 		XC:<sub>FillValue</sub> = NaNf ;
> 		XC:coordinate = "YC XC" ;
> 		XC:units = "degrees<sub>east</sub>" ;
> 		XC:standard<sub>name</sub> = "longitude" ;
> 		XC:long<sub>name</sub> = "longitude" ;
> 		XC:axis = "X" ;
> 	float YC(YC) ;
> 		YC:<sub>FillValue</sub> = NaNf ;
> 		YC:coordinate = "YC XC" ;
> 		YC:units = "degrees<sub>north</sub>" ;
> 		YC:standard<sub>name</sub> = "latitude" ;
> 		YC:long<sub>name</sub> = "latitude" ;
> 		YC:axis = "Y" ;
> 	float Depth(YC, XC) ;
> 		Depth:<sub>FillValue</sub> = NaNf ;
> 		Depth:coordinate = "XC YC" ;
> 		Depth:units = "m" ;
> 		Depth:standard<sub>name</sub> = "ocean<sub>depth</sub>" ;
> 		Depth:long<sub>name</sub> = "ocean depth" ;
> 	float rA(YC, XC) ;
> 		rA:<sub>FillValue</sub> = NaNf ;
> 		rA:coordinate = "YC XC" ;
> 		rA:units = "m2" ;
> 		rA:standard<sub>name</sub> = "cell<sub>area</sub>" ;
> 		rA:long<sub>name</sub> = "cell area" ;
> 	float SItices(time, YC, XC) ;
> 		SItices:<sub>FillValue</sub> = NaNf ;
> 		SItices:units = "K" ;
> 		SItices:long<sub>name</sub> = "Surface Temperature over Sea-Ice (area weighted)" ;
> 		SItices:standard<sub>name</sub> = "SItices" ;
> 		SItices:mate = "SIarea" ;
> 		SItices:coordinates = "Depth rA iter" ;
> }


<a id="orgec146f6"></a>

### Sea ice zonal transport effective thickness

> netcdf bsose<sub>i122</sub><sub>2013to2017</sub><sub>1day</sub><sub>SeaIceUHeff</sub> {
> dimensions:
> 	time = 1826 ;
> 	YC = 588 ;
> 	XG = 2160 ;
> variables:
> 	int64 iter(time) ;
> 		iter:long<sub>name</sub> = "model timestep number" ;
> 		iter:standard<sub>name</sub> = "timestep" ;
> 	int64 time(time) ;
> 		time:long<sub>name</sub> = "Time" ;
> 		time:standard<sub>name</sub> = "time" ;
> 		time:axis = "T" ;
> 		time:units = "seconds since 2012-12-01" ;
> 		time:calendar = "proleptic<sub>gregorian</sub>" ;
> 	float YC(YC) ;
> 		YC:<sub>FillValue</sub> = NaNf ;
> 		YC:coordinate = "YC XC" ;
> 		YC:units = "degrees<sub>north</sub>" ;
> 		YC:standard<sub>name</sub> = "latitude" ;
> 		YC:long<sub>name</sub> = "latitude" ;
> 		YC:axis = "Y" ;
> 	float XG(XG) ;
> 		XG:<sub>FillValue</sub> = NaNf ;
> 		XG:coordinate = "YG XG" ;
> 		XG:units = "degrees<sub>east</sub>" ;
> 		XG:standard<sub>name</sub> = "longitude<sub>at</sub><sub>f</sub><sub>location</sub>" ;
> 		XG:long<sub>name</sub> = "longitude" ;
> 		XG:axis = "X" ;
> 		XG:c<sub>grid</sub><sub>axis</sub><sub>shift</sub> = -0.5 ;
> 	float dxC(YC, XG) ;
> 		dxC:<sub>FillValue</sub> = NaNf ;
> 		dxC:coordinate = "YC XG" ;
> 		dxC:units = "m" ;
> 		dxC:standard<sub>name</sub> = "cell<sub>x</sub><sub>size</sub><sub>at</sub><sub>u</sub><sub>location</sub>" ;
> 		dxC:long<sub>name</sub> = "cell x size" ;
> 	float rAw(YC, XG) ;
> 		rAw:<sub>FillValue</sub> = NaNf ;
> 		rAw:coordinate = "YG XC" ;
> 		rAw:units = "m2" ;
> 		rAw:standard<sub>name</sub> = "cell<sub>area</sub><sub>at</sub><sub>u</sub><sub>location</sub>" ;
> 		rAw:long<sub>name</sub> = "cell area" ;
> 	float dyG(YC, XG) ;
> 		dyG:<sub>FillValue</sub> = NaNf ;
> 		dyG:coordinate = "YC XG" ;
> 		dyG:units = "m" ;
> 		dyG:standard<sub>name</sub> = "cell<sub>y</sub><sub>size</sub><sub>at</sub><sub>u</sub><sub>location</sub>" ;
> 		dyG:long<sub>name</sub> = "cell y size" ;
> 	float SIuheff(time, YC, XG) ;
> 		SIuheff:<sub>FillValue</sub> = NaNf ;
> 		SIuheff:units = "m<sup>2</sup>/s" ;
> 		SIuheff:long<sub>name</sub> = "Zonal      Transport of eff ice thickn (centered)" ;
> 		SIuheff:standard<sub>name</sub> = "SIuheff" ;
> 		SIuheff:mate = "SIvheff" ;
> 		SIuheff:coordinates = "rAw dyG dxC iter" ;
> }


<a id="orgf8a0379"></a>

### Sea ice zonal velocity

> netcdf bsose<sub>i122</sub><sub>2013to2017</sub><sub>1day</sub><sub>SeaIceUvel</sub> {
> dimensions:
> 	time = 1826 ;
> 	YC = 588 ;
> 	XG = 2160 ;
> variables:
> 	int64 iter(time) ;
> 		iter:long<sub>name</sub> = "model timestep number" ;
> 		iter:standard<sub>name</sub> = "timestep" ;
> 	int64 time(time) ;
> 		time:long<sub>name</sub> = "Time" ;
> 		time:standard<sub>name</sub> = "time" ;
> 		time:axis = "T" ;
> 		time:units = "seconds since 2012-12-01" ;
> 		time:calendar = "proleptic<sub>gregorian</sub>" ;
> 	float YC(YC) ;
> 		YC:<sub>FillValue</sub> = NaNf ;
> 		YC:coordinate = "YC XC" ;
> 		YC:units = "degrees<sub>north</sub>" ;
> 		YC:standard<sub>name</sub> = "latitude" ;
> 		YC:long<sub>name</sub> = "latitude" ;
> 		YC:axis = "Y" ;
> 	float XG(XG) ;
> 		XG:<sub>FillValue</sub> = NaNf ;
> 		XG:coordinate = "YG XG" ;
> 		XG:units = "degrees<sub>east</sub>" ;
> 		XG:standard<sub>name</sub> = "longitude<sub>at</sub><sub>f</sub><sub>location</sub>" ;
> 		XG:long<sub>name</sub> = "longitude" ;
> 		XG:axis = "X" ;
> 		XG:c<sub>grid</sub><sub>axis</sub><sub>shift</sub> = -0.5 ;
> 	float dxC(YC, XG) ;
> 		dxC:<sub>FillValue</sub> = NaNf ;
> 		dxC:coordinate = "YC XG" ;
> 		dxC:units = "m" ;
> 		dxC:standard<sub>name</sub> = "cell<sub>x</sub><sub>size</sub><sub>at</sub><sub>u</sub><sub>location</sub>" ;
> 		dxC:long<sub>name</sub> = "cell x size" ;
> 	float rAw(YC, XG) ;
> 		rAw:<sub>FillValue</sub> = NaNf ;
> 		rAw:coordinate = "YG XC" ;
> 		rAw:units = "m2" ;
> 		rAw:standard<sub>name</sub> = "cell<sub>area</sub><sub>at</sub><sub>u</sub><sub>location</sub>" ;
> 		rAw:long<sub>name</sub> = "cell area" ;
> 	float dyG(YC, XG) ;
> 		dyG:<sub>FillValue</sub> = NaNf ;
> 		dyG:coordinate = "YC XG" ;
> 		dyG:units = "m" ;
> 		dyG:standard<sub>name</sub> = "cell<sub>y</sub><sub>size</sub><sub>at</sub><sub>u</sub><sub>location</sub>" ;
> 		dyG:long<sub>name</sub> = "cell y size" ;
> 	float SIuice(time, YC, XG) ;
> 		SIuice:<sub>FillValue</sub> = NaNf ;
> 		SIuice:units = "m/s" ;
> 		SIuice:long<sub>name</sub> = "SEAICE zonal ice velocity, >0 from West to East" ;
> 		SIuice:standard<sub>name</sub> = "SIuice" ;
> 		SIuice:mate = "SIvice" ;
> 		SIuice:coordinates = "rAw dyG dxC iter" ;
> }


<a id="orgedaedf7"></a>

### Sea ice meridional transport effective thickness

> netcdf bsose<sub>i122</sub><sub>2013to2017</sub><sub>1day</sub><sub>SeaIceVHeff</sub> {
> dimensions:
> 	time = 1826 ;
> 	XC = 2160 ;
> 	YG = 588 ;
> variables:
> 	int64 iter(time) ;
> 		iter:long<sub>name</sub> = "model timestep number" ;
> 		iter:standard<sub>name</sub> = "timestep" ;
> 	int64 time(time) ;
> 		time:long<sub>name</sub> = "Time" ;
> 		time:standard<sub>name</sub> = "time" ;
> 		time:axis = "T" ;
> 		time:units = "seconds since 2012-12-01" ;
> 		time:calendar = "proleptic<sub>gregorian</sub>" ;
> 	float XC(XC) ;
> 		XC:<sub>FillValue</sub> = NaNf ;
> 		XC:coordinate = "YC XC" ;
> 		XC:units = "degrees<sub>east</sub>" ;
> 		XC:standard<sub>name</sub> = "longitude" ;
> 		XC:long<sub>name</sub> = "longitude" ;
> 		XC:axis = "X" ;
> 	float YG(YG) ;
> 		YG:<sub>FillValue</sub> = NaNf ;
> 		YG:units = "degrees<sub>north</sub>" ;
> 		YG:long<sub>name</sub> = "latitude" ;
> 		YG:standard<sub>name</sub> = "latitude<sub>at</sub><sub>f</sub><sub>location</sub>" ;
> 		YG:axis = "Y" ;
> 		YG:c<sub>grid</sub><sub>axis</sub><sub>shift</sub> = -0.5 ;
> 	float rAs(YG, XC) ;
> 		rAs:<sub>FillValue</sub> = NaNf ;
> 		rAs:units = "m2" ;
> 		rAs:long<sub>name</sub> = "cell area" ;
> 		rAs:standard<sub>name</sub> = "cell<sub>area</sub><sub>at</sub><sub>v</sub><sub>location</sub>" ;
> 	float dxG(YG, XC) ;
> 		dxG:<sub>FillValue</sub> = NaNf ;
> 		dxG:coordinate = "YG XC" ;
> 		dxG:units = "m" ;
> 		dxG:standard<sub>name</sub> = "cell<sub>x</sub><sub>size</sub><sub>at</sub><sub>v</sub><sub>location</sub>" ;
> 		dxG:long<sub>name</sub> = "cell x size" ;
> 	float dyC(YG, XC) ;
> 		dyC:<sub>FillValue</sub> = NaNf ;
> 		dyC:coordinate = "YG XC" ;
> 		dyC:units = "m" ;
> 		dyC:standard<sub>name</sub> = "cell<sub>y</sub><sub>size</sub><sub>at</sub><sub>v</sub><sub>location</sub>" ;
> 		dyC:long<sub>name</sub> = "cell y size" ;
> 	float SIvheff(time, YG, XC) ;
> 		SIvheff:<sub>FillValue</sub> = NaNf ;
> 		SIvheff:units = "m<sup>2</sup>/s" ;
> 		SIvheff:long<sub>name</sub> = "Meridional Transport of eff ice thickn (centered)" ;
> 		SIvheff:standard<sub>name</sub> = "SIvheff" ;
> 		SIvheff:mate = "SIuheff" ;
> 		SIvheff:coordinates = "dxG rAs iter dyC" ;
> }


<a id="orgc0d5441"></a>

### Sea ice meridional velocity

> netcdf bsose<sub>i122</sub><sub>2013to2017</sub><sub>1day</sub><sub>SeaIceVvel</sub> {
> dimensions:
> 	time = 1826 ;
> 	XC = 2160 ;
> 	YG = 588 ;
> variables:
> 	int64 iter(time) ;
> 		iter:long<sub>name</sub> = "model timestep number" ;
> 		iter:standard<sub>name</sub> = "timestep" ;
> 	int64 time(time) ;
> 		time:long<sub>name</sub> = "Time" ;
> 		time:standard<sub>name</sub> = "time" ;
> 		time:axis = "T" ;
> 		time:units = "seconds since 2012-12-01" ;
> 		time:calendar = "proleptic<sub>gregorian</sub>" ;
> 	float XC(XC) ;
> 		XC:<sub>FillValue</sub> = NaNf ;
> 		XC:coordinate = "YC XC" ;
> 		XC:units = "degrees<sub>east</sub>" ;
> 		XC:standard<sub>name</sub> = "longitude" ;
> 		XC:long<sub>name</sub> = "longitude" ;
> 		XC:axis = "X" ;
> 	float YG(YG) ;
> 		YG:<sub>FillValue</sub> = NaNf ;
> 		YG:units = "degrees<sub>north</sub>" ;
> 		YG:long<sub>name</sub> = "latitude" ;
> 		YG:standard<sub>name</sub> = "latitude<sub>at</sub><sub>f</sub><sub>location</sub>" ;
> 		YG:axis = "Y" ;
> 		YG:c<sub>grid</sub><sub>axis</sub><sub>shift</sub> = -0.5 ;
> 	float rAs(YG, XC) ;
> 		rAs:<sub>FillValue</sub> = NaNf ;
> 		rAs:units = "m2" ;
> 		rAs:long<sub>name</sub> = "cell area" ;
> 		rAs:standard<sub>name</sub> = "cell<sub>area</sub><sub>at</sub><sub>v</sub><sub>location</sub>" ;
> 	float dxG(YG, XC) ;
> 		dxG:<sub>FillValue</sub> = NaNf ;
> 		dxG:coordinate = "YG XC" ;
> 		dxG:units = "m" ;
> 		dxG:standard<sub>name</sub> = "cell<sub>x</sub><sub>size</sub><sub>at</sub><sub>v</sub><sub>location</sub>" ;
> 		dxG:long<sub>name</sub> = "cell x size" ;
> 	float dyC(YG, XC) ;
> 		dyC:<sub>FillValue</sub> = NaNf ;
> 		dyC:coordinate = "YG XC" ;
> 		dyC:units = "m" ;
> 		dyC:standard<sub>name</sub> = "cell<sub>y</sub><sub>size</sub><sub>at</sub><sub>v</sub><sub>location</sub>" ;
> 		dyC:long<sub>name</sub> = "cell y size" ;
> 	float SIvice(time, YG, XC) ;
> 		SIvice:<sub>FillValue</sub> = NaNf ;
> 		SIvice:units = "m/s" ;
> 		SIvice:long<sub>name</sub> = "SEAICE merid. ice velocity, >0 from South to North" ;
> 		SIvice:standard<sub>name</sub> = "SIvice" ;
> 		SIvice:mate = "SIuice" ;
> 		SIvice:coordinates = "dxG rAs iter dyC" ;
> }


<a id="orgb39e948"></a>

### Sea ice effective snow thickness

> netcdf bsose<sub>i122</sub><sub>2013to2017</sub><sub>1day</sub><sub>SnowHeff</sub> {
> dimensions:
> 	time = 1826 ;
> 	XC = 2160 ;
> 	YC = 588 ;
> variables:
> 	int64 iter(time) ;
> 		iter:long<sub>name</sub> = "model timestep number" ;
> 		iter:standard<sub>name</sub> = "timestep" ;
> 	int64 time(time) ;
> 		time:long<sub>name</sub> = "Time" ;
> 		time:standard<sub>name</sub> = "time" ;
> 		time:axis = "T" ;
> 		time:units = "seconds since 2012-12-01" ;
> 		time:calendar = "proleptic<sub>gregorian</sub>" ;
> 	float XC(XC) ;
> 		XC:<sub>FillValue</sub> = NaNf ;
> 		XC:coordinate = "YC XC" ;
> 		XC:units = "degrees<sub>east</sub>" ;
> 		XC:standard<sub>name</sub> = "longitude" ;
> 		XC:long<sub>name</sub> = "longitude" ;
> 		XC:axis = "X" ;
> 	float YC(YC) ;
> 		YC:<sub>FillValue</sub> = NaNf ;
> 		YC:coordinate = "YC XC" ;
> 		YC:units = "degrees<sub>north</sub>" ;
> 		YC:standard<sub>name</sub> = "latitude" ;
> 		YC:long<sub>name</sub> = "latitude" ;
> 		YC:axis = "Y" ;
> 	float Depth(YC, XC) ;
> 		Depth:<sub>FillValue</sub> = NaNf ;
> 		Depth:coordinate = "XC YC" ;
> 		Depth:units = "m" ;
> 		Depth:standard<sub>name</sub> = "ocean<sub>depth</sub>" ;
> 		Depth:long<sub>name</sub> = "ocean depth" ;
> 	float rA(YC, XC) ;
> 		rA:<sub>FillValue</sub> = NaNf ;
> 		rA:coordinate = "YC XC" ;
> 		rA:units = "m2" ;
> 		rA:standard<sub>name</sub> = "cell<sub>area</sub>" ;
> 		rA:long<sub>name</sub> = "cell area" ;
> 	float SIhsnow(time, YC, XC) ;
> 		SIhsnow:<sub>FillValue</sub> = NaNf ;
> 		SIhsnow:units = "m" ;
> 		SIhsnow:long<sub>name</sub> = "SEAICE effective snow thickness" ;
> 		SIhsnow:standard<sub>name</sub> = "SIhsnow" ;
> 		SIhsnow:coordinates = "Depth rA iter" ;
> }


<a id="org98336b2"></a>

### Sea ice snow precipitation over sea ice

> netcdf bsose<sub>i122</sub><sub>2013to2017</sub><sub>1day</sub><sub>SnowPrecip</sub> {
> dimensions:
> 	time = 1826 ;
> 	XC = 2160 ;
> 	YC = 588 ;
> variables:
> 	int64 iter(time) ;
> 		iter:long<sub>name</sub> = "model timestep number" ;
> 		iter:standard<sub>name</sub> = "timestep" ;
> 	int64 time(time) ;
> 		time:long<sub>name</sub> = "Time" ;
> 		time:standard<sub>name</sub> = "time" ;
> 		time:axis = "T" ;
> 		time:units = "seconds since 2012-12-01" ;
> 		time:calendar = "proleptic<sub>gregorian</sub>" ;
> 	float XC(XC) ;
> 		XC:<sub>FillValue</sub> = NaNf ;
> 		XC:coordinate = "YC XC" ;
> 		XC:units = "degrees<sub>east</sub>" ;
> 		XC:standard<sub>name</sub> = "longitude" ;
> 		XC:long<sub>name</sub> = "longitude" ;
> 		XC:axis = "X" ;
> 	float YC(YC) ;
> 		YC:<sub>FillValue</sub> = NaNf ;
> 		YC:coordinate = "YC XC" ;
> 		YC:units = "degrees<sub>north</sub>" ;
> 		YC:standard<sub>name</sub> = "latitude" ;
> 		YC:long<sub>name</sub> = "latitude" ;
> 		YC:axis = "Y" ;
> 	float Depth(YC, XC) ;
> 		Depth:<sub>FillValue</sub> = NaNf ;
> 		Depth:coordinate = "XC YC" ;
> 		Depth:units = "m" ;
> 		Depth:standard<sub>name</sub> = "ocean<sub>depth</sub>" ;
> 		Depth:long<sub>name</sub> = "ocean depth" ;
> 	float rA(YC, XC) ;
> 		rA:<sub>FillValue</sub> = NaNf ;
> 		rA:coordinate = "YC XC" ;
> 		rA:units = "m2" ;
> 		rA:standard<sub>name</sub> = "cell<sub>area</sub>" ;
> 		rA:long<sub>name</sub> = "cell area" ;
> 	float SIsnPrcp(time, YC, XC) ;
> 		SIsnPrcp:<sub>FillValue</sub> = NaNf ;
> 		SIsnPrcp:units = "kg/m<sup>2</sup>/s" ;
> 		SIsnPrcp:long<sub>name</sub> = "Snow precip. (+=dw) over Sea-Ice (area weighted)" ;
> 		SIsnPrcp:standard<sub>name</sub> = "SIsnPrcp" ;
> 		SIsnPrcp:coordinates = "Depth rA iter" ;
> }


<a id="org6a6e72b"></a>

### Zonal wind stress at surface

> netcdf bsose<sub>i122</sub><sub>2013to2017</sub><sub>1day</sub><sub>oceTAUX</sub> {
> dimensions:
> 	time = 1826 ;
> 	YC = 588 ;
> 	XG = 2160 ;
> variables:
> 	int64 iter(time) ;
> 		iter:long<sub>name</sub> = "model timestep number" ;
> 		iter:standard<sub>name</sub> = "timestep" ;
> 	int64 time(time) ;
> 		time:long<sub>name</sub> = "Time" ;
> 		time:standard<sub>name</sub> = "time" ;
> 		time:axis = "T" ;
> 		time:units = "seconds since 2012-12-01" ;
> 		time:calendar = "proleptic<sub>gregorian</sub>" ;
> 	float YC(YC) ;
> 		YC:<sub>FillValue</sub> = NaNf ;
> 		YC:coordinate = "YC XC" ;
> 		YC:units = "degrees<sub>north</sub>" ;
> 		YC:standard<sub>name</sub> = "latitude" ;
> 		YC:long<sub>name</sub> = "latitude" ;
> 		YC:axis = "Y" ;
> 	float XG(XG) ;
> 		XG:<sub>FillValue</sub> = NaNf ;
> 		XG:coordinate = "YG XG" ;
> 		XG:units = "degrees<sub>east</sub>" ;
> 		XG:standard<sub>name</sub> = "longitude<sub>at</sub><sub>f</sub><sub>location</sub>" ;
> 		XG:long<sub>name</sub> = "longitude" ;
> 		XG:axis = "X" ;
> 		XG:c<sub>grid</sub><sub>axis</sub><sub>shift</sub> = -0.5 ;
> 	float dxC(YC, XG) ;
> 		dxC:<sub>FillValue</sub> = NaNf ;
> 		dxC:coordinate = "YC XG" ;
> 		dxC:units = "m" ;
> 		dxC:standard<sub>name</sub> = "cell<sub>x</sub><sub>size</sub><sub>at</sub><sub>u</sub><sub>location</sub>" ;
> 		dxC:long<sub>name</sub> = "cell x size" ;
> 	float rAw(YC, XG) ;
> 		rAw:<sub>FillValue</sub> = NaNf ;
> 		rAw:coordinate = "YG XC" ;
> 		rAw:units = "m2" ;
> 		rAw:standard<sub>name</sub> = "cell<sub>area</sub><sub>at</sub><sub>u</sub><sub>location</sub>" ;
> 		rAw:long<sub>name</sub> = "cell area" ;
> 	float dyG(YC, XG) ;
> 		dyG:<sub>FillValue</sub> = NaNf ;
> 		dyG:coordinate = "YC XG" ;
> 		dyG:units = "m" ;
> 		dyG:standard<sub>name</sub> = "cell<sub>y</sub><sub>size</sub><sub>at</sub><sub>u</sub><sub>location</sub>" ;
> 		dyG:long<sub>name</sub> = "cell y size" ;
> 	float oceTAUX(time, YC, XG) ;
> 		oceTAUX:<sub>FillValue</sub> = NaNf ;
> 		oceTAUX:units = "N/m<sup>2</sup>" ;
> 		oceTAUX:long<sub>name</sub> = "zonal surface wind stress, >0 increases uVel" ;
> 		oceTAUX:standard<sub>name</sub> = "oceTAUX" ;
> 		oceTAUX:mate = "oceTAUY" ;
> 		oceTAUX:coordinates = "rAw dyG dxC iter" ;
> }


<a id="orge50e9a6"></a>

### Meridional wind stress at surface

> netcdf bsose<sub>i122</sub><sub>2013to2017</sub><sub>1day</sub><sub>oceTAUY</sub> {
> dimensions:
> 	time = 1826 ;
> 	XC = 2160 ;
> 	YG = 588 ;
> variables:
> 	int64 iter(time) ;
> 		iter:long<sub>name</sub> = "model timestep number" ;
> 		iter:standard<sub>name</sub> = "timestep" ;
> 	int64 time(time) ;
> 		time:long<sub>name</sub> = "Time" ;
> 		time:standard<sub>name</sub> = "time" ;
> 		time:axis = "T" ;
> 		time:units = "seconds since 2012-12-01" ;
> 		time:calendar = "proleptic<sub>gregorian</sub>" ;
> 	float XC(XC) ;
> 		XC:<sub>FillValue</sub> = NaNf ;
> 		XC:coordinate = "YC XC" ;
> 		XC:units = "degrees<sub>east</sub>" ;
> 		XC:standard<sub>name</sub> = "longitude" ;
> 		XC:long<sub>name</sub> = "longitude" ;
> 		XC:axis = "X" ;
> 	float YG(YG) ;
> 		YG:<sub>FillValue</sub> = NaNf ;
> 		YG:units = "degrees<sub>north</sub>" ;
> 		YG:long<sub>name</sub> = "latitude" ;
> 		YG:standard<sub>name</sub> = "latitude<sub>at</sub><sub>f</sub><sub>location</sub>" ;
> 		YG:axis = "Y" ;
> 		YG:c<sub>grid</sub><sub>axis</sub><sub>shift</sub> = -0.5 ;
> 	float rAs(YG, XC) ;
> 		rAs:<sub>FillValue</sub> = NaNf ;
> 		rAs:units = "m2" ;
> 		rAs:long<sub>name</sub> = "cell area" ;
> 		rAs:standard<sub>name</sub> = "cell<sub>area</sub><sub>at</sub><sub>v</sub><sub>location</sub>" ;
> 	float dxG(YG, XC) ;
> 		dxG:<sub>FillValue</sub> = NaNf ;
> 		dxG:coordinate = "YG XC" ;
> 		dxG:units = "m" ;
> 		dxG:standard<sub>name</sub> = "cell<sub>x</sub><sub>size</sub><sub>at</sub><sub>v</sub><sub>location</sub>" ;
> 		dxG:long<sub>name</sub> = "cell x size" ;
> 	float dyC(YG, XC) ;
> 		dyC:<sub>FillValue</sub> = NaNf ;
> 		dyC:coordinate = "YG XC" ;
> 		dyC:units = "m" ;
> 		dyC:standard<sub>name</sub> = "cell<sub>y</sub><sub>size</sub><sub>at</sub><sub>v</sub><sub>location</sub>" ;
> 		dyC:long<sub>name</sub> = "cell y size" ;
> 	float oceTAUY(time, YG, XC) ;
> 		oceTAUY:<sub>FillValue</sub> = NaNf ;
> 		oceTAUY:units = "N/m<sup>2</sup>" ;
> 		oceTAUY:long<sub>name</sub> = "meridional surf. wind stress, >0 increases vVel" ;
> 		oceTAUY:standard<sub>name</sub> = "oceTAUY" ;
> 		oceTAUY:mate = "oceTAUX" ;
> 		oceTAUY:coordinates = "dxG rAs iter dyC" ;
> }


<a id="orgd114522"></a>

### Total salt flux

> netcdf bsose<sub>i122</sub><sub>2013to2017</sub><sub>1day</sub><sub>surfSflx</sub> {
> dimensions:
> 	time = 1826 ;
> 	XC = 2160 ;
> 	YC = 588 ;
> variables:
> 	int64 iter(time) ;
> 		iter:long<sub>name</sub> = "model timestep number" ;
> 		iter:standard<sub>name</sub> = "timestep" ;
> 	int64 time(time) ;
> 		time:long<sub>name</sub> = "Time" ;
> 		time:standard<sub>name</sub> = "time" ;
> 		time:axis = "T" ;
> 		time:units = "seconds since 2012-12-01" ;
> 		time:calendar = "proleptic<sub>gregorian</sub>" ;
> 	float XC(XC) ;
> 		XC:<sub>FillValue</sub> = NaNf ;
> 		XC:coordinate = "YC XC" ;
> 		XC:units = "degrees<sub>east</sub>" ;
> 		XC:standard<sub>name</sub> = "longitude" ;
> 		XC:long<sub>name</sub> = "longitude" ;
> 		XC:axis = "X" ;
> 	float YC(YC) ;
> 		YC:<sub>FillValue</sub> = NaNf ;
> 		YC:coordinate = "YC XC" ;
> 		YC:units = "degrees<sub>north</sub>" ;
> 		YC:standard<sub>name</sub> = "latitude" ;
> 		YC:long<sub>name</sub> = "latitude" ;
> 		YC:axis = "Y" ;
> 	float Depth(YC, XC) ;
> 		Depth:<sub>FillValue</sub> = NaNf ;
> 		Depth:coordinate = "XC YC" ;
> 		Depth:units = "m" ;
> 		Depth:standard<sub>name</sub> = "ocean<sub>depth</sub>" ;
> 		Depth:long<sub>name</sub> = "ocean depth" ;
> 	float rA(YC, XC) ;
> 		rA:<sub>FillValue</sub> = NaNf ;
> 		rA:coordinate = "YC XC" ;
> 		rA:units = "m2" ;
> 		rA:standard<sub>name</sub> = "cell<sub>area</sub>" ;
> 		rA:long<sub>name</sub> = "cell area" ;
> 	float SFLUX(time, YC, XC) ;
> 		SFLUX:<sub>FillValue</sub> = NaNf ;
> 		SFLUX:units = "g/m<sup>2</sup>/s" ;
> 		SFLUX:long<sub>name</sub> = "total salt flux (match salt-content variations), >0 increases salt" ;
> 		SFLUX:standard<sub>name</sub> = "SFLUX" ;
> 		SFLUX:coordinates = "Depth rA iter" ;
> }


<a id="orgcad33e4"></a>

### Total heat flux

> netcdf bsose<sub>i122</sub><sub>2013to2017</sub><sub>1day</sub><sub>surfTflx</sub> {
> dimensions:
> 	time = 1826 ;
> 	XC = 2160 ;
> 	YC = 588 ;
> variables:
> 	int64 iter(time) ;
> 		iter:long<sub>name</sub> = "model timestep number" ;
> 		iter:standard<sub>name</sub> = "timestep" ;
> 	int64 time(time) ;
> 		time:long<sub>name</sub> = "Time" ;
> 		time:standard<sub>name</sub> = "time" ;
> 		time:axis = "T" ;
> 		time:units = "seconds since 2012-12-01" ;
> 		time:calendar = "proleptic<sub>gregorian</sub>" ;
> 	float XC(XC) ;
> 		XC:<sub>FillValue</sub> = NaNf ;
> 		XC:coordinate = "YC XC" ;
> 		XC:units = "degrees<sub>east</sub>" ;
> 		XC:standard<sub>name</sub> = "longitude" ;
> 		XC:long<sub>name</sub> = "longitude" ;
> 		XC:axis = "X" ;
> 	float YC(YC) ;
> 		YC:<sub>FillValue</sub> = NaNf ;
> 		YC:coordinate = "YC XC" ;
> 		YC:units = "degrees<sub>north</sub>" ;
> 		YC:standard<sub>name</sub> = "latitude" ;
> 		YC:long<sub>name</sub> = "latitude" ;
> 		YC:axis = "Y" ;
> 	float Depth(YC, XC) ;
> 		Depth:<sub>FillValue</sub> = NaNf ;
> 		Depth:coordinate = "XC YC" ;
> 		Depth:units = "m" ;
> 		Depth:standard<sub>name</sub> = "ocean<sub>depth</sub>" ;
> 		Depth:long<sub>name</sub> = "ocean depth" ;
> 	float rA(YC, XC) ;
> 		rA:<sub>FillValue</sub> = NaNf ;
> 		rA:coordinate = "YC XC" ;
> 		rA:units = "m2" ;
> 		rA:standard<sub>name</sub> = "cell<sub>area</sub>" ;
> 		rA:long<sub>name</sub> = "cell area" ;
> 	float TFLUX(time, YC, XC) ;
> 		TFLUX:<sub>FillValue</sub> = NaNf ;
> 		TFLUX:units = "W/m<sup>2</sup>" ;
> 		TFLUX:long<sub>name</sub> = "total heat flux (match heat-content variations), >0 increases theta" ;
> 		TFLUX:standard<sub>name</sub> = "TFLUX" ;
> 		TFLUX:coordinates = "Depth rA iter" ;
> }


<a id="org68eff42"></a>

### Atmospheric surface pressure

> netcdf bsose<sub>i122</sub><sub>2013to2017</sub><sub>daily</sub><sub>EXFpress</sub> {
> dimensions:
> 	time = 1826 ;
> 	XC = 2160 ;
> 	YC = 588 ;
> variables:
> 	int64 iter(time) ;
> 		iter:long<sub>name</sub> = "model timestep number" ;
> 		iter:standard<sub>name</sub> = "timestep" ;
> 	int64 time(time) ;
> 		time:long<sub>name</sub> = "Time" ;
> 		time:standard<sub>name</sub> = "time" ;
> 		time:axis = "T" ;
> 		time:units = "seconds since 2012-12-01" ;
> 		time:calendar = "proleptic<sub>gregorian</sub>" ;
> 	float XC(XC) ;
> 		XC:<sub>FillValue</sub> = NaNf ;
> 		XC:coordinate = "YC XC" ;
> 		XC:units = "degrees<sub>east</sub>" ;
> 		XC:standard<sub>name</sub> = "longitude" ;
> 		XC:long<sub>name</sub> = "longitude" ;
> 		XC:axis = "X" ;
> 	float YC(YC) ;
> 		YC:<sub>FillValue</sub> = NaNf ;
> 		YC:coordinate = "YC XC" ;
> 		YC:units = "degrees<sub>north</sub>" ;
> 		YC:standard<sub>name</sub> = "latitude" ;
> 		YC:long<sub>name</sub> = "latitude" ;
> 		YC:axis = "Y" ;
> 	float Depth(YC, XC) ;
> 		Depth:<sub>FillValue</sub> = NaNf ;
> 		Depth:coordinate = "XC YC" ;
> 		Depth:units = "m" ;
> 		Depth:standard<sub>name</sub> = "ocean<sub>depth</sub>" ;
> 		Depth:long<sub>name</sub> = "ocean depth" ;
> 	float rA(YC, XC) ;
> 		rA:<sub>FillValue</sub> = NaNf ;
> 		rA:coordinate = "YC XC" ;
> 		rA:units = "m2" ;
> 		rA:standard<sub>name</sub> = "cell<sub>area</sub>" ;
> 		rA:long<sub>name</sub> = "cell area" ;
> 	float EXFpress(time, YC, XC) ;
> 		EXFpress:<sub>FillValue</sub> = NaNf ;
> 		EXFpress:units = "N/m<sup>2</sup>" ;
> 		EXFpress:long<sub>name</sub> = "atmospheric pressure field" ;
> 		EXFpress:standard<sub>name</sub> = "EXFpress" ;
> 		EXFpress:coordinates = "Depth rA iter" ;
> }


<a id="orgea1d349"></a>

## BRAN

[CSIRO maintained dataset](https://geonetwork.nci.org.au/geonetwork/srv/eng/catalog.search#/metadata/f9372_7752_2015_3718), but available on gadi via the [BRAN data](https://dapds00.nci.org.au/thredds/catalog/gb6/BRAN/BRAN2020/catalog.html) repository.


<a id="orgea7c781"></a>

## CAWCR

The [CICE GitHub online setup guide](https://github.com/CICE-Consortium/CICE/wiki/CICE-Input-Data#descriptions-of-tar-files) explains that ocean surface gravity wave
information is necessary for floe size distribution (FSD). Explanation is given
that illustrates from just one day with the only variable present from [WaveWatch
III,](https://github.com/NOAA-EMC/WW3) `efreq` ($m^2/s$) the importance of including wave information in sea ice
modelling. The [CICE Documentation](https://cice-consortium-cice.readthedocs.io/en/main/#) states the basic problem in sea ice modelling,
the one CICE aims to model, is ice thickness distribution (ITD). However, CICE
can optionally attempt to model FSD. It does this via a probability distribution
function that charcacterises the variability of the horizontal floe size through
vertical and lateral growth, new freezing, breaking waves and welding
of floes. This work is based on \citep:horvatPrognosticModelSeaice2015

-   Available to download from a [THREDDS](https://data-cbr.csiro.au/thredds/catalog/catch_all/CMAR_CAWCR-Wave_archive/CAWCR_Wave_Hindcast_aggregate/gridded/catalog.html)


<a id="org057ba5d"></a>

# Pre-processing Forcing and Initial Conditions


<a id="orgcfbc076"></a>

## Mimicing CICE's current options for forcing

In the [CICE6 documentation for the forcing namelist fields](https://cice-consortium-cice.readthedocs.io/en/main/user_guide/ug_case_settings.html#forcing-nml) the `atm_data_type` and
the `ocn_data_type` appear to be the only location from which to control the type
of atmosphere and ocean forcing datasets, respectively, that one can provide
CICE6. It's noted that neither [2.3](#org150c1e8) or [2.7](#orgea1d349) are offered as supported options.
Hence if I am to get CICE to be forced with [2.3](#org150c1e8) and [2.7](#orgea1d349) then I've got
essentially two options:

1.  write my own Fortran sub-routine that can be used directly by CICE6, or
2.  re-configure [2.3](#org150c1e8) and [2.7](#orgea1d349) to \`look' like an existing supported CICE6 forcing
    option.

Do to my lack of skill in writing Fortran and my initial thought that getting
[2.3](#org150c1e8) and [2.7](#orgea1d349) to \`look' like an existing supported CICE6 forcing option, I chose
the second option. This turned out to be less than trivial.


<a id="org1b51396"></a>

### Finding which dataset to mimic

Fortunately all the Fortran sub-routines that perform the reading-in of physical
forcing datasets are found in the [ice<sub>forcing.F90</sub>](file:///Users/dpath2o/src/CICE/cicecore/cicedynB/general/ice_forcing.F90) file. The sub-routine `JRA55_tx1_files`
performs the initial reading-in of JRA55 files. This seemed like the most
logical place to start mimicing the atmosphere. Hence using the [1.4.4](#org0eecc9b)
(from above) I began in earnest to re-grid two years worth of [2.3](#org150c1e8) dataset.
Similarily, for the ocean forcing, `ocn_data_ncar_init_3D` performs the reading-in
of 3D ocean data. The beginning comments in that sub-routine provide some useful
information and just a few lines below those comments the variable names are
given. Using this information I began in earnest to re-grid two years worth of
[2.7](#orgea1d349) dataset.


<a id="org8a1d213"></a>

### Creating Python module [afim.py](file:///Users/dpath2o/PHD/src/python/afim.py)

After beginning to do the above re-gridding/mimicing in a Python notebook it
became apparent that the resources used by iPython/Jupyter were too taxing and
hence a script-based (command-line) interface is best to utilise [afim.py](file:///Users/dpath2o/PHD/src/python/afim.py). At
present, the class `cice_prep` performs the re-gridding of [2.3](#org150c1e8) and [2.7](#orgea1d349) via a
required JSON, here's the [current example](file:///Users/dpath2o/PHD/src/python/afim.json). So performing a re-gridding is simply
done via the following commands, which are [here in this example](file:///Users/dpath2o/PHD/src/python/CICE/prepare_cice_forcing_script.py):

    import os
    import afim
    cp = afim.cice_prep('/Users/dpath2o/PHD/src/python/afim.json')
    # do ERA5 regridding and derivations for forcing
    cp.regrid_wrapper(ds_name='ERA5')
    # do BRAN regridding and derivations for forcing
    cp.regrid_wrapper(ds_name='BRAN')


<a id="org08accaa"></a>

## One-to-One Forcing Discusion

[As discussed earlier](#org28fbdde), CICE6, \`out-of-the-box' comes with three example
or test runs based on the underlying spatial grid: a [three-degree global
cartesian](#orge8e3b2d), a [one-degree global cartesian](#org05ab192), and a [one-degree global tri-polar](#org79cf8e6).
[The CICE6 documentation](https://cice-consortium-cice.readthedocs.io/en/main/user_guide/ug_case_settings.html#forcing-nml) indicates a very narrow range of supported forcing
datasets, however, the documentation does describe what variables are
nice-to-haves for running CICE6 in stand-alone. However, users are cautioned
that stand-alone is not the intended \`scientific' mode of CICE6, and that
coupling with an ocean and/or atmosphere is preferred for obtaining more
scientific results. From that [discussion in the documentation](https://cice-consortium-cice.readthedocs.io/en/main/developer_guide/dg_forcing.html), I created a table
for future use in mapping [2.3](#org150c1e8) and [2.7](#orgea1d349) data to in case it becomes
useful/beneficial to force CICE6 writing my own Fortran sub-routine.


<a id="orga0f92be"></a>

### Table of Stand-alone CICE6 Forcing Fields

<table border="2" cellspacing="0" cellpadding="6" rules="groups" frame="hsides">


<colgroup>
<col  class="org-left" />

<col  class="org-left" />

<col  class="org-left" />

<col  class="org-left" />

<col  class="org-left" />

<col  class="org-left" />

<col  class="org-left" />

<col  class="org-left" />
</colgroup>
<thead>
<tr>
<th scope="col" class="org-left">CICE6 short</th>
<th scope="col" class="org-left">CICE6 description</th>
<th scope="col" class="org-left">Units</th>
<th scope="col" class="org-left">Dataset</th>
<th scope="col" class="org-left">Short Name</th>
<th scope="col" class="org-left">Long Name</th>
<th scope="col" class="org-left">Derivation</th>
<th scope="col" class="org-left">Notes</th>
</tr>
</thead>

<tbody>
<tr>
<td class="org-left">zlvl</td>
<td class="org-left">atmosphere level height</td>
<td class="org-left">m</td>
<td class="org-left">-</td>
<td class="org-left">-</td>
<td class="org-left">-</td>
<td class="org-left">-</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">uatm</td>
<td class="org-left">model grid i-direction wind velocity component</td>
<td class="org-left">m/s</td>
<td class="org-left">ERA5</td>
<td class="org-left">10u</td>
<td class="org-left">10 metre U wind component</td>
<td class="org-left">no</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">vatm</td>
<td class="org-left">model grid j-direction wind velocity component</td>
<td class="org-left">m/s</td>
<td class="org-left">ERA5</td>
<td class="org-left">10v</td>
<td class="org-left">10 metre V wind component</td>
<td class="org-left">no</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">strax</td>
<td class="org-left">model grid i-direction wind stress</td>
<td class="org-left">N/m<sup>2</sup></td>
<td class="org-left">-</td>
<td class="org-left">-</td>
<td class="org-left">Mean eastward turbulent surface stress</td>
<td class="org-left">no</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">stray</td>
<td class="org-left">model grid j-direction wind stress</td>
<td class="org-left">N/m<sup>2</sup></td>
<td class="org-left">-</td>
<td class="org-left">-</td>
<td class="org-left">Mean northward turbulent surface stress</td>
<td class="org-left">no</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">potT</td>
<td class="org-left">air potential temperature</td>
<td class="org-left">K</td>
<td class="org-left">-</td>
<td class="org-left">-</td>
<td class="org-left">-</td>
<td class="org-left">yes</td>
<td class="org-left"><a href="https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.potential_temperature.html">derivation</a></td>
</tr>


<tr>
<td class="org-left">Tair</td>
<td class="org-left">air termperature</td>
<td class="org-left">K</td>
<td class="org-left">ERA5</td>
<td class="org-left">2t</td>
<td class="org-left">2 metre temperature</td>
<td class="org-left">no</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">Qa</td>
<td class="org-left">specific humidity</td>
<td class="org-left">kg/kg</td>
<td class="org-left">ERA5</td>
<td class="org-left">-</td>
<td class="org-left">-</td>
<td class="org-left">yes</td>
<td class="org-left"><a href="https://confluence.ecmwf.int/pages/viewpage.action?pageId=171411214">derivation</a></td>
</tr>


<tr>
<td class="org-left">rhoa</td>
<td class="org-left">air density</td>
<td class="org-left">kg/m<sup>3</sup></td>
<td class="org-left">ERA5</td>
<td class="org-left">-</td>
<td class="org-left">-</td>
<td class="org-left">yes</td>
<td class="org-left">compute <a href="https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.relative_humidity_from_specific_humidity.html">RH</a> then <a href="https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.mixing_ratio_from_relative_humidity.html?highlight=mixing%20ratio#metpy.calc.mixing_ratio_from_relative_humidity">mixing ratio</a> then <a href="https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.density.html?highlight=density">rho</a></td>
</tr>


<tr>
<td class="org-left">flw</td>
<td class="org-left">incoming longwave radiation</td>
<td class="org-left">W/m<sup>2</sup></td>
<td class="org-left">ERA5</td>
<td class="org-left">msdwlwrf</td>
<td class="org-left">Mean surface downward long-wave radiation flux</td>
<td class="org-left">no</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">fsw</td>
<td class="org-left">incoming shortwave radiation</td>
<td class="org-left">W/m<sup>2</sup></td>
<td class="org-left">ERA5</td>
<td class="org-left">msdwswrf</td>
<td class="org-left">Mean surface downward short-wave radiation flux</td>
<td class="org-left">no</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">swvdr</td>
<td class="org-left">sw down, visible, direct</td>
<td class="org-left">W/m<sup>2</sup></td>
<td class="org-left">-</td>
<td class="org-left">-</td>
<td class="org-left">-</td>
<td class="org-left">no</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">swvdf</td>
<td class="org-left">sw down, visible, diffuse</td>
<td class="org-left">W/m<sup>2</sup></td>
<td class="org-left">-</td>
<td class="org-left">-</td>
<td class="org-left">-</td>
<td class="org-left">yes</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">swidr</td>
<td class="org-left">sw down, near IR, direct</td>
<td class="org-left">W/m<sup>2</sup></td>
<td class="org-left">-</td>
<td class="org-left">-</td>
<td class="org-left">-</td>
<td class="org-left">yes</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">swidf</td>
<td class="org-left">sw down, near IR, diffuse</td>
<td class="org-left">W/m<sup>2</sup></td>
<td class="org-left">-</td>
<td class="org-left">-</td>
<td class="org-left">-</td>
<td class="org-left">yes</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">frain</td>
<td class="org-left">rainfall rate</td>
<td class="org-left">kg/(m<sup>2</sup> s)</td>
<td class="org-left">ERA5</td>
<td class="org-left">mtpr</td>
<td class="org-left">Mean total precipitation rate</td>
<td class="org-left">no</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">fsnow</td>
<td class="org-left">snowfall rate</td>
<td class="org-left">kg/(m<sup>2</sup> s)</td>
<td class="org-left">-</td>
<td class="org-left">-</td>
<td class="org-left">-</td>
<td class="org-left">no</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">uocn</td>
<td class="org-left">ocean current, x-direction</td>
<td class="org-left">m/s</td>
<td class="org-left">BRAN</td>
<td class="org-left">u</td>
<td class="org-left">eastward velocity</td>
<td class="org-left">no</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">vocn</td>
<td class="org-left">ocean current, y-direction</td>
<td class="org-left">m/s</td>
<td class="org-left">BRAN</td>
<td class="org-left">v</td>
<td class="org-left">northward velocity</td>
<td class="org-left">no</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">ss<sub>tltx</sub></td>
<td class="org-left">sea surface slope, x-direction</td>
<td class="org-left">m</td>
<td class="org-left">BRAN</td>
<td class="org-left">-</td>
<td class="org-left">-</td>
<td class="org-left">yes</td>
<td class="org-left">see \eqref{orgfddb70b} ; use <code>eta</code></td>
</tr>


<tr>
<td class="org-left">ss<sub>tlty</sub></td>
<td class="org-left">sea surface slope, y-direction</td>
<td class="org-left">m</td>
<td class="org-left">BRAN</td>
<td class="org-left">-</td>
<td class="org-left">-</td>
<td class="org-left">yes</td>
<td class="org-left">see \eqref{orgfddb70b} ; use <code>eta</code></td>
</tr>


<tr>
<td class="org-left">sss</td>
<td class="org-left">sea surface salinity</td>
<td class="org-left">ppt</td>
<td class="org-left">BRAN</td>
<td class="org-left">salt</td>
<td class="org-left">salinity</td>
<td class="org-left">no</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">sst</td>
<td class="org-left">sea surface temperature</td>
<td class="org-left">C</td>
<td class="org-left">BRAN</td>
<td class="org-left">temp</td>
<td class="org-left">temperature</td>
<td class="org-left">no</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">frzmlt</td>
<td class="org-left">freezing/melting potential</td>
<td class="org-left">W/m<sup>2</sup></td>
<td class="org-left">-</td>
<td class="org-left">-</td>
<td class="org-left">-</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">frzmlt<sub>init</sub></td>
<td class="org-left">frzmlt used in current time step</td>
<td class="org-left">W/m<sup>2</sup></td>
<td class="org-left">-</td>
<td class="org-left">-</td>
<td class="org-left">-</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">Tf</td>
<td class="org-left">freezing temperature</td>
<td class="org-left">C</td>
<td class="org-left">-</td>
<td class="org-left">-</td>
<td class="org-left">-</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">qdp</td>
<td class="org-left">deep ocean heat flux</td>
<td class="org-left">W/m<sup>2</sup></td>
<td class="org-left">-</td>
<td class="org-left">-</td>
<td class="org-left">-</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">\eqref{org03a0abf}</td>
</tr>


<tr>
<td class="org-left">hmix</td>
<td class="org-left">mixed layer depth</td>
<td class="org-left">m</td>
<td class="org-left">BRAN</td>
<td class="org-left">mld</td>
<td class="org-left">mixed layer depth</td>
<td class="org-left">no</td>
<td class="org-left">&#xa0;</td>
</tr>
</tbody>
</table>

1.  Computing Short-wave downward visibile diffuse

    \begin{equation}
    \label{org42131a3}
    SW_{diff} = 0.952 - 1.041\exp{-\exp{(2.3-4.702k_t)}}
    \end{equation}
    
    where,
    
    \begin{equation}
    
    k_t = SW / SW_{TOA}
    \end{equation}
    
    and $SW_{TOA}$ is the short-wave at the top-of-the-atmosphere

2.  Computing two-metre wind speed components

    \begin{equation}
    \label{orgf066b2d}
    u2 = (u10 * 4.87) / ln ((67.8 * 10) − 5.42)
    \end{equation}

3.  Computing sea surface slope

    \begin{equation}
    \label{orgfddb70b}
    ss_tltx = eta/dx
    ss_tlty = eta/dy
    \end{equation}
    
    where $dx$ and $dy$ is the grid spacing in meters in the longitudinal and latitudinal, respectively

4.  Computing deep ocean heat flux

    To compute the ocean heat flux at the mixed layer depth the following reference
    was used \citep{fofonoffAlgorithmsComputationFundamental1983}
    
    \begin{equation}
    \label{org03a0abf}
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
    
    where $T_{mld}$, $S_{mld}$ and $P_{mld}$ are the temperature (C), salinity (psu)
    and pressure (dbar) at the mixed layer depth, and $sst$ is the sea surface
    temperature (C).


<a id="org373a9d3"></a>

# External Drives


<a id="org6c2b476"></a>

## ioa01 (5000GB capacity &#x2013; 4940GB used)


<a id="org4552813"></a>

### Datasets

1.  BRAN (3000GB)

2.  JRA55do (140GB)

3.  MERRAv2 (1800GB &#x2013; incomplete)

    -   currently only contains surface fluxes ('flx') and radiative ('rad') components of the dataset
    -   only partial 'rad' upto mid-2015 (add another 400GB)
    -   also 'ocn' and 'slv' (see [2.5](#org6bedea4) above) need to be downloaded, with estimated
        600GB and 1400GB space required, respectively
    -   estimated total storage space for one decade is 4300GB, hence will require own
        5TB disk


<a id="org0a77467"></a>

## ioa02 (5000GB capacity &#x2013; 4170GB used)


<a id="org8c6f352"></a>

### Datasets

1.  ERA5 (2300GB)

2.  CAWCR (925GB)

3.  SOSE (675GB &#x2013; incomplete)

    missing U and V (both roughly 500GB each)

4.  C-GLORS (7GB)

5.  COREv2 (30GB)


<a id="orgc15bedc"></a>

## ioa03 (5000GB capacity &#x2013; varying used)

1.  Fast Ice (8GB)

2.  model<sub>input</sub>

3.  model<sub>output</sub>

4.  observations

5.  orography

6.  satellite


<a id="orgfd9826b"></a>

## Tracking Data Downloads

<table border="2" cellspacing="0" cellpadding="6" rules="groups" frame="hsides">


<colgroup>
<col  class="org-left" />

<col  class="org-left" />

<col  class="org-left" />

<col  class="org-left" />

<col  class="org-left" />

<col  class="org-left" />

<col  class="org-left" />

<col  class="org-left" />

<col  class="org-left" />

<col  class="org-left" />
</colgroup>
<thead>
<tr>
<th scope="col" class="org-left">Downloaded</th>
<th scope="col" class="org-left">Mandatory</th>
<th scope="col" class="org-left">Might</th>
<th scope="col" class="org-left">Reanalysis</th>
<th scope="col" class="org-left">&#xa0;</th>
<th scope="col" class="org-left">&#xa0;</th>
<th scope="col" class="org-left">&#xa0;</th>
<th scope="col" class="org-left">Download</th>
<th scope="col" class="org-left">&#xa0;</th>
<th scope="col" class="org-left">&#xa0;</th>
</tr>


<tr>
<th scope="col" class="org-left">&#xa0;</th>
<th scope="col" class="org-left">CICE</th>
<th scope="col" class="org-left">Be</th>
<th scope="col" class="org-left">or</th>
<th scope="col" class="org-left">Field</th>
<th scope="col" class="org-left">Time</th>
<th scope="col" class="org-left">Time</th>
<th scope="col" class="org-left">Location</th>
<th scope="col" class="org-left">Field</th>
<th scope="col" class="org-left">&#xa0;</th>
</tr>


<tr>
<th scope="col" class="org-left">&#xa0;</th>
<th scope="col" class="org-left">Field</th>
<th scope="col" class="org-left">Useful</th>
<th scope="col" class="org-left">Gridded</th>
<th scope="col" class="org-left">Short</th>
<th scope="col" class="org-left">Step</th>
<th scope="col" class="org-left">Coveage</th>
<th scope="col" class="org-left">&#xa0;</th>
<th scope="col" class="org-left">Long</th>
<th scope="col" class="org-left">Units</th>
</tr>


<tr>
<th scope="col" class="org-left">&#xa0;</th>
<th scope="col" class="org-left">&#xa0;</th>
<th scope="col" class="org-left">&#xa0;</th>
<th scope="col" class="org-left">Model</th>
<th scope="col" class="org-left">Name</th>
<th scope="col" class="org-left">&#xa0;</th>
<th scope="col" class="org-left">&#xa0;</th>
<th scope="col" class="org-left">&#xa0;</th>
<th scope="col" class="org-left">Name</th>
<th scope="col" class="org-left">&#xa0;</th>
</tr>
</thead>

<tbody>
<tr>
<td class="org-left">X</td>
<td class="org-left">x</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">ERA5</td>
<td class="org-left">2d</td>
<td class="org-left">1h</td>
<td class="org-left"><span class="timestamp-wrapper"><span class="timestamp">[2010-01-01 Fri 00:00]</span></span>-<span class="timestamp-wrapper"><span class="timestamp">[2019-12-31 Tue 23:00]</span></span></td>
<td class="org-left">gadi</td>
<td class="org-left">2 metre dewpoint temperature</td>
<td class="org-left">K</td>
</tr>


<tr>
<td class="org-left">X</td>
<td class="org-left">x</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">ERA5</td>
<td class="org-left">2t</td>
<td class="org-left">1h</td>
<td class="org-left"><span class="timestamp-wrapper"><span class="timestamp">[2010-01-01 Fri 00:00]</span></span>-<span class="timestamp-wrapper"><span class="timestamp">[2019-12-31 Tue 23:00]</span></span></td>
<td class="org-left">gadi</td>
<td class="org-left">2 metre temperature</td>
<td class="org-left">K</td>
</tr>


<tr>
<td class="org-left">x</td>
<td class="org-left">x</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">ERA5</td>
<td class="org-left">10u</td>
<td class="org-left">1h</td>
<td class="org-left"><span class="timestamp-wrapper"><span class="timestamp">[2010-01-01 Fri 00:00]</span></span>-<span class="timestamp-wrapper"><span class="timestamp">[2019-12-31 Tue 23:00]</span></span></td>
<td class="org-left">gadi</td>
<td class="org-left">10 metre u wind component</td>
<td class="org-left">m/s</td>
</tr>


<tr>
<td class="org-left">x</td>
<td class="org-left">x</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">ERA5</td>
<td class="org-left">10v</td>
<td class="org-left">1h</td>
<td class="org-left"><span class="timestamp-wrapper"><span class="timestamp">[2010-01-01 Fri 00:00]</span></span>-<span class="timestamp-wrapper"><span class="timestamp">[2019-12-31 Tue 23:00]</span></span></td>
<td class="org-left">gadi</td>
<td class="org-left">10 metre v wind component</td>
<td class="org-left">m/s</td>
</tr>


<tr>
<td class="org-left">x</td>
<td class="org-left">x</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">ERA5</td>
<td class="org-left">metss</td>
<td class="org-left">1h</td>
<td class="org-left"><span class="timestamp-wrapper"><span class="timestamp">[2010-01-01 Fri 00:00]</span></span>-<span class="timestamp-wrapper"><span class="timestamp">[2019-12-31 Tue 23:00]</span></span></td>
<td class="org-left">gadi</td>
<td class="org-left">Mean eastward turbulent surface stress</td>
<td class="org-left">N/m<sup>2</sup></td>
</tr>


<tr>
<td class="org-left">x</td>
<td class="org-left">x</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">ERA5</td>
<td class="org-left">mntss</td>
<td class="org-left">1h</td>
<td class="org-left"><span class="timestamp-wrapper"><span class="timestamp">[2010-01-01 Fri 00:00]</span></span>-<span class="timestamp-wrapper"><span class="timestamp">[2019-12-31 Tue 23:00]</span></span></td>
<td class="org-left">gadi</td>
<td class="org-left">Mean northward turbulent surface stress</td>
<td class="org-left">N/m<sup>2</sup></td>
</tr>


<tr>
<td class="org-left">x</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">x</td>
<td class="org-left">ERA5</td>
<td class="org-left">mror</td>
<td class="org-left">1h</td>
<td class="org-left"><span class="timestamp-wrapper"><span class="timestamp">[2010-01-01 Fri 00:00]</span></span>-<span class="timestamp-wrapper"><span class="timestamp">[2019-12-31 Tue 23:00]</span></span></td>
<td class="org-left">gadi</td>
<td class="org-left">Mean runoff rate</td>
<td class="org-left">kg/(m<sup>2</sup> s)</td>
</tr>


<tr>
<td class="org-left">x</td>
<td class="org-left">x</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">ERA5</td>
<td class="org-left">msdrswrf</td>
<td class="org-left">1h</td>
<td class="org-left"><span class="timestamp-wrapper"><span class="timestamp">[2010-01-01 Fri 00:00]</span></span>-<span class="timestamp-wrapper"><span class="timestamp">[2019-12-31 Tue 23:00]</span></span></td>
<td class="org-left">gadi</td>
<td class="org-left">Mean surface direct short-wave radiation flux</td>
<td class="org-left">W/m<sup>2</sup></td>
</tr>


<tr>
<td class="org-left">x</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">x</td>
<td class="org-left">ERA5</td>
<td class="org-left">msdrswrfcs</td>
<td class="org-left">1h</td>
<td class="org-left"><span class="timestamp-wrapper"><span class="timestamp">[2010-01-01 Fri 00:00]</span></span>-<span class="timestamp-wrapper"><span class="timestamp">[2019-12-31 Tue 23:00]</span></span></td>
<td class="org-left">gadi</td>
<td class="org-left">Mean surface direct short-wave radiation flux, clear sky</td>
<td class="org-left">W/m<sup>2</sup></td>
</tr>


<tr>
<td class="org-left">x</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">x</td>
<td class="org-left">ERA5</td>
<td class="org-left">msdwlwrf</td>
<td class="org-left">1h</td>
<td class="org-left"><span class="timestamp-wrapper"><span class="timestamp">[2010-01-01 Fri 00:00]</span></span>-<span class="timestamp-wrapper"><span class="timestamp">[2019-12-31 Tue 23:00]</span></span></td>
<td class="org-left">gadi</td>
<td class="org-left">Mean surface downward long-wave radiation flux</td>
<td class="org-left">W/m<sup>2</sup></td>
</tr>


<tr>
<td class="org-left">x</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">x</td>
<td class="org-left">ERA5</td>
<td class="org-left">msdwlwrfcs</td>
<td class="org-left">1h</td>
<td class="org-left"><span class="timestamp-wrapper"><span class="timestamp">[2010-01-01 Fri 00:00]</span></span>-<span class="timestamp-wrapper"><span class="timestamp">[2019-12-31 Tue 23:00]</span></span></td>
<td class="org-left">gadi</td>
<td class="org-left">Mean surface downward long-wave radiation flux, clear sky</td>
<td class="org-left">W/m<sup>2</sup></td>
</tr>


<tr>
<td class="org-left">x</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">x</td>
<td class="org-left">ERA5</td>
<td class="org-left">msdwswrf</td>
<td class="org-left">1h</td>
<td class="org-left"><span class="timestamp-wrapper"><span class="timestamp">[2010-01-01 Fri 00:00]</span></span>-<span class="timestamp-wrapper"><span class="timestamp">[2019-12-31 Tue 23:00]</span></span></td>
<td class="org-left">gadi</td>
<td class="org-left">Mean surface downward short-wave radiation flux</td>
<td class="org-left">W/m<sup>2</sup></td>
</tr>


<tr>
<td class="org-left">x</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">x</td>
<td class="org-left">ERA5</td>
<td class="org-left">msdwswrfcs</td>
<td class="org-left">1h</td>
<td class="org-left"><span class="timestamp-wrapper"><span class="timestamp">[2010-01-01 Fri 00:00]</span></span>-<span class="timestamp-wrapper"><span class="timestamp">[2019-12-31 Tue 23:00]</span></span></td>
<td class="org-left">gadi</td>
<td class="org-left">Mean surface downward short-wave radiation flux, clear sky</td>
<td class="org-left">W/m<sup>2</sup></td>
</tr>


<tr>
<td class="org-left">x</td>
<td class="org-left">x</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">ERA5</td>
<td class="org-left">msl</td>
<td class="org-left">1h</td>
<td class="org-left"><span class="timestamp-wrapper"><span class="timestamp">[2010-01-01 Fri 00:00]</span></span>-<span class="timestamp-wrapper"><span class="timestamp">[2019-12-31 Tue 23:00]</span></span></td>
<td class="org-left">gadi</td>
<td class="org-left">Mean sea level pressure</td>
<td class="org-left">Pa</td>
</tr>


<tr>
<td class="org-left">x</td>
<td class="org-left">x</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">ERA5</td>
<td class="org-left">msnlwrf</td>
<td class="org-left">1h</td>
<td class="org-left"><span class="timestamp-wrapper"><span class="timestamp">[2010-01-01 Fri 00:00]</span></span>-<span class="timestamp-wrapper"><span class="timestamp">[2019-12-31 Tue 23:00]</span></span></td>
<td class="org-left">gadi</td>
<td class="org-left">Mean surface net long-wave radiation flux</td>
<td class="org-left">W/m<sup>2</sup></td>
</tr>


<tr>
<td class="org-left">x</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">x</td>
<td class="org-left">ERA5</td>
<td class="org-left">msnlwrfcs</td>
<td class="org-left">1h</td>
<td class="org-left"><span class="timestamp-wrapper"><span class="timestamp">[2010-01-01 Fri 00:00]</span></span>-<span class="timestamp-wrapper"><span class="timestamp">[2019-12-31 Tue 23:00]</span></span></td>
<td class="org-left">gadi</td>
<td class="org-left">Mean surface net long-wave radiation flux, clear sky</td>
<td class="org-left">W/m<sup>2</sup></td>
</tr>


<tr>
<td class="org-left">x</td>
<td class="org-left">x</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">ERA5</td>
<td class="org-left">msnswrf</td>
<td class="org-left">1h</td>
<td class="org-left"><span class="timestamp-wrapper"><span class="timestamp">[2010-01-01 Fri 00:00]</span></span>-<span class="timestamp-wrapper"><span class="timestamp">[2019-12-31 Tue 23:00]</span></span></td>
<td class="org-left">gadi</td>
<td class="org-left">Mean surface net short-wave radiation flux</td>
<td class="org-left">W/m<sup>2</sup></td>
</tr>


<tr>
<td class="org-left">x</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">x</td>
<td class="org-left">ERA5</td>
<td class="org-left">msnswrfcs</td>
<td class="org-left">1h</td>
<td class="org-left"><span class="timestamp-wrapper"><span class="timestamp">[2010-01-01 Fri 00:00]</span></span>-<span class="timestamp-wrapper"><span class="timestamp">[2019-12-31 Tue 23:00]</span></span></td>
<td class="org-left">gadi</td>
<td class="org-left">Mean surface net short-wave radiation flux, clear sky</td>
<td class="org-left">W/m<sup>2</sup></td>
</tr>


<tr>
<td class="org-left">x</td>
<td class="org-left">x</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">ERA5</td>
<td class="org-left">msr</td>
<td class="org-left">1h</td>
<td class="org-left"><span class="timestamp-wrapper"><span class="timestamp">[2010-01-01 Fri 00:00]</span></span>-<span class="timestamp-wrapper"><span class="timestamp">[2019-12-31 Tue 23:00]</span></span></td>
<td class="org-left">gadi</td>
<td class="org-left">Mean snowfall rate</td>
<td class="org-left">kg/(m<sup>2</sup> s)</td>
</tr>


<tr>
<td class="org-left">x</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">x</td>
<td class="org-left">ERA5</td>
<td class="org-left">msror</td>
<td class="org-left">1h</td>
<td class="org-left"><span class="timestamp-wrapper"><span class="timestamp">[2010-01-01 Fri 00:00]</span></span>-<span class="timestamp-wrapper"><span class="timestamp">[2019-12-31 Tue 23:00]</span></span></td>
<td class="org-left">gadi</td>
<td class="org-left">Mean surface runoff rate</td>
<td class="org-left">kg/(m<sup>2</sup> s)</td>
</tr>


<tr>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">x</td>
<td class="org-left">ERA5</td>
<td class="org-left">mtdwswrf</td>
<td class="org-left">1h</td>
<td class="org-left"><span class="timestamp-wrapper"><span class="timestamp">[2010-01-01 Fri 00:00]</span></span>-<span class="timestamp-wrapper"><span class="timestamp">[2019-12-31 Tue 23:00]</span></span></td>
<td class="org-left">gadi</td>
<td class="org-left">Mean top downward short-wave radiation flux</td>
<td class="org-left">W/m<sup>2</sup></td>
</tr>


<tr>
<td class="org-left">x</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">x</td>
<td class="org-left">ERA5</td>
<td class="org-left">mtnlwrf</td>
<td class="org-left">1h</td>
<td class="org-left"><span class="timestamp-wrapper"><span class="timestamp">[2010-01-01 Fri 00:00]</span></span>-<span class="timestamp-wrapper"><span class="timestamp">[2019-12-31 Tue 23:00]</span></span></td>
<td class="org-left">gadi</td>
<td class="org-left">Mean top net long-wave radiation flux</td>
<td class="org-left">W/m<sup>2</sup></td>
</tr>


<tr>
<td class="org-left">x</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">x</td>
<td class="org-left">ERA5</td>
<td class="org-left">mtnlwrfcs</td>
<td class="org-left">1h</td>
<td class="org-left"><span class="timestamp-wrapper"><span class="timestamp">[2010-01-01 Fri 00:00]</span></span>-<span class="timestamp-wrapper"><span class="timestamp">[2019-12-31 Tue 23:00]</span></span></td>
<td class="org-left">gadi</td>
<td class="org-left">Mean top net long-wave radiation flux, clear sky</td>
<td class="org-left">W/m<sup>2</sup></td>
</tr>


<tr>
<td class="org-left">x</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">x</td>
<td class="org-left">ERA5</td>
<td class="org-left">mtnswrf</td>
<td class="org-left">1h</td>
<td class="org-left"><span class="timestamp-wrapper"><span class="timestamp">[2010-01-01 Fri 00:00]</span></span>-<span class="timestamp-wrapper"><span class="timestamp">[2019-12-31 Tue 23:00]</span></span></td>
<td class="org-left">gadi</td>
<td class="org-left">Mean top net short-wave radiation flux</td>
<td class="org-left">W/m<sup>2</sup></td>
</tr>


<tr>
<td class="org-left">x</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">x</td>
<td class="org-left">ERA5</td>
<td class="org-left">mtnswrfcs</td>
<td class="org-left">1h</td>
<td class="org-left"><span class="timestamp-wrapper"><span class="timestamp">[2010-01-01 Fri 00:00]</span></span>-<span class="timestamp-wrapper"><span class="timestamp">[2019-12-31 Tue 23:00]</span></span></td>
<td class="org-left">gadi</td>
<td class="org-left">Mean top net short-wave radiation flux, clear sky</td>
<td class="org-left">W/m<sup>2</sup></td>
</tr>


<tr>
<td class="org-left">x</td>
<td class="org-left">x</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">ERA5</td>
<td class="org-left">mtpr</td>
<td class="org-left">1h</td>
<td class="org-left"><span class="timestamp-wrapper"><span class="timestamp">[2010-01-01 Fri 00:00]</span></span>-<span class="timestamp-wrapper"><span class="timestamp">[2019-12-31 Tue 23:00]</span></span></td>
<td class="org-left">gadi</td>
<td class="org-left">Mean total precipitation rate</td>
<td class="org-left">kg/(m<sup>2</sup> s)</td>
</tr>


<tr>
<td class="org-left">x</td>
<td class="org-left">x</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">ERA5</td>
<td class="org-left">sp</td>
<td class="org-left">1h</td>
<td class="org-left"><span class="timestamp-wrapper"><span class="timestamp">[2010-01-01 Fri 00:00]</span></span>-<span class="timestamp-wrapper"><span class="timestamp">[2019-12-31 Tue 23:00]</span></span></td>
<td class="org-left">gadi</td>
<td class="org-left">Surface pressure</td>
<td class="org-left">Pa</td>
</tr>


<tr>
<td class="org-left">x</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">x</td>
<td class="org-left">CAWCR</td>
<td class="org-left">dataset</td>
<td class="org-left">1h</td>
<td class="org-left"><span class="timestamp-wrapper"><span class="timestamp">[2010-01-01 Fri 00:00]</span></span>-<span class="timestamp-wrapper"><span class="timestamp">[2021-12-31 Fri 23:00]</span></span></td>
<td class="org-left"><a href="https://data-cbr.csiro.au/thredds/catalog/catch_all/CMAR_CAWCR-Wave_archive/CAWCR_Wave_Hindcast_aggregate/gridded/catalog.html">CAWCR</a></td>
<td class="org-left">complete dataset</td>
<td class="org-left">various</td>
</tr>


<tr>
<td class="org-left">x</td>
<td class="org-left">x</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">BRAN</td>
<td class="org-left">u</td>
<td class="org-left">1h</td>
<td class="org-left"><span class="timestamp-wrapper"><span class="timestamp">[2010-01-01 Fri 00:00]</span></span>-<span class="timestamp-wrapper"><span class="timestamp">[2019-12-31 Tue 23:00]</span></span></td>
<td class="org-left">gadi</td>
<td class="org-left">eastward velocity</td>
<td class="org-left">m/s</td>
</tr>


<tr>
<td class="org-left">x</td>
<td class="org-left">x</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">BRAN</td>
<td class="org-left">v</td>
<td class="org-left">1h</td>
<td class="org-left"><span class="timestamp-wrapper"><span class="timestamp">[2010-01-01 Fri 00:00]</span></span>-<span class="timestamp-wrapper"><span class="timestamp">[2019-12-31 Tue 23:00]</span></span></td>
<td class="org-left">gadi</td>
<td class="org-left">northward velocity</td>
<td class="org-left">m/s</td>
</tr>


<tr>
<td class="org-left">x</td>
<td class="org-left">x</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">BRAN</td>
<td class="org-left">eta</td>
<td class="org-left">1h</td>
<td class="org-left"><span class="timestamp-wrapper"><span class="timestamp">[2010-01-01 Fri 00:00]</span></span>-<span class="timestamp-wrapper"><span class="timestamp">[2019-12-31 Tue 23:00]</span></span></td>
<td class="org-left">gadi</td>
<td class="org-left">surface height on T cells [Boussinesq (volume conserving) model]</td>
<td class="org-left">m</td>
</tr>


<tr>
<td class="org-left">x</td>
<td class="org-left">x</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">BRAN</td>
<td class="org-left">salt</td>
<td class="org-left">1h</td>
<td class="org-left"><span class="timestamp-wrapper"><span class="timestamp">[2010-01-01 Fri 00:00]</span></span>-<span class="timestamp-wrapper"><span class="timestamp">[2019-12-31 Tue 23:00]</span></span></td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">x</td>
<td class="org-left">x</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">BRAN</td>
<td class="org-left">temp</td>
<td class="org-left">1h</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">x</td>
<td class="org-left">x</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">BRAN</td>
<td class="org-left">mld</td>
<td class="org-left">1h</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">mixed layer depth determined by density criteria</td>
<td class="org-left">m</td>
</tr>


<tr>
<td class="org-left">x</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">x</td>
<td class="org-left">BRAN</td>
<td class="org-left">ice<sub>force</sub></td>
<td class="org-left">1h</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">contains two variables of short-wave that might be useful</td>
<td class="org-left">W/m<sup>2</sup></td>
</tr>
</tbody>
</table>


<a id="org29b629f"></a>

# Documenting Runs


<a id="org8fe0643"></a>

## Initial Trials


<a id="orgfe74489"></a>

### CICE Installation

[1.1](#org5179e48) was cloned onto local machine


<a id="org49c392a"></a>

### Grid, Forcing and Initial Conditions

Downloaded [1](#org28fbdde) forcing data, unzipped and copied data into directory
structure that matched [CICE github online setup guide](https://github.com/CICE-Consortium/CICE/wiki/CICE-Input-Data#directory-structure)


<a id="org389d754"></a>

### Test Run \`\`Vanilla''

A successful test run of a \`\`vanilla'' setup of CICE. Test directory
setup is preserved in [case1](file:///Users/dpath2o/PHD/MODELS/src/CICE/case1) and output is preserved in [case1](file:///Users/dpath2o/cice-dirs/runs/case1). What I
mean by \`\`vanilla'' is that I made no changes to any internal files
and simply ran the following (as per the instructions):

    cd ~/PHD/MODELS/src/CICE/
    ./cice.setup --case case1 --mach conda --env macos
    cd case1
    ./cice.build
    ./cice.run


<a id="org8150af8"></a>

## Sandbox


<a id="orgd05f616"></a>

### GX1 fail <span class="timestamp-wrapper"><span class="timestamp">[2022-05-23 Mon 06:00]</span></span>

Edited `ice_in` file replacing `gx3` enteries with `gx1` . While this ran
through to competion the results were garbarge


<a id="org1e539ea"></a>

### GX1 success <span class="timestamp-wrapper"><span class="timestamp">[2022-05-23 Mon 14:00]</span></span>

    ./cice.setup -m conda -e macos -c sandbox_gx1 -g gx1
    cd sandbox_gx1
    ./cice.build
    ./cice.run

![img](./CICE_runs/sandbox_gx1_bin/cice6-gx1_sice.png "This is the 'sice' field for GX1 vanilla")


<a id="orgc03d9b5"></a>

### GX3 success <span class="timestamp-wrapper"><span class="timestamp">[2022-05-23 Mon 14:30]</span></span>

    ./cice.setup -m conda -e macos -c sandbox_gx3 -g gx3
    cd sandbox_gx3
    ./cice.build
    ./cice.run

![img](./CICE_runs/sandbox_gx3_bin/cice6-gx3_sice.png "This is the 'sice' field for GX3 vanilla")


<a id="org40b9209"></a>

### TX1 (Binary Grid) partial success <span class="timestamp-wrapper"><span class="timestamp">[2022-06-06 Mon]</span></span>

This run through to completion with the following commands.

    ./cice.setup -m conda -e macos -c sandbox_tx1 -g tx1
    cd sandbox_tx1
    ./cice.build
    ./cice.run

![img](./CICE_runs/sandbox_tx1_bin/cice6-tx1_sice_BINgrid.png "This is the 'sice' field for TX1 vanilla")


<a id="org21a9a22"></a>

### TX1 (NetCDF Grid &#x2013; ACCESS-OM) partial success <span class="timestamp-wrapper"><span class="timestamp">[2022-06-06 Mon]</span></span>

The grid file can be viewed on the jupyter notebook [here](../src/py/cice_analysis.ipynb)

    ./cice.setup -m conda -e macos -c sandbox_tx1_nc -g tx1
    cd sandbox_tx1_nc

After initial setup I edited <./CICE_runs/sandbox_tx1_nc/ice_in> file
changing the following lines

    grid_format  = 'nc'
    grid_type    = 'tripole'
    grid_file    = '/Users/dpath2o/cice-dirs/input/CICE_data/grid/tx1/grid_tx1.nc'
    kmt_file     = 'unknown_kmt_file'

This runs through to completion with the following commands but the
results (see [264](#org7326c6f), below)

    ./cice.build
    ./cice.run

![img](./CICE_runs/sandbox_tx1_nc/cice6-tx1_sice_NCgrid.png "This is the 'sice' field for ACCESS-OM2 grid")


<a id="org5f83a3b"></a>

### TXp25 (NetCDF Grid) from ACCESS-OM2 forced by ERA5 and BRAN

Here is the cice setup

    ./cice.setup -m conda -e macos -c sandbox_tx1_nc -g tx1
    cd sandbox_tx1_nc

After initial setup I edited <./CICE_runs/sandbox_tx1_nc/ice_in> file
changing the following lines

    grid_format  = 'nc'
    grid_type    = 'tripole'
    grid_file    = '/Users/dpath2o/cice-dirs/input/CICE_data/grid/tx1/grid_tx1.nc'
    kmt_file     = 'unknown_kmt_file'

This runs through to completion with the following commands but the
results (see [264](#org7326c6f), below)

    ./cice.build
    ./cice.run

1.  First considered effort

    <span class="timestamp-wrapper"><span class="timestamp">[2022-06-07 Tue 09:00]</span></span>
    See post on [CICE Consortium Forum](https://bb.cgd.ucar.edu/cesm/threads/cice6-sandboxing.7412/post-45021) for details

2.  CICE setup alternate grid

    Create a new [configuration optional script file](file:///Users/dpath2o/src/CICE/configuration/scripts/options/set_nml.afim_0p25) that has the name-list parameters that you want start with, something like this:
    
    > dt                     = 3600.0
    > runtype                = 'initial'
    > year<sub>init</sub>              = 2010
    > use<sub>leap</sub><sub>years</sub>         = .false.
    > use<sub>restart</sub><sub>time</sub>       = .false.
    > ice<sub>ic</sub>                 = 'none'
    > grid<sub>format</sub>            = 'nc'
    > grid<sub>type</sub>              = 'tripole'
    > grid<sub>file</sub>              = '/Volumes/ioa03/cice-dirs/input/AFIM/grid/0p25/grid.nc'
    > kmt<sub>file</sub>               = '/Volumes/ioa03/cice-dirs/input/AFIM/grid/0p25/kmt.nc'
    > bathymetry<sub>file</sub>        = '/Volumes/ioa03/cice-dirs/input/AFIM/grid/0p25/topog.nc'
    > maskhalo<sub>dyn</sub>           = .true.
    > maskhalo<sub>remap</sub>         = .true.
    > maskhalo<sub>bound</sub>         = .true.
    > fyear<sub>init</sub>             = 2010
    > atm<sub>data</sub><sub>format</sub>        = 'nc'
    > atm<sub>data</sub><sub>type</sub>          = 'JRA55<sub>tx1</sub>'
    > atm<sub>data</sub><sub>dir</sub>           = '*Volumes/ioa03/cice-dirs/input/AFIM/forcing/0p25*'
    > precip<sub>units</sub>           = 'mm<sub>per</sub><sub>month</sub>'
    > ocn<sub>data</sub><sub>dir</sub>           = '/Volumes/ioa03/cice-dirs/input/AFIM/forcing/0p25/daily'
    > bgc<sub>data</sub><sub>dir</sub>           = ''
    > distribution<sub>wght</sub><sub>file</sub> = ''
    
    Then run CICE setup script with this new unique namelist file
    
        ./cice.setup -m conda -e macos -c sandbox_gx0p25 -s sandbox_gx0p25


<a id="org5eb5303"></a>

## AFIM stand-alone


<a id="org6fc8564"></a>

# Python Notebook and Scripts


<a id="org1f7b356"></a>

## [AFIM CICE Analysis](../src/python/CICE/cice_analysis.ipynb)

