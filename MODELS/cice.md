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
### GX3 Grid

Header from the [GX3 grid file](file:///Volumes/ioa03/cice-dirs/input/CICE_data/grid/gx3/grid_gx3.nc):

> netcdf grid~gx3~ { dimensions: y = 116 ; x = 100 ; variables: double
> ulat(y, x) ; ulat:units = \"rad\" ; double ulon(y, x) ; ulon:units =
> \"rad\" ; double htn(y, x) ; htn:units = \"cm\" ; double hte(y, x) ;
> hte:units = \"cm\" ; double hus(y, x) ; hus:units = \"cm\" ; double
> huw(y, x) ; huw:units = \"cm\" ; double angle(y, x) ; angle:units =
> \"rad\" ;
>
> // global attributes: :history = \"Tue Apr 17 16:11:59 2007: ncatted
> -astandard~name~,,d,, global~gx3~.grid.nc`\n`{=latex}\", \"Tue Apr 17
> 16:11:38 2007: ncatted -aunits,huw,c,c,cm
> global~gx3~.grid.nc`\n`{=latex}\", \"Tue Apr 17 16:11:34 2007: ncatted
> -aunits,hus,c,c,cm global~gx3~.grid.nc`\n`{=latex}\", \"Tue Apr 17
> 16:11:30 2007: ncatted -aunits,hte,c,c,cm
> global~gx3~.grid.nc`\n`{=latex}\", \"Tue Apr 17 16:11:27 2007: ncatted
> -aunits,htn,c,c,cm global~gx3~.grid.nc`\n`{=latex}\", \"Tue Apr 17
> 16:09:35 2007: ncatted -aunits,ulat,o,c,rad
> global~gx3~.grid.nc`\n`{=latex}\", \"Tue Apr 17 16:09:32 2007: ncatted
> -aunits,ulon,o,c,rad global~gx3~.grid.nc`\n`{=latex}\", \"Tue Apr 17
> 16:09:19 2007: ncatted -aunits,angle,c,c,rad
> global~gx3~.grid.nc`\n`{=latex}\", \"Tue Apr 17 16:09:04 2007: ncatted
> -aunits,ulon,c,c,rad global~gx3~.grid.nc`\n`{=latex}\", \"Tue Apr 17
> 16:08:59 2007: ncatted -aunits,ulat,c,c,rad
> global~gx3~.grid.nc`\n`{=latex}\", \"Tue Apr 17 16:06:06 2007: ncatted
> -aunits,ulon,c,c,degree~east~ global~gx3~.grid.nc`\n`{=latex}\", \"Tue
> Apr 17 16:05:53 2007: ncatted -aunits,ulat,c,c,degree~north~
> global~gx3~.grid.nc`\n`{=latex}\", \"Tue Apr 17 16:05:15 2007:
> ncrename -aunits,standard~name~ global~gx3~.grid.nc`\n`{=latex}\",
> \"Tue Apr 17 16:04:40 2007: ncatted -aunits,ulon,c,c,longitude
> global~gx3~.grid.nc`\n`{=latex}\", \"Tue Apr 17 16:04:31 2007: ncatted
> -aunits,ulat,c,c,latitude global~gx3~.grid.nc`\n`{=latex}\", \"Tue Apr
> 17 16:02:44 2007: ncatted -aConventions,global,c,c,CF-1.0
> global~gx3~.grid.nc`\n`{=latex}\", \"Tue Apr 17 15:55:18 2007: ncap -O
> -sulat=ulat\*0.0174533;ulon=ulon\*0.0174533;htn=htn\*100.0;hte=hte\*100.0;hus=hus\*100.0;huw=huw\*100.0;angle=angle\*0.0174533
> global~gx3~.grid.nc temp.nc`\n`{=latex}\", \"Tue Apr 17 15:34:39 2007:
> ncrename -vunspecified,ulat -vunspecified~2~,ulon -vunspecified~3~,htn
> -vunspecified~4~,hte -vunspecified~5~,hus -vunspecified~6~,huw
> -vunspecified~7~,angle temp.nc`\n`{=latex}\", \"Tue Apr 17 15:33:13
> 2007: ncks -x -a -vlongitude,latitude,unspecified~1~,t
> global~gx3~.grid.nc temp.nc`\n`{=latex}\", \"Tue Apr 17 15:32:26 2007:
> ncrename -dlongitude,x -dlatitude,y global~gx3~.grid.nc`\n`{=latex}\",
> \"Tue Apr 17 15:31:50 2007: ncwa -O -aunspecified~1~
> global~gx3~.grid.nc global~gx3~.grid.nc`\n`{=latex}\", \"Tue Apr 17
> 15:31:34 2007: ncwa -at global~gx3~.grid.nc
> global~gx3~.grid.nc`\n`{=latex}\", \"Tue Apr 17 15:31:16 2007: ncatted
> -a,,d,, global~gx3~.grid.nc`\n`{=latex}\", \"Tue Apr 17 15:30:50 BST
> 2007 - XCONV V1.91 Development\" ; :Conventions = \"CF-1.0\" ; }

### GX3 Bathymetry

Header from the [GX3 bathymetry
file](file:///Volumes/ioa03/cice-dirs/input/CICE_data/grid/gx3/global_gx3.bathy.nc):

> netcdf global~gx3~.bathy { dimensions: nj = 116 ; ni = 100 ;
> variables: float Bathymetry(nj, ni) ; Bathymetry:long~name~ = \"ocean
> bathymetry for grounding scheme\" ; Bathymetry:units = \"m\" ;
> Bathymetry:coordinates = \"TLON TLAT\" ; Bathymetry:missing~value~ =
> 1.e+30f ; Bathymetry:~FillValue~ = 1.e+30f ; float TLAT(nj, ni) ;
> TLAT:long~name~ = \"T grid center latitude\" ; TLAT:units =
> \"degrees~north~\" ; TLAT:missing~value~ = 1.e+30f ; TLAT:~FillValue~
> = 1.e+30f ; float TLON(nj, ni) ; TLON:long~name~ = \"T grid center
> longitude\" ; TLON:units = \"degrees~east~\" ; TLON:missing~value~ =
> 1.e+30f ; TLON:~FillValue~ = 1.e+30f ;
>
> // global attributes: :title = \"ocean bathymetry for grounding
> scheme\" ; :source = \"created from bathy~ORCA025LIM~.std\" ; :history
> = \"Mon Oct 29 15:57:35 2018: ncatted -a conventions,global,d,,
> gx3~bathy~.nc gx3~bathyTP12~.nc`\n`{=latex}\", \"Mon Oct 29 15:56:16
> 2018: ncatted -a source,global,o,c,created from bathy~ORCA025LIM~.std
> gx3~bathyTP10~.nc gx3~bathyTP11~.nc`\n`{=latex}\", \"Mon Oct 29
> 15:54:19 2018: ncatted -a contents,global,d,, gx3~bathyTP9~.nc
> gx3~bathyTP10~.nc`\n`{=latex}\", \"Mon Oct 29 15:53:44 2018: ncatted
> -a comment,global,d,, gx3~bathyTP8~.nc gx3~bathyTP9~.nc`\n`{=latex}\",
> \"Mon Oct 29 15:53:36 2018: ncatted -a comment2,global,d,,
> gx3~bathyTP7~.nc gx3~bathyTP8~.nc`\n`{=latex}\", \"Mon Oct 29 15:53:10
> 2018: ncatted -a comment3,global,d,, gx3~bathyTP6~.nc
> gx3~bathyTP7~.nc`\n`{=latex}\", \"Mon Oct 29 15:27:25 2018: ncatted -a
> title,global,o,c,ocean bathymetry for grounding scheme
> gx3~bathyTP5~.nc gx3~bathyTP6~.nc`\n`{=latex}\", \"Mon Oct 29 15:11:32
> 2018: ncks -x -v time,time~bounds~ gx3~bathyTP4~.nc
> gx3~bathyTP5~.nc`\n`{=latex}\", \"Mon Oct 29 15:10:09 2018: ncatted -a
> units,Bathymetry,o,c,m gx3~bathyTP3~.nc
> gx3~bathyTP4~.nc`\n`{=latex}\", \"Mon Oct 29 15:06:47 2018: ncatted -a
> long~name~,Bathymetry,o,c,ocean bathymetry for grounding scheme
> gx3~bathyTP2~.nc gx3~bathyTP3~.nc`\n`{=latex}\", \"Mon Oct 29 14:57:31
> 2018: ncks -x -v time~bounds~ gx3~bathyTP~.nc
> gx3~bathyTP2~.nc`\n`{=latex}\", \"Mon Oct 29 14:57:06 2018: ncks -x -v
> time gx3~bathyTP~.nc gx3~bathyTP2~.nc`\n`{=latex}\", \"Mon Oct 29
> 14:55:17 2018: ncks -x -v ULAT,ULON gx3~bathyTP~.nc
> gx3~bathyTP2~.nc`\n`{=latex}\", \"Thu Oct 25 18:26:01 2018: ncrename
> -v tarea,Bathymetry tempgx3.nc`\n`{=latex}\", \"Thu Oct 25 17:26:26
> 2018: ncks -x -v hi tempgx3.nc tempgx3out.nc`\n`{=latex}\", \"This
> dataset was created on 2018-10-25 at 17:24:50.3\" ; :io~flavor~ =
> \"io~netcdf~\" ; :NCO = \"4.4.2\" ; }

### GX3 Land Mask

Header from the [GX3 land mask
file](file:///Volumes/ioa03/cice-dirs/input/CICE_data/grid/gx3/kmt_gx3.nc):

> netcdf kmt~gx3~ { dimensions: y = 116 ; x = 100 ; variables: float
> kmt(y, x) ; kmt:units = \"1\" ;
>
> // global attributes: :history = \"Tue Apr 17 16:12:58 2007: ncatted
> -astandard~name~,,d,, global~gx3~.kmt.nc`\n`{=latex}\", \"Tue Apr 17
> 16:01:24 2007: ncatted -astandard~name~,,c,c,model~levelnumber~
> global~gx3~.kmt.nc`\n`{=latex}\", \"Tue Apr 17 16:00:39 2007: ncatted
> -aunits,,c,c,1 global~gx3~.kmt.nc`\n`{=latex}\", \"Tue Apr 17 15:58:38
> 2007: ncatted -aConventions,global,c,c,CF-1.0
> global~gx3~.kmt.nc`\n`{=latex}\", \"Tue Apr 17 15:27:57 2007: ncrename
> -dlongitude,x -dlatitude,y temp.nc`\n`{=latex}\", \"Tue Apr 17
> 15:27:39 2007: ncwa -aunspecified~1~ temp.nc temp.nc`\n`{=latex}\",
> \"Tue Apr 17 15:27:32 2007: ncwa -at temp.nc temp.nc`\n`{=latex}\",
> \"Tue Apr 17 15:26:58 2007: ncrename -vunspecified,kmt
> temp.nc`\n`{=latex}\", \"Tue Apr 17 15:26:21 2007: ncks -C
> -vunspecified global~gx3~.kmt.nc temp.nc`\n`{=latex}\", \"Tue Apr 17
> 15:25:37 2007: ncatted -a,,d,, global~gx3~.kmt.nc`\n`{=latex}\", \"Tue
> Apr 17 15:25:13 BST 2007 - XCONV V1.91 Development\" ; :Conventions =
> \"CF-1.0\" ; }

### GX3 Forcing

Three datasets are used for forcing: [*GX3 JRA55*]{.spurious-link
target="GX3 JRA55"}, [*GX3 NCAR*]{.spurious-link target="GX3 NCAR"}, and
[*GX3 WW3*]{.spurious-link target="GX3 WW3"}

1.  GX3 JRA55

    Separated into five yearly NetCDF files that have three-hourly data
    -- i.e.
    [8XDAILY](file:///Volumes/ioa03/cice-dirs/input/CICE_data/forcing/gx3/JRA55/8XDAILY).

    Dictories contents:

    > -rw-r-xr-- 1 dpath2o staff 905M 22 Feb 2020
    > JRA55~gx303hrforcing2005~.nc -rw-r-xr-- 1 dpath2o staff 905M 22
    > Feb 2020 JRA55~gx303hrforcing2006~.nc -rw-r-xr-- 1 dpath2o staff
    > 905M 22 Feb 2020 JRA55~gx303hrforcing2007~.nc -rw-r-xr-- 1 dpath2o
    > staff 907M 22 Feb 2020 JRA55~gx303hrforcing2008~.nc -rw-r-xr-- 1
    > dpath2o staff 905M 22 Feb 2020 JRA55~gx303hrforcing2009~.nc

    Here is the header information of the first file

    > netcdf JRA55~gx303hrforcing2005~ { dimensions: time = UNLIMITED ;
    > // (2920 currently) nj = 116 ; ni = 100 ; variables: float
    > time(time) ; time:long~name~ = \"model time\" ; time:units =
    > \"days since 1900-12-31 00:00:00\" ; time:calendar = \"standard\"
    > ; float LON(nj, ni) ; LON:long~name~ = \"Longitude\" ; LON:units =
    > \"degrees~east~\" ; float LAT(nj, ni) ; LAT:long~name~ =
    > \"Latitude\" ; LAT:units = \"degrees~north~\" ; float glbrad(time,
    > nj, ni) ; glbrad:long~name~ = \"downward surface shortwave\" ;
    > glbrad:units = \"W/m^2^\" ; glbrad:coordinates = \"LON LAT\" ;
    > glbrad:~FillValue~ = 1.e+30f ; float dlwsfc(time, nj, ni) ;
    > dlwsfc:long~name~ = \"downward surface longwave\" ; dlwsfc:units =
    > \"W/m^2^\" ; dlwsfc:coordinates = \"LON LAT\" ; dlwsfc:~FillValue~
    > = 1.e+30f ; float wndewd(time, nj, ni) ; wndewd:long~name~ =
    > \"x-ward winds\" ; wndewd:units = \"m/s\" ; wndewd:coordinates =
    > \"LON LAT\" ; wndewd:~FillValue~ = 1.e+30f ; float wndnwd(time,
    > nj, ni) ; wndnwd:long~name~ = \"y-ward winds\" ; wndnwd:units =
    > \"m/s\" ; wndnwd:coordinates = \"LON LAT\" ; wndnwd:~FillValue~ =
    > 1.e+30f ; float airtmp(time, nj, ni) ; airtmp:long~name~ = \"Air
    > temperature\" ; airtmp:units = \"Kelvin\" ; airtmp:coordinates =
    > \"LON LAT\" ; airtmp:~FillValue~ = 1.e+30f ; float spchmd(time,
    > nj, ni) ; spchmd:long~name~ = \"Specific Humidity\" ; spchmd:units
    > = \"kg/kg\" ; spchmd:coordinates = \"LON LAT\" ;
    > spchmd:~FillValue~ = 1.e+30f ; float ttlpcp(time, nj, ni) ;
    > ttlpcp:long~name~ = \"Precipitation\" ; ttlpcp:units = \"kg m-2
    > s-1\" ; ttlpcp:coordinates = \"LON LAT\" ; ttlpcp:~FillValue~ =
    > 1.e+30f ; }

2.  GX3 NCAR

    Separate into
    [4XDAILY](file:///Volumes/ioa03/cice-dirs/input/CICE_data/forcing/gx3/NCAR_bulk/4XDAILY)
    (i.e. six-hourly) and
    [monthly](file:///Volumes/ioa03/cice-dirs/input/CICE_data/forcing/gx3/NCAR_bulk/MONTHLY)
    files, which are all binary, and hence NIL header information is
    listed here. The contents of the files can be dedecued from the file
    names:

      filename   my guess
      ---------- ------------------------------
      cldf       cloud
      prec       precipitation
      swdn       net downward shortwave
      dn10       dew-point temperature @ 10 m
      q~10~      specific humidity @ 10 m
      t~10~      air temperature @ 10 m
      u~10~      zonal wind speed @ 10 m
      v~10~      meridional wind speed @ 10 m

    The directory contents are provided as:

    1.  GX3 NCAR 4XDAILY

        ```{=org}
        #+BEGIN_QUOTE:
        ```
        -rw-r--r-- 1 dpath2o staff 1.1M 12 Dec 2003 cldf.1997.dat
        -rw-r--r-- 1 dpath2o staff 1.1M 12 Dec 2003 prec.1997.dat
        -rw-r--r-- 1 dpath2o staff 1.1M 12 Dec 2003 swdn.1997.dat

        ```{=org}
        #+END_QUOTE
        ```

    2.  GX3 NCAR MONTHLY

        > -rw-r--r-- 1 dpath2o staff 129M 12 Dec 2003 dn10.1997.dat
        > -rw-r--r-- 1 dpath2o staff 129M 12 Dec 2003 q~10~.1997.dat
        > -rw-r--r-- 1 dpath2o staff 129M 12 Dec 2003 t~10~.1997.dat
        > -rw-r--r-- 1 dpath2o staff 129M 12 Dec 2003 u~10~.1997.dat
        > -rw-r--r-- 1 dpath2o staff 129M 12 Dec 2003 v~10~.1997.dat

3.  GX3 WW3

    [One wave spectral
    file](file:///Volumes/ioa03/cice-dirs/input/CICE_data/forcing/gx3/WW3/ww3.20100101_efreq_remapgx3.nc)
    is given with a total of four time steps and 25 frequencies. The
    header is provided as:

    > time = UNLIMITED ; // (4 currently) ni = 100 ; nj = 116 ; f = 25 ;
    > variables: double time(time) ; time:standard~name~ = \"time\" ;
    > time:long~name~ = \"julian day (UT)\" ; time:units = \"days since
    > 1990-01-01 00:00:00\" ; time:calendar = \"standard\" ; time:axis =
    > \"T\" ; float TLON(nj, ni) ; TLON:standard~name~ = \"longitude\" ;
    > TLON:long~name~ = \"longitude\" ; TLON:units = \"degrees~east~\" ;
    > TLON:~CoordinateAxisType~ = \"Lon\" ; float TLAT(nj, ni) ;
    > TLAT:standard~name~ = \"latitude\" ; TLAT:long~name~ =
    > \"latitude\" ; TLAT:units = \"degrees~north~\" ;
    > TLAT:~CoordinateAxisType~ = \"Lat\" ; float f(f) ; f:long~name~ =
    > \"wave~frequency~\" ; f:units = \"s-1\" ; f:axis = \"Hz\" ;
    > f:standard~name~ = \"wave~frequency~\" ; float efreq(time, f, nj,
    > ni) ; efreq:standard~name~ =
    > \"power~spectraldensityofsurfaceelevation~\" ; efreq:long~name~ =
    > \"wave~elevationspectrum~\" ; efreq:units = \"m2 s\" ;
    > efreq:coordinates = \"TLAT TLON\" ; efreq:~FillValue~ =
    > 9.96921e+36f ; efreq:missing~value~ = 9.96921e+36f ;
    > efreq:globwave~name~ =
    > \"power~spectraldensityofsurfaceelevation~\" ;
    >
    > // global attributes: :CDI = \"Climate Data Interface version
    > 1.9.7.1 (<http://mpimet.mpg.de/cdi>)\" ; :history = \"Wed Sep 04
    > 15:05:27 2019: cdo chname,ef,efreq ww3.20100101~efremapgx3~.nc
    > ww3.20100101~efreqremapgx3~.nc`\nWed`{=latex} Sep 04 13:37:05
    > 2019: cdo remapnn,../gx3~tgrid~\_.nc ww3.20100101~ef~\_.nc
    > ww3.20100101~efremapgx3~.nc`\nWed`{=latex} Sep 4 13:36:54 2019:
    > ncks -x -v MAPSTA ww3.20100101~ef~.nc
    > ww3.20100101~ef~\_.nc`\nWed`{=latex} Sep 4 13:28:07 2019: ncatted
    > -a coordinates,ef,c,c,longitude latitude ww3.20100101~ef~.nc\" ;
    > :Conventions = \"CF-1.6\" ; :WAVEWATCH~IIIversionnumber~ =
    > \"5.16\" ; :WAVEWATCH~IIIswitches~ = \"F90 NOGRB NOPA LRB4 NC4 PR3
    > UQ FLX0 LN1 ST4 BT1 DB1 MLIM NL1 TR0 BS0 REF0 IS0 IC4 XX0 WNT1
    > WNX1 CRT1 CRX1 O0 O1 O2 O3 O4 O5 O6 O7 O11 SHRD\" ; :product~name~
    > = \"ww3.20100101~ef~.nc\" ; :area = \"POP 1 degree grid (gx1v6b)\"
    > ; :latitude~resolution~ = \"n/a\" ; :longitude~resolution~ =
    > \"n/a\" ; :southernmost~latitude~ = \"-79.22052\" ;
    > :northernmost~latitude~ = \"89.70641\" ; :westernmost~longitude~ =
    > \"1.4731102E-02\" ; :easternmost~longitude~ = \"359.9960\" ;
    > :minimum~altitude~ = \"-12000 m\" ; :maximum~altitude~ = \"9000
    > m\" ; :altitude~resolution~ = \"n/a\" ; :start~date~ =
    > \"2010-01-01 00:00:00\" ; :stop~date~ = \"2010-01-01 18:00:00\" ;
    > :NCO = \"netCDF Operators version 4.7.9 (Homepage =
    > <http://nco.sf.net>, Code = <http://github.com/nco/nco>)\" ; :CDO
    > = \"Climate Data Operators version 1.9.7.1
    > (<http://mpimet.mpg.de/cdo>)\" ; }

### GX3 Initial Conditions

There are a total of 12 initial condition files provided. With each file
containing NIL time step across the This is interesting to me as I
thought only one would be required. Nonetheless, the header of the first
file is given below.

> netcdf iced~gx3v6~.2005-01-01 { dimensions: ni = 100 ; nj = 116 ; ncat
> = 5 ; variables: double uvel(nj, ni) ; double vvel(nj, ni) ; double
> scale~factor~(nj, ni) ; double swvdr(nj, ni) ; double swvdf(nj, ni) ;
> double swidr(nj, ni) ; double swidf(nj, ni) ; double strocnxT(nj, ni)
> ; double strocnyT(nj, ni) ; double stressp~1~(nj, ni) ; double
> stressp~2~(nj, ni) ; double stressp~3~(nj, ni) ; double stressp~4~(nj,
> ni) ; double stressm~1~(nj, ni) ; double stressm~2~(nj, ni) ; double
> stressm~3~(nj, ni) ; double stressm~4~(nj, ni) ; double
> stress12~1~(nj, ni) ; double stress12~2~(nj, ni) ; double
> stress12~3~(nj, ni) ; double stress12~4~(nj, ni) ; double iceumask(nj,
> ni) ; double sst(nj, ni) ; double frzmlt(nj, ni) ; double
> frz~onset~(nj, ni) ; double fsnow(nj, ni) ; double aicen(ncat, nj, ni)
> ; double vicen(ncat, nj, ni) ; double vsnon(ncat, nj, ni) ; double
> Tsfcn(ncat, nj, ni) ; double iage(ncat, nj, ni) ; double FY(ncat, nj,
> ni) ; double alvl(ncat, nj, ni) ; double vlvl(ncat, nj, ni) ; double
> apnd(ncat, nj, ni) ; double hpnd(ncat, nj, ni) ; double ipnd(ncat, nj,
> ni) ; double dhs(ncat, nj, ni) ; double ffrac(ncat, nj, ni) ; double
> fbrn(ncat, nj, ni) ; double first~ice~(ncat, nj, ni) ; double
> sice001(ncat, nj, ni) ; double qice001(ncat, nj, ni) ; double
> sice002(ncat, nj, ni) ; double qice002(ncat, nj, ni) ; double
> sice003(ncat, nj, ni) ; double qice003(ncat, nj, ni) ; double
> sice004(ncat, nj, ni) ; double qice004(ncat, nj, ni) ; double
> sice005(ncat, nj, ni) ; double qice005(ncat, nj, ni) ; double
> sice006(ncat, nj, ni) ; double qice006(ncat, nj, ni) ; double
> sice007(ncat, nj, ni) ; double qice007(ncat, nj, ni) ; double
> qsno001(ncat, nj, ni) ;
>
> // global attributes: :istep1 = 87672 ; :time = 315619200. ;
> :time~forc~ = 315619200. ; :nyr = 11 ; :month = 1 ; :mday = 1 ; :sec =
> 0 ; :created = \"2021-03-19\" ; :history = \"Wed Apr 7 12:29:43 2021:
> ncatted --attribute created,global,c,c,2021-03-19
> iced~gx3v6~.2005-01-01.nc\" ; :NCO = \"netCDF Operators version 4.9.5
> (Homepage = <http://nco.sf.net>, Code = <http://github.com/nco/nco>)\"
> ; }

## GX1 Description

Global One-Degree Cartesian

### GX1 Grid

For some reason the [GX1 grid
file](file:///Volumes/ioa03/cice-dirs/input/CICE_data/grid/gx1/grid_gx1.bin)
is in binary format and the header information is not easily accessible.

### GX1 Bathymetry

Header from the [GX1 bathymetry
file](file:///Volumes/ioa03/cice-dirs/input/CICE_data/grid/gx1/grid_gx1.bin):

> netcdf global~gx1~.bathy { dimensions: nj = 384 ; ni = 320 ;
> variables: float Bathymetry(nj, ni) ; Bathymetry:long~name~ = \"ocean
> bathymetry for grounding scheme\" ; Bathymetry:units = \"m\" ;
> Bathymetry:coordinates = \"TLON TLAT\" ; Bathymetry:missing~value~ =
> 1.e+30f ; Bathymetry:~FillValue~ = 1.e+30f ; float TLAT(nj, ni) ;
> TLAT:long~name~ = \"T grid center latitude\" ; TLAT:units =
> \"degrees~north~\" ; TLAT:missing~value~ = 1.e+30f ; TLAT:~FillValue~
> = 1.e+30f ; float TLON(nj, ni) ; TLON:long~name~ = \"T grid center
> longitude\" ; TLON:units = \"degrees~east~\" ; TLON:missing~value~ =
> 1.e+30f ; TLON:~FillValue~ = 1.e+30f ;
>
> // global attributes: :title = \"ocean bathymetry for grounding
> scheme\" ; :source = \"created from bathy~ORCA025LIM~.std\" ; :history
> = \"Mon Oct 29 16:01:01 2018: ncatted -a source,global,o,c,created
> from bathy~ORCA025LIM~.std gx1~bathyTP10~.nc
> gx1~bathyTP11~.nc`\n`{=latex}\", \"Mon Oct 29 16:00:02 2018: ncatted
> -a contents,global,d,, gx1~bathyTP9~.nc
> gx1~bathyTP10~.nc`\n`{=latex}\", \"Mon Oct 29 15:59:40 2018: ncatted
> -a comment,global,d,, gx1~bathyTP8~.nc gx1~bathyTP9~.nc`\n`{=latex}\",
> \"Mon Oct 29 15:59:33 2018: ncatted -a comment2,global,d,,
> gx1~bathyTP7~.nc gx1~bathyTP8~.nc`\n`{=latex}\", \"Mon Oct 29 15:59:26
> 2018: ncatted -a comment3,global,d,, gx1~bathyTP6~.nc
> gx1~bathyTP7~.nc`\n`{=latex}\", \"Mon Oct 29 15:59:13 2018: ncatted -a
> conventions,global,d,, gx1~bathyTP5~.nc
> gx1~bathyTP6~.nc`\n`{=latex}\", \"Mon Oct 29 15:27:15 2018: ncatted -a
> title,global,o,c,ocean bathymetry for grounding scheme
> gx1~bathyTP4~.nc gx1~bathyTP5~.nc`\n`{=latex}\", \"Mon Oct 29 15:13:20
> 2018: ncks -x -v time,time~bounds~ gx1~bathyTP3~.nc
> gx1~bathyTP4~.nc`\n`{=latex}\", \"Mon Oct 29 15:13:04 2018: ncatted -a
> units,Bathymetry,o,c,m gx1~bathyTP2~.nc
> gx1~bathyTP3~.nc`\n`{=latex}\", \"Mon Oct 29 15:12:46 2018: ncatted -a
> long~name~,Bathymetry,o,c,ocean bathymetry for grounding scheme
> gx1~bathy~.nc gx1~bathyTP2~.nc`\n`{=latex}\", \"Fri Oct 26 13:28:47
> 2018: ncrename -v tarea,Bathymetry tempgx1.nc`\n`{=latex}\", \"Thu Oct
> 25 17:17:42 2018: ncks -x -v
> ANGLE,ANGLET,NCAT,ULAT,ULON,VGRDa,aice,blkmask,hi,taubx,tauby,time~bounds~,tmask,uarea,uvel,vvel
> tempgx1.nc tempgx1out.nc`\n`{=latex}\", \"This dataset was created on
> 2018-10-22 at 18:12:11.5\" ; :io~flavor~ = \"io~netcdf~\" ; :NCO =
> \"4.4.2\" ; }

### GX1 Land Mask

For some reason the [GX1 land mask
file](file:///Volumes/ioa03/cice-dirs/input/CICE_data/grid/gx1/kmt_gx1.bin)
is in binary format and the header information is not easily accessible.

### GX1 Forcing

Four datasets are used for forcing: [*GX1 CESM*]{.spurious-link
target="GX1 CESM"}, [*GX1 COREII*]{.spurious-link target="GX1 COREII"},
[*GX1 JRA55*]{.spurious-link target="GX1 JRA55"}, and \[\[GX1 WOA\]\]. I
find it interesting to note that WaveWatch III data is not provided.

1.  GX1 CESM

    One file provides monthly summarised oceanographic data in the
    following format. Note that the time dimension is only 12 increments
    in length.

    > netcdf ocean~forcingclim2Dgx1~ { dimensions: time = 12 ; nj = 384
    > ; ni = 320 ; variables: double area(nj, ni) ; area:long~name~ =
    > \"area of grid cell in radians squared\" ; area:units = \"area\" ;
    > int mask(nj, ni) ; mask:units = \"unitless\" ; mask:long~name~ =
    > \"domain maskr\" ; mask:~FillValue~ = -2147483647 ;
    > mask:missing~value~ = -2147483647 ; double xc(nj, ni) ;
    > xc:missing~value~ = 9.96920996838687e+36 ; xc:~FillValue~ =
    > 9.96920996838687e+36 ; xc:units = \"degrees east\" ; xc:long~name~
    > = \"longitude of grid cell center\" ; double yc(nj, ni) ;
    > yc:missing~value~ = 9.96920996838687e+36 ; yc:~FillValue~ =
    > 9.96920996838687e+36 ; yc:units = \"degrees north\" ;
    > yc:long~name~ = \"latitude of grid cell center\" ; float
    > time(time) ; time:calendar = \"noleap\" ; time:long~name~ =
    > \"observation time\" ; time:units = \"days since 0001-01-01
    > 00:00:00\" ; float S(time, nj, ni) ; S:units = \"ppt\" ;
    > S:long~name~ = \"salinity\" ; S:~FillValue~ = 9.96921e+36f ; float
    > T(time, nj, ni) ; T:units = \"degC\" ; T:long~name~ =
    > \"temperature\" ; T:~FillValue~ = 9.96921e+36f ; float U(time, nj,
    > ni) ; U:units = \"m/s\" ; U:long~name~ = \"u ocean current\" ;
    > U:~FillValue~ = 9.96921e+36f ; float V(time, nj, ni) ; V:units =
    > \"m/s\" ; V:long~name~ = \"v ocean current\" ; V:~FillValue~ =
    > 9.96921e+36f ; float dhdx(time, nj, ni) ; dhdx:units = \"m/m\" ;
    > dhdx:long~name~ = \"ocean surface slope: zonal\" ;
    > dhdx:~FillValue~ = 9.96921e+36f ; float dhdy(time, nj, ni) ;
    > dhdy:units = \"m/m\" ; dhdy:long~name~ = \"ocean surface slope:
    > meridional\" ; dhdy:~FillValue~ = 9.96921e+36f ; float hblt(time,
    > nj, ni) ; hblt:units = \"m\" ; hblt:long~name~ = \"boundary layer
    > depth\" ; hblt:~FillValue~ = 9.96921e+36f ; float qdp(time, nj,
    > ni) ; qdp:units = \"W/m^2^\" ; qdp:long~name~ = \"ocean heat flux
    > convergence\" ; qdp:~FillValue~ = 9.96921e+36f ;
    >
    > // global attributes: :creation~date~ = \"Tue Feb 17 09:45:48 MST
    > 2015\" ; :comment = \"This data is on the displaced pole grid
    > gx1v5\" ; :calendar = \"standard\" ; :author = \"D. Bailey\" ;
    > :note3 = \"qdp is computed from depth summed ocean column\" ;
    > :note2 = \"all fields interpolated to T-grid\" ; :note1 = \"fields
    > computed from years 402 to 1510 monthly means from pop\" ;
    > :description = \"Input data for DOCN7 mixed layer model from
    > b.e11.B1850C5CN.f09~g16~.005\" ; :source = \"pop~frc~.ncl\" ;
    > :conventions = \"CCSM data model domain description\" ; :title =
    > \"Monthly averaged ocean forcing from POP output\" ; :history =
    > \"Wed Mar 24 11:42:29 2021: ncatted --glb~attadd~
    > description2=From the CESM-LE control run as in Kay et al. 2015
    > pop~frc~.gx1.nc\" ; :description2 = \"From the CESM-LE control run
    > as in Kay et al. 2015\" ; :NCO = \"netCDF Operators version 4.9.5
    > (Homepage = <http://nco.sf.net>, Code =
    > <http://github.com/nco/nco>)\" ; }

2.  GX1 COREII

    These files are given as binary with two sub-directories
    [4XDAILY](file:///Volumes/ioa03/cice-dirs/input/CICE_data/forcing/gx1/COREII/4XDAILY)
    and
    [MONTHLY](file:///Volumes/ioa03/cice-dirs/input/CICE_data/forcing/gx1/COREII/MONTHLY).
    The file names indicate the representative variable contained
    within.

    Here is the directory listing of the
    [4XDAILY](file:///Volumes/ioa03/cice-dirs/input/CICE_data/forcing/gx1/COREII/4XDAILY/):

    > -rwxr-xr-x 1 dpath2o staff 1.3G 30 Aug 2017 q~10~.2005.dat
    > -rwxr-xr-x 1 dpath2o staff 1.3G 30 Aug 2017 q~10~.2006.dat
    > -rwxr-xr-x 1 dpath2o staff 1.3G 30 Aug 2017 q~10~.2007.dat
    > -rwxr-xr-x 1 dpath2o staff 1.3G 30 Aug 2017 q~10~.2008.dat
    > -rwxr-xr-x 1 dpath2o staff 1.3G 30 Aug 2017 q~10~.2009.dat
    > -rwxr-xr-x 1 dpath2o staff 1.3G 30 Aug 2017 t~10~.2005.dat
    > -rwxr-xr-x 1 dpath2o staff 1.3G 30 Aug 2017 t~10~.2006.dat
    > -rwxr-xr-x 1 dpath2o staff 1.3G 30 Aug 2017 t~10~.2007.dat
    > -rwxr-xr-x 1 dpath2o staff 1.3G 30 Aug 2017 t~10~.2008.dat
    > -rwxr-xr-x 1 dpath2o staff 1.3G 30 Aug 2017 t~10~.2009.dat
    > -rwxr-xr-x 1 dpath2o staff 1.3G 30 Aug 2017 u~10~.2005.dat
    > -rwxr-xr-x 1 dpath2o staff 1.3G 30 Aug 2017 u~10~.2006.dat
    > -rwxr-xr-x 1 dpath2o staff 1.3G 30 Aug 2017 u~10~.2007.dat
    > -rwxr-xr-x 1 dpath2o staff 1.3G 30 Aug 2017 u~10~.2008.dat
    > -rwxr-xr-x 1 dpath2o staff 1.3G 31 Aug 2017 u~10~.2009.dat
    > -rwxr-xr-x 1 dpath2o staff 1.3G 30 Aug 2017 v~10~.2005.dat
    > -rwxr-xr-x 1 dpath2o staff 1.3G 30 Aug 2017 v~10~.2006.dat
    > -rwxr-xr-x 1 dpath2o staff 1.3G 30 Aug 2017 v~10~.2007.dat
    > -rwxr-xr-x 1 dpath2o staff 1.3G 30 Aug 2017 v~10~.2008.dat
    > -rwxr-xr-x 1 dpath2o staff 1.3G 31 Aug 2017 v~10~.2009.dat

    Here is the listing of the
    [MONTHLY](file:///Volumes/ioa03/cice-dirs/input/CICE_data/forcing/gx1/COREII/MONTHLY):

    > -rw-r--r-- 1 dpath2o staff 11M 12 Jul 2012 cldf.omip.dat
    > -rw-r--r-- 1 dpath2o staff 11M 12 Jul 2012 prec.nmyr.dat

3.  GX1 JRA55

    The forcing files are nearly identical to
    [*GX3 JRA55*]{.spurious-link target="GX3 JRA55"} with the exception
    of the change in grid. This is most notable in comparing the file
    sizes -- [*GX3 JRA55*]{.spurious-link target="GX3 JRA55"} files are
    0.9GB and these GX1 JRA55 files are 9.5GB. Here is the directory
    listing of the
    [8XDAILY](file:///Volumes/ioa03/cice-dirs/input/CICE_data/forcing/gx1/JRA55/8XDAILY/):

    > -rw-r--r--@ 1 dpath2o staff 9.4G 22 Aug 2019
    > JRA55~03hrforcing2005~.nc -rw-r--r-- 1 dpath2o staff 9.4G 22 Aug
    > 2019 JRA55~03hrforcing2006~.nc -rw-r--r-- 1 dpath2o staff 9.4G 22
    > Aug 2019 JRA55~03hrforcing2007~.nc -rw-r--r-- 1 dpath2o staff 9.4G
    > 22 Aug 2019 JRA55~03hrforcing2008~.nc -rw-r--r-- 1 dpath2o staff
    > 9.4G 22 Aug 2019 JRA55~03hrforcing2009~.nc

    Here is the header information of the first file:

    > netcdf JRA55~03hrforcing2005~ { dimensions: nj = 384 ; ni = 320 ;
    > time = UNLIMITED ; // (2920 currently) variables: float LAT(nj,
    > ni) ; LAT:long~name~ = \"Latitude\" ; LAT:units =
    > \"degrees~north~\" ; float LON(nj, ni) ; LON:long~name~ =
    > \"Longitude\" ; LON:units = \"degrees~east~\" ; float airtmp(time,
    > nj, ni) ; airtmp:long~name~ = \"Air temperature\" ; airtmp:units =
    > \"Kelvin\" ; airtmp:coordinates = \"LON LAT\" ; airtmp:~FillValue~
    > = 1.e+30f ; float dlwsfc(time, nj, ni) ; dlwsfc:long~name~ =
    > \"downward surface longwave\" ; dlwsfc:units = \"W/m^2^\" ;
    > dlwsfc:coordinates = \"LON LAT\" ; dlwsfc:~FillValue~ = 1.e+30f ;
    > float glbrad(time, nj, ni) ; glbrad:long~name~ = \"downward
    > surface shortwave\" ; glbrad:units = \"W/m^2^\" ;
    > glbrad:coordinates = \"LON LAT\" ; glbrad:~FillValue~ = 1.e+30f ;
    > float spchmd(time, nj, ni) ; spchmd:long~name~ = \"Specific
    > Humidity\" ; spchmd:units = \"kg/kg\" ; spchmd:coordinates = \"LON
    > LAT\" ; spchmd:~FillValue~ = 1.e+30f ; float time(time) ;
    > time:long~name~ = \"model time\" ; time:units = \"days since
    > 1900-12-31 00:00:00\" ; time:calendar = \"standard\" ; float
    > ttlpcp(time, nj, ni) ; ttlpcp:long~name~ = \"Total precipiation\"
    > ; ttlpcp:units = \"kg m-2 s-1\" ; ttlpcp:coordinates = \"LON LAT\"
    > ; ttlpcp:~FillValue~ = 1.e+30f ; ttlpcp:comment = \"derived from
    > bucket dumps\" ; float wndewd(time, nj, ni) ; wndewd:long~name~ =
    > \"x-ward winds\" ; wndewd:units = \"m/s\" ; wndewd:coordinates =
    > \"LON LAT\" ; wndewd:~FillValue~ = 1.e+30f ; float wndnwd(time,
    > nj, ni) ; wndnwd:long~name~ = \"y-ward winds\" ; wndnwd:units =
    > \"m/s\" ; wndnwd:coordinates = \"LON LAT\" ; wndnwd:~FillValue~ =
    > 1.e+30f ;
    >
    > // global attributes: :history = \"Mon Aug 5 15:31:32 2019: ncks
    > -v time,LON,LAT,glbrad,dlwsfc,wndewd,wndnwd,airtmp,spchmd,ttlpcp
    > JRA55~03hrforcing20050101~.00.netrad.nc
    > JRA55~03hrforcing20050101~.00.nc\" ; :NCO = \"netCDF Operators
    > version 4.7.5 (Homepage = <http://nco.sf.net>, Code =
    > <http://github.com/nco/nco>)\" ; }

4.  GX1 WOA

    Monthly statistics for silicate and nitrate are provided across two
    files. The header information for these two files is provided.

    [Nitrate](file:///Volumes/ioa03/cice-dirs/input/CICE_data/forcing/gx1/WOA/MONTHLY/nitrate_climatologyWOA_gx1v6f_20150107.nc)
    file:

    > netcdf nitrate~climatologyWOAgx1v6f20150107~ { dimensions: time =
    > 12 ; nj = 384 ; ni = 320 ; variables: float time(time) ; int
    > nj(nj) ; int ni(ni) ; float TLAT(nj, ni) ; TLAT:~FillValue~ =
    > 1.e+30f ; float TLON(nj, ni) ; TLON:~FillValue~ = 1.e+30f ; double
    > nitrate(time, nj, ni) ; nitrate:~FillValue~ = 1.00000001504747e+30
    > ; }

    [Silicate](file:///Volumes/ioa03/cice-dirs/input/CICE_data/forcing/gx1/WOA/MONTHLY/silicate_climatologyWOA_gx1v6f_20150107.nc)
    file:

    > netcdf silicate~climatologyWOAgx1v6f20150107~ { dimensions: time =
    > 12 ; nj = 384 ; ni = 320 ; variables: float time(time) ;
    > time:~FillValue~ = 9.96921e+36f ; int nj(nj) ; int ni(ni) ; float
    > TLAT(nj, ni) ; TLAT:~FillValue~ = 1.e+30f ; float TLON(nj, ni) ;
    > TLON:~FillValue~ = 1.e+30f ; float silicate(time, nj, ni) ;
    > silicate:~FillValue~ = 1.e+30f ; }

### GX1 Initial Conditions

These files are nearly identical to the
[*GX3 Initial Conditions*]{.spurious-link
target="GX3 Initial Conditions"} **other than the underlying grid**. So
there are a total of 12 initial condition files provided. With each file
containing NIL time step across the This is interesting to me as I
thought only one would be required. Nonetheless, the header of the first
file is given below.

> netcdf iced~gx1v6~.2005-01-01 { dimensions: ni = 320 ; nj = 384 ; ncat
> = 5 ; variables: double uvel(nj, ni) ; double vvel(nj, ni) ; double
> scale~factor~(nj, ni) ; double swvdr(nj, ni) ; double swvdf(nj, ni) ;
> double swidr(nj, ni) ; double swidf(nj, ni) ; double strocnxT(nj, ni)
> ; double strocnyT(nj, ni) ; double stressp~1~(nj, ni) ; double
> stressp~2~(nj, ni) ; double stressp~3~(nj, ni) ; double stressp~4~(nj,
> ni) ; double stressm~1~(nj, ni) ; double stressm~2~(nj, ni) ; double
> stressm~3~(nj, ni) ; double stressm~4~(nj, ni) ; double
> stress12~1~(nj, ni) ; double stress12~2~(nj, ni) ; double
> stress12~3~(nj, ni) ; double stress12~4~(nj, ni) ; double iceumask(nj,
> ni) ; double sst(nj, ni) ; double frzmlt(nj, ni) ; double
> frz~onset~(nj, ni) ; double fsnow(nj, ni) ; double aicen(ncat, nj, ni)
> ; double vicen(ncat, nj, ni) ; double vsnon(ncat, nj, ni) ; double
> Tsfcn(ncat, nj, ni) ; double iage(ncat, nj, ni) ; double FY(ncat, nj,
> ni) ; double alvl(ncat, nj, ni) ; double vlvl(ncat, nj, ni) ; double
> apnd(ncat, nj, ni) ; double hpnd(ncat, nj, ni) ; double ipnd(ncat, nj,
> ni) ; double dhs(ncat, nj, ni) ; double ffrac(ncat, nj, ni) ; double
> fbrn(ncat, nj, ni) ; double first~ice~(ncat, nj, ni) ; double
> sice001(ncat, nj, ni) ; double qice001(ncat, nj, ni) ; double
> sice002(ncat, nj, ni) ; double qice002(ncat, nj, ni) ; double
> sice003(ncat, nj, ni) ; double qice003(ncat, nj, ni) ; double
> sice004(ncat, nj, ni) ; double qice004(ncat, nj, ni) ; double
> sice005(ncat, nj, ni) ; double qice005(ncat, nj, ni) ; double
> sice006(ncat, nj, ni) ; double qice006(ncat, nj, ni) ; double
> sice007(ncat, nj, ni) ; double qice007(ncat, nj, ni) ; double
> qsno001(ncat, nj, ni) ;
>
> // global attributes: :istep1 = 87672 ; :time = 315619200. ;
> :time~forc~ = 315619200. ; :nyr = 11 ; :month = 1 ; :mday = 1 ; :sec =
> 0 ; :created = \"2021-03-30\" ; :history = \"Fri Apr 2 08:48:44 2021:
> ncatted --attribute created,global,c,c,2021-03-30
> iced~gx1v6~.2005-01-01.nc\" ; :NCO = \"netCDF Operators version 4.9.5
> (Homepage = <http://nco.sf.net>, Code = <http://github.com/nco/nco>)\"
> ; }

## TX1 Description

Global One-Degree Tri-polar

### TX1 Grid

Here this [grid
file](file:///Volumes/ioa03/cice-dirs/input/CICE_data/grid/tx1/grid_tx1.nc)
is given in both NetCDF and binary format, hence here is the header
contents of this file. This header contents indicates are fair amount
more information than the [*GX3 Grid*]{.spurious-link target="GX3 Grid"}
file.

> netcdf grid~tx1~ { dimensions: nx = 360 ; ny = 300 ; nc = 4 ;
> variables: double ulat(ny, nx) ; ulat:units = \"radians\" ; ulat:title
> = \"Latitude of U points\" ; double ulon(ny, nx) ; ulon:units =
> \"radians\" ; ulon:title = \"Longitude of U points\" ; double tlat(ny,
> nx) ; tlat:units = \"radians\" ; tlat:title = \"Latitude of T points\"
> ; double tlon(ny, nx) ; tlon:units = \"radians\" ; tlon:title =
> \"Longitude of T points\" ; double htn(ny, nx) ; htn:units = \"cm\" ;
> htn:title = \"Width of T cells on N side\" ; double hte(ny, nx) ;
> hte:units = \"cm\" ; hte:title = \"Width of T cells on E side\" ;
> double hun(ny, nx) ; hun:units = \"cm\" ; hun:title = \"Width of U
> cells on N side\" ; double hue(ny, nx) ; hue:units = \"cm\" ;
> hue:title = \"Width of U cells on E side\" ; double angle(ny, nx) ;
> angle:units = \"radians\" ; angle:title = \"Rotation Angle of U
> cells\" ; double angleT(ny, nx) ; angleT:units = \"radians\" ;
> angleT:title = \"Rotation Angle of T cells\" ; double tarea(ny, nx) ;
> tarea:units = \"m^2^\" ; tarea:title = \"area in m^2^ for T cells\" ;
> double uarea(ny, nx) ; uarea:units = \"m^2^\" ; uarea:title = \"area
> in m^2^ for U cells\" ; double lont~bonds~(nc, ny, nx) ;
> lont~bonds~:units = \"radians\" ; lont~bonds~:title = \"long. of T
> cell vertices\" ; double latt~bonds~(nc, ny, nx) ; latt~bonds~:units =
> \"radians\" ; latt~bonds~:title = \"lat. of T cell vertices\" ; double
> lonu~bonds~(nc, ny, nx) ; lonu~bonds~:units = \"radians\" ;
> lonu~bonds~:title = \"long. of U cell vertices\" ; double
> latu~bonds~(nc, ny, nx) ; latu~bonds~:units = \"radians\" ;
> latu~bonds~:title = \"lat. of U cell vertices\" ; }

### TX1 Bathymetry

The [bathymetry
file](file:///Volumes/ioa03/cice-dirs/input/CICE_data/grid/tx1/tx1_bathy.nc)
for the tri-polar grid is pretty straight forward:

> netcdf tx1~bathy~ { dimensions: dim~j~ = 240 ; dim~i~ = 360 ;
> variables: double Bathymetry(dim~j~, dim~i~) ; Bathymetry:long~name~ =
> \"ocean bathymetry for grounding scheme\" ; Bathymetry:units = \"m\" ;
> Bathymetry:~FillValue~ = 1.e+30 ; double lat(dim~j~, dim~i~) ;
> lat:long~name~ = \"Latitude\" ; lat:units = \"degrees~north~\" ;
> lat:~FillValue~ = 1.e+30 ; double lon(dim~j~, dim~i~) ; lon:long~name~
> = \"Longitude\" ; lon:units = \"degrees~east~\" ; lon:~FillValue~ =
> 1.e+30 ;
>
> // global attributes: :history = \"Wed Feb 10 21:25:30 2021: ncatted
> -a units,Bathymetry,o,c,m temptx1.nc tx1~bathy~.nc`\n`{=latex}\",
> \"Wed Feb 10 21:24:26 2021: ncatted -a long~name~,Bathymetry,o,c,ocean
> bathymetry for grounding scheme temptx1.nc temptx1.nc`\n`{=latex}\",
> \"Wed Feb 10 21:23:52 2021: ncks -x -v dx,dy,mask temptx1.nc
> temptx1.nc`\n`{=latex}\", \"Wed Feb 10 21:22:19 2021: ncrename -v
> ang,Bathymetry temptx1.nc\" ; :NCO = \"4.4.2\" ; }

### TX1 Land Mask

The [land mask
file](file:///Volumes/ioa03/cice-dirs/input/CICE_data/grid/tx1/kmt_tx1.bin)
is binary and cannot be easily shown here.

### TX1 Forcing

The only dataset utilised for the tri-polar forcing is JRA55. It is
equivalent in it\'s structure and content to the other two JRA55
([*GX3 JRA55*]{.spurious-link target="GX3 JRA55"} and
[*GX1 JRA55*]{.spurious-link target="GX1 JRA55"}), but just re-gridded
onto the [*TX1 Grid*]{.spurious-link target="TX1 Grid"}.

### TX1 Initial Conditions

There is no version 6 files provided [*CICE Consortium*]{.spurious-link
target="CICE Consortium"}

## Initial Conditions Discussion

All three (GX3, GX1 and TX1) cases of
[*CICE Consortium Test Data*]{.spurious-link
target="CICE Consortium Test Data"} have the same structure of the
initial condition files. Each case provides twelve files named for each
month. The contents of those files show 56 gridded parameters over two
and three dimensions (one categorical dimension and two spatial
dimensions). There is not a time dimension. Hence time-step is implied
as monthly with 12 files named for each month -- it is not known wether
this is a snap snot or an average over the month. The table below is
made in order to consolidate the initial condition parameters should at
some point it be determined that initial conditions should be provided.
The short name and description (the first two columns) are obtained from
correlating the [*CICE Consortium Test Data*]{.spurious-link
target="CICE Consortium Test Data"} information and the information
found in the documentation, online
[here](https://cice-consortium-cice.readthedocs.io/en/main/cice_index.html).

  CICE6 short         CICE6 description                                    Dims           Units
  ------------------- ---------------------------------------------------- -------------- --------
  uvel                sea ice velocity, eastward                           ni, nj         m/s
  vvel                sea ice velocity, northward                          ni, nj         m/s
  scale~factor~       timestep offsets in shortwave radiation components   ni, nj         \-
  swvdr               incoming shortwave radiation, visible direct         ni, nj         W/m^2^
  swvdf               incoming shortwave radiation, visible diffuse        ni, nj         W/m^2^
  swidr               incoming shortwave radiation, infrared direct        ni, nj         W/m^2^
  swidf               incoming shortwave radiation, infrared diffuse       ni, nj         W/m^2^
  strocnxT            ice-ocean stress                                     ni, nj         N/m
  strocnyT            ice-ocean stress                                     ni, nj         N/m
  stressp\_\[1-4\]    internal strain tensor                               ni, nj         N/m
  stressm\_\[1-4\]    internal strain tensor                               ni, nj         N/m
  stress12\_\[1-4\]   internal strain tensor                               ni, nj         N/m
  iceumask            the mask where there is ice in a gridcell            ni, nj         \-
  sst                 sea surface temperature                              ni, nj         C
  frzmlt              ice ocean heat exhange                               ni, nj         W/m^2^
  frz~onset~          Julian day in which the onset of freezing occurs     ni, nj         \-
  ulat                u-grid latitudes                                     ni, nj         \-
  ulon                u-grid longitudes                                    ni, nj         \-
  tlat                t-grid latitudes                                     ni, nj         \-
  tlon                t-grid latitudes                                     ni, nj         \-
  aicen               ice fraction per category                            ncat, ni, nj   \-
  vicen               ice volume per unit area per category                ncat, ni, nj   \-
  vsnon               snow volume per unit area per category               ncat, ni, nj   \-
  Tsfcn               surface snow/ice temperature per category            ncat, ni, nj   C
  iage                ice age                                              ncat, ni, nj   \-
  FY                  first year ice area                                  ncat, ni, nj   \-
  alvl                fraction of level ice area                           ncat, ni, nj   \-
  vlvl                volume of the level ice area                         ncat, ni, nj   \-
  apnd                fraction of ponds per cell                           ncat, ni, nj   \-
  hpnd                depth of ponds per cell                              ncat, ni, nj   \-
  ipnd                other pond characteristics                           ncat, ni, nj   \-
  dhs                 \-                                                   ncat, ni, nj   \-
  ffrac               \-                                                   ncat, ni, nj   \-
  fbrs                \-                                                   ncat, ni, nj   \-
  sice00\[1-5\]       \-                                                   ncat, ni, nj   \-
  qice00\[1-5\]       \-                                                   ncat, ni, nj   \-
  qsno00\[1-5\]       \-                                                   ncat, ni, nj   \-

# Climatological Datasets

The following external datasets are relevant to stand-alone
[*CICE6*]{.spurious-link target="CICE Consortium"}. Below is given
information that is helpful to reference for this purpose.

## COREv2

[COREv2](https://data1.gfdl.noaa.gov/nomads/forms/core/COREv2/CIAF_v2.html)
is a high quality, but now aging, global one degree and one-hour
atmospheric climate dataset with from years 1948-2009. I have downloaded
the following files:

> -rw-r--r-- 1 dpath2o staff 102M 24 Oct 2012
> ncar~precip~.1948-2009.23OCT2012.nc -rw-r--r-- 1 dpath2o staff 3.0G 24
> Oct 2012 ncar~rad~.1948-2009.23OCT2012.nc -rw-r--r-- 1 dpath2o staff
> 6.1G 24 Oct 2012 q~10~.1948-2009.23OCT2012.nc -rw-r--r-- 1 dpath2o
> staff 509K 18 Jun 2009 runoff.15JUNE2009.nc -rw-r--r-- 1 dpath2o staff
> 386M 26 Jul 2016 runoff.daitren.iaf.20120419.nc -rw-r--r-- 1 dpath2o
> staff 6.1G 24 Oct 2012 slp.1948-2009.23OCT2012.nc -rw-r--r-- 1 dpath2o
> staff 6.1G 24 Oct 2012 u~10~.1948-2009.23OCT2012.nc -rw-r--r-- 1
> dpath2o staff 6.1G 24 Oct 2012 v~10~.1948-2009.23OCT2012.nc

Header of u10 file:

> netcdf u~10~.1948-2009.23OCT2012 { dimensions: LAT = 94 ; LON = 192 ;
> TIME = UNLIMITED ; // (90520 currently) bnds = 2 ; variables: double
> LAT(LAT) ; LAT:units = \"degrees~north~\" ; LAT:point~spacing~ =
> \"uneven\" ; LAT:axis = \"Y\" ; double LON(LON) ; LON:units =
> \"degrees~east~\" ; LON:modulo = 360. ; LON:point~spacing~ = \"even\"
> ; LON:axis = \"X\" ; double TIME(TIME) ; TIME:units = \"days since
> 1948-01-01 00:00:00\" ; TIME:axis = \"T\" ; TIME:bounds =
> \"TIME~bnds~\" ; TIME:time~origin~ = \"1-JAN-1948\" ; TIME:calendar =
> \"NOLEAP\" ; double TIME~bnds~(TIME, bnds) ; float U~10MOD~(TIME, LAT,
> LON) ; U~10MOD~:missing~value~ = -1.e+34f ; U~10MOD~:~FillValue~ =
> -1.e+34f ; U~10MOD~:long~name~ = \"10m U Wind\" ; U~10MOD~:units =
> \"m/s\" ;
>
> // global attributes: :history = \"Thu Sep 27 07:10:36 2012: ncatted
> -O -a bounds,LAT,d,c,LAT~bnds~
> u~10mod~.1948-2009.nomads.nc`\n`{=latex}\", \"Thu Sep 27 07:10:08
> 2012: ncks -x -v LAT~bnds~ u~10mod~.1948-2009.new.nc
> u~10mod~.1948-2009.nomads.nc`\n`{=latex}\", \"FERRET V6.725
> 18-Sep-12\" ; :Conventions = \"CF-1.0\" ; :NCO = \"4.0.3\" ; }

## C-GLORS

The CMCC Global Ocean Physical Reanalysis System
([C-GLORS](http://c-glors.cmcc.it/index/index.html)) is used to simulate
the state of the ocean in the last decades. It consists of a variational
data assimilation system (OceanVar), capable of assimilating all in-situ
observations along with altimetry data, and a forecast step performed by
the ocean model NEMO coupled with the LIM2 sea-ice model.

I have downloaded the MLD files for 2010-2020 and here is the header
information:

> netcdf CGv7~2010MLD~ { dimensions: y = 1050 ; x = 1442 ; time~counter~
> = UNLIMITED ; // (365 currently) axis~nbounds~ = 2 ; variables: float
> nav~lat~(y, x) ; nav~lat~:standard~name~ = \"latitude\" ;
> nav~lat~:long~name~ = \"Latitude\" ; nav~lat~:units =
> \"degrees~north~\" ; nav~lat~:nav~model~ = \"grid~T~\" ; float
> nav~lon~(y, x) ; nav~lon~:standard~name~ = \"longitude\" ;
> nav~lon~:long~name~ = \"Longitude\" ; nav~lon~:units =
> \"degrees~east~\" ; nav~lon~:nav~model~ = \"grid~T~\" ; float
> somxl010(time~counter~, y, x) ; somxl010:standard~name~ =
> \"ocean~mixedlayerthicknessdefinedbysigmatheta~\" ;
> somxl010:long~name~ = \"Mixed Layer Depth (dsigma = 0.01 wrt 10m)\" ;
> somxl010:units = \"m\" ; somxl010:online~operation~ = \"average\" ;
> somxl010:interval~operation~ = \"1080 s\" ; somxl010:interval~write~ =
> \"1 d\" ; somxl010:cell~methods~ = \"time: mean (interval: 1080 s)\" ;
> somxl010:~FillValue~ = 1.e+20f ; somxl010:missing~value~ = 1.e+20f ;
> somxl010:coordinates = \"time~centered~ nav~lon~ nav~lat~\" ; double
> time~centered~(time~counter~) ; time~centered~:standard~name~ =
> \"time\" ; time~centered~:long~name~ = \"Time axis\" ;
> time~centered~:calendar = \"gregorian\" ; time~centered~:units =
> \"seconds since 1950-01-01 00:00:00\" ; time~centered~:time~origin~ =
> \"1950-01-01 00:00:00\" ; time~centered~:bounds =
> \"time~centeredbounds~\" ; double time~centeredbounds~(time~counter~,
> axis~nbounds~) ; double time~counter~(time~counter~) ;
> time~counter~:axis = \"T\" ; time~counter~:standard~name~ = \"time\" ;
> time~counter~:long~name~ = \"Time axis\" ; time~counter~:calendar =
> \"gregorian\" ; time~counter~:units = \"seconds since 1950-01-01
> 00:00:00\" ; time~counter~:time~origin~ = \"1950-01-01 00:00:00\" ;
> time~counter~:bounds = \"time~counterbounds~\" ; double
> time~counterbounds~(time~counter~, axis~nbounds~) ;
>
> // global attributes: :name = \"NEMO~1d2009122720100106~\" ;
> :description = \"ocean T grid variables\" ; :title = \"ocean T grid
> variables\" ; :Conventions = \"CF-1.5\" ; :production = \"An IPSL
> model\" ; :timeStamp = \"2016-Aug-14 13:15:38 CEST\" ; }

## ERA5

$1/4$ degree spatial and 1 hour temporal resolutions. A thorough
description of all the fields in ERA5 is given
[here](https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview).
[This is
also](https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation#ERA5:datadocumentation-Table3),
another good location to resource for information on the ERA5 fields.

[COSIMA GitHub Working
Group/Issue](https://github.com/COSIMA/access-om2/issues/242). For a
comprehensive listing on gadi of mapped ERA5 fields [see Andrew
Kiss\'s](https://github.com/COSIMA/access-om2/issues/242#issuecomment-913908910)
helpful table -- that thread also contains other helpful for information
re. ERA5 setup for CICE6.

-   [Latitudes and Longitudes (grid) in text (ascii)
    format](https://rda.ucar.edu/datasets/ds633.0/docs/ERA5.025deg_lats_lons.txt)
-   [Pressure Levels (levels) in text (ascii)
    format](https://rda.ucar.edu/datasets/ds633.0/docs/ERA5.025deg_std_pres_levels.txt)
-   [ERA5 Atmospheric Surface
    Analysis](https://rda.ucar.edu/datasets/ds633.0/docs/ds633.0.e5.oper.an.sfc.grib1.table.web.txt)
-   [ERA5 Atmospheric Pressure Level
    Analysis](https://rda.ucar.edu/datasets/ds633.0/docs/ds633.0.e5.oper.an.pl.grib1.table.web.txt)
-   [ERA5 Accummulated
    Fields](https://rda.ucar.edu/datasets/ds633.0/docs/ds633.0.e5.oper.fc.sfc.accumu.grib1.table.web.txt)
-   [ERA5 Instantaneous
    Fields](https://rda.ucar.edu/datasets/ds633.0/docs/ds633.0.e5.oper.fc.sfc.instan.grib1.table.web.txt)

The following ERA5 fields have been downloaded from gadi to my [local
repository](file:///Volumes/ioa02/reanalysis/ERA5) and organised in same
year sub-directory structure for the years 2010-2019.

### [10u](file:///Volumes/ioa02/reanalysis/ERA5/10u)

Ten-metre zonal surface winds

> netcdf \\10u~era5opersfc20100101~-20100131 { dimensions: longitude =
> 1440 ; latitude = 721 ; time = 744 ; variables: float
> longitude(longitude) ; longitude:units = \"degrees~east~\" ;
> longitude:long~name~ = \"longitude\" ; float latitude(latitude) ;
> latitude:units = \"degrees~north~\" ; latitude:long~name~ =
> \"latitude\" ; int time(time) ; time:units = \"hours since 1900-01-01
> 00:00:00.0\" ; time:long~name~ = \"time\" ; time:calendar =
> \"gregorian\" ; short u10(time, latitude, longitude) ;
> u10:scale~factor~ = 0.00113594702823864 ; u10:add~offset~ =
> 7.62408107433744 ; u10:~FillValue~ = -32767s ; u10:missing~value~ =
> -32767s ; u10:units = \"m s\*\*-1\" ; u10:long~name~ = \"10 metre U
> wind component\" ;
>
> // global attributes: :Conventions = \"CF-1.6\" ; :history =
> \"2020-09-28 12:42:00 UTC+1000 by era5~replicationtools~-1.2.1: mv
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/10u/2010*.10u~era5opersfc20100101~-20100131.tmp
> /g/data/rt52/era5/single-levels/reanalysis/10u/2010/10u~era5opersfc20100101~-20100131.nc`\n`{=latex}\",
> \"2020-09-28 12:40:58 UTC+1000 by era5~replicationtools~-1.2.1: nccopy
> -k4 -d5 -s
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/10u/2010/10u~era5opersfc20100101~-20100131.nc
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/10u/2010*.10u~era5opersfc20100101~-20100131.tmp`\n`{=latex}\",
> \"2020-09-28 12:36:04 UTC+1000 by era5~replicationtools~-1.2.1: curl
> --connect-timeout 20 --show-error --silent --max-time 36000 -o
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/10u/2010/10u~era5opersfc20100101~-20100131.nc
> [http://136.156.133.25/cache-compute-0008/cache/data0/adaptor.mars.internal-1601260358.5388222-19803-23-f152e830-9d4e-468b-a146-f2132486120c.nc\\n](http://136.156.133.25/cache-compute-0008/cache/data0/adaptor.mars.internal-1601260358.5388222-19803-23-f152e830-9d4e-468b-a146-f2132486120c.nc\n)\",
> \"2020-09-28 02:33:37 GMT by grib~tonetcdf~-2.16.0:
> /opt/ecmwf/eccodes/bin/grib~tonetcdf~ -S param -o
> /cache/data0/adaptor.mars.internal-1601260358.5388222-19803-23-f152e830-9d4e-468b-a146-f2132486120c.nc
> /cache/tmp/f152e830-9d4e-468b-a146-f2132486120c-adaptor.mars.internal-1601260358.5393667-19803-6-tmp.grib\"
> ; :license = \"Licence to use Copernicus Products:
> <https://apps.ecmwf.int/datasets/licences/copernicus/>\" ; :summary =
> \"ERA5 is the fifth generation ECMWF atmospheric reanalysis of the
> global climate. This file is part of the ERA5 replica hosted at NCI
> Australia. For more information please see
> <http://dx.doi.org/10.25914/5f48874388857>\" ; :title = \"ERA5
> single-levels reanalysis 10m~ucomponentofwind~ 20100101-20100131\" ; }

### [10v](file:///Volumes/ioa02/reanalysis/ERA5/10v)

Ten-metre meridional surface winds

> netcdf \\10v~era5opersfc20100101~-20100131 { dimensions: longitude =
> 1440 ; latitude = 721 ; time = 744 ; variables: float
> longitude(longitude) ; longitude:units = \"degrees~east~\" ;
> longitude:long~name~ = \"longitude\" ; float latitude(latitude) ;
> latitude:units = \"degrees~north~\" ; latitude:long~name~ =
> \"latitude\" ; int time(time) ; time:units = \"hours since 1900-01-01
> 00:00:00.0\" ; time:long~name~ = \"time\" ; time:calendar =
> \"gregorian\" ; short v10(time, latitude, longitude) ;
> v10:scale~factor~ = 0.00113620990606891 ; v10:add~offset~ =
> -6.70077104684757 ; v10:~FillValue~ = -32767s ; v10:missing~value~ =
> -32767s ; v10:units = \"m s\*\*-1\" ; v10:long~name~ = \"10 metre V
> wind component\" ;
>
> // global attributes: :Conventions = \"CF-1.6\" ; :history =
> \"2020-09-28 12:40:15 UTC+1000 by era5~replicationtools~-1.2.1: mv
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/10v/2010*.10v~era5opersfc20100101~-20100131.tmp
> /g/data/rt52/era5/single-levels/reanalysis/10v/2010/10v~era5opersfc20100101~-20100131.nc`\n`{=latex}\",
> \"2020-09-28 12:39:12 UTC+1000 by era5~replicationtools~-1.2.1: nccopy
> -k4 -d5 -s
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/10v/2010/10v~era5opersfc20100101~-20100131.nc
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/10v/2010*.10v~era5opersfc20100101~-20100131.tmp`\n`{=latex}\",
> \"2020-09-28 12:35:24 UTC+1000 by era5~replicationtools~-1.2.1: curl
> --connect-timeout 20 --show-error --silent --max-time 36000 -o
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/10v/2010/10v~era5opersfc20100101~-20100131.nc
> [http://136.156.132.153/cache-compute-0002/cache/data1/adaptor.mars.internal-1601260363.3011663-8591-19-7ca36757-cc7c-4890-9fc2-43150cac7864.nc\\n](http://136.156.132.153/cache-compute-0002/cache/data1/adaptor.mars.internal-1601260363.3011663-8591-19-7ca36757-cc7c-4890-9fc2-43150cac7864.nc\n)\",
> \"2020-09-28 02:33:37 GMT by grib~tonetcdf~-2.16.0:
> /opt/ecmwf/eccodes/bin/grib~tonetcdf~ -S param -o
> /cache/data1/adaptor.mars.internal-1601260363.3011663-8591-19-7ca36757-cc7c-4890-9fc2-43150cac7864.nc
> /cache/tmp/7ca36757-cc7c-4890-9fc2-43150cac7864-adaptor.mars.internal-1601260363.3017175-8591-7-tmp.grib\"
> ; :license = \"Licence to use Copernicus Products:
> <https://apps.ecmwf.int/datasets/licences/copernicus/>\" ; :summary =
> \"ERA5 is the fifth generation ECMWF atmospheric reanalysis of the
> global climate. This file is part of the ERA5 replica hosted at NCI
> Australia. For more information please see
> <http://dx.doi.org/10.25914/5f48874388857>\" ; :title = \"ERA5
> single-levels reanalysis 10m~vcomponentofwind~ 20100101-20100131\" ; }

### [2d](file:///Volumes/ioa02/reanalysis/ERA5/2d)

Two-metre dew-point temperature

> netcdf \\2d~era5opersfc20100101~-20100131 { dimensions: longitude =
> 1440 ; latitude = 721 ; time = 744 ; variables: float
> longitude(longitude) ; longitude:units = \"degrees~east~\" ;
> longitude:long~name~ = \"longitude\" ; float latitude(latitude) ;
> latitude:units = \"degrees~north~\" ; latitude:long~name~ =
> \"latitude\" ; int time(time) ; time:units = \"hours since 1900-01-01
> 00:00:00.0\" ; time:long~name~ = \"time\" ; time:calendar =
> \"gregorian\" ; short d2m(time, latitude, longitude) ;
> d2m:scale~factor~ = 0.00140529551830638 ; d2m:add~offset~ =
> 257.339492054389 ; d2m:~FillValue~ = -32767s ; d2m:missing~value~ =
> -32767s ; d2m:units = \"K\" ; d2m:long~name~ = \"2 metre dewpoint
> temperature\" ;
>
> // global attributes: :Conventions = \"CF-1.6\" ; :history =
> \"2020-09-28 12:43:43 UTC+1000 by era5~replicationtools~-1.2.1: mv
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/2d/2010*.2d~era5opersfc20100101~-20100131.tmp
> /g/data/rt52/era5/single-levels/reanalysis/2d/2010/2d~era5opersfc20100101~-20100131.nc`\n`{=latex}\",
> \"2020-09-28 12:42:44 UTC+1000 by era5~replicationtools~-1.2.1: nccopy
> -k4 -d5 -s
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/2d/2010/2d~era5opersfc20100101~-20100131.nc
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/2d/2010*.2d~era5opersfc20100101~-20100131.tmp`\n`{=latex}\",
> \"2020-09-28 12:37:40 UTC+1000 by era5~replicationtools~-1.2.1: curl
> --connect-timeout 20 --show-error --silent --max-time 36000 -o
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/2d/2010/2d~era5opersfc20100101~-20100131.nc
> [http://136.156.132.201/cache-compute-0004/cache/data5/adaptor.mars.internal-1601260378.023689-11063-21-0f2ccb24-24bb-4785-bfe1-620949bb0f34.nc\\n](http://136.156.132.201/cache-compute-0004/cache/data5/adaptor.mars.internal-1601260378.023689-11063-21-0f2ccb24-24bb-4785-bfe1-620949bb0f34.nc\n)\",
> \"2020-09-28 02:33:43 GMT by grib~tonetcdf~-2.16.0:
> /opt/ecmwf/eccodes/bin/grib~tonetcdf~ -S param -o
> /cache/data5/adaptor.mars.internal-1601260378.023689-11063-21-0f2ccb24-24bb-4785-bfe1-620949bb0f34.nc
> /cache/tmp/0f2ccb24-24bb-4785-bfe1-620949bb0f34-adaptor.mars.internal-1601260378.0242786-11063-3-tmp.grib\"
> ; :license = \"Licence to use Copernicus Products:
> <https://apps.ecmwf.int/datasets/licences/copernicus/>\" ; :summary =
> \"ERA5 is the fifth generation ECMWF atmospheric reanalysis of the
> global climate. This file is part of the ERA5 replica hosted at NCI
> Australia. For more information please see
> <http://dx.doi.org/10.25914/5f48874388857>\" ; :title = \"ERA5
> single-levels reanalysis 2m~dewpointtemperature~ 20100101-20100131\" ;
> }

### [2t](file:///Volumes/ioa02/reanalysis/ERA5/2t)

Two-metre air temperature

> netcdf \\2t~era5opersfc20100101~-20100131 { dimensions: longitude =
> 1440 ; latitude = 721 ; time = 744 ; variables: float
> longitude(longitude) ; longitude:units = \"degrees~east~\" ;
> longitude:long~name~ = \"longitude\" ; float latitude(latitude) ;
> latitude:units = \"degrees~north~\" ; latitude:long~name~ =
> \"latitude\" ; int time(time) ; time:units = \"hours since 1900-01-01
> 00:00:00.0\" ; time:long~name~ = \"time\" ; time:calendar =
> \"gregorian\" ; short t2m(time, latitude, longitude) ;
> t2m:scale~factor~ = 0.00163849236925778 ; t2m:add~offset~ =
> 268.255307157624 ; t2m:~FillValue~ = -32767s ; t2m:missing~value~ =
> -32767s ; t2m:units = \"K\" ; t2m:long~name~ = \"2 metre temperature\"
> ;
>
> // global attributes: :Conventions = \"CF-1.6\" ; :history =
> \"2020-09-28 12:42:00 UTC+1000 by era5~replicationtools~-1.2.1: mv
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/2t/2010*.2t~era5opersfc20100101~-20100131.tmp
> /g/data/rt52/era5/single-levels/reanalysis/2t/2010/2t~era5opersfc20100101~-20100131.nc`\n`{=latex}\",
> \"2020-09-28 12:41:06 UTC+1000 by era5~replicationtools~-1.2.1: nccopy
> -k4 -d5 -s
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/2t/2010/2t~era5opersfc20100101~-20100131.nc
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/2t/2010*.2t~era5opersfc20100101~-20100131.tmp`\n`{=latex}\",
> \"2020-09-28 12:37:30 UTC+1000 by era5~replicationtools~-1.2.1: curl
> --connect-timeout 20 --show-error --silent --max-time 36000 -o
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/2t/2010/2t~era5opersfc20100101~-20100131.nc
> [http://136.156.133.36/cache-compute-0010/cache/data6/adaptor.mars.internal-1601260376.2580519-20638-15-e6b74d84-3516-47c1-900e-48a2e567b012.nc\\n](http://136.156.133.36/cache-compute-0010/cache/data6/adaptor.mars.internal-1601260376.2580519-20638-15-e6b74d84-3516-47c1-900e-48a2e567b012.nc\n)\",
> \"2020-09-28 02:33:41 GMT by grib~tonetcdf~-2.16.0:
> /opt/ecmwf/eccodes/bin/grib~tonetcdf~ -S param -o
> /cache/data6/adaptor.mars.internal-1601260376.2580519-20638-15-e6b74d84-3516-47c1-900e-48a2e567b012.nc
> /cache/tmp/e6b74d84-3516-47c1-900e-48a2e567b012-adaptor.mars.internal-1601260376.258729-20638-3-tmp.grib\"
> ; :license = \"Licence to use Copernicus Products:
> <https://apps.ecmwf.int/datasets/licences/copernicus/>\" ; :summary =
> \"ERA5 is the fifth generation ECMWF atmospheric reanalysis of the
> global climate. This file is part of the ERA5 replica hosted at NCI
> Australia. For more information please see
> <http://dx.doi.org/10.25914/5f48874388857>\" ; :title = \"ERA5
> single-levels reanalysis 2m~temperature~ 20100101-20100131\" ; }

### [metss](file:///Volumes/ioa02/reanalysis/ERA5/metss)

Mean eastward turbulent surface stress

> netcdf metss~era5opersfc20100101~-20100131 { dimensions: longitude =
> 1440 ; latitude = 721 ; time = 744 ; variables: float
> longitude(longitude) ; longitude:units = \"degrees~east~\" ;
> longitude:long~name~ = \"longitude\" ; float latitude(latitude) ;
> latitude:units = \"degrees~north~\" ; latitude:long~name~ =
> \"latitude\" ; int time(time) ; time:units = \"hours since 1900-01-01
> 00:00:00.0\" ; time:long~name~ = \"time\" ; time:calendar =
> \"gregorian\" ; short metss(time, latitude, longitude) ;
> metss:scale~factor~ = 0.000366911687951802 ; metss:add~offset~ =
> -0.80675398578172 ; metss:~FillValue~ = -32767s ; metss:missing~value~
> = -32767s ; metss:units = \"N m\*\*-2\" ; metss:long~name~ = \"Mean
> eastward turbulent surface stress\" ;
>
> // global attributes: :Conventions = \"CF-1.6\" ; :history =
> \"2020-10-06 09:10:44 UTC+1100 by era5~replicationtools~-1.2.1: mv
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/metss/2010*.metss~era5opersfc20100101~-20100131.tmp
> /g/data/rt52/era5/single-levels/reanalysis/metss/2010/metss~era5opersfc20100101~-20100131.nc`\n`{=latex}\",
> \"2020-10-06 09:09:59 UTC+1100 by era5~replicationtools~-1.2.1: nccopy
> -k4 -d5 -s
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/metss/2010/metss~era5opersfc20100101~-20100131.nc
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/metss/2010*.metss~era5opersfc20100101~-20100131.tmp`\n`{=latex}\",
> \"2020-10-06 09:04:24 UTC+1100 by era5~replicationtools~-1.2.1: curl
> --connect-timeout 20 --show-error --silent --max-time 36000 -o
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/metss/2010/metss~era5opersfc20100101~-20100131.nc
> [http://136.156.132.110/cache-compute-0001/cache/data2/adaptor.mars.internal-1601935170.1456351-18982-33-a6b9dca0-74b0-4ce5-8017-3660c3fbe2d5.nc\\n](http://136.156.132.110/cache-compute-0001/cache/data2/adaptor.mars.internal-1601935170.1456351-18982-33-a6b9dca0-74b0-4ce5-8017-3660c3fbe2d5.nc\n)\",
> \"2020-10-05 22:00:31 GMT by grib~tonetcdf~-2.16.0:
> /opt/ecmwf/eccodes/bin/grib~tonetcdf~ -S param -o
> /cache/data2/adaptor.mars.internal-1601935170.1456351-18982-33-a6b9dca0-74b0-4ce5-8017-3660c3fbe2d5.nc
> /cache/tmp/a6b9dca0-74b0-4ce5-8017-3660c3fbe2d5-adaptor.mars.internal-1601935170.1462164-18982-8-tmp.grib\"
> ; :license = \"Licence to use Copernicus Products:
> <https://apps.ecmwf.int/datasets/licences/copernicus/>\" ; :summary =
> \"ERA5 is the fifth generation ECMWF atmospheric reanalysis of the
> global climate. This file is part of the ERA5 replica hosted at NCI
> Australia. For more information please see
> <http://dx.doi.org/10.25914/5f48874388857>\" ; :title = \"ERA5
> single-levels reanalysis mean~eastwardturbulentsurfacestress~
> 20100101-20100131\" ; }

### [mntss](file:///Volumes/ioa02/reanalysis/ERA5/mntss)

Mean northward turbulent surface stress

> netcdf mntss~era5opersfc20100101~-20100131 { dimensions: longitude =
> 1440 ; latitude = 721 ; time = 744 ; variables: float
> longitude(longitude) ; longitude:units = \"degrees~east~\" ;
> longitude:long~name~ = \"longitude\" ; float latitude(latitude) ;
> latitude:units = \"degrees~north~\" ; latitude:long~name~ =
> \"latitude\" ; int time(time) ; time:units = \"hours since 1900-01-01
> 00:00:00.0\" ; time:long~name~ = \"time\" ; time:calendar =
> \"gregorian\" ; short mntss(time, latitude, longitude) ;
> mntss:scale~factor~ = 0.000288605498999811 ; mntss:add~offset~ =
> 0.0627023578645138 ; mntss:~FillValue~ = -32767s ;
> mntss:missing~value~ = -32767s ; mntss:units = \"N m\*\*-2\" ;
> mntss:long~name~ = \"Mean northward turbulent surface stress\" ;
>
> // global attributes: :Conventions = \"CF-1.6\" ; :history =
> \"2020-10-06 09:07:45 UTC+1100 by era5~replicationtools~-1.2.1: mv
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/mntss/2010*.mntss~era5opersfc20100101~-20100131.tmp
> /g/data/rt52/era5/single-levels/reanalysis/mntss/2010/mntss~era5opersfc20100101~-20100131.nc`\n`{=latex}\",
> \"2020-10-06 09:06:59 UTC+1100 by era5~replicationtools~-1.2.1: nccopy
> -k4 -d5 -s
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/mntss/2010/mntss~era5opersfc20100101~-20100131.nc
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/mntss/2010*.mntss~era5opersfc20100101~-20100131.tmp`\n`{=latex}\",
> \"2020-10-06 09:03:25 UTC+1100 by era5~replicationtools~-1.2.1: curl
> --connect-timeout 20 --show-error --silent --max-time 36000 -o
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/mntss/2010/mntss~era5opersfc20100101~-20100131.nc
> [http://136.156.133.37/cache-compute-0011/cache/data9/adaptor.mars.internal-1601935174.829866-29967-28-db022499-88a8-4c1e-9b82-746968321737.nc\\n](http://136.156.133.37/cache-compute-0011/cache/data9/adaptor.mars.internal-1601935174.829866-29967-28-db022499-88a8-4c1e-9b82-746968321737.nc\n)\",
> \"2020-10-05 22:00:30 GMT by grib~tonetcdf~-2.16.0:
> /opt/ecmwf/eccodes/bin/grib~tonetcdf~ -S param -o
> /cache/data9/adaptor.mars.internal-1601935174.829866-29967-28-db022499-88a8-4c1e-9b82-746968321737.nc
> /cache/tmp/db022499-88a8-4c1e-9b82-746968321737-adaptor.mars.internal-1601935174.8304174-29967-7-tmp.grib\"
> ; :license = \"Licence to use Copernicus Products:
> <https://apps.ecmwf.int/datasets/licences/copernicus/>\" ; :summary =
> \"ERA5 is the fifth generation ECMWF atmospheric reanalysis of the
> global climate. This file is part of the ERA5 replica hosted at NCI
> Australia. For more information please see
> <http://dx.doi.org/10.25914/5f48874388857>\" ; :title = \"ERA5
> single-levels reanalysis mean~northwardturbulentsurfacestress~
> 20100101-20100131\" ; }

### [mror](file:///Volumes/ioa02/reanalysis/ERA5/mror)

Mean runoff rate

> netcdf mror~era5opersfc20100101~-20100131 { dimensions: longitude =
> 1440 ; latitude = 721 ; time = 744 ; variables: float
> longitude(longitude) ; longitude:units = \"degrees~east~\" ;
> longitude:long~name~ = \"longitude\" ; float latitude(latitude) ;
> latitude:units = \"degrees~north~\" ; latitude:long~name~ =
> \"latitude\" ; int time(time) ; time:units = \"hours since 1900-01-01
> 00:00:00.0\" ; time:long~name~ = \"time\" ; time:calendar =
> \"gregorian\" ; short mror(time, latitude, longitude) ;
> mror:scale~factor~ = 3.5196145727257e-07 ; mror:add~offset~ =
> 0.011532369108993 ; mror:~FillValue~ = -32767s ; mror:missing~value~ =
> -32767s ; mror:units = \"kg m\*\*-2 s\*\*-1\" ; mror:long~name~ =
> \"Mean runoff rate\" ;
>
> // global attributes: :Conventions = \"CF-1.6\" ; :history =
> \"2020-10-06 09:09:27 UTC+1100 by era5~replicationtools~-1.2.1: mv
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/mror/2010*.mror~era5opersfc20100101~-20100131.tmp
> /g/data/rt52/era5/single-levels/reanalysis/mror/2010/mror~era5opersfc20100101~-20100131.nc`\n`{=latex}\",
> \"2020-10-06 09:09:11 UTC+1100 by era5~replicationtools~-1.2.1: nccopy
> -k4 -d5 -s
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/mror/2010/mror~era5opersfc20100101~-20100131.nc
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/mror/2010*.mror~era5opersfc20100101~-20100131.tmp`\n`{=latex}\",
> \"2020-10-06 09:04:54 UTC+1100 by era5~replicationtools~-1.2.1: curl
> --connect-timeout 20 --show-error --silent --max-time 36000 -o
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/mror/2010/mror~era5opersfc20100101~-20100131.nc
> [http://136.156.132.198/cache-compute-0003/cache/data2/adaptor.mars.internal-1601935221.2261062-30150-20-349424a8-a7ab-4eae-8d46-b300f17ce397.nc\\n](http://136.156.132.198/cache-compute-0003/cache/data2/adaptor.mars.internal-1601935221.2261062-30150-20-349424a8-a7ab-4eae-8d46-b300f17ce397.nc\n)\",
> \"2020-10-05 22:01:06 GMT by grib~tonetcdf~-2.16.0:
> /opt/ecmwf/eccodes/bin/grib~tonetcdf~ -S param -o
> /cache/data2/adaptor.mars.internal-1601935221.2261062-30150-20-349424a8-a7ab-4eae-8d46-b300f17ce397.nc
> /cache/tmp/349424a8-a7ab-4eae-8d46-b300f17ce397-adaptor.mars.internal-1601935221.2266717-30150-7-tmp.grib\"
> ; :license = \"Licence to use Copernicus Products:
> <https://apps.ecmwf.int/datasets/licences/copernicus/>\" ; :summary =
> \"ERA5 is the fifth generation ECMWF atmospheric reanalysis of the
> global climate. This file is part of the ERA5 replica hosted at NCI
> Australia. For more information please see
> <http://dx.doi.org/10.25914/5f48874388857>\" ; :title = \"ERA5
> single-levels reanalysis mean~runoffrate~ 20100101-20100131\" ; }

### [msdrswrf](file:///Volumes/ioa02/reanalysis/ERA5/msdrswrf)

Mean surface direct short-wave radiation flux

> netcdf msdrswrf~era5opersfc20100101~-20100131 { dimensions: longitude
> = 1440 ; latitude = 721 ; time = 744 ; variables: float
> longitude(longitude) ; longitude:units = \"degrees~east~\" ;
> longitude:long~name~ = \"longitude\" ; float latitude(latitude) ;
> latitude:units = \"degrees~north~\" ; latitude:long~name~ =
> \"latitude\" ; int time(time) ; time:units = \"hours since 1900-01-01
> 00:00:00.0\" ; time:long~name~ = \"time\" ; time:calendar =
> \"gregorian\" ; short msdrswrf(time, latitude, longitude) ;
> msdrswrf:scale~factor~ = 0.0188225779378328 ; msdrswrf:add~offset~ =
> 616.740588711031 ; msdrswrf:~FillValue~ = -32767s ;
> msdrswrf:missing~value~ = -32767s ; msdrswrf:units = \"W m\*\*-2\" ;
> msdrswrf:long~name~ = \"Mean surface direct short-wave radiation
> flux\" ;
>
> // global attributes: :Conventions = \"CF-1.6\" ; :history =
> \"2020-10-06 09:14:19 UTC+1100 by era5~replicationtools~-1.2.1: mv
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msdrswrf/2010*.msdrswrf~era5opersfc20100101~-20100131.tmp
> /g/data/rt52/era5/single-levels/reanalysis/msdrswrf/2010/msdrswrf~era5opersfc20100101~-20100131.nc`\n`{=latex}\",
> \"2020-10-06 09:13:44 UTC+1100 by era5~replicationtools~-1.2.1: nccopy
> -k4 -d5 -s
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msdrswrf/2010/msdrswrf~era5opersfc20100101~-20100131.nc
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msdrswrf/2010*.msdrswrf~era5opersfc20100101~-20100131.tmp`\n`{=latex}\",
> \"2020-10-06 09:08:14 UTC+1100 by era5~replicationtools~-1.2.1: curl
> --connect-timeout 20 --show-error --silent --max-time 36000 -o
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msdrswrf/2010/msdrswrf~era5opersfc20100101~-20100131.nc
> [http://136.156.133.36/cache-compute-0010/cache/data1/adaptor.mars.internal-1601935408.957139-15597-14-b57f2409-d32e-448d-bf0e-ee6d5e9e4503.nc\\n](http://136.156.133.36/cache-compute-0010/cache/data1/adaptor.mars.internal-1601935408.957139-15597-14-b57f2409-d32e-448d-bf0e-ee6d5e9e4503.nc\n)\",
> \"2020-10-05 22:04:10 GMT by grib~tonetcdf~-2.16.0:
> /opt/ecmwf/eccodes/bin/grib~tonetcdf~ -S param -o
> /cache/data1/adaptor.mars.internal-1601935408.957139-15597-14-b57f2409-d32e-448d-bf0e-ee6d5e9e4503.nc
> /cache/tmp/b57f2409-d32e-448d-bf0e-ee6d5e9e4503-adaptor.mars.internal-1601935408.9578357-15597-6-tmp.grib\"
> ; :license = \"Licence to use Copernicus Products:
> <https://apps.ecmwf.int/datasets/licences/copernicus/>\" ; :summary =
> \"ERA5 is the fifth generation ECMWF atmospheric reanalysis of the
> global climate. This file is part of the ERA5 replica hosted at NCI
> Australia. For more information please see
> <http://dx.doi.org/10.25914/5f48874388857>\" ; :title = \"ERA5
> single-levels reanalysis mean~surfacedirectshortwaveradiationflux~
> 20100101-20100131\" ; }

### [msdrswrfcs](file:///Volumes/ioa02/reanalysis/ERA5/msdrswrfcs)

Mean surface direct short-wave radiation flux, clear sky

> netcdf msdrswrfcs~era5opersfc20100101~-20100131 { dimensions:
> longitude = 1440 ; latitude = 721 ; time = 744 ; variables: float
> longitude(longitude) ; longitude:units = \"degrees~east~\" ;
> longitude:long~name~ = \"longitude\" ; float latitude(latitude) ;
> latitude:units = \"degrees~north~\" ; latitude:long~name~ =
> \"latitude\" ; int time(time) ; time:units = \"hours since 1900-01-01
> 00:00:00.0\" ; time:long~name~ = \"time\" ; time:calendar =
> \"gregorian\" ; short msdrswrfcs(time, latitude, longitude) ;
> msdrswrfcs:scale~factor~ = 0.0182069529855187 ; msdrswrfcs:add~offset~
> = 596.569021523507 ; msdrswrfcs:~FillValue~ = -32767s ;
> msdrswrfcs:missing~value~ = -32767s ; msdrswrfcs:units = \"W m\*\*-2\"
> ; msdrswrfcs:long~name~ = \"Mean surface direct short-wave radiation
> flux, clear sky\" ;
>
> // global attributes: :Conventions = \"CF-1.6\" ; :history =
> \"2020-10-06 09:12:53 UTC+1100 by era5~replicationtools~-1.2.1: mv
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msdrswrfcs/2010*.msdrswrfcs~era5opersfc20100101~-20100131.tmp
> /g/data/rt52/era5/single-levels/reanalysis/msdrswrfcs/2010/msdrswrfcs~era5opersfc20100101~-20100131.nc`\n`{=latex}\",
> \"2020-10-06 09:12:21 UTC+1100 by era5~replicationtools~-1.2.1: nccopy
> -k4 -d5 -s
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msdrswrfcs/2010/msdrswrfcs~era5opersfc20100101~-20100131.nc
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msdrswrfcs/2010*.msdrswrfcs~era5opersfc20100101~-20100131.tmp`\n`{=latex}\",
> \"2020-10-06 09:08:21 UTC+1100 by era5~replicationtools~-1.2.1: curl
> --connect-timeout 20 --show-error --silent --max-time 36000 -o
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msdrswrfcs/2010/msdrswrfcs~era5opersfc20100101~-20100131.nc
> [http://136.156.132.153/cache-compute-0002/cache/data3/adaptor.mars.internal-1601935426.7869933-12103-12-b44cb567-fe87-4047-bde0-d85c747fde81.nc\\n](http://136.156.132.153/cache-compute-0002/cache/data3/adaptor.mars.internal-1601935426.7869933-12103-12-b44cb567-fe87-4047-bde0-d85c747fde81.nc\n)\",
> \"2020-10-05 22:04:43 GMT by grib~tonetcdf~-2.16.0:
> /opt/ecmwf/eccodes/bin/grib~tonetcdf~ -S param -o
> /cache/data3/adaptor.mars.internal-1601935426.7869933-12103-12-b44cb567-fe87-4047-bde0-d85c747fde81.nc
> /cache/tmp/b44cb567-fe87-4047-bde0-d85c747fde81-adaptor.mars.internal-1601935426.787668-12103-5-tmp.grib\"
> ; :license = \"Licence to use Copernicus Products:
> <https://apps.ecmwf.int/datasets/licences/copernicus/>\" ; :summary =
> \"ERA5 is the fifth generation ECMWF atmospheric reanalysis of the
> global climate. This file is part of the ERA5 replica hosted at NCI
> Australia. For more information please see
> <http://dx.doi.org/10.25914/5f48874388857>\" ; :title = \"ERA5
> single-levels reanalysis
> mean~surfacedirectshortwaveradiationfluxclearsky~ 20100101-20100131\"
> ; }

### [msdwlwrf](file:///Volumes/ioa02/reanalysis/ERA5/msdwlwrf)

Mean surface downward long-wave radiation flux

> netcdf msdwlwrf~era5opersfc20100101~-20100131 { dimensions: longitude
> = 1440 ; latitude = 721 ; time = 744 ; variables: float
> longitude(longitude) ; longitude:units = \"degrees~east~\" ;
> longitude:long~name~ = \"longitude\" ; float latitude(latitude) ;
> latitude:units = \"degrees~north~\" ; latitude:long~name~ =
> \"latitude\" ; int time(time) ; time:units = \"hours since 1900-01-01
> 00:00:00.0\" ; time:long~name~ = \"time\" ; time:calendar =
> \"gregorian\" ; short msdwlwrf(time, latitude, longitude) ;
> msdwlwrf:scale~factor~ = 0.00672298478721479 ; msdwlwrf:add~offset~ =
> 289.4882766912 ; msdwlwrf:~FillValue~ = -32767s ;
> msdwlwrf:missing~value~ = -32767s ; msdwlwrf:units = \"W m\*\*-2\" ;
> msdwlwrf:long~name~ = \"Mean surface downward long-wave radiation
> flux\" ;
>
> // global attributes: :Conventions = \"CF-1.6\" ; :history =
> \"2020-10-06 09:09:32 UTC+1100 by era5~replicationtools~-1.2.1: mv
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msdwlwrf/2010*.msdwlwrf~era5opersfc20100101~-20100131.tmp
> /g/data/rt52/era5/single-levels/reanalysis/msdwlwrf/2010/msdwlwrf~era5opersfc20100101~-20100131.nc`\n`{=latex}\",
> \"2020-10-06 09:08:27 UTC+1100 by era5~replicationtools~-1.2.1: nccopy
> -k4 -d5 -s
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msdwlwrf/2010/msdwlwrf~era5opersfc20100101~-20100131.nc
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msdwlwrf/2010*.msdwlwrf~era5opersfc20100101~-20100131.tmp`\n`{=latex}\",
> \"2020-10-06 09:01:48 UTC+1100 by era5~replicationtools~-1.2.1: curl
> --connect-timeout 20 --show-error --silent --max-time 36000 -o
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msdwlwrf/2010/msdwlwrf~era5opersfc20100101~-20100131.nc
> [http://136.156.133.46/cache-compute-0015/cache/data3/adaptor.mars.internal-1601935022.9836493-31252-13-a31306dc-69ab-4a06-bc9e-ad88c60a42ca.nc\\n](http://136.156.133.46/cache-compute-0015/cache/data3/adaptor.mars.internal-1601935022.9836493-31252-13-a31306dc-69ab-4a06-bc9e-ad88c60a42ca.nc\n)\",
> \"2020-10-05 21:57:45 GMT by grib~tonetcdf~-2.16.0:
> /opt/ecmwf/eccodes/bin/grib~tonetcdf~ -S param -o
> /cache/data3/adaptor.mars.internal-1601935022.9836493-31252-13-a31306dc-69ab-4a06-bc9e-ad88c60a42ca.nc
> /cache/tmp/a31306dc-69ab-4a06-bc9e-ad88c60a42ca-adaptor.mars.internal-1601935022.9842257-31252-5-tmp.grib\"
> ; :license = \"Licence to use Copernicus Products:
> <https://apps.ecmwf.int/datasets/licences/copernicus/>\" ; :summary =
> \"ERA5 is the fifth generation ECMWF atmospheric reanalysis of the
> global climate. This file is part of the ERA5 replica hosted at NCI
> Australia. For more information please see
> <http://dx.doi.org/10.25914/5f48874388857>\" ; :title = \"ERA5
> single-levels reanalysis mean~surfacedownwardlongwaveradiationflux~
> 20100101-20100131\" ; }

### [msdwlwrfcs](file:///Volumes/ioa02/reanalysis/ERA5/msdwlwrfcs)

Mean surface downward long-wave radiation flux, clear sky

> netcdf msdwlwrfcs~era5opersfc20100101~-20100131 { dimensions:
> longitude = 1440 ; latitude = 721 ; time = 744 ; variables: float
> longitude(longitude) ; longitude:units = \"degrees~east~\" ;
> longitude:long~name~ = \"longitude\" ; float latitude(latitude) ;
> latitude:units = \"degrees~north~\" ; latitude:long~name~ =
> \"latitude\" ; int time(time) ; time:units = \"hours since 1900-01-01
> 00:00:00.0\" ; time:long~name~ = \"time\" ; time:calendar =
> \"gregorian\" ; short msdwlwrfcs(time, latitude, longitude) ;
> msdwlwrfcs:scale~factor~ = 0.00657778122565039 ;
> msdwlwrfcs:add~offset~ = 287.159476612317 ; msdwlwrfcs:~FillValue~ =
> -32767s ; msdwlwrfcs:missing~value~ = -32767s ; msdwlwrfcs:units = \"W
> m\*\*-2\" ; msdwlwrfcs:long~name~ = \"Mean surface downward long-wave
> radiation flux, clear sky\" ;
>
> // global attributes: :Conventions = \"CF-1.6\" ; :history =
> \"2020-10-06 09:16:35 UTC+1100 by era5~replicationtools~-1.2.1: mv
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msdwlwrfcs/2010*.msdwlwrfcs~era5opersfc20100101~-20100131.tmp
> /g/data/rt52/era5/single-levels/reanalysis/msdwlwrfcs/2010/msdwlwrfcs~era5opersfc20100101~-20100131.nc`\n`{=latex}\",
> \"2020-10-06 09:15:39 UTC+1100 by era5~replicationtools~-1.2.1: nccopy
> -k4 -d5 -s
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msdwlwrfcs/2010/msdwlwrfcs~era5opersfc20100101~-20100131.nc
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msdwlwrfcs/2010*.msdwlwrfcs~era5opersfc20100101~-20100131.tmp`\n`{=latex}\",
> \"2020-10-06 09:09:11 UTC+1100 by era5~replicationtools~-1.2.1: curl
> --connect-timeout 20 --show-error --silent --max-time 36000 -o
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msdwlwrfcs/2010/msdwlwrfcs~era5opersfc20100101~-20100131.nc
> [http://136.156.132.105/cache-compute-0000/cache/data4/adaptor.mars.internal-1601935467.8587792-11179-31-6f536992-db19-476c-8205-d5ebf4a7e4c2.nc\\n](http://136.156.132.105/cache-compute-0000/cache/data4/adaptor.mars.internal-1601935467.8587792-11179-31-6f536992-db19-476c-8205-d5ebf4a7e4c2.nc\n)\",
> \"2020-10-05 22:05:11 GMT by grib~tonetcdf~-2.16.0:
> /opt/ecmwf/eccodes/bin/grib~tonetcdf~ -S param -o
> /cache/data4/adaptor.mars.internal-1601935467.8587792-11179-31-6f536992-db19-476c-8205-d5ebf4a7e4c2.nc
> /cache/tmp/6f536992-db19-476c-8205-d5ebf4a7e4c2-adaptor.mars.internal-1601935467.859394-11179-9-tmp.grib\"
> ; :license = \"Licence to use Copernicus Products:
> <https://apps.ecmwf.int/datasets/licences/copernicus/>\" ; :summary =
> \"ERA5 is the fifth generation ECMWF atmospheric reanalysis of the
> global climate. This file is part of the ERA5 replica hosted at NCI
> Australia. For more information please see
> <http://dx.doi.org/10.25914/5f48874388857>\" ; :title = \"ERA5
> single-levels reanalysis
> mean~surfacedownwardlongwaveradiationfluxclearsky~ 20100101-20100131\"
> ; }

### [msdwswrf](file:///Volumes/ioa02/reanalysis/ERA5/msdwswrf)

Mean surface downward short-wave radiation flux

> netcdf msdwswrf~era5opersfc20100101~-20100131 { dimensions: longitude
> = 1440 ; latitude = 721 ; time = 744 ; variables: float
> longitude(longitude) ; longitude:units = \"degrees~east~\" ;
> longitude:long~name~ = \"longitude\" ; float latitude(latitude) ;
> latitude:units = \"degrees~north~\" ; latitude:long~name~ =
> \"latitude\" ; int time(time) ; time:units = \"hours since 1900-01-01
> 00:00:00.0\" ; time:long~name~ = \"time\" ; time:calendar =
> \"gregorian\" ; short msdwswrf(time, latitude, longitude) ;
> msdwswrf:scale~factor~ = 0.0201353707292509 ; msdwswrf:add~offset~ =
> 659.755557314635 ; msdwswrf:~FillValue~ = -32767s ;
> msdwswrf:missing~value~ = -32767s ; msdwswrf:units = \"W m\*\*-2\" ;
> msdwswrf:long~name~ = \"Mean surface downward short-wave radiation
> flux\" ;
>
> // global attributes: :Conventions = \"CF-1.6\" ; :history =
> \"2020-10-06 09:05:31 UTC+1100 by era5~replicationtools~-1.2.1: mv
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msdwswrf/2010*.msdwswrf~era5opersfc20100101~-20100131.tmp
> /g/data/rt52/era5/single-levels/reanalysis/msdwswrf/2010/msdwswrf~era5opersfc20100101~-20100131.nc`\n`{=latex}\",
> \"2020-10-06 09:04:54 UTC+1100 by era5~replicationtools~-1.2.1: nccopy
> -k4 -d5 -s
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msdwswrf/2010/msdwswrf~era5opersfc20100101~-20100131.nc
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msdwswrf/2010*.msdwswrf~era5opersfc20100101~-20100131.tmp`\n`{=latex}\",
> \"2020-10-06 09:01:28 UTC+1100 by era5~replicationtools~-1.2.1: curl
> --connect-timeout 20 --show-error --silent --max-time 36000 -o
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msdwswrf/2010/msdwswrf~era5opersfc20100101~-20100131.nc
> [http://136.156.132.105/cache-compute-0000/cache/data1/adaptor.mars.internal-1601935012.7223494-5877-21-fad3c536-1dd2-422a-9b5d-47a7db1e9840.nc\\n](http://136.156.132.105/cache-compute-0000/cache/data1/adaptor.mars.internal-1601935012.7223494-5877-21-fad3c536-1dd2-422a-9b5d-47a7db1e9840.nc\n)\",
> \"2020-10-05 21:57:36 GMT by grib~tonetcdf~-2.16.0:
> /opt/ecmwf/eccodes/bin/grib~tonetcdf~ -S param -o
> /cache/data1/adaptor.mars.internal-1601935012.7223494-5877-21-fad3c536-1dd2-422a-9b5d-47a7db1e9840.nc
> /cache/tmp/fad3c536-1dd2-422a-9b5d-47a7db1e9840-adaptor.mars.internal-1601935012.7229025-5877-7-tmp.grib\"
> ; :license = \"Licence to use Copernicus Products:
> <https://apps.ecmwf.int/datasets/licences/copernicus/>\" ; :summary =
> \"ERA5 is the fifth generation ECMWF atmospheric reanalysis of the
> global climate. This file is part of the ERA5 replica hosted at NCI
> Australia. For more information please see
> <http://dx.doi.org/10.25914/5f48874388857>\" ; :title = \"ERA5
> single-levels reanalysis mean~surfacedownwardshortwaveradiationflux~
> 20100101-20100131\" ; }

### [msdwswrfcs](file:///Volumes/ioa02/reanalysis/ERA5/msdwswrfcs)

Mean surface downward short-wave radiation flux, clear sky

> netcdf msdwswrfcs~era5opersfc20100101~-20100131 { dimensions:
> longitude = 1440 ; latitude = 721 ; time = 744 ; variables: float
> longitude(longitude) ; longitude:units = \"degrees~east~\" ;
> longitude:long~name~ = \"longitude\" ; float latitude(latitude) ;
> latitude:units = \"degrees~north~\" ; latitude:long~name~ =
> \"latitude\" ; int time(time) ; time:units = \"hours since 1900-01-01
> 00:00:00.0\" ; time:long~name~ = \"time\" ; time:calendar =
> \"gregorian\" ; short msdwswrfcs(time, latitude, longitude) ;
> msdwswrfcs:scale~factor~ = 0.0194467863519143 ; msdwswrfcs:add~offset~
> = 637.193401606824 ; msdwswrfcs:~FillValue~ = -32767s ;
> msdwswrfcs:missing~value~ = -32767s ; msdwswrfcs:units = \"W m\*\*-2\"
> ; msdwswrfcs:long~name~ = \"Mean surface downward short-wave radiation
> flux, clear sky\" ;
>
> // global attributes: :Conventions = \"CF-1.6\" ; :history =
> \"2020-10-06 09:12:58 UTC+1100 by era5~replicationtools~-1.2.1: mv
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msdwswrfcs/2010*.msdwswrfcs~era5opersfc20100101~-20100131.tmp
> /g/data/rt52/era5/single-levels/reanalysis/msdwswrfcs/2010/msdwswrfcs~era5opersfc20100101~-20100131.nc`\n`{=latex}\",
> \"2020-10-06 09:12:25 UTC+1100 by era5~replicationtools~-1.2.1: nccopy
> -k4 -d5 -s
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msdwswrfcs/2010/msdwswrfcs~era5opersfc20100101~-20100131.nc
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msdwswrfcs/2010*.msdwswrfcs~era5opersfc20100101~-20100131.tmp`\n`{=latex}\",
> \"2020-10-06 09:08:27 UTC+1100 by era5~replicationtools~-1.2.1: curl
> --connect-timeout 20 --show-error --silent --max-time 36000 -o
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msdwswrfcs/2010/msdwswrfcs~era5opersfc20100101~-20100131.nc
> [http://136.156.132.235/cache-compute-0006/cache/data8/adaptor.mars.internal-1601935429.6491008-27899-16-e96963ed-62ea-4560-b279-6b1ed587aa92.nc\\n](http://136.156.132.235/cache-compute-0006/cache/data8/adaptor.mars.internal-1601935429.6491008-27899-16-e96963ed-62ea-4560-b279-6b1ed587aa92.nc\n)\",
> \"2020-10-05 22:04:51 GMT by grib~tonetcdf~-2.16.0:
> /opt/ecmwf/eccodes/bin/grib~tonetcdf~ -S param -o
> /cache/data8/adaptor.mars.internal-1601935429.6491008-27899-16-e96963ed-62ea-4560-b279-6b1ed587aa92.nc
> /cache/tmp/e96963ed-62ea-4560-b279-6b1ed587aa92-adaptor.mars.internal-1601935429.6498413-27899-8-tmp.grib\"
> ; :license = \"Licence to use Copernicus Products:
> <https://apps.ecmwf.int/datasets/licences/copernicus/>\" ; :summary =
> \"ERA5 is the fifth generation ECMWF atmospheric reanalysis of the
> global climate. This file is part of the ERA5 replica hosted at NCI
> Australia. For more information please see
> <http://dx.doi.org/10.25914/5f48874388857>\" ; :title = \"ERA5
> single-levels reanalysis
> mean~surfacedownwardshortwaveradiationfluxclearsky~
> 20100101-20100131\" ; }

### [msl](file:///Volumes/ioa02/reanalysis/ERA5/msl)

Mean sea level pressure

> netcdf msl~era5opersfc20100101~-20100131 { dimensions: longitude =
> 1440 ; latitude = 721 ; time = 744 ; variables: float
> longitude(longitude) ; longitude:units = \"degrees~east~\" ;
> longitude:long~name~ = \"longitude\" ; float latitude(latitude) ;
> latitude:units = \"degrees~north~\" ; latitude:long~name~ =
> \"latitude\" ; int time(time) ; time:units = \"hours since 1900-01-01
> 00:00:00.0\" ; time:long~name~ = \"time\" ; time:calendar =
> \"gregorian\" ; short msl(time, latitude, longitude) ;
> msl:scale~factor~ = 0.196378160621366 ; msl:add~offset~ =
> 100651.52681092 ; msl:~FillValue~ = -32767s ; msl:missing~value~ =
> -32767s ; msl:units = \"Pa\" ; msl:long~name~ = \"Mean sea level
> pressure\" ; msl:standard~name~ = \"air~pressureatmeansealevel~\" ;
>
> // global attributes: :Conventions = \"CF-1.6\" ; :history =
> \"2020-09-28 12:39:13 UTC+1000 by era5~replicationtools~-1.2.1: mv
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msl/2010*.msl~era5opersfc20100101~-20100131.tmp
> /g/data/rt52/era5/single-levels/reanalysis/msl/2010/msl~era5opersfc20100101~-20100131.nc`\n`{=latex}\",
> \"2020-09-28 12:38:23 UTC+1000 by era5~replicationtools~-1.2.1: nccopy
> -k4 -d5 -s
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msl/2010/msl~era5opersfc20100101~-20100131.nc
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msl/2010*.msl~era5opersfc20100101~-20100131.tmp`\n`{=latex}\",
> \"2020-09-28 12:35:03 UTC+1000 by era5~replicationtools~-1.2.1: curl
> --connect-timeout 20 --show-error --silent --max-time 36000 -o
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msl/2010/msl~era5opersfc20100101~-20100131.nc
> [http://136.156.133.39/cache-compute-0012/cache/data5/adaptor.mars.internal-1601260296.1220453-31017-5-1e4f2d61-ed30-44f1-bdfb-012d72c51c74.nc\\n](http://136.156.133.39/cache-compute-0012/cache/data5/adaptor.mars.internal-1601260296.1220453-31017-5-1e4f2d61-ed30-44f1-bdfb-012d72c51c74.nc\n)\",
> \"2020-09-28 02:32:49 GMT by grib~tonetcdf~-2.16.0:
> /opt/ecmwf/eccodes/bin/grib~tonetcdf~ -S param -o
> /cache/data5/adaptor.mars.internal-1601260296.1220453-31017-5-1e4f2d61-ed30-44f1-bdfb-012d72c51c74.nc
> /cache/tmp/1e4f2d61-ed30-44f1-bdfb-012d72c51c74-adaptor.mars.internal-1601260296.12292-31017-3-tmp.grib\"
> ; :license = \"Licence to use Copernicus Products:
> <https://apps.ecmwf.int/datasets/licences/copernicus/>\" ; :summary =
> \"ERA5 is the fifth generation ECMWF atmospheric reanalysis of the
> global climate. This file is part of the ERA5 replica hosted at NCI
> Australia. For more information please see
> <http://dx.doi.org/10.25914/5f48874388857>\" ; :title = \"ERA5
> single-levels reanalysis mean~sealevelpressure~ 20100101-20100131\" ;
> }

### [msnlwrf](file:///Volumes/ioa02/reanalysis/ERA5/msnlwrf)

Mean surface net long-wave radiation flux

> netcdf msnlwrf~era5opersfc20100101~-20100131 { dimensions: longitude =
> 1440 ; latitude = 721 ; time = 744 ; variables: float
> longitude(longitude) ; longitude:units = \"degrees~east~\" ;
> longitude:long~name~ = \"longitude\" ; float latitude(latitude) ;
> latitude:units = \"degrees~north~\" ; latitude:long~name~ =
> \"latitude\" ; int time(time) ; time:units = \"hours since 1900-01-01
> 00:00:00.0\" ; time:long~name~ = \"time\" ; time:calendar =
> \"gregorian\" ; short msnlwrf(time, latitude, longitude) ;
> msnlwrf:scale~factor~ = 0.00765515516156936 ; msnlwrf:add~offset~ =
> -185.368207460393 ; msnlwrf:~FillValue~ = -32767s ;
> msnlwrf:missing~value~ = -32767s ; msnlwrf:units = \"W m\*\*-2\" ;
> msnlwrf:long~name~ = \"Mean surface net long-wave radiation flux\" ;
>
> // global attributes: :Conventions = \"CF-1.6\" ; :history =
> \"2020-10-06 09:09:18 UTC+1100 by era5~replicationtools~-1.2.1: mv
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msnlwrf/2010*.msnlwrf~era5opersfc20100101~-20100131.tmp
> /g/data/rt52/era5/single-levels/reanalysis/msnlwrf/2010/msnlwrf~era5opersfc20100101~-20100131.nc`\n`{=latex}\",
> \"2020-10-06 09:08:14 UTC+1100 by era5~replicationtools~-1.2.1: nccopy
> -k4 -d5 -s
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msnlwrf/2010/msnlwrf~era5opersfc20100101~-20100131.nc
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msnlwrf/2010*.msnlwrf~era5opersfc20100101~-20100131.tmp`\n`{=latex}\",
> \"2020-10-06 09:02:18 UTC+1100 by era5~replicationtools~-1.2.1: curl
> --connect-timeout 20 --show-error --silent --max-time 36000 -o
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msnlwrf/2010/msnlwrf~era5opersfc20100101~-20100131.nc
> [http://136.156.132.153/cache-compute-0002/cache/data6/adaptor.mars.internal-1601935105.2208984-18608-5-8a196cc4-a7d0-47dd-b326-528e0e8ffb31.nc\\n](http://136.156.132.153/cache-compute-0002/cache/data6/adaptor.mars.internal-1601935105.2208984-18608-5-8a196cc4-a7d0-47dd-b326-528e0e8ffb31.nc\n)\",
> \"2020-10-05 21:59:08 GMT by grib~tonetcdf~-2.16.0:
> /opt/ecmwf/eccodes/bin/grib~tonetcdf~ -S param -o
> /cache/data6/adaptor.mars.internal-1601935105.2208984-18608-5-8a196cc4-a7d0-47dd-b326-528e0e8ffb31.nc
> /cache/tmp/8a196cc4-a7d0-47dd-b326-528e0e8ffb31-adaptor.mars.internal-1601935105.2214694-18608-1-tmp.grib\"
> ; :license = \"Licence to use Copernicus Products:
> <https://apps.ecmwf.int/datasets/licences/copernicus/>\" ; :summary =
> \"ERA5 is the fifth generation ECMWF atmospheric reanalysis of the
> global climate. This file is part of the ERA5 replica hosted at NCI
> Australia. For more information please see
> <http://dx.doi.org/10.25914/5f48874388857>\" ; :title = \"ERA5
> single-levels reanalysis mean~surfacenetlongwaveradiationflux~
> 20100101-20100131\" ; }

### [msnlwrfcs](file:///Volumes/ioa02/reanalysis/ERA5/msnlwrfcs)

Mean surface net long-wave radiation flux, clear sky

> netcdf msnlwrfcs~era5opersfc20100101~-20100131 { dimensions: longitude
> = 1440 ; latitude = 721 ; time = 744 ; variables: float
> longitude(longitude) ; longitude:units = \"degrees~east~\" ;
> longitude:long~name~ = \"longitude\" ; float latitude(latitude) ;
> latitude:units = \"degrees~north~\" ; latitude:long~name~ =
> \"latitude\" ; int time(time) ; time:units = \"hours since 1900-01-01
> 00:00:00.0\" ; time:long~name~ = \"time\" ; time:calendar =
> \"gregorian\" ; short msnlwrfcs(time, latitude, longitude) ;
> msnlwrfcs:scale~factor~ = 0.00686127063368837 ; msnlwrfcs:add~offset~
> = -207.683118135317 ; msnlwrfcs:~FillValue~ = -32767s ;
> msnlwrfcs:missing~value~ = -32767s ; msnlwrfcs:units = \"W m\*\*-2\" ;
> msnlwrfcs:long~name~ = \"Mean surface net long-wave radiation flux,
> clear sky\" ;
>
> // global attributes: :Conventions = \"CF-1.6\" ; :history =
> \"2020-10-06 09:11:31 UTC+1100 by era5~replicationtools~-1.2.1: mv
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msnlwrfcs/2010*.msnlwrfcs~era5opersfc20100101~-20100131.tmp
> /g/data/rt52/era5/single-levels/reanalysis/msnlwrfcs/2010/msnlwrfcs~era5opersfc20100101~-20100131.nc`\n`{=latex}\",
> \"2020-10-06 09:10:36 UTC+1100 by era5~replicationtools~-1.2.1: nccopy
> -k4 -d5 -s
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msnlwrfcs/2010/msnlwrfcs~era5opersfc20100101~-20100131.nc
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msnlwrfcs/2010*.msnlwrfcs~era5opersfc20100101~-20100131.tmp`\n`{=latex}\",
> \"2020-10-06 09:06:58 UTC+1100 by era5~replicationtools~-1.2.1: curl
> --connect-timeout 20 --show-error --silent --max-time 36000 -o
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msnlwrfcs/2010/msnlwrfcs~era5opersfc20100101~-20100131.nc
> [http://136.156.133.25/cache-compute-0008/cache/data7/adaptor.mars.internal-1601935291.563169-20885-28-16218919-c599-4310-b83c-d3e5b2d7ef7d.nc\\n](http://136.156.133.25/cache-compute-0008/cache/data7/adaptor.mars.internal-1601935291.563169-20885-28-16218919-c599-4310-b83c-d3e5b2d7ef7d.nc\n)\",
> \"2020-10-05 22:02:20 GMT by grib~tonetcdf~-2.16.0:
> /opt/ecmwf/eccodes/bin/grib~tonetcdf~ -S param -o
> /cache/data7/adaptor.mars.internal-1601935291.563169-20885-28-16218919-c599-4310-b83c-d3e5b2d7ef7d.nc
> /cache/tmp/16218919-c599-4310-b83c-d3e5b2d7ef7d-adaptor.mars.internal-1601935291.5818343-20885-10-tmp.grib\"
> ; :license = \"Licence to use Copernicus Products:
> <https://apps.ecmwf.int/datasets/licences/copernicus/>\" ; :summary =
> \"ERA5 is the fifth generation ECMWF atmospheric reanalysis of the
> global climate. This file is part of the ERA5 replica hosted at NCI
> Australia. For more information please see
> <http://dx.doi.org/10.25914/5f48874388857>\" ; :title = \"ERA5
> single-levels reanalysis mean~surfacenetlongwaveradiationfluxclearsky~
> 20100101-20100131\" ; }

### [msnswrf](file:///Volumes/ioa02/reanalysis/ERA5/msnswrf)

Mean surface net short-wave radiation flux

> netcdf msnswrf~era5opersfc20100101~-20100131 { dimensions: longitude =
> 1440 ; latitude = 721 ; time = 744 ; variables: float
> longitude(longitude) ; longitude:units = \"degrees~east~\" ;
> longitude:long~name~ = \"longitude\" ; float latitude(latitude) ;
> latitude:units = \"degrees~north~\" ; latitude:long~name~ =
> \"latitude\" ; int time(time) ; time:units = \"hours since 1900-01-01
> 00:00:00.0\" ; time:long~name~ = \"time\" ; time:calendar =
> \"gregorian\" ; short msnswrf(time, latitude, longitude) ;
> msnswrf:scale~factor~ = 0.0171468954572505 ; msnswrf:add~offset~ =
> 561.835176552271 ; msnswrf:~FillValue~ = -32767s ;
> msnswrf:missing~value~ = -32767s ; msnswrf:units = \"W m\*\*-2\" ;
> msnswrf:long~name~ = \"Mean surface net short-wave radiation flux\" ;
>
> // global attributes: :Conventions = \"CF-1.6\" ; :history =
> \"2020-10-06 09:07:34 UTC+1100 by era5~replicationtools~-1.2.1: mv
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msnswrf/2010*.msnswrf~era5opersfc20100101~-20100131.tmp
> /g/data/rt52/era5/single-levels/reanalysis/msnswrf/2010/msnswrf~era5opersfc20100101~-20100131.nc`\n`{=latex}\",
> \"2020-10-06 09:06:58 UTC+1100 by era5~replicationtools~-1.2.1: nccopy
> -k4 -d5 -s
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msnswrf/2010/msnswrf~era5opersfc20100101~-20100131.nc
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msnswrf/2010*.msnswrf~era5opersfc20100101~-20100131.tmp`\n`{=latex}\",
> \"2020-10-06 09:02:13 UTC+1100 by era5~replicationtools~-1.2.1: curl
> --connect-timeout 20 --show-error --silent --max-time 36000 -o
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msnswrf/2010/msnswrf~era5opersfc20100101~-20100131.nc
> [http://136.156.133.39/cache-compute-0012/cache/data8/adaptor.mars.internal-1601935069.7217205-6425-1-6f222a72-e8c4-4c60-acb7-d3b89ffc8e8f.nc\\n](http://136.156.133.39/cache-compute-0012/cache/data8/adaptor.mars.internal-1601935069.7217205-6425-1-6f222a72-e8c4-4c60-acb7-d3b89ffc8e8f.nc\n)\",
> \"2020-10-05 21:58:59 GMT by grib~tonetcdf~-2.16.0:
> /opt/ecmwf/eccodes/bin/grib~tonetcdf~ -S param -o
> /cache/data8/adaptor.mars.internal-1601935069.7217205-6425-1-6f222a72-e8c4-4c60-acb7-d3b89ffc8e8f.nc
> /cache/tmp/6f222a72-e8c4-4c60-acb7-d3b89ffc8e8f-adaptor.mars.internal-1601935069.723566-6425-1-tmp.grib\"
> ; :license = \"Licence to use Copernicus Products:
> <https://apps.ecmwf.int/datasets/licences/copernicus/>\" ; :summary =
> \"ERA5 is the fifth generation ECMWF atmospheric reanalysis of the
> global climate. This file is part of the ERA5 replica hosted at NCI
> Australia. For more information please see
> <http://dx.doi.org/10.25914/5f48874388857>\" ; :title = \"ERA5
> single-levels reanalysis mean~surfacenetshortwaveradiationflux~
> 20100101-20100131\" ; }

### [msnswrfcs](file:///Volumes/ioa02/reanalysis/ERA5/msnswrfcs)

Mean surface net short-wave radiation flux, clear sky

> netcdf msnswrfcs~era5opersfc20100101~-20100131 { dimensions: longitude
> = 1440 ; latitude = 721 ; time = 744 ; variables: float
> longitude(longitude) ; longitude:units = \"degrees~east~\" ;
> longitude:long~name~ = \"longitude\" ; float latitude(latitude) ;
> latitude:units = \"degrees~north~\" ; latitude:long~name~ =
> \"latitude\" ; int time(time) ; time:units = \"hours since 1900-01-01
> 00:00:00.0\" ; time:long~name~ = \"time\" ; time:calendar =
> \"gregorian\" ; short msnswrfcs(time, latitude, longitude) ;
> msnswrfcs:scale~factor~ = 0.0169985923122702 ; msnswrfcs:add~offset~ =
> 556.975875703844 ; msnswrfcs:~FillValue~ = -32767s ;
> msnswrfcs:missing~value~ = -32767s ; msnswrfcs:units = \"W m\*\*-2\" ;
> msnswrfcs:long~name~ = \"Mean surface net short-wave radiation flux,
> clear sky\" ;
>
> // global attributes: :Conventions = \"CF-1.6\" ; :history =
> \"2020-10-06 09:12:36 UTC+1100 by era5~replicationtools~-1.2.1: mv
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msnswrfcs/2010*.msnswrfcs~era5opersfc20100101~-20100131.tmp
> /g/data/rt52/era5/single-levels/reanalysis/msnswrfcs/2010/msnswrfcs~era5opersfc20100101~-20100131.nc`\n`{=latex}\",
> \"2020-10-06 09:12:01 UTC+1100 by era5~replicationtools~-1.2.1: nccopy
> -k4 -d5 -s
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msnswrfcs/2010/msnswrfcs~era5opersfc20100101~-20100131.nc
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msnswrfcs/2010*.msnswrfcs~era5opersfc20100101~-20100131.tmp`\n`{=latex}\",
> \"2020-10-06 09:06:10 UTC+1100 by era5~replicationtools~-1.2.1: curl
> --connect-timeout 20 --show-error --silent --max-time 36000 -o
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msnswrfcs/2010/msnswrfcs~era5opersfc20100101~-20100131.nc
> [http://136.156.133.39/cache-compute-0012/cache/data6/adaptor.mars.internal-1601935253.1640558-4989-18-bababa4c-ec9e-4ba6-98cd-9ad90e6be27b.nc\\n](http://136.156.133.39/cache-compute-0012/cache/data6/adaptor.mars.internal-1601935253.1640558-4989-18-bababa4c-ec9e-4ba6-98cd-9ad90e6be27b.nc\n)\",
> \"2020-10-05 22:01:38 GMT by grib~tonetcdf~-2.16.0:
> /opt/ecmwf/eccodes/bin/grib~tonetcdf~ -S param -o
> /cache/data6/adaptor.mars.internal-1601935253.1640558-4989-18-bababa4c-ec9e-4ba6-98cd-9ad90e6be27b.nc
> /cache/tmp/bababa4c-ec9e-4ba6-98cd-9ad90e6be27b-adaptor.mars.internal-1601935253.164695-4989-6-tmp.grib\"
> ; :license = \"Licence to use Copernicus Products:
> <https://apps.ecmwf.int/datasets/licences/copernicus/>\" ; :summary =
> \"ERA5 is the fifth generation ECMWF atmospheric reanalysis of the
> global climate. This file is part of the ERA5 replica hosted at NCI
> Australia. For more information please see
> <http://dx.doi.org/10.25914/5f48874388857>\" ; :title = \"ERA5
> single-levels reanalysis
> mean~surfacenetshortwaveradiationfluxclearsky~ 20100101-20100131\" ; }

### [msr](file:///Volumes/ioa02/reanalysis/ERA5/msr)

Mean snowfall rate

> netcdf msr~era5opersfc20100101~-20100131 { dimensions: longitude =
> 1440 ; latitude = 721 ; time = 744 ; variables: float
> longitude(longitude) ; longitude:units = \"degrees~east~\" ;
> longitude:long~name~ = \"longitude\" ; float latitude(latitude) ;
> latitude:units = \"degrees~north~\" ; latitude:long~name~ =
> \"latitude\" ; int time(time) ; time:units = \"hours since 1900-01-01
> 00:00:00.0\" ; time:long~name~ = \"time\" ; time:calendar =
> \"gregorian\" ; short msr(time, latitude, longitude) ;
> msr:scale~factor~ = 3.16955222732488e-08 ; msr:add~offset~ =
> 0.00103853548280527 ; msr:~FillValue~ = -32767s ; msr:missing~value~ =
> -32767s ; msr:units = \"kg m\*\*-2 s\*\*-1\" ; msr:long~name~ = \"Mean
> snowfall rate\" ;
>
> // global attributes: :Conventions = \"CF-1.6\" ; :history =
> \"2020-10-06 09:04:04 UTC+1100 by era5~replicationtools~-1.2.1: mv
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msr/2010*.msr~era5opersfc20100101~-20100131.tmp
> /g/data/rt52/era5/single-levels/reanalysis/msr/2010/msr~era5opersfc20100101~-20100131.nc`\n`{=latex}\",
> \"2020-10-06 09:03:42 UTC+1100 by era5~replicationtools~-1.2.1: nccopy
> -k4 -d5 -s
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msr/2010/msr~era5opersfc20100101~-20100131.nc
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msr/2010*.msr~era5opersfc20100101~-20100131.tmp`\n`{=latex}\",
> \"2020-10-06 09:00:18 UTC+1100 by era5~replicationtools~-1.2.1: curl
> --connect-timeout 20 --show-error --silent --max-time 36000 -o
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msr/2010/msr~era5opersfc20100101~-20100131.nc
> [http://136.156.133.25/cache-compute-0008/cache/data8/adaptor.mars.internal-1601934958.8138766-19540-22-829a84f4-b413-4b71-b15a-99bb8c146d46.nc\\n](http://136.156.133.25/cache-compute-0008/cache/data8/adaptor.mars.internal-1601934958.8138766-19540-22-829a84f4-b413-4b71-b15a-99bb8c146d46.nc\n)\",
> \"2020-10-05 21:56:43 GMT by grib~tonetcdf~-2.16.0:
> /opt/ecmwf/eccodes/bin/grib~tonetcdf~ -S param -o
> /cache/data8/adaptor.mars.internal-1601934958.8138766-19540-22-829a84f4-b413-4b71-b15a-99bb8c146d46.nc
> /cache/tmp/829a84f4-b413-4b71-b15a-99bb8c146d46-adaptor.mars.internal-1601934958.8144617-19540-8-tmp.grib\"
> ; :license = \"Licence to use Copernicus Products:
> <https://apps.ecmwf.int/datasets/licences/copernicus/>\" ; :summary =
> \"ERA5 is the fifth generation ECMWF atmospheric reanalysis of the
> global climate. This file is part of the ERA5 replica hosted at NCI
> Australia. For more information please see
> <http://dx.doi.org/10.25914/5f48874388857>\" ; :title = \"ERA5
> single-levels reanalysis mean~snowfallrate~ 20100101-20100131\" ; }

### [msror](file:///Volumes/ioa02/reanalysis/ERA5/msror)

Mean surface runoff rate

> netcdf msror~era5opersfc20100101~-20100131 { dimensions: longitude =
> 1440 ; latitude = 721 ; time = 744 ; variables: float
> longitude(longitude) ; longitude:units = \"degrees~east~\" ;
> longitude:long~name~ = \"longitude\" ; float latitude(latitude) ;
> latitude:units = \"degrees~north~\" ; latitude:long~name~ =
> \"latitude\" ; int time(time) ; time:units = \"hours since 1900-01-01
> 00:00:00.0\" ; time:long~name~ = \"time\" ; time:calendar =
> \"gregorian\" ; short msror(time, latitude, longitude) ;
> msror:scale~factor~ = 3.49596662796261e-07 ; msror:add~offset~ =
> 0.0114548842531823 ; msror:~FillValue~ = -32767s ;
> msror:missing~value~ = -32767s ; msror:units = \"kg m\*\*-2 s\*\*-1\"
> ; msror:long~name~ = \"Mean surface runoff rate\" ;
>
> // global attributes: :Conventions = \"CF-1.6\" ; :history =
> \"2020-10-06 09:02:38 UTC+1100 by era5~replicationtools~-1.2.1: mv
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msror/2010*.msror~era5opersfc20100101~-20100131.tmp
> /g/data/rt52/era5/single-levels/reanalysis/msror/2010/msror~era5opersfc20100101~-20100131.nc`\n`{=latex}\",
> \"2020-10-06 09:02:24 UTC+1100 by era5~replicationtools~-1.2.1: nccopy
> -k4 -d5 -s
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msror/2010/msror~era5opersfc20100101~-20100131.nc
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/msror/2010*.msror~era5opersfc20100101~-20100131.tmp`\n`{=latex}\",
> \"2020-10-06 08:58:36 UTC+1100 by era5~replicationtools~-1.2.1: curl
> --connect-timeout 20 --show-error --silent --max-time 36000 -o
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/msror/2010/msror~era5opersfc20100101~-20100131.nc
> [http://136.156.133.46/cache-compute-0015/cache/data2/adaptor.mars.internal-1601934851.1730347-28836-22-a7717ea4-1d51-4b28-a876-87ce3356d63c.nc\\n](http://136.156.133.46/cache-compute-0015/cache/data2/adaptor.mars.internal-1601934851.1730347-28836-22-a7717ea4-1d51-4b28-a876-87ce3356d63c.nc\n)\",
> \"2020-10-05 21:55:17 GMT by grib~tonetcdf~-2.16.0:
> /opt/ecmwf/eccodes/bin/grib~tonetcdf~ -S param -o
> /cache/data2/adaptor.mars.internal-1601934851.1730347-28836-22-a7717ea4-1d51-4b28-a876-87ce3356d63c.nc
> /cache/tmp/a7717ea4-1d51-4b28-a876-87ce3356d63c-adaptor.mars.internal-1601934851.173735-28836-8-tmp.grib\"
> ; :license = \"Licence to use Copernicus Products:
> <https://apps.ecmwf.int/datasets/licences/copernicus/>\" ; :summary =
> \"ERA5 is the fifth generation ECMWF atmospheric reanalysis of the
> global climate. This file is part of the ERA5 replica hosted at NCI
> Australia. For more information please see
> <http://dx.doi.org/10.25914/5f48874388857>\" ; :title = \"ERA5
> single-levels reanalysis mean~surfacerunoffrate~ 20100101-20100131\" ;
> }

### [mtdwswrf](file:///Volumes/ioa02/reanalysis/ERA5/mtdwswrf)

Mean top downward short-wave radiation flux

> netcdf mtdwswrf~era5opersfc20100101~-20100131 { dimensions: longitude
> = 1440 ; latitude = 721 ; time = 744 ; variables: float
> longitude(longitude) ; longitude:units = \"degrees~east~\" ;
> longitude:long~name~ = \"longitude\" ; float latitude(latitude) ;
> latitude:units = \"degrees~north~\" ; latitude:long~name~ =
> \"latitude\" ; int time(time) ; time:units = \"hours since 1900-01-01
> 00:00:00.0\" ; time:long~name~ = \"time\" ; time:calendar =
> \"gregorian\" ; short mtdwswrf(time, latitude, longitude) ;
> mtdwswrf:scale~factor~ = 0.021422889994354 ; mtdwswrf:add~offset~ =
> 701.942413555003 ; mtdwswrf:~FillValue~ = -32767s ;
> mtdwswrf:missing~value~ = -32767s ; mtdwswrf:units = \"W m\*\*-2\" ;
> mtdwswrf:long~name~ = \"Mean top downward short-wave radiation flux\"
> ;
>
> // global attributes: :Conventions = \"CF-1.6\" ; :history =
> \"2020-10-06 09:11:01 UTC+1100 by era5~replicationtools~-1.2.1: mv
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/mtdwswrf/2010*.mtdwswrf~era5opersfc20100101~-20100131.tmp
> /g/data/rt52/era5/single-levels/reanalysis/mtdwswrf/2010/mtdwswrf~era5opersfc20100101~-20100131.nc`\n`{=latex}\",
> \"2020-10-06 09:10:32 UTC+1100 by era5~replicationtools~-1.2.1: nccopy
> -k4 -d5 -s
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/mtdwswrf/2010/mtdwswrf~era5opersfc20100101~-20100131.nc
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/mtdwswrf/2010*.mtdwswrf~era5opersfc20100101~-20100131.tmp`\n`{=latex}\",
> \"2020-10-06 09:06:38 UTC+1100 by era5~replicationtools~-1.2.1: curl
> --connect-timeout 20 --show-error --silent --max-time 36000 -o
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/mtdwswrf/2010/mtdwswrf~era5opersfc20100101~-20100131.nc
> [http://136.156.132.236/cache-compute-0007/cache/data7/adaptor.mars.internal-1601935311.62434-12055-32-fb7e93c7-424f-4538-bf13-ba57fce74df4.nc\\n](http://136.156.132.236/cache-compute-0007/cache/data7/adaptor.mars.internal-1601935311.62434-12055-32-fb7e93c7-424f-4538-bf13-ba57fce74df4.nc\n)\",
> \"2020-10-05 22:02:35 GMT by grib~tonetcdf~-2.16.0:
> /opt/ecmwf/eccodes/bin/grib~tonetcdf~ -S param -o
> /cache/data7/adaptor.mars.internal-1601935311.62434-12055-32-fb7e93c7-424f-4538-bf13-ba57fce74df4.nc
> /cache/tmp/fb7e93c7-424f-4538-bf13-ba57fce74df4-adaptor.mars.internal-1601935311.6249025-12055-14-tmp.grib\"
> ; :license = \"Licence to use Copernicus Products:
> <https://apps.ecmwf.int/datasets/licences/copernicus/>\" ; :summary =
> \"ERA5 is the fifth generation ECMWF atmospheric reanalysis of the
> global climate. This file is part of the ERA5 replica hosted at NCI
> Australia. For more information please see
> <http://dx.doi.org/10.25914/5f48874388857>\" ; :title = \"ERA5
> single-levels reanalysis mean~topdownwardshortwaveradiationflux~
> 20100101-20100131\" ; }

### [mtnlwrf](file:///Volumes/ioa02/reanalysis/ERA5/mtnlwrf)

Mean top net long-wave radiation flux

> netcdf mtnlwrf~era5opersfc20100101~-20100131 { dimensions: longitude =
> 1440 ; latitude = 721 ; time = 744 ; variables: float
> longitude(longitude) ; longitude:units = \"degrees~east~\" ;
> longitude:long~name~ = \"longitude\" ; float latitude(latitude) ;
> latitude:units = \"degrees~north~\" ; latitude:long~name~ =
> \"latitude\" ; int time(time) ; time:units = \"hours since 1900-01-01
> 00:00:00.0\" ; time:long~name~ = \"time\" ; time:calendar =
> \"gregorian\" ; short mtnlwrf(time, latitude, longitude) ;
> mtnlwrf:scale~factor~ = 0.00547009042080898 ; mtnlwrf:add~offset~ =
> -239.618579771773 ; mtnlwrf:~FillValue~ = -32767s ;
> mtnlwrf:missing~value~ = -32767s ; mtnlwrf:units = \"W m\*\*-2\" ;
> mtnlwrf:long~name~ = \"Mean top net long-wave radiation flux\" ;
>
> // global attributes: :Conventions = \"CF-1.6\" ; :history =
> \"2020-10-06 09:09:23 UTC+1100 by era5~replicationtools~-1.2.1: mv
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/mtnlwrf/2010*.mtnlwrf~era5opersfc20100101~-20100131.tmp
> /g/data/rt52/era5/single-levels/reanalysis/mtnlwrf/2010/mtnlwrf~era5opersfc20100101~-20100131.nc`\n`{=latex}\",
> \"2020-10-06 09:08:21 UTC+1100 by era5~replicationtools~-1.2.1: nccopy
> -k4 -d5 -s
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/mtnlwrf/2010/mtnlwrf~era5opersfc20100101~-20100131.nc
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/mtnlwrf/2010*.mtnlwrf~era5opersfc20100101~-20100131.tmp`\n`{=latex}\",
> \"2020-10-06 09:02:24 UTC+1100 by era5~replicationtools~-1.2.1: curl
> --connect-timeout 20 --show-error --silent --max-time 36000 -o
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/mtnlwrf/2010/mtnlwrf~era5opersfc20100101~-20100131.nc
> [http://136.156.133.41/cache-compute-0013/cache/data1/adaptor.mars.internal-1601935118.9578016-15806-13-c7cd5c57-4899-4528-8267-f01896333fba.nc\\n](http://136.156.133.41/cache-compute-0013/cache/data1/adaptor.mars.internal-1601935118.9578016-15806-13-c7cd5c57-4899-4528-8267-f01896333fba.nc\n)\",
> \"2020-10-05 21:59:22 GMT by grib~tonetcdf~-2.16.0:
> /opt/ecmwf/eccodes/bin/grib~tonetcdf~ -S param -o
> /cache/data1/adaptor.mars.internal-1601935118.9578016-15806-13-c7cd5c57-4899-4528-8267-f01896333fba.nc
> /cache/tmp/c7cd5c57-4899-4528-8267-f01896333fba-adaptor.mars.internal-1601935118.9584107-15806-4-tmp.grib\"
> ; :license = \"Licence to use Copernicus Products:
> <https://apps.ecmwf.int/datasets/licences/copernicus/>\" ; :summary =
> \"ERA5 is the fifth generation ECMWF atmospheric reanalysis of the
> global climate. This file is part of the ERA5 replica hosted at NCI
> Australia. For more information please see
> <http://dx.doi.org/10.25914/5f48874388857>\" ; :title = \"ERA5
> single-levels reanalysis mean~topnetlongwaveradiationflux~
> 20100101-20100131\" ; }

### [mtnlwrfcs](file:///Volumes/ioa02/reanalysis/ERA5/mtnlwrfcs)

Mean top net long-wave radiation flux, clear sky

> netcdf mtnlwrfcs~era5opersfc20100101~-20100131 { dimensions: longitude
> = 1440 ; latitude = 721 ; time = 744 ; variables: float
> longitude(longitude) ; longitude:units = \"degrees~east~\" ;
> longitude:long~name~ = \"longitude\" ; float latitude(latitude) ;
> latitude:units = \"degrees~north~\" ; latitude:long~name~ =
> \"latitude\" ; int time(time) ; time:units = \"hours since 1900-01-01
> 00:00:00.0\" ; time:long~name~ = \"time\" ; time:calendar =
> \"gregorian\" ; short mtnlwrfcs(time, latitude, longitude) ;
> mtnlwrfcs:scale~factor~ = 0.00447323906109327 ; mtnlwrfcs:add~offset~
> = -272.350747361718 ; mtnlwrfcs:~FillValue~ = -32767s ;
> mtnlwrfcs:missing~value~ = -32767s ; mtnlwrfcs:units = \"W m\*\*-2\" ;
> mtnlwrfcs:long~name~ = \"Mean top net long-wave radiation flux, clear
> sky\" ;
>
> // global attributes: :Conventions = \"CF-1.6\" ; :history =
> \"2020-10-06 09:12:24 UTC+1100 by era5~replicationtools~-1.2.1: mv
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/mtnlwrfcs/2010*.mtnlwrfcs~era5opersfc20100101~-20100131.tmp
> /g/data/rt52/era5/single-levels/reanalysis/mtnlwrfcs/2010/mtnlwrfcs~era5opersfc20100101~-20100131.nc`\n`{=latex}\",
> \"2020-10-06 09:11:30 UTC+1100 by era5~replicationtools~-1.2.1: nccopy
> -k4 -d5 -s
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/mtnlwrfcs/2010/mtnlwrfcs~era5opersfc20100101~-20100131.nc
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/mtnlwrfcs/2010*.mtnlwrfcs~era5opersfc20100101~-20100131.tmp`\n`{=latex}\",
> \"2020-10-06 09:06:09 UTC+1100 by era5~replicationtools~-1.2.1: curl
> --connect-timeout 20 --show-error --silent --max-time 36000 -o
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/mtnlwrfcs/2010/mtnlwrfcs~era5opersfc20100101~-20100131.nc
> [http://136.156.132.110/cache-compute-0001/cache/data6/adaptor.mars.internal-1601935245.2007327-31712-13-9fe894bb-6f03-4d13-b576-5f138e35f11b.nc\\n](http://136.156.132.110/cache-compute-0001/cache/data6/adaptor.mars.internal-1601935245.2007327-31712-13-9fe894bb-6f03-4d13-b576-5f138e35f11b.nc\n)\",
> \"2020-10-05 22:01:29 GMT by grib~tonetcdf~-2.16.0:
> /opt/ecmwf/eccodes/bin/grib~tonetcdf~ -S param -o
> /cache/data6/adaptor.mars.internal-1601935245.2007327-31712-13-9fe894bb-6f03-4d13-b576-5f138e35f11b.nc
> /cache/tmp/9fe894bb-6f03-4d13-b576-5f138e35f11b-adaptor.mars.internal-1601935245.2013144-31712-6-tmp.grib\"
> ; :license = \"Licence to use Copernicus Products:
> <https://apps.ecmwf.int/datasets/licences/copernicus/>\" ; :summary =
> \"ERA5 is the fifth generation ECMWF atmospheric reanalysis of the
> global climate. This file is part of the ERA5 replica hosted at NCI
> Australia. For more information please see
> <http://dx.doi.org/10.25914/5f48874388857>\" ; :title = \"ERA5
> single-levels reanalysis mean~topnetlongwaveradiationfluxclearsky~
> 20100101-20100131\" ; }

### [mtnswrf](file:///Volumes/ioa02/reanalysis/ERA5/mtnswrf)

Mean top net short-wave radiation flux

> netcdf mtnswrf~era5opersfc20100101~-20100131 { dimensions: longitude =
> 1440 ; latitude = 721 ; time = 744 ; variables: float
> longitude(longitude) ; longitude:units = \"degrees~east~\" ;
> longitude:long~name~ = \"longitude\" ; float latitude(latitude) ;
> latitude:units = \"degrees~north~\" ; latitude:long~name~ =
> \"latitude\" ; int time(time) ; time:units = \"hours since 1900-01-01
> 00:00:00.0\" ; time:long~name~ = \"time\" ; time:calendar =
> \"gregorian\" ; short mtnswrf(time, latitude, longitude) ;
> mtnswrf:scale~factor~ = 0.0207843758106603 ; mtnswrf:add~offset~ =
> 681.020857812095 ; mtnswrf:~FillValue~ = -32767s ;
> mtnswrf:missing~value~ = -32767s ; mtnswrf:units = \"W m\*\*-2\" ;
> mtnswrf:long~name~ = \"Mean top net short-wave radiation flux\" ;
>
> // global attributes: :Conventions = \"CF-1.6\" ; :history =
> \"2020-10-06 09:08:09 UTC+1100 by era5~replicationtools~-1.2.1: mv
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/mtnswrf/2010*.mtnswrf~era5opersfc20100101~-20100131.tmp
> /g/data/rt52/era5/single-levels/reanalysis/mtnswrf/2010/mtnswrf~era5opersfc20100101~-20100131.nc`\n`{=latex}\",
> \"2020-10-06 09:07:33 UTC+1100 by era5~replicationtools~-1.2.1: nccopy
> -k4 -d5 -s
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/mtnswrf/2010/mtnswrf~era5opersfc20100101~-20100131.nc
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/mtnswrf/2010*.mtnswrf~era5opersfc20100101~-20100131.tmp`\n`{=latex}\",
> \"2020-10-06 09:03:06 UTC+1100 by era5~replicationtools~-1.2.1: curl
> --connect-timeout 20 --show-error --silent --max-time 36000 -o
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/mtnswrf/2010/mtnswrf~era5opersfc20100101~-20100131.nc
> [http://136.156.132.201/cache-compute-0004/cache/data9/adaptor.mars.internal-1601935110.2694445-4501-11-b5809659-6840-445d-a99a-233947a79209.nc\\n](http://136.156.132.201/cache-compute-0004/cache/data9/adaptor.mars.internal-1601935110.2694445-4501-11-b5809659-6840-445d-a99a-233947a79209.nc\n)\",
> \"2020-10-05 21:59:15 GMT by grib~tonetcdf~-2.16.0:
> /opt/ecmwf/eccodes/bin/grib~tonetcdf~ -S param -o
> /cache/data9/adaptor.mars.internal-1601935110.2694445-4501-11-b5809659-6840-445d-a99a-233947a79209.nc
> /cache/tmp/b5809659-6840-445d-a99a-233947a79209-adaptor.mars.internal-1601935110.270101-4501-4-tmp.grib\"
> ; :license = \"Licence to use Copernicus Products:
> <https://apps.ecmwf.int/datasets/licences/copernicus/>\" ; :summary =
> \"ERA5 is the fifth generation ECMWF atmospheric reanalysis of the
> global climate. This file is part of the ERA5 replica hosted at NCI
> Australia. For more information please see
> <http://dx.doi.org/10.25914/5f48874388857>\" ; :title = \"ERA5
> single-levels reanalysis mean~topnetshortwaveradiationflux~
> 20100101-20100131\" ; }

### [mtnswrfcs](file:///Volumes/ioa02/reanalysis/ERA5/mtnswrfcs)

Mean top net short-wave radiation flux, clear sky

> netcdf mtnswrfcs~era5opersfc20100101~-20100131 { dimensions: longitude
> = 1440 ; latitude = 721 ; time = 744 ; variables: float
> longitude(longitude) ; longitude:units = \"degrees~east~\" ;
> longitude:long~name~ = \"longitude\" ; float latitude(latitude) ;
> latitude:units = \"degrees~north~\" ; latitude:long~name~ =
> \"latitude\" ; int time(time) ; time:units = \"hours since 1900-01-01
> 00:00:00.0\" ; time:long~name~ = \"time\" ; time:calendar =
> \"gregorian\" ; short mtnswrfcs(time, latitude, longitude) ;
> mtnswrfcs:scale~factor~ = 0.0203156234263653 ; mtnswrfcs:add~offset~ =
> 665.661717188287 ; mtnswrfcs:~FillValue~ = -32767s ;
> mtnswrfcs:missing~value~ = -32767s ; mtnswrfcs:units = \"W m\*\*-2\" ;
> mtnswrfcs:long~name~ = \"Mean top net short-wave radiation flux, clear
> sky\" ;
>
> // global attributes: :Conventions = \"CF-1.6\" ; :history =
> \"2020-10-06 09:10:27 UTC+1100 by era5~replicationtools~-1.2.1: mv
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/mtnswrfcs/2010*.mtnswrfcs~era5opersfc20100101~-20100131.tmp
> /g/data/rt52/era5/single-levels/reanalysis/mtnswrfcs/2010/mtnswrfcs~era5opersfc20100101~-20100131.nc`\n`{=latex}\",
> \"2020-10-06 09:09:55 UTC+1100 by era5~replicationtools~-1.2.1: nccopy
> -k4 -d5 -s
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/mtnswrfcs/2010/mtnswrfcs~era5opersfc20100101~-20100131.nc
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/mtnswrfcs/2010*.mtnswrfcs~era5opersfc20100101~-20100131.tmp`\n`{=latex}\",
> \"2020-10-06 09:06:03 UTC+1100 by era5~replicationtools~-1.2.1: curl
> --connect-timeout 20 --show-error --silent --max-time 36000 -o
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/mtnswrfcs/2010/mtnswrfcs~era5opersfc20100101~-20100131.nc
> [http://136.156.133.37/cache-compute-0011/cache/data8/adaptor.mars.internal-1601935241.8915317-1334-3-4c3516f9-0313-410d-9c58-6f6f823a7e3f.nc\\n](http://136.156.133.37/cache-compute-0011/cache/data8/adaptor.mars.internal-1601935241.8915317-1334-3-4c3516f9-0313-410d-9c58-6f6f823a7e3f.nc\n)\",
> \"2020-10-05 22:01:23 GMT by grib~tonetcdf~-2.16.0:
> /opt/ecmwf/eccodes/bin/grib~tonetcdf~ -S param -o
> /cache/data8/adaptor.mars.internal-1601935241.8915317-1334-3-4c3516f9-0313-410d-9c58-6f6f823a7e3f.nc
> /cache/tmp/4c3516f9-0313-410d-9c58-6f6f823a7e3f-adaptor.mars.internal-1601935241.8921776-1334-1-tmp.grib\"
> ; :license = \"Licence to use Copernicus Products:
> <https://apps.ecmwf.int/datasets/licences/copernicus/>\" ; :summary =
> \"ERA5 is the fifth generation ECMWF atmospheric reanalysis of the
> global climate. This file is part of the ERA5 replica hosted at NCI
> Australia. For more information please see
> <http://dx.doi.org/10.25914/5f48874388857>\" ; :title = \"ERA5
> single-levels reanalysis mean~topnetshortwaveradiationfluxclearsky~
> 20100101-20100131\" ; }

### [mtpr](file:///Volumes/ioa02/reanalysis/ERA5/mtpr)

Mean total precipitation rate

> netcdf mtpr~era5opersfc20100101~-20100131 { dimensions: longitude =
> 1440 ; latitude = 721 ; time = 744 ; variables: float
> longitude(longitude) ; longitude:units = \"degrees~east~\" ;
> longitude:long~name~ = \"longitude\" ; float latitude(latitude) ;
> latitude:units = \"degrees~north~\" ; latitude:long~name~ =
> \"latitude\" ; int time(time) ; time:units = \"hours since 1900-01-01
> 00:00:00.0\" ; time:long~name~ = \"time\" ; time:calendar =
> \"gregorian\" ; short mtpr(time, latitude, longitude) ;
> mtpr:scale~factor~ = 3.71643823606165e-07 ; mtpr:add~offset~ =
> 0.0121772815242796 ; mtpr:~FillValue~ = -32767s ; mtpr:missing~value~
> = -32767s ; mtpr:units = \"kg m\*\*-2 s\*\*-1\" ; mtpr:long~name~ =
> \"Mean total precipitation rate\" ;
>
> // global attributes: :Conventions = \"CF-1.6\" ; :history =
> \"2020-10-06 09:14:14 UTC+1100 by era5~replicationtools~-1.2.1: mv
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/mtpr/2010*.mtpr~era5opersfc20100101~-20100131.tmp
> /g/data/rt52/era5/single-levels/reanalysis/mtpr/2010/mtpr~era5opersfc20100101~-20100131.nc`\n`{=latex}\",
> \"2020-10-06 09:13:40 UTC+1100 by era5~replicationtools~-1.2.1: nccopy
> -k4 -d5 -s
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/mtpr/2010/mtpr~era5opersfc20100101~-20100131.nc
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/mtpr/2010*.mtpr~era5opersfc20100101~-20100131.tmp`\n`{=latex}\",
> \"2020-10-06 09:07:33 UTC+1100 by era5~replicationtools~-1.2.1: curl
> --connect-timeout 20 --show-error --silent --max-time 36000 -o
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/mtpr/2010/mtpr~era5opersfc20100101~-20100131.nc
> [http://136.156.133.46/cache-compute-0015/cache/data8/adaptor.mars.internal-1601935342.505529-1432-7-5cc4d452-44f9-435e-8afc-b8d12664c4a7.nc\\n](http://136.156.133.46/cache-compute-0015/cache/data8/adaptor.mars.internal-1601935342.505529-1432-7-5cc4d452-44f9-435e-8afc-b8d12664c4a7.nc\n)\",
> \"2020-10-05 22:03:07 GMT by grib~tonetcdf~-2.16.0:
> /opt/ecmwf/eccodes/bin/grib~tonetcdf~ -S param -o
> /cache/data8/adaptor.mars.internal-1601935342.505529-1432-7-5cc4d452-44f9-435e-8afc-b8d12664c4a7.nc
> /cache/tmp/5cc4d452-44f9-435e-8afc-b8d12664c4a7-adaptor.mars.internal-1601935342.5060759-1432-3-tmp.grib\"
> ; :license = \"Licence to use Copernicus Products:
> <https://apps.ecmwf.int/datasets/licences/copernicus/>\" ; :summary =
> \"ERA5 is the fifth generation ECMWF atmospheric reanalysis of the
> global climate. This file is part of the ERA5 replica hosted at NCI
> Australia. For more information please see
> <http://dx.doi.org/10.25914/5f48874388857>\" ; :title = \"ERA5
> single-levels reanalysis mean~totalprecipitationrate~
> 20100101-20100131\" ; }

### [sp](file:///Volumes/ioa02/reanalysis/ERA5/sp)

Surface pressure

> netcdf sp~era5opersfc20100101~-20100131 { dimensions: longitude = 1440
> ; latitude = 721 ; time = 744 ; variables: float longitude(longitude)
> ; longitude:units = \"degrees~east~\" ; longitude:long~name~ =
> \"longitude\" ; float latitude(latitude) ; latitude:units =
> \"degrees~north~\" ; latitude:long~name~ = \"latitude\" ; int
> time(time) ; time:units = \"hours since 1900-01-01 00:00:00.0\" ;
> time:long~name~ = \"time\" ; time:calendar = \"gregorian\" ; short
> sp(time, latitude, longitude) ; sp:scale~factor~ = 0.882383148089512 ;
> sp:add~offset~ = 77345.396699051 ; sp:~FillValue~ = -32767s ;
> sp:missing~value~ = -32767s ; sp:units = \"Pa\" ; sp:long~name~ =
> \"Surface pressure\" ; sp:standard~name~ = \"surface~airpressure~\" ;
>
> // global attributes: :Conventions = \"CF-1.6\" ; :history =
> \"2020-09-28 12:39:13 UTC+1000 by era5~replicationtools~-1.2.1: mv
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/sp/2010*.sp~era5opersfc20100101~-20100131.tmp
> /g/data/rt52/era5/single-levels/reanalysis/sp/2010/sp~era5opersfc20100101~-20100131.nc`\n`{=latex}\",
> \"2020-09-28 12:38:26 UTC+1000 by era5~replicationtools~-1.2.1: nccopy
> -k4 -d5 -s
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/sp/2010/sp~era5opersfc20100101~-20100131.nc
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/sp/2010*.sp~era5opersfc20100101~-20100131.tmp`\n`{=latex}\",
> \"2020-09-28 12:34:44 UTC+1000 by era5~replicationtools~-1.2.1: curl
> --connect-timeout 20 --show-error --silent --max-time 36000 -o
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/sp/2010/sp~era5opersfc20100101~-20100131.nc
> [http://136.156.132.198/cache-compute-0003/cache/data4/adaptor.mars.internal-1601260248.6544905-25194-9-c80c2295-8d73-4f03-8c30-c28361fb17ff.nc\\n](http://136.156.132.198/cache-compute-0003/cache/data4/adaptor.mars.internal-1601260248.6544905-25194-9-c80c2295-8d73-4f03-8c30-c28361fb17ff.nc\n)\",
> \"2020-09-28 02:32:17 GMT by grib~tonetcdf~-2.16.0:
> /opt/ecmwf/eccodes/bin/grib~tonetcdf~ -S param -o
> /cache/data4/adaptor.mars.internal-1601260248.6544905-25194-9-c80c2295-8d73-4f03-8c30-c28361fb17ff.nc
> /cache/tmp/c80c2295-8d73-4f03-8c30-c28361fb17ff-adaptor.mars.internal-1601260248.6552908-25194-4-tmp.grib\"
> ; :license = \"Licence to use Copernicus Products:
> <https://apps.ecmwf.int/datasets/licences/copernicus/>\" ; :summary =
> \"ERA5 is the fifth generation ECMWF atmospheric reanalysis of the
> global climate. This file is part of the ERA5 replica hosted at NCI
> Australia. For more information please see
> <http://dx.doi.org/10.25914/5f48874388857>\" ; :title = \"ERA5
> single-levels reanalysis surface~pressure~ 20100101-20100131\" ; }

### [z](file:///Volumes/ioa02/reanalysis/ERA5/z)

Geopotential

> netcdf z~era5opersfc20100101~-20100131 { dimensions: longitude = 1440
> ; latitude = 721 ; time = 744 ; variables: float longitude(longitude)
> ; longitude:units = \"degrees~east~\" ; longitude:long~name~ =
> \"longitude\" ; float latitude(latitude) ; latitude:units =
> \"degrees~north~\" ; latitude:long~name~ = \"latitude\" ; int
> time(time) ; time:units = \"hours since 1900-01-01 00:00:00.0\" ;
> time:long~name~ = \"time\" ; time:calendar = \"gregorian\" ; short
> z(time, latitude, longitude) ; z:scale~factor~ = 0.897939420788 ;
> z:add~offset~ = 28035.3894091959 ; z:~FillValue~ = -32767s ;
> z:missing~value~ = -32767s ; z:units = \"m\*\*2 s\*\*-2\" ;
> z:long~name~ = \"Geopotential\" ; z:standard~name~ = \"geopotential\"
> ;
>
> // global attributes: :Conventions = \"CF-1.6\" ; :history =
> \"2020-09-24 12:45:20 UTC+1000 by era5~replicationtools~-1.2.1: mv
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/z/2010*.z~era5opersfc20100101~-20100131.tmp
> /g/data/rt52/era5/single-levels/reanalysis/z/2010/z~era5opersfc20100101~-20100131.nc`\n`{=latex}\",
> \"2020-09-24 12:45:06 UTC+1000 by era5~replicationtools~-1.2.1: nccopy
> -k4 -d5 -s
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/z/2010/z~era5opersfc20100101~-20100131.nc
> *g/data/rt52/admin/incoming/era5/single-levels/reanalysis/z/2010*.z~era5opersfc20100101~-20100131.tmp`\n`{=latex}\",
> \"2020-09-24 12:41:04 UTC+1000 by era5~replicationtools~-1.2.1: curl
> --connect-timeout 20 --show-error --silent --max-time 36000 -o
> /g/data/zt77/admin/incoming/era5/single-levels/reanalysis/z/2010/z~era5opersfc20100101~-20100131.nc
> [http://136.156.133.41/cache-compute-0013/cache/data7/adaptor.mars.internal-1600915045.6589952-3753-13-24833832-41d7-4def-96fb-c55f6af6e05a.nc\\n](http://136.156.133.41/cache-compute-0013/cache/data7/adaptor.mars.internal-1600915045.6589952-3753-13-24833832-41d7-4def-96fb-c55f6af6e05a.nc\n)\",
> \"2020-09-24 02:38:37 GMT by grib~tonetcdf~-2.16.0:
> /opt/ecmwf/eccodes/bin/grib~tonetcdf~ -S param -o
> /cache/data7/adaptor.mars.internal-1600915045.6589952-3753-13-24833832-41d7-4def-96fb-c55f6af6e05a.nc
> /cache/tmp/24833832-41d7-4def-96fb-c55f6af6e05a-adaptor.mars.internal-1600915045.660192-3753-4-tmp.grib\"
> ; :license = \"Licence to use Copernicus Products:
> <https://apps.ecmwf.int/datasets/licences/copernicus/>\" ; :summary =
> \"ERA5 is the fifth generation ECMWF atmospheric reanalysis of the
> global climate. This file is part of the ERA5 replica hosted at NCI
> Australia. For more information please see
> <http://dx.doi.org/10.25914/5f48874388857>\" ; :title = \"ERA5
> single-levels reanalysis geopotential 20100101-20100131\" ; }

## JRA55do

[JRA55do](https://climate.mri-jma.go.jp/pub/ocean/JRA55-do/) is a $1/2$
degree spatial and three-hour temporal resolution atmospheric climate
dataset. Version 1.5 has been downloaded to my [local
repository](file:///Volumes/ioa01/reanalysis/JRA55do/) for the following
parameters.

### [huss](file:///Volumes/ioa01/reanalysis/JRA55do/huss)

Near-surface specific humidity

> netcdf
> huss~input4MIPsatmosphericStateOMIPMRI~-JRA55-do-1-5-0~gr201001010000~-201012312100
> { dimensions: time = UNLIMITED ; // (2920 currently) lat = 320 ; lon =
> 640 ; bnds = 2 ; variables: double time(time) ; time:bounds =
> \"time~bnds~\" ; time:units = \"days since 1900-01-01 00:00:00\" ;
> time:calendar = \"gregorian\" ; time:axis = \"T\" ; time:long~name~ =
> \"time\" ; time:standard~name~ = \"time\" ; double time~bnds~(time,
> bnds) ; double lat(lat) ; lat:bounds = \"lat~bnds~\" ; lat:units =
> \"degrees~north~\" ; lat:axis = \"Y\" ; lat:long~name~ = \"Latitude\"
> ; lat:standard~name~ = \"latitude\" ; double lat~bnds~(lat, bnds) ;
> double lon(lon) ; lon:bounds = \"lon~bnds~\" ; lon:units =
> \"degrees~east~\" ; lon:axis = \"X\" ; lon:long~name~ = \"Longitude\"
> ; lon:standard~name~ = \"longitude\" ; double lon~bnds~(lon, bnds) ;
> double height ; height:units = \"m\" ; height:axis = \"Z\" ;
> height:positive = \"up\" ; height:long~name~ = \"height\" ;
> height:standard~name~ = \"height\" ; float huss(time, lat, lon) ;
> huss:standard~name~ = \"specific~humidity~\" ; huss:long~name~ =
> \"Near-Surface Specific Humidity\" ; huss:comment = \"Near-surface
> (usually, 2 meter) specific humidity\" ; huss:units = \"1\" ;
> huss:original~units~ = \"1.0\" ; huss:history = \"2020-09-15T15:01:10Z
> altered by CMOR: Converted units from \\\'1.0\\\' to \\\'1\\\'.
> 2020-09-15T15:01:10Z altered by CMOR: Treated scalar dimension:
> \\\'height\\\'.\" ; huss:cell~methods~ = \"area: mean time: point\" ;
> huss:cell~measures~ = \"area: areacella\" ; huss:coordinates =
> \"height\" ; huss:missing~value~ = 1.e+20f ; huss:~FillValue~ =
> 1.e+20f ;
>
> // global attributes: :Conventions = \"CF-1.7 CMIP-6.2\" ;
> :activity~id~ = \"input4MIPs\" ; :cell~measures~ = \"area: areacella\"
> ; :comment = \"Based on JRA-55 reanalysis (1958-01 to 2020-07)\" ;
> :contact = \"Hiroyuki Tsujino (htsujino@mri-jma.go.jp)\" ;
> :creation~date~ = \"2020-09-15T15:01:18Z\" ; :data~specsversion~ =
> \"01.00.32\" ; :dataset~category~ = \"atmosphericState\" ;
> :external~variables~ = \"areacella\" ; :frequency = \"3hrPt\" ;
> :further~infourl~ =
> \"<http://climate.mri-jma.go.jp/~htsujino/jra55do.html>\" ; :grid =
> \"data regridded to the normal atmosphere TL319 gaussian grid (320x640
> latxlon) from a reduced TL319 gaussian grid\" ; :grid~label~ = \"gr\"
> ; :history = \"2020-09-15T15:01:18Z; CMOR rewrote data to be
> consistent with input4MIPs, CMIP6, and CF-1.7 CMIP-6.2 standards\" ;
> :institution = \"Meteorological Research Institute, Tsukuba, Ibaraki
> 305-0052, Japan\" ; :institution~id~ = \"MRI\" ; :mip~era~ = \"CMIP6\"
> ; :nominal~resolution~ = \"50 km\" ; :product = \"reanalysis\" ;
> :realm = \"atmos\" ; :references = \"Tsujino et al., 2018: JRA-55
> based surface dataset for driving ocean-sea-ice models (JRA55-do),
> Ocean Modelling, 130(1), pp 79-139.
> <https://doi.org/10.1016/j.ocemod.2018.07.002>\" ; :region =
> \"global~ocean~\" ; :release~year~ = \"2020\" ; :source = \"MRI
> JRA55-do 1.5.0: Atmospheric state generated for OMIP based on the
> JRA-55 reanalysis\" ; :source~description~ = \"Atmospheric state and
> terrestrial runoff datasets produced by MRI for the OMIP experiment of
> CMIP6\" ; :source~id~ = \"MRI-JRA55-do-1-5-0\" ; :source~type~ =
> \"satellite~blended~\" ; :source~version~ = \"1.5.0\" ; :table~id~ =
> \"input4MIPs~A3hrPt~\" ; :table~info~ = \"Creation Date:(14 September
> 2020) MD5:d10deb7a45fa93934ec0a16c9d534e45\" ; :target~mip~ = \"OMIP\"
> ; :title = \"MRI JRA55-do 1.5.0 dataset prepared for input4MIPs\" ;
> :tracking~id~ = \"hdl:21.14100/369cad98-ed94-4040-a47d-72258fd68584\"
> ; :variable~id~ = \"huss\" ; :license = \"OMIP boundary condition data
> produced by MRI is licensed under a Creative Commons
> Attribution-\[NonCommercial-\]ShareAlike 4.0 International License
> (<https://creativecommons.org/licenses>). Consult
> <https://pcmdi.llnl.gov/CMIP6/TermsOfUse> for terms of use governing
> input4MIPs output, including citation requirements and proper
> acknowledgment. Further information about this data, including some
> limitations, can be found via the further~infourl~ (recorded as a
> global attribute in this file). The data producers and data providers
> make no warranty, either express or implied, including, but not
> limited to, warranties of merchantability and fitness for a particular
> purpose. All liabilities arising from the supply of the information
> (including any liability arising in negligence) are excluded to the
> fullest extent permitted by law.\" ; :cmor~version~ = \"3.6.0\" ; }

### [prra](file:///Volumes/ioa01/reanalysis/JRA55do/prra)

Rainfall flux

> netcdf
> prra~input4MIPsatmosphericStateOMIPMRI~-JRA55-do-1-5-0~gr201001010130~-201012312230
> { dimensions: time = UNLIMITED ; // (2920 currently) lat = 320 ; lon =
> 640 ; bnds = 2 ; variables: double time(time) ; time:bounds =
> \"time~bnds~\" ; time:units = \"days since 1900-01-01 00:00:00\" ;
> time:calendar = \"gregorian\" ; time:axis = \"T\" ; time:long~name~ =
> \"time\" ; time:standard~name~ = \"time\" ; double time~bnds~(time,
> bnds) ; double lat(lat) ; lat:bounds = \"lat~bnds~\" ; lat:units =
> \"degrees~north~\" ; lat:axis = \"Y\" ; lat:long~name~ = \"Latitude\"
> ; lat:standard~name~ = \"latitude\" ; double lat~bnds~(lat, bnds) ;
> double lon(lon) ; lon:bounds = \"lon~bnds~\" ; lon:units =
> \"degrees~east~\" ; lon:axis = \"X\" ; lon:long~name~ = \"Longitude\"
> ; lon:standard~name~ = \"longitude\" ; double lon~bnds~(lon, bnds) ;
> float prra(time, lat, lon) ; prra:standard~name~ = \"rainfall~flux~\"
> ; prra:long~name~ = \"Rainfall Flux\" ; prra:comment = \"In accordance
> with common usage in geophysical disciplines, \\\'flux\\\' implies per
> unit area, called \\\'flux density\\\' in physics\" ; prra:units =
> \"kg m-2 s-1\" ; prra:cell~methods~ = \"area: time: mean\" ;
> prra:cell~measures~ = \"area: areacella\" ; prra:missing~value~ =
> 1.e+20f ; prra:~FillValue~ = 1.e+20f ;
>
> // global attributes: :Conventions = \"CF-1.7 CMIP-6.2\" ;
> :activity~id~ = \"input4MIPs\" ; :cell~measures~ = \"area: areacella\"
> ; :comment = \"Based on JRA-55 reanalysis (1958-01 to 2020-07)\" ;
> :contact = \"Hiroyuki Tsujino (htsujino@mri-jma.go.jp)\" ;
> :creation~date~ = \"2020-09-15T06:03:44Z\" ; :data~specsversion~ =
> \"01.00.32\" ; :dataset~category~ = \"atmosphericState\" ;
> :external~variables~ = \"areacella\" ; :frequency = \"3hr\" ;
> :further~infourl~ =
> \"<http://climate.mri-jma.go.jp/~htsujino/jra55do.html>\" ; :grid =
> \"data regridded to the normal atmosphere TL319 gaussian grid (320x640
> latxlon) from a reduced TL319 gaussian grid\" ; :grid~label~ = \"gr\"
> ; :history = \"2020-09-15T06:03:44Z; CMOR rewrote data to be
> consistent with input4MIPs, CMIP6, and CF-1.7 CMIP-6.2 standards\" ;
> :institution = \"Meteorological Research Institute, Tsukuba, Ibaraki
> 305-0052, Japan\" ; :institution~id~ = \"MRI\" ; :mip~era~ = \"CMIP6\"
> ; :nominal~resolution~ = \"50 km\" ; :product = \"reanalysis\" ;
> :realm = \"atmos\" ; :references = \"Tsujino et al., 2018: JRA-55
> based surface dataset for driving ocean-sea-ice models (JRA55-do),
> Ocean Modelling, 130(1), pp 79-139.
> <https://doi.org/10.1016/j.ocemod.2018.07.002>\" ; :region =
> \"global~ocean~\" ; :release~year~ = \"2020\" ; :source = \"MRI
> JRA55-do 1.5.0: Atmospheric state generated for OMIP based on the
> JRA-55 reanalysis\" ; :source~description~ = \"Atmospheric state and
> terrestrial runoff datasets produced by MRI for the OMIP experiment of
> CMIP6\" ; :source~id~ = \"MRI-JRA55-do-1-5-0\" ; :source~type~ =
> \"satellite~blended~\" ; :source~version~ = \"1.5.0\" ; :table~id~ =
> \"input4MIPs~A3hr~\" ; :table~info~ = \"Creation Date:(14 September
> 2020) MD5:d10deb7a45fa93934ec0a16c9d534e45\" ; :target~mip~ = \"OMIP\"
> ; :title = \"MRI JRA55-do 1.5.0 dataset prepared for input4MIPs\" ;
> :tracking~id~ = \"hdl:21.14100/0ebc0818-9a24-47ad-a6d5-112506908d87\"
> ; :variable~id~ = \"prra\" ; :license = \"OMIP boundary condition data
> produced by MRI is licensed under a Creative Commons
> Attribution-\[NonCommercial-\]ShareAlike 4.0 International License
> (<https://creativecommons.org/licenses>). Consult
> <https://pcmdi.llnl.gov/CMIP6/TermsOfUse> for terms of use governing
> input4MIPs output, including citation requirements and proper
> acknowledgment. Further information about this data, including some
> limitations, can be found via the further~infourl~ (recorded as a
> global attribute in this file). The data producers and data providers
> make no warranty, either express or implied, including, but not
> limited to, warranties of merchantability and fitness for a particular
> purpose. All liabilities arising from the supply of the information
> (including any liability arising in negligence) are excluded to the
> fullest extent permitted by law.\" ; :cmor~version~ = \"3.6.0\" ; }

### [prsn](file:///Volumes/ioa01/reanalysis/JRA55do/prsn/)

Snowfall flux

> netcdf
> prsn~input4MIPsatmosphericStateOMIPMRI~-JRA55-do-1-5-0~gr201001010130~-201012312230
> { dimensions: time = UNLIMITED ; // (2920 currently) lat = 320 ; lon =
> 640 ; bnds = 2 ; variables: double time(time) ; time:bounds =
> \"time~bnds~\" ; time:units = \"days since 1900-01-01 00:00:00\" ;
> time:calendar = \"gregorian\" ; time:axis = \"T\" ; time:long~name~ =
> \"time\" ; time:standard~name~ = \"time\" ; double time~bnds~(time,
> bnds) ; double lat(lat) ; lat:bounds = \"lat~bnds~\" ; lat:units =
> \"degrees~north~\" ; lat:axis = \"Y\" ; lat:long~name~ = \"Latitude\"
> ; lat:standard~name~ = \"latitude\" ; double lat~bnds~(lat, bnds) ;
> double lon(lon) ; lon:bounds = \"lon~bnds~\" ; lon:units =
> \"degrees~east~\" ; lon:axis = \"X\" ; lon:long~name~ = \"Longitude\"
> ; lon:standard~name~ = \"longitude\" ; double lon~bnds~(lon, bnds) ;
> float prsn(time, lat, lon) ; prsn:standard~name~ = \"snowfall~flux~\"
> ; prsn:long~name~ = \"Snowfall Flux\" ; prsn:comment = \"At surface;
> includes precipitation of all forms of water in the solid phase\" ;
> prsn:units = \"kg m-2 s-1\" ; prsn:cell~methods~ = \"area: time:
> mean\" ; prsn:cell~measures~ = \"area: areacella\" ;
> prsn:missing~value~ = 1.e+20f ; prsn:~FillValue~ = 1.e+20f ;
>
> // global attributes: :Conventions = \"CF-1.7 CMIP-6.2\" ;
> :activity~id~ = \"input4MIPs\" ; :cell~measures~ = \"area: areacella\"
> ; :comment = \"Based on JRA-55 reanalysis (1958-01 to 2020-07)\" ;
> :contact = \"Hiroyuki Tsujino (htsujino@mri-jma.go.jp)\" ;
> :creation~date~ = \"2020-09-15T07:53:10Z\" ; :data~specsversion~ =
> \"01.00.32\" ; :dataset~category~ = \"atmosphericState\" ;
> :external~variables~ = \"areacella\" ; :frequency = \"3hr\" ;
> :further~infourl~ =
> \"<http://climate.mri-jma.go.jp/~htsujino/jra55do.html>\" ; :grid =
> \"data regridded to the normal atmosphere TL319 gaussian grid (320x640
> latxlon) from a reduced TL319 gaussian grid\" ; :grid~label~ = \"gr\"
> ; :history = \"2020-09-15T07:53:10Z; CMOR rewrote data to be
> consistent with input4MIPs, CMIP6, and CF-1.7 CMIP-6.2 standards\" ;
> :institution = \"Meteorological Research Institute, Tsukuba, Ibaraki
> 305-0052, Japan\" ; :institution~id~ = \"MRI\" ; :mip~era~ = \"CMIP6\"
> ; :nominal~resolution~ = \"50 km\" ; :product = \"reanalysis\" ;
> :realm = \"atmos\" ; :references = \"Tsujino et al., 2018: JRA-55
> based surface dataset for driving ocean-sea-ice models (JRA55-do),
> Ocean Modelling, 130(1), pp 79-139.
> <https://doi.org/10.1016/j.ocemod.2018.07.002>\" ; :region =
> \"global~ocean~\" ; :release~year~ = \"2020\" ; :source = \"MRI
> JRA55-do 1.5.0: Atmospheric state generated for OMIP based on the
> JRA-55 reanalysis\" ; :source~description~ = \"Atmospheric state and
> terrestrial runoff datasets produced by MRI for the OMIP experiment of
> CMIP6\" ; :source~id~ = \"MRI-JRA55-do-1-5-0\" ; :source~type~ =
> \"satellite~blended~\" ; :source~version~ = \"1.5.0\" ; :table~id~ =
> \"input4MIPs~A3hr~\" ; :table~info~ = \"Creation Date:(14 September
> 2020) MD5:d10deb7a45fa93934ec0a16c9d534e45\" ; :target~mip~ = \"OMIP\"
> ; :title = \"MRI JRA55-do 1.5.0 dataset prepared for input4MIPs\" ;
> :tracking~id~ = \"hdl:21.14100/11edbd94-09f9-4e92-b307-c69ee1238c4e\"
> ; :variable~id~ = \"prsn\" ; :license = \"OMIP boundary condition data
> produced by MRI is licensed under a Creative Commons
> Attribution-\[NonCommercial-\]ShareAlike 4.0 International License
> (<https://creativecommons.org/licenses>). Consult
> <https://pcmdi.llnl.gov/CMIP6/TermsOfUse> for terms of use governing
> input4MIPs output, including citation requirements and proper
> acknowledgment. Further information about this data, including some
> limitations, can be found via the further~infourl~ (recorded as a
> global attribute in this file). The data producers and data providers
> make no warranty, either express or implied, including, but not
> limited to, warranties of merchantability and fitness for a particular
> purpose. All liabilities arising from the supply of the information
> (including any liability arising in negligence) are excluded to the
> fullest extent permitted by law.\" ; :cmor~version~ = \"3.6.0\" ; }

### [psl](file:///Volumes/ioa01/reanalysis/JRA55do/psl)

Mean surface pressure

> netcdf
> psl~input4MIPsatmosphericStateOMIPMRI~-JRA55-do-1-5-0~gr201001010000~-201012312100
> { dimensions: time = UNLIMITED ; // (2920 currently) lat = 320 ; lon =
> 640 ; bnds = 2 ; variables: double time(time) ; time:bounds =
> \"time~bnds~\" ; time:units = \"days since 1900-01-01 00:00:00\" ;
> time:calendar = \"gregorian\" ; time:axis = \"T\" ; time:long~name~ =
> \"time\" ; time:standard~name~ = \"time\" ; double time~bnds~(time,
> bnds) ; double lat(lat) ; lat:bounds = \"lat~bnds~\" ; lat:units =
> \"degrees~north~\" ; lat:axis = \"Y\" ; lat:long~name~ = \"Latitude\"
> ; lat:standard~name~ = \"latitude\" ; double lat~bnds~(lat, bnds) ;
> double lon(lon) ; lon:bounds = \"lon~bnds~\" ; lon:units =
> \"degrees~east~\" ; lon:axis = \"X\" ; lon:long~name~ = \"Longitude\"
> ; lon:standard~name~ = \"longitude\" ; double lon~bnds~(lon, bnds) ;
> float psl(time, lat, lon) ; psl:standard~name~ =
> \"air~pressureatmeansealevel~\" ; psl:long~name~ = \"Sea Level
> Pressure\" ; psl:comment = \"Sea Level Pressure\" ; psl:units = \"Pa\"
> ; psl:cell~methods~ = \"area: mean time: point\" ; psl:cell~measures~
> = \"area: areacella\" ; psl:missing~value~ = 1.e+20f ; psl:~FillValue~
> = 1.e+20f ;
>
> // global attributes: :Conventions = \"CF-1.7 CMIP-6.2\" ;
> :activity~id~ = \"input4MIPs\" ; :cell~measures~ = \"area: areacella\"
> ; :comment = \"Based on JRA-55 reanalysis (1958-01 to 2020-07)\" ;
> :contact = \"Hiroyuki Tsujino (htsujino@mri-jma.go.jp)\" ;
> :creation~date~ = \"2020-09-15T17:20:36Z\" ; :data~specsversion~ =
> \"01.00.32\" ; :dataset~category~ = \"atmosphericState\" ;
> :external~variables~ = \"areacella\" ; :frequency = \"3hrPt\" ;
> :further~infourl~ =
> \"<http://climate.mri-jma.go.jp/~htsujino/jra55do.html>\" ; :grid =
> \"data regridded to the normal atmosphere TL319 gaussian grid (320x640
> latxlon) from a reduced TL319 gaussian grid\" ; :grid~label~ = \"gr\"
> ; :history = \"2020-09-15T17:20:36Z; CMOR rewrote data to be
> consistent with input4MIPs, CMIP6, and CF-1.7 CMIP-6.2 standards\" ;
> :institution = \"Meteorological Research Institute, Tsukuba, Ibaraki
> 305-0052, Japan\" ; :institution~id~ = \"MRI\" ; :mip~era~ = \"CMIP6\"
> ; :nominal~resolution~ = \"50 km\" ; :product = \"reanalysis\" ;
> :realm = \"atmos\" ; :references = \"Tsujino et al., 2018: JRA-55
> based surface dataset for driving ocean-sea-ice models (JRA55-do),
> Ocean Modelling, 130(1), pp 79-139.
> <https://doi.org/10.1016/j.ocemod.2018.07.002>\" ; :region =
> \"global~ocean~\" ; :release~year~ = \"2020\" ; :source = \"MRI
> JRA55-do 1.5.0: Atmospheric state generated for OMIP based on the
> JRA-55 reanalysis\" ; :source~description~ = \"Atmospheric state and
> terrestrial runoff datasets produced by MRI for the OMIP experiment of
> CMIP6\" ; :source~id~ = \"MRI-JRA55-do-1-5-0\" ; :source~type~ =
> \"satellite~blended~\" ; :source~version~ = \"1.5.0\" ; :table~id~ =
> \"input4MIPs~A3hrPt~\" ; :table~info~ = \"Creation Date:(14 September
> 2020) MD5:d10deb7a45fa93934ec0a16c9d534e45\" ; :target~mip~ = \"OMIP\"
> ; :title = \"MRI JRA55-do 1.5.0 dataset prepared for input4MIPs\" ;
> :tracking~id~ = \"hdl:21.14100/4e5646ac-45d5-40cd-a999-1554bb4cd75a\"
> ; :variable~id~ = \"psl\" ; :license = \"OMIP boundary condition data
> produced by MRI is licensed under a Creative Commons
> Attribution-\[NonCommercial-\]ShareAlike 4.0 International License
> (<https://creativecommons.org/licenses>). Consult
> <https://pcmdi.llnl.gov/CMIP6/TermsOfUse> for terms of use governing
> input4MIPs output, including citation requirements and proper
> acknowledgment. Further information about this data, including some
> limitations, can be found via the further~infourl~ (recorded as a
> global attribute in this file). The data producers and data providers
> make no warranty, either express or implied, including, but not
> limited to, warranties of merchantability and fitness for a particular
> purpose. All liabilities arising from the supply of the information
> (including any liability arising in negligence) are excluded to the
> fullest extent permitted by law.\" ; :cmor~version~ = \"3.6.0\" ; }

### [rlds](file:///Volumes/ioa01/reanalysis/JRA55do/rlds)

Surface downwelling longwave flux

> netcdf
> rlds~input4MIPsatmosphericStateOMIPMRI~-JRA55-do-1-5-0~gr201001010130~-201012312230
> { dimensions: time = UNLIMITED ; // (2920 currently) lat = 320 ; lon =
> 640 ; bnds = 2 ; variables: double time(time) ; time:bounds =
> \"time~bnds~\" ; time:units = \"days since 1900-01-01 00:00:00\" ;
> time:calendar = \"gregorian\" ; time:axis = \"T\" ; time:long~name~ =
> \"time\" ; time:standard~name~ = \"time\" ; double time~bnds~(time,
> bnds) ; double lat(lat) ; lat:bounds = \"lat~bnds~\" ; lat:units =
> \"degrees~north~\" ; lat:axis = \"Y\" ; lat:long~name~ = \"Latitude\"
> ; lat:standard~name~ = \"latitude\" ; double lat~bnds~(lat, bnds) ;
> double lon(lon) ; lon:bounds = \"lon~bnds~\" ; lon:units =
> \"degrees~east~\" ; lon:axis = \"X\" ; lon:long~name~ = \"Longitude\"
> ; lon:standard~name~ = \"longitude\" ; double lon~bnds~(lon, bnds) ;
> float rlds(time, lat, lon) ; rlds:standard~name~ =
> \"surface~downwellinglongwavefluxinair~\" ; rlds:long~name~ =
> \"Surface Downwelling Longwave Radiation\" ; rlds:comment = \"The
> surface called \\\'surface\\\' means the lower boundary of the
> atmosphere. \\\'longwave\\\' means longwave radiation. Downwelling
> radiation is radiation from above. It does not mean \\\'net
> downward\\\'. When thought of as being incident on a surface, a
> radiative flux is sometimes called \\\'irradiance\\\'. In addition, it
> is identical with the quantity measured by a cosine-collector
> light-meter and sometimes called \\\'vector irradiance\\\'. In
> accordance with common usage in geophysical disciplines, \\\'flux\\\'
> implies per unit area, called \\\'flux density\\\' in physics.\" ;
> rlds:units = \"W m-2\" ; rlds:cell~methods~ = \"area: time: mean\" ;
> rlds:cell~measures~ = \"area: areacella\" ; rlds:missing~value~ =
> 1.e+20f ; rlds:~FillValue~ = 1.e+20f ;
>
> // global attributes: :Conventions = \"CF-1.7 CMIP-6.2\" ;
> :activity~id~ = \"input4MIPs\" ; :cell~measures~ = \"area: areacella\"
> ; :comment = \"Based on JRA-55 reanalysis (1958-01 to 2020-07)\" ;
> :contact = \"Hiroyuki Tsujino (htsujino@mri-jma.go.jp)\" ;
> :creation~date~ = \"2020-09-15T10:15:08Z\" ; :data~specsversion~ =
> \"01.00.32\" ; :dataset~category~ = \"atmosphericState\" ;
> :external~variables~ = \"areacella\" ; :frequency = \"3hr\" ;
> :further~infourl~ =
> \"<http://climate.mri-jma.go.jp/~htsujino/jra55do.html>\" ; :grid =
> \"data regridded to the normal atmosphere TL319 gaussian grid (320x640
> latxlon) from a reduced TL319 gaussian grid\" ; :grid~label~ = \"gr\"
> ; :history = \"2020-09-15T10:15:08Z; CMOR rewrote data to be
> consistent with input4MIPs, CMIP6, and CF-1.7 CMIP-6.2 standards\" ;
> :institution = \"Meteorological Research Institute, Tsukuba, Ibaraki
> 305-0052, Japan\" ; :institution~id~ = \"MRI\" ; :mip~era~ = \"CMIP6\"
> ; :nominal~resolution~ = \"50 km\" ; :product = \"reanalysis\" ;
> :realm = \"atmos\" ; :references = \"Tsujino et al., 2018: JRA-55
> based surface dataset for driving ocean-sea-ice models (JRA55-do),
> Ocean Modelling, 130(1), pp 79-139.
> <https://doi.org/10.1016/j.ocemod.2018.07.002>\" ; :region =
> \"global~ocean~\" ; :release~year~ = \"2020\" ; :source = \"MRI
> JRA55-do 1.5.0: Atmospheric state generated for OMIP based on the
> JRA-55 reanalysis\" ; :source~description~ = \"Atmospheric state and
> terrestrial runoff datasets produced by MRI for the OMIP experiment of
> CMIP6\" ; :source~id~ = \"MRI-JRA55-do-1-5-0\" ; :source~type~ =
> \"satellite~blended~\" ; :source~version~ = \"1.5.0\" ; :table~id~ =
> \"input4MIPs~A3hr~\" ; :table~info~ = \"Creation Date:(14 September
> 2020) MD5:d10deb7a45fa93934ec0a16c9d534e45\" ; :target~mip~ = \"OMIP\"
> ; :title = \"MRI JRA55-do 1.5.0 dataset prepared for input4MIPs\" ;
> :tracking~id~ = \"hdl:21.14100/9a1230a5-bb91-4cb5-9c99-179ba031866e\"
> ; :variable~id~ = \"rlds\" ; :license = \"OMIP boundary condition data
> produced by MRI is licensed under a Creative Commons
> Attribution-\[NonCommercial-\]ShareAlike 4.0 International License
> (<https://creativecommons.org/licenses>). Consult
> <https://pcmdi.llnl.gov/CMIP6/TermsOfUse> for terms of use governing
> input4MIPs output, including citation requirements and proper
> acknowledgment. Further information about this data, including some
> limitations, can be found via the further~infourl~ (recorded as a
> global attribute in this file). The data producers and data providers
> make no warranty, either express or implied, including, but not
> limited to, warranties of merchantability and fitness for a particular
> purpose. All liabilities arising from the supply of the information
> (including any liability arising in negligence) are excluded to the
> fullest extent permitted by law.\" ; :cmor~version~ = \"3.6.0\" ; }

### [rsds](file:///Volumes/ioa01/reanalysis/JRA55do/rsds)

Surface downwelling shortwave flux

> netcdf
> rsds~input4MIPsatmosphericStateOMIPMRI~-JRA55-do-1-5-0~gr201001010130~-201012312230
> { dimensions: time = UNLIMITED ; // (2920 currently) lat = 320 ; lon =
> 640 ; bnds = 2 ; variables: double time(time) ; time:bounds =
> \"time~bnds~\" ; time:units = \"days since 1900-01-01 00:00:00\" ;
> time:calendar = \"gregorian\" ; time:axis = \"T\" ; time:long~name~ =
> \"time\" ; time:standard~name~ = \"time\" ; double time~bnds~(time,
> bnds) ; double lat(lat) ; lat:bounds = \"lat~bnds~\" ; lat:units =
> \"degrees~north~\" ; lat:axis = \"Y\" ; lat:long~name~ = \"Latitude\"
> ; lat:standard~name~ = \"latitude\" ; double lat~bnds~(lat, bnds) ;
> double lon(lon) ; lon:bounds = \"lon~bnds~\" ; lon:units =
> \"degrees~east~\" ; lon:axis = \"X\" ; lon:long~name~ = \"Longitude\"
> ; lon:standard~name~ = \"longitude\" ; double lon~bnds~(lon, bnds) ;
> float rsds(time, lat, lon) ; rsds:standard~name~ =
> \"surface~downwellingshortwavefluxinair~\" ; rsds:long~name~ =
> \"Surface Downwelling Shortwave Radiation\" ; rsds:comment = \"Surface
> solar irradiance for UV calculations.\" ; rsds:units = \"W m-2\" ;
> rsds:cell~methods~ = \"area: time: mean\" ; rsds:cell~measures~ =
> \"area: areacella\" ; rsds:missing~value~ = 1.e+20f ; rsds:~FillValue~
> = 1.e+20f ;
>
> // global attributes: :Conventions = \"CF-1.7 CMIP-6.2\" ;
> :activity~id~ = \"input4MIPs\" ; :cell~measures~ = \"area: areacella\"
> ; :comment = \"Based on JRA-55 reanalysis (1958-01 to 2020-07)\" ;
> :contact = \"Hiroyuki Tsujino (htsujino@mri-jma.go.jp)\" ;
> :creation~date~ = \"2020-09-15T12:28:57Z\" ; :data~specsversion~ =
> \"01.00.32\" ; :dataset~category~ = \"atmosphericState\" ;
> :external~variables~ = \"areacella\" ; :frequency = \"3hr\" ;
> :further~infourl~ =
> \"<http://climate.mri-jma.go.jp/~htsujino/jra55do.html>\" ; :grid =
> \"data regridded to the normal atmosphere TL319 gaussian grid (320x640
> latxlon) from a reduced TL319 gaussian grid\" ; :grid~label~ = \"gr\"
> ; :history = \"2020-09-15T12:28:57Z; CMOR rewrote data to be
> consistent with input4MIPs, CMIP6, and CF-1.7 CMIP-6.2 standards\" ;
> :institution = \"Meteorological Research Institute, Tsukuba, Ibaraki
> 305-0052, Japan\" ; :institution~id~ = \"MRI\" ; :mip~era~ = \"CMIP6\"
> ; :nominal~resolution~ = \"50 km\" ; :product = \"reanalysis\" ;
> :realm = \"atmos\" ; :references = \"Tsujino et al., 2018: JRA-55
> based surface dataset for driving ocean-sea-ice models (JRA55-do),
> Ocean Modelling, 130(1), pp 79-139.
> <https://doi.org/10.1016/j.ocemod.2018.07.002>\" ; :region =
> \"global~ocean~\" ; :release~year~ = \"2020\" ; :source = \"MRI
> JRA55-do 1.5.0: Atmospheric state generated for OMIP based on the
> JRA-55 reanalysis\" ; :source~description~ = \"Atmospheric state and
> terrestrial runoff datasets produced by MRI for the OMIP experiment of
> CMIP6\" ; :source~id~ = \"MRI-JRA55-do-1-5-0\" ; :source~type~ =
> \"satellite~blended~\" ; :source~version~ = \"1.5.0\" ; :table~id~ =
> \"input4MIPs~A3hr~\" ; :table~info~ = \"Creation Date:(14 September
> 2020) MD5:d10deb7a45fa93934ec0a16c9d534e45\" ; :target~mip~ = \"OMIP\"
> ; :title = \"MRI JRA55-do 1.5.0 dataset prepared for input4MIPs\" ;
> :tracking~id~ = \"hdl:21.14100/f376d4e4-6396-44b4-ab8b-c7220a1afcb3\"
> ; :variable~id~ = \"rsds\" ; :license = \"OMIP boundary condition data
> produced by MRI is licensed under a Creative Commons
> Attribution-\[NonCommercial-\]ShareAlike 4.0 International License
> (<https://creativecommons.org/licenses>). Consult
> <https://pcmdi.llnl.gov/CMIP6/TermsOfUse> for terms of use governing
> input4MIPs output, including citation requirements and proper
> acknowledgment. Further information about this data, including some
> limitations, can be found via the further~infourl~ (recorded as a
> global attribute in this file). The data producers and data providers
> make no warranty, either express or implied, including, but not
> limited to, warranties of merchantability and fitness for a particular
> purpose. All liabilities arising from the supply of the information
> (including any liability arising in negligence) are excluded to the
> fullest extent permitted by law.\" ; :cmor~version~ = \"3.6.0\" ; }

### [tas](file:///Volumes/ioa01/reanalysis/JRA55do/tas)

Near-surface air temperature

> netcdf
> tas~input4MIPsatmosphericStateOMIPMRI~-JRA55-do-1-5-0~gr201001010000~-201012312100
> { dimensions: time = UNLIMITED ; // (2920 currently) lat = 320 ; lon =
> 640 ; bnds = 2 ; variables: double time(time) ; time:bounds =
> \"time~bnds~\" ; time:units = \"days since 1900-01-01 00:00:00\" ;
> time:calendar = \"gregorian\" ; time:axis = \"T\" ; time:long~name~ =
> \"time\" ; time:standard~name~ = \"time\" ; double time~bnds~(time,
> bnds) ; double lat(lat) ; lat:bounds = \"lat~bnds~\" ; lat:units =
> \"degrees~north~\" ; lat:axis = \"Y\" ; lat:long~name~ = \"Latitude\"
> ; lat:standard~name~ = \"latitude\" ; double lat~bnds~(lat, bnds) ;
> double lon(lon) ; lon:bounds = \"lon~bnds~\" ; lon:units =
> \"degrees~east~\" ; lon:axis = \"X\" ; lon:long~name~ = \"Longitude\"
> ; lon:standard~name~ = \"longitude\" ; double lon~bnds~(lon, bnds) ;
> double height ; height:units = \"m\" ; height:axis = \"Z\" ;
> height:positive = \"up\" ; height:long~name~ = \"height\" ;
> height:standard~name~ = \"height\" ; float tas(time, lat, lon) ;
> tas:standard~name~ = \"air~temperature~\" ; tas:long~name~ =
> \"Near-Surface Air Temperature\" ; tas:comment = \"near-surface
> (usually, 2 meter) air temperature\" ; tas:units = \"K\" ;
> tas:cell~methods~ = \"area: mean time: point\" ; tas:cell~measures~ =
> \"area: areacella\" ; tas:history = \"2020-09-15T19:35:02Z altered by
> CMOR: Treated scalar dimension: \\\'height\\\'.\" ; tas:coordinates =
> \"height\" ; tas:missing~value~ = 1.e+20f ; tas:~FillValue~ = 1.e+20f
> ;
>
> // global attributes: :Conventions = \"CF-1.7 CMIP-6.2\" ;
> :activity~id~ = \"input4MIPs\" ; :cell~measures~ = \"area: areacella\"
> ; :comment = \"Based on JRA-55 reanalysis (1958-01 to 2020-07)\" ;
> :contact = \"Hiroyuki Tsujino (htsujino@mri-jma.go.jp)\" ;
> :creation~date~ = \"2020-09-15T19:35:07Z\" ; :data~specsversion~ =
> \"01.00.32\" ; :dataset~category~ = \"atmosphericState\" ;
> :external~variables~ = \"areacella\" ; :frequency = \"3hrPt\" ;
> :further~infourl~ =
> \"<http://climate.mri-jma.go.jp/~htsujino/jra55do.html>\" ; :grid =
> \"data regridded to the normal atmosphere TL319 gaussian grid (320x640
> latxlon) from a reduced TL319 gaussian grid\" ; :grid~label~ = \"gr\"
> ; :history = \"2020-09-15T19:35:07Z; CMOR rewrote data to be
> consistent with input4MIPs, CMIP6, and CF-1.7 CMIP-6.2 standards\" ;
> :institution = \"Meteorological Research Institute, Tsukuba, Ibaraki
> 305-0052, Japan\" ; :institution~id~ = \"MRI\" ; :mip~era~ = \"CMIP6\"
> ; :nominal~resolution~ = \"50 km\" ; :product = \"reanalysis\" ;
> :realm = \"atmos\" ; :references = \"Tsujino et al., 2018: JRA-55
> based surface dataset for driving ocean-sea-ice models (JRA55-do),
> Ocean Modelling, 130(1), pp 79-139.
> <https://doi.org/10.1016/j.ocemod.2018.07.002>\" ; :region =
> \"global~ocean~\" ; :release~year~ = \"2020\" ; :source = \"MRI
> JRA55-do 1.5.0: Atmospheric state generated for OMIP based on the
> JRA-55 reanalysis\" ; :source~description~ = \"Atmospheric state and
> terrestrial runoff datasets produced by MRI for the OMIP experiment of
> CMIP6\" ; :source~id~ = \"MRI-JRA55-do-1-5-0\" ; :source~type~ =
> \"satellite~blended~\" ; :source~version~ = \"1.5.0\" ; :table~id~ =
> \"input4MIPs~A3hrPt~\" ; :table~info~ = \"Creation Date:(14 September
> 2020) MD5:d10deb7a45fa93934ec0a16c9d534e45\" ; :target~mip~ = \"OMIP\"
> ; :title = \"MRI JRA55-do 1.5.0 dataset prepared for input4MIPs\" ;
> :tracking~id~ = \"hdl:21.14100/a0582304-8b18-4297-b548-f33f0b40eccb\"
> ; :variable~id~ = \"tas\" ; :license = \"OMIP boundary condition data
> produced by MRI is licensed under a Creative Commons
> Attribution-\[NonCommercial-\]ShareAlike 4.0 International License
> (<https://creativecommons.org/licenses>). Consult
> <https://pcmdi.llnl.gov/CMIP6/TermsOfUse> for terms of use governing
> input4MIPs output, including citation requirements and proper
> acknowledgment. Further information about this data, including some
> limitations, can be found via the further~infourl~ (recorded as a
> global attribute in this file). The data producers and data providers
> make no warranty, either express or implied, including, but not
> limited to, warranties of merchantability and fitness for a particular
> purpose. All liabilities arising from the supply of the information
> (including any liability arising in negligence) are excluded to the
> fullest extent permitted by law.\" ; :cmor~version~ = \"3.6.0\" ; }

### [ts](file:///Volumes/ioa01/reanalysis/JRA55do/ts)

Surface temperature

> netcdf
> ts~input4MIPsatmosphericStateOMIPMRI~-JRA55-do-1-5-0~gr201001010000~-201012312100
> { dimensions: time = UNLIMITED ; // (2920 currently) lat = 320 ; lon =
> 640 ; bnds = 2 ; variables: double time(time) ; time:bounds =
> \"time~bnds~\" ; time:units = \"days since 1900-01-01 00:00:00\" ;
> time:calendar = \"gregorian\" ; time:axis = \"T\" ; time:long~name~ =
> \"time\" ; time:standard~name~ = \"time\" ; double time~bnds~(time,
> bnds) ; double lat(lat) ; lat:bounds = \"lat~bnds~\" ; lat:units =
> \"degrees~north~\" ; lat:axis = \"Y\" ; lat:long~name~ = \"Latitude\"
> ; lat:standard~name~ = \"latitude\" ; double lat~bnds~(lat, bnds) ;
> double lon(lon) ; lon:bounds = \"lon~bnds~\" ; lon:units =
> \"degrees~east~\" ; lon:axis = \"X\" ; lon:long~name~ = \"Longitude\"
> ; lon:standard~name~ = \"longitude\" ; double lon~bnds~(lon, bnds) ;
> float ts(time, lat, lon) ; ts:standard~name~ =
> \"surface~temperature~\" ; ts:long~name~ = \"Surface Temperature\" ;
> ts:comment = \"Temperature of the lower boundary of the atmosphere\" ;
> ts:units = \"K\" ; ts:cell~methods~ = \"area: mean time: point\" ;
> ts:cell~measures~ = \"area: areacella\" ; ts:missing~value~ = 1.e+20f
> ; ts:~FillValue~ = 1.e+20f ;
>
> // global attributes: :Conventions = \"CF-1.7 CMIP-6.2\" ;
> :activity~id~ = \"input4MIPs\" ; :cell~measures~ = \"area: areacella\"
> ; :comment = \"Based on JRA-55 reanalysis (1958-01 to 2020-07)\" ;
> :contact = \"Hiroyuki Tsujino (htsujino@mri-jma.go.jp)\" ;
> :creation~date~ = \"2020-09-15T21:48:05Z\" ; :data~specsversion~ =
> \"01.00.32\" ; :dataset~category~ = \"atmosphericState\" ;
> :external~variables~ = \"areacella\" ; :frequency = \"3hrPt\" ;
> :further~infourl~ =
> \"<http://climate.mri-jma.go.jp/~htsujino/jra55do.html>\" ; :grid =
> \"data regridded to the normal atmosphere TL319 gaussian grid (320x640
> latxlon) from a reduced TL319 gaussian grid\" ; :grid~label~ = \"gr\"
> ; :history = \"2020-09-15T21:48:05Z; CMOR rewrote data to be
> consistent with input4MIPs, CMIP6, and CF-1.7 CMIP-6.2 standards\" ;
> :institution = \"Meteorological Research Institute, Tsukuba, Ibaraki
> 305-0052, Japan\" ; :institution~id~ = \"MRI\" ; :mip~era~ = \"CMIP6\"
> ; :nominal~resolution~ = \"50 km\" ; :product = \"reanalysis\" ;
> :realm = \"atmos\" ; :references = \"Tsujino et al., 2018: JRA-55
> based surface dataset for driving ocean-sea-ice models (JRA55-do),
> Ocean Modelling, 130(1), pp 79-139.
> <https://doi.org/10.1016/j.ocemod.2018.07.002>\" ; :region =
> \"global~ocean~\" ; :release~year~ = \"2020\" ; :source = \"MRI
> JRA55-do 1.5.0: Atmospheric state generated for OMIP based on the
> JRA-55 reanalysis\" ; :source~description~ = \"Atmospheric state and
> terrestrial runoff datasets produced by MRI for the OMIP experiment of
> CMIP6\" ; :source~id~ = \"MRI-JRA55-do-1-5-0\" ; :source~type~ =
> \"satellite~blended~\" ; :source~version~ = \"1.5.0\" ; :table~id~ =
> \"input4MIPs~A3hrPt~\" ; :table~info~ = \"Creation Date:(14 September
> 2020) MD5:d10deb7a45fa93934ec0a16c9d534e45\" ; :target~mip~ = \"OMIP\"
> ; :title = \"MRI JRA55-do 1.5.0 dataset prepared for input4MIPs\" ;
> :tracking~id~ = \"hdl:21.14100/d725c7b9-ce5a-4fe7-8f65-3003acf15323\"
> ; :variable~id~ = \"ts\" ; :license = \"OMIP boundary condition data
> produced by MRI is licensed under a Creative Commons
> Attribution-\[NonCommercial-\]ShareAlike 4.0 International License
> (<https://creativecommons.org/licenses>). Consult
> <https://pcmdi.llnl.gov/CMIP6/TermsOfUse> for terms of use governing
> input4MIPs output, including citation requirements and proper
> acknowledgment. Further information about this data, including some
> limitations, can be found via the further~infourl~ (recorded as a
> global attribute in this file). The data producers and data providers
> make no warranty, either express or implied, including, but not
> limited to, warranties of merchantability and fitness for a particular
> purpose. All liabilities arising from the supply of the information
> (including any liability arising in negligence) are excluded to the
> fullest extent permitted by law.\" ; :cmor~version~ = \"3.6.0\" ; }

### [uas](file:///Volumes/ioa01/reanalysis/JRA55do/uas)

Eastward near-surface wind

> netcdf
> uas~input4MIPsatmosphericStateOMIPMRI~-JRA55-do-1-5-0~gr201001010000~-201012312100
> { dimensions: time = UNLIMITED ; // (2920 currently) lat = 320 ; lon =
> 640 ; bnds = 2 ; variables: double time(time) ; time:bounds =
> \"time~bnds~\" ; time:units = \"days since 1900-01-01 00:00:00\" ;
> time:calendar = \"gregorian\" ; time:axis = \"T\" ; time:long~name~ =
> \"time\" ; time:standard~name~ = \"time\" ; double time~bnds~(time,
> bnds) ; double lat(lat) ; lat:bounds = \"lat~bnds~\" ; lat:units =
> \"degrees~north~\" ; lat:axis = \"Y\" ; lat:long~name~ = \"Latitude\"
> ; lat:standard~name~ = \"latitude\" ; double lat~bnds~(lat, bnds) ;
> double lon(lon) ; lon:bounds = \"lon~bnds~\" ; lon:units =
> \"degrees~east~\" ; lon:axis = \"X\" ; lon:long~name~ = \"Longitude\"
> ; lon:standard~name~ = \"longitude\" ; double lon~bnds~(lon, bnds) ;
> double height ; height:units = \"m\" ; height:axis = \"Z\" ;
> height:positive = \"up\" ; height:long~name~ = \"height\" ;
> height:standard~name~ = \"height\" ; float uas(time, lat, lon) ;
> uas:standard~name~ = \"eastward~wind~\" ; uas:long~name~ = \"Eastward
> Near-Surface Wind\" ; uas:comment = \"Eastward component of the
> near-surface wind\" ; uas:units = \"m s-1\" ; uas:cell~methods~ =
> \"area: mean time: point\" ; uas:cell~measures~ = \"area: areacella\"
> ; uas:history = \"2020-09-16T01:30:40Z altered by CMOR: Treated scalar
> dimension: \\\'height\\\'.\" ; uas:coordinates = \"height\" ;
> uas:missing~value~ = 1.e+20f ; uas:~FillValue~ = 1.e+20f ;
>
> // global attributes: :Conventions = \"CF-1.7 CMIP-6.2\" ;
> :activity~id~ = \"input4MIPs\" ; :cell~measures~ = \"area: areacella\"
> ; :comment = \"Based on JRA-55 reanalysis (1958-01 to 2020-07)\" ;
> :contact = \"Hiroyuki Tsujino (htsujino@mri-jma.go.jp)\" ;
> :creation~date~ = \"2020-09-16T01:30:45Z\" ; :data~specsversion~ =
> \"01.00.32\" ; :dataset~category~ = \"atmosphericState\" ;
> :external~variables~ = \"areacella\" ; :frequency = \"3hrPt\" ;
> :further~infourl~ =
> \"<http://climate.mri-jma.go.jp/~htsujino/jra55do.html>\" ; :grid =
> \"data regridded to the normal atmosphere TL319 gaussian grid (320x640
> latxlon) from a reduced TL319 gaussian grid\" ; :grid~label~ = \"gr\"
> ; :history = \"2020-09-16T01:30:45Z; CMOR rewrote data to be
> consistent with input4MIPs, CMIP6, and CF-1.7 CMIP-6.2 standards\" ;
> :institution = \"Meteorological Research Institute, Tsukuba, Ibaraki
> 305-0052, Japan\" ; :institution~id~ = \"MRI\" ; :mip~era~ = \"CMIP6\"
> ; :nominal~resolution~ = \"50 km\" ; :product = \"reanalysis\" ;
> :realm = \"atmos\" ; :references = \"Tsujino et al., 2018: JRA-55
> based surface dataset for driving ocean-sea-ice models (JRA55-do),
> Ocean Modelling, 130(1), pp 79-139.
> <https://doi.org/10.1016/j.ocemod.2018.07.002>\" ; :region =
> \"global~ocean~\" ; :release~year~ = \"2020\" ; :source = \"MRI
> JRA55-do 1.5.0: Atmospheric state generated for OMIP based on the
> JRA-55 reanalysis\" ; :source~description~ = \"Atmospheric state and
> terrestrial runoff datasets produced by MRI for the OMIP experiment of
> CMIP6\" ; :source~id~ = \"MRI-JRA55-do-1-5-0\" ; :source~type~ =
> \"satellite~blended~\" ; :source~version~ = \"1.5.0\" ; :table~id~ =
> \"input4MIPs~A3hrPt~\" ; :table~info~ = \"Creation Date:(14 September
> 2020) MD5:d10deb7a45fa93934ec0a16c9d534e45\" ; :target~mip~ = \"OMIP\"
> ; :title = \"MRI JRA55-do 1.5.0 dataset prepared for input4MIPs\" ;
> :tracking~id~ = \"hdl:21.14100/9e2256e3-29b7-4e8a-906a-47a98a256308\"
> ; :variable~id~ = \"uas\" ; :license = \"OMIP boundary condition data
> produced by MRI is licensed under a Creative Commons
> Attribution-\[NonCommercial-\]ShareAlike 4.0 International License
> (<https://creativecommons.org/licenses>). Consult
> <https://pcmdi.llnl.gov/CMIP6/TermsOfUse> for terms of use governing
> input4MIPs output, including citation requirements and proper
> acknowledgment. Further information about this data, including some
> limitations, can be found via the further~infourl~ (recorded as a
> global attribute in this file). The data producers and data providers
> make no warranty, either express or implied, including, but not
> limited to, warranties of merchantability and fitness for a particular
> purpose. All liabilities arising from the supply of the information
> (including any liability arising in negligence) are excluded to the
> fullest extent permitted by law.\" ; :cmor~version~ = \"3.6.0\" ; }

### [vas](file:///Volumes/ioa01/reanalysis/JRA55do/vas)

Northward near-surface wind

> netcdf
> vas~input4MIPsatmosphericStateOMIPMRI~-JRA55-do-1-5-0~gr201001010000~-201012312100
> { dimensions: time = UNLIMITED ; // (2920 currently) lat = 320 ; lon =
> 640 ; bnds = 2 ; variables: double time(time) ; time:bounds =
> \"time~bnds~\" ; time:units = \"days since 1900-01-01 00:00:00\" ;
> time:calendar = \"gregorian\" ; time:axis = \"T\" ; time:long~name~ =
> \"time\" ; time:standard~name~ = \"time\" ; double time~bnds~(time,
> bnds) ; double lat(lat) ; lat:bounds = \"lat~bnds~\" ; lat:units =
> \"degrees~north~\" ; lat:axis = \"Y\" ; lat:long~name~ = \"Latitude\"
> ; lat:standard~name~ = \"latitude\" ; double lat~bnds~(lat, bnds) ;
> double lon(lon) ; lon:bounds = \"lon~bnds~\" ; lon:units =
> \"degrees~east~\" ; lon:axis = \"X\" ; lon:long~name~ = \"Longitude\"
> ; lon:standard~name~ = \"longitude\" ; double lon~bnds~(lon, bnds) ;
> double height ; height:units = \"m\" ; height:axis = \"Z\" ;
> height:positive = \"up\" ; height:long~name~ = \"height\" ;
> height:standard~name~ = \"height\" ; float vas(time, lat, lon) ;
> vas:standard~name~ = \"northward~wind~\" ; vas:long~name~ =
> \"Northward Near-Surface Wind\" ; vas:comment = \"Northward component
> of the near surface wind\" ; vas:units = \"m s-1\" ; vas:cell~methods~
> = \"area: mean time: point\" ; vas:cell~measures~ = \"area:
> areacella\" ; vas:history = \"2020-09-16T04:10:06Z altered by CMOR:
> Treated scalar dimension: \\\'height\\\'.\" ; vas:coordinates =
> \"height\" ; vas:missing~value~ = 1.e+20f ; vas:~FillValue~ = 1.e+20f
> ;
>
> // global attributes: :Conventions = \"CF-1.7 CMIP-6.2\" ;
> :activity~id~ = \"input4MIPs\" ; :cell~measures~ = \"area: areacella\"
> ; :comment = \"Based on JRA-55 reanalysis (1958-01 to 2020-07)\" ;
> :contact = \"Hiroyuki Tsujino (htsujino@mri-jma.go.jp)\" ;
> :creation~date~ = \"2020-09-16T04:10:11Z\" ; :data~specsversion~ =
> \"01.00.32\" ; :dataset~category~ = \"atmosphericState\" ;
> :external~variables~ = \"areacella\" ; :frequency = \"3hrPt\" ;
> :further~infourl~ =
> \"<http://climate.mri-jma.go.jp/~htsujino/jra55do.html>\" ; :grid =
> \"data regridded to the normal atmosphere TL319 gaussian grid (320x640
> latxlon) from a reduced TL319 gaussian grid\" ; :grid~label~ = \"gr\"
> ; :history = \"2020-09-16T04:10:11Z; CMOR rewrote data to be
> consistent with input4MIPs, CMIP6, and CF-1.7 CMIP-6.2 standards\" ;
> :institution = \"Meteorological Research Institute, Tsukuba, Ibaraki
> 305-0052, Japan\" ; :institution~id~ = \"MRI\" ; :mip~era~ = \"CMIP6\"
> ; :nominal~resolution~ = \"50 km\" ; :product = \"reanalysis\" ;
> :realm = \"atmos\" ; :references = \"Tsujino et al., 2018: JRA-55
> based surface dataset for driving ocean-sea-ice models (JRA55-do),
> Ocean Modelling, 130(1), pp 79-139.
> <https://doi.org/10.1016/j.ocemod.2018.07.002>\" ; :region =
> \"global~ocean~\" ; :release~year~ = \"2020\" ; :source = \"MRI
> JRA55-do 1.5.0: Atmospheric state generated for OMIP based on the
> JRA-55 reanalysis\" ; :source~description~ = \"Atmospheric state and
> terrestrial runoff datasets produced by MRI for the OMIP experiment of
> CMIP6\" ; :source~id~ = \"MRI-JRA55-do-1-5-0\" ; :source~type~ =
> \"satellite~blended~\" ; :source~version~ = \"1.5.0\" ; :table~id~ =
> \"input4MIPs~A3hrPt~\" ; :table~info~ = \"Creation Date:(14 September
> 2020) MD5:d10deb7a45fa93934ec0a16c9d534e45\" ; :target~mip~ = \"OMIP\"
> ; :title = \"MRI JRA55-do 1.5.0 dataset prepared for input4MIPs\" ;
> :tracking~id~ = \"hdl:21.14100/b80b6346-812b-4171-9be7-75718695dd0e\"
> ; :variable~id~ = \"vas\" ; :license = \"OMIP boundary condition data
> produced by MRI is licensed under a Creative Commons
> Attribution-\[NonCommercial-\]ShareAlike 4.0 International License
> (<https://creativecommons.org/licenses>). Consult
> <https://pcmdi.llnl.gov/CMIP6/TermsOfUse> for terms of use governing
> input4MIPs output, including citation requirements and proper
> acknowledgment. Further information about this data, including some
> limitations, can be found via the further~infourl~ (recorded as a
> global attribute in this file). The data producers and data providers
> make no warranty, either express or implied, including, but not
> limited to, warranties of merchantability and fitness for a particular
> purpose. All liabilities arising from the supply of the information
> (including any liability arising in negligence) are excluded to the
> fullest extent permitted by law.\" ; :cmor~version~ = \"3.6.0\" ; }

## MERRAv2

[MERRAv2](https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/) is a $1/2$
degree spatial and 1 hour temporal resolution atmospheric climate
dataset maintained by NASA. gadi has a repository described
[here](http://climate-cms.wikis.unsw.edu.au/MERRA2) and an [excel
spreadsheet](file:///Users/dpath2o/PHD/references/data/MERRA2_variables.xlsx)
has also been provided to further describe this repository.

### Surface Flux Diagnostics (\'flx\')

I have downloaded years 2010-2019 files from gadi to my [local
repository](file:///Volumes/ioa01/reanalysis/MERRAv2/flx). Here is the
header of the first file:

> netcdf MERRA2~300~.tavg1~2dflxNx~.20100101 { dimensions: lon = 576 ;
> lat = 361 ; time = UNLIMITED ; // (24 currently) variables: double
> lon(lon) ; lon:long~name~ = \"longitude\" ; lon:units =
> \"degrees~east~\" ; lon:vmax = 1.e+15f ; lon:vmin = -1.e+15f ;
> lon:valid~range~ = -1.e+15f, 1.e+15f ; double lat(lat) ;
> lat:long~name~ = \"latitude\" ; lat:units = \"degrees~north~\" ;
> lat:vmax = 1.e+15f ; lat:vmin = -1.e+15f ; lat:valid~range~ =
> -1.e+15f, 1.e+15f ; int time(time) ; time:long~name~ = \"time\" ;
> time:units = \"minutes since 2010-01-01 00:30:00\" ;
> time:time~increment~ = 10000 ; time:begin~date~ = 20100101 ;
> time:begin~time~ = 3000 ; time:vmax = 1.e+15f ; time:vmin = -1.e+15f ;
> time:valid~range~ = -1.e+15f, 1.e+15f ; float BSTAR(time, lat, lon) ;
> BSTAR:long~name~ = \"surface~bouyancyscale~\" ; BSTAR:units = \"m
> s-2\" ; BSTAR:~FillValue~ = 1.e+15f ; BSTAR:missing~value~ = 1.e+15f ;
> BSTAR:fmissing~value~ = 1.e+15f ; BSTAR:scale~factor~ = 1.f ;
> BSTAR:add~offset~ = 0.f ; BSTAR:standard~name~ =
> \"surface~bouyancyscale~\" ; BSTAR:vmax = 1.e+15f ; BSTAR:vmin =
> -1.e+15f ; BSTAR:valid~range~ = -1.e+15f, 1.e+15f ; float CDH(time,
> lat, lon) ; CDH:long~name~ = \"surface~exchangecoefficientforheat~\" ;
> CDH:units = \"kg m-2 s-1\" ; CDH:~FillValue~ = 1.e+15f ;
> CDH:missing~value~ = 1.e+15f ; CDH:fmissing~value~ = 1.e+15f ;
> CDH:scale~factor~ = 1.f ; CDH:add~offset~ = 0.f ; CDH:standard~name~ =
> \"surface~exchangecoefficientforheat~\" ; CDH:vmax = 1.e+15f ;
> CDH:vmin = -1.e+15f ; CDH:valid~range~ = -1.e+15f, 1.e+15f ; float
> CDM(time, lat, lon) ; CDM:long~name~ =
> \"surface~exchangecoefficientformomentum~\" ; CDM:units = \"kg m-2
> s-1\" ; CDM:~FillValue~ = 1.e+15f ; CDM:missing~value~ = 1.e+15f ;
> CDM:fmissing~value~ = 1.e+15f ; CDM:scale~factor~ = 1.f ;
> CDM:add~offset~ = 0.f ; CDM:standard~name~ =
> \"surface~exchangecoefficientformomentum~\" ; CDM:vmax = 1.e+15f ;
> CDM:vmin = -1.e+15f ; CDM:valid~range~ = -1.e+15f, 1.e+15f ; float
> CDQ(time, lat, lon) ; CDQ:long~name~ =
> \"surface~exchangecoefficientformoisture~\" ; CDQ:units = \"kg m-2
> s-1\" ; CDQ:~FillValue~ = 1.e+15f ; CDQ:missing~value~ = 1.e+15f ;
> CDQ:fmissing~value~ = 1.e+15f ; CDQ:scale~factor~ = 1.f ;
> CDQ:add~offset~ = 0.f ; CDQ:standard~name~ =
> \"surface~exchangecoefficientformoisture~\" ; CDQ:vmax = 1.e+15f ;
> CDQ:vmin = -1.e+15f ; CDQ:valid~range~ = -1.e+15f, 1.e+15f ; float
> CN(time, lat, lon) ; CN:long~name~ =
> \"surface~neutraldragcoefficient~\" ; CN:units = \"1\" ;
> CN:~FillValue~ = 1.e+15f ; CN:missing~value~ = 1.e+15f ;
> CN:fmissing~value~ = 1.e+15f ; CN:scale~factor~ = 1.f ; CN:add~offset~
> = 0.f ; CN:standard~name~ = \"surface~neutraldragcoefficient~\" ;
> CN:vmax = 1.e+15f ; CN:vmin = -1.e+15f ; CN:valid~range~ = -1.e+15f,
> 1.e+15f ; float DISPH(time, lat, lon) ; DISPH:long~name~ =
> \"zero~planedisplacementheight~\" ; DISPH:units = \"m\" ;
> DISPH:~FillValue~ = 1.e+15f ; DISPH:missing~value~ = 1.e+15f ;
> DISPH:fmissing~value~ = 1.e+15f ; DISPH:scale~factor~ = 1.f ;
> DISPH:add~offset~ = 0.f ; DISPH:standard~name~ =
> \"zero~planedisplacementheight~\" ; DISPH:vmax = 1.e+15f ; DISPH:vmin
> = -1.e+15f ; DISPH:valid~range~ = -1.e+15f, 1.e+15f ; float
> EFLUX(time, lat, lon) ; EFLUX:long~name~ = \"total~latentenergyflux~\"
> ; EFLUX:units = \"W m-2\" ; EFLUX:~FillValue~ = 1.e+15f ;
> EFLUX:missing~value~ = 1.e+15f ; EFLUX:fmissing~value~ = 1.e+15f ;
> EFLUX:scale~factor~ = 1.f ; EFLUX:add~offset~ = 0.f ;
> EFLUX:standard~name~ = \"total~latentenergyflux~\" ; EFLUX:vmax =
> 1.e+15f ; EFLUX:vmin = -1.e+15f ; EFLUX:valid~range~ = -1.e+15f,
> 1.e+15f ; float EVAP(time, lat, lon) ; EVAP:long~name~ =
> \"evaporation~fromturbulence~\" ; EVAP:units = \"kg m-2 s-1\" ;
> EVAP:~FillValue~ = 1.e+15f ; EVAP:missing~value~ = 1.e+15f ;
> EVAP:fmissing~value~ = 1.e+15f ; EVAP:scale~factor~ = 1.f ;
> EVAP:add~offset~ = 0.f ; EVAP:standard~name~ =
> \"evaporation~fromturbulence~\" ; EVAP:vmax = 1.e+15f ; EVAP:vmin =
> -1.e+15f ; EVAP:valid~range~ = -1.e+15f, 1.e+15f ; float FRCAN(time,
> lat, lon) ; FRCAN:long~name~ = \"areal~fractionofanvilshowers~\" ;
> FRCAN:units = \"1\" ; FRCAN:~FillValue~ = 1.e+15f ;
> FRCAN:missing~value~ = 1.e+15f ; FRCAN:fmissing~value~ = 1.e+15f ;
> FRCAN:scale~factor~ = 1.f ; FRCAN:add~offset~ = 0.f ;
> FRCAN:standard~name~ = \"areal~fractionofanvilshowers~\" ; FRCAN:vmax
> = 1.e+15f ; FRCAN:vmin = -1.e+15f ; FRCAN:valid~range~ = -1.e+15f,
> 1.e+15f ; float FRCCN(time, lat, lon) ; FRCCN:long~name~ =
> \"areal~fractionofconvectiveshowers~\" ; FRCCN:units = \"1\" ;
> FRCCN:~FillValue~ = 1.e+15f ; FRCCN:missing~value~ = 1.e+15f ;
> FRCCN:fmissing~value~ = 1.e+15f ; FRCCN:scale~factor~ = 1.f ;
> FRCCN:add~offset~ = 0.f ; FRCCN:standard~name~ =
> \"areal~fractionofconvectiveshowers~\" ; FRCCN:vmax = 1.e+15f ;
> FRCCN:vmin = -1.e+15f ; FRCCN:valid~range~ = -1.e+15f, 1.e+15f ; float
> FRCLS(time, lat, lon) ; FRCLS:long~name~ =
> \"areal~fractionofnonanvillargescaleshowers~\" ; FRCLS:units = \"1\" ;
> FRCLS:~FillValue~ = 1.e+15f ; FRCLS:missing~value~ = 1.e+15f ;
> FRCLS:fmissing~value~ = 1.e+15f ; FRCLS:scale~factor~ = 1.f ;
> FRCLS:add~offset~ = 0.f ; FRCLS:standard~name~ =
> \"areal~fractionofnonanvillargescaleshowers~\" ; FRCLS:vmax = 1.e+15f
> ; FRCLS:vmin = -1.e+15f ; FRCLS:valid~range~ = -1.e+15f, 1.e+15f ;
> float FRSEAICE(time, lat, lon) ; FRSEAICE:long~name~ =
> \"ice~coveredfractionoftile~\" ; FRSEAICE:units = \"1\" ;
> FRSEAICE:~FillValue~ = 1.e+15f ; FRSEAICE:missing~value~ = 1.e+15f ;
> FRSEAICE:fmissing~value~ = 1.e+15f ; FRSEAICE:scale~factor~ = 1.f ;
> FRSEAICE:add~offset~ = 0.f ; FRSEAICE:standard~name~ =
> \"ice~coveredfractionoftile~\" ; FRSEAICE:vmax = 1.e+15f ;
> FRSEAICE:vmin = -1.e+15f ; FRSEAICE:valid~range~ = -1.e+15f, 1.e+15f ;
> float GHTSKIN(time, lat, lon) ; GHTSKIN:long~name~ =
> \"Ground~heatingforskintemp~\" ; GHTSKIN:units = \"W m-2\" ;
> GHTSKIN:~FillValue~ = 1.e+15f ; GHTSKIN:missing~value~ = 1.e+15f ;
> GHTSKIN:fmissing~value~ = 1.e+15f ; GHTSKIN:scale~factor~ = 1.f ;
> GHTSKIN:add~offset~ = 0.f ; GHTSKIN:standard~name~ =
> \"Ground~heatingforskintemp~\" ; GHTSKIN:vmax = 1.e+15f ; GHTSKIN:vmin
> = -1.e+15f ; GHTSKIN:valid~range~ = -1.e+15f, 1.e+15f ; float
> HFLUX(time, lat, lon) ; HFLUX:long~name~ =
> \"sensible~heatfluxfromturbulence~\" ; HFLUX:units = \"W m-2\" ;
> HFLUX:~FillValue~ = 1.e+15f ; HFLUX:missing~value~ = 1.e+15f ;
> HFLUX:fmissing~value~ = 1.e+15f ; HFLUX:scale~factor~ = 1.f ;
> HFLUX:add~offset~ = 0.f ; HFLUX:standard~name~ =
> \"sensible~heatfluxfromturbulence~\" ; HFLUX:vmax = 1.e+15f ;
> HFLUX:vmin = -1.e+15f ; HFLUX:valid~range~ = -1.e+15f, 1.e+15f ; float
> HLML(time, lat, lon) ; HLML:long~name~ = \"surface~layerheight~\" ;
> HLML:units = \"m\" ; HLML:~FillValue~ = 1.e+15f ; HLML:missing~value~
> = 1.e+15f ; HLML:fmissing~value~ = 1.e+15f ; HLML:scale~factor~ = 1.f
> ; HLML:add~offset~ = 0.f ; HLML:standard~name~ =
> \"surface~layerheight~\" ; HLML:vmax = 1.e+15f ; HLML:vmin = -1.e+15f
> ; HLML:valid~range~ = -1.e+15f, 1.e+15f ; float NIRDF(time, lat, lon)
> ; NIRDF:long~name~ = \"surface~downwellingnearinfrareddiffuseflux~\" ;
> NIRDF:units = \"W m-2\" ; NIRDF:~FillValue~ = 1.e+15f ;
> NIRDF:missing~value~ = 1.e+15f ; NIRDF:fmissing~value~ = 1.e+15f ;
> NIRDF:scale~factor~ = 1.f ; NIRDF:add~offset~ = 0.f ;
> NIRDF:standard~name~ = \"surface~downwellingnearinfrareddiffuseflux~\"
> ; NIRDF:vmax = 1.e+15f ; NIRDF:vmin = -1.e+15f ; NIRDF:valid~range~ =
> -1.e+15f, 1.e+15f ; float NIRDR(time, lat, lon) ; NIRDR:long~name~ =
> \"surface~downwellingnearinfraredbeamflux~\" ; NIRDR:units = \"W m-2\"
> ; NIRDR:~FillValue~ = 1.e+15f ; NIRDR:missing~value~ = 1.e+15f ;
> NIRDR:fmissing~value~ = 1.e+15f ; NIRDR:scale~factor~ = 1.f ;
> NIRDR:add~offset~ = 0.f ; NIRDR:standard~name~ =
> \"surface~downwellingnearinfraredbeamflux~\" ; NIRDR:vmax = 1.e+15f ;
> NIRDR:vmin = -1.e+15f ; NIRDR:valid~range~ = -1.e+15f, 1.e+15f ; float
> PBLH(time, lat, lon) ; PBLH:long~name~ =
> \"planetary~boundarylayerheight~\" ; PBLH:units = \"m\" ;
> PBLH:~FillValue~ = 1.e+15f ; PBLH:missing~value~ = 1.e+15f ;
> PBLH:fmissing~value~ = 1.e+15f ; PBLH:scale~factor~ = 1.f ;
> PBLH:add~offset~ = 0.f ; PBLH:standard~name~ =
> \"planetary~boundarylayerheight~\" ; PBLH:vmax = 1.e+15f ; PBLH:vmin =
> -1.e+15f ; PBLH:valid~range~ = -1.e+15f, 1.e+15f ; float PGENTOT(time,
> lat, lon) ; PGENTOT:long~name~ =
> \"Total~columnproductionofprecipitation~\" ; PGENTOT:units = \"kg m-2
> s-1\" ; PGENTOT:~FillValue~ = 1.e+15f ; PGENTOT:missing~value~ =
> 1.e+15f ; PGENTOT:fmissing~value~ = 1.e+15f ; PGENTOT:scale~factor~ =
> 1.f ; PGENTOT:add~offset~ = 0.f ; PGENTOT:standard~name~ =
> \"Total~columnproductionofprecipitation~\" ; PGENTOT:vmax = 1.e+15f ;
> PGENTOT:vmin = -1.e+15f ; PGENTOT:valid~range~ = -1.e+15f, 1.e+15f ;
> float PRECANV(time, lat, lon) ; PRECANV:long~name~ =
> \"anvil~precipitation~\" ; PRECANV:units = \"kg m-2 s-1\" ;
> PRECANV:~FillValue~ = 1.e+15f ; PRECANV:missing~value~ = 1.e+15f ;
> PRECANV:fmissing~value~ = 1.e+15f ; PRECANV:scale~factor~ = 1.f ;
> PRECANV:add~offset~ = 0.f ; PRECANV:standard~name~ =
> \"anvil~precipitation~\" ; PRECANV:vmax = 1.e+15f ; PRECANV:vmin =
> -1.e+15f ; PRECANV:valid~range~ = -1.e+15f, 1.e+15f ; float
> PRECCON(time, lat, lon) ; PRECCON:long~name~ =
> \"convective~precipitation~\" ; PRECCON:units = \"kg m-2 s-1\" ;
> PRECCON:~FillValue~ = 1.e+15f ; PRECCON:missing~value~ = 1.e+15f ;
> PRECCON:fmissing~value~ = 1.e+15f ; PRECCON:scale~factor~ = 1.f ;
> PRECCON:add~offset~ = 0.f ; PRECCON:standard~name~ =
> \"convective~precipitation~\" ; PRECCON:vmax = 1.e+15f ; PRECCON:vmin
> = -1.e+15f ; PRECCON:valid~range~ = -1.e+15f, 1.e+15f ; float
> PRECLSC(time, lat, lon) ; PRECLSC:long~name~ =
> \"nonanvil~largescaleprecipitation~\" ; PRECLSC:units = \"kg m-2 s-1\"
> ; PRECLSC:~FillValue~ = 1.e+15f ; PRECLSC:missing~value~ = 1.e+15f ;
> PRECLSC:fmissing~value~ = 1.e+15f ; PRECLSC:scale~factor~ = 1.f ;
> PRECLSC:add~offset~ = 0.f ; PRECLSC:standard~name~ =
> \"nonanvil~largescaleprecipitation~\" ; PRECLSC:vmax = 1.e+15f ;
> PRECLSC:vmin = -1.e+15f ; PRECLSC:valid~range~ = -1.e+15f, 1.e+15f ;
> float PRECSNO(time, lat, lon) ; PRECSNO:long~name~ = \"snowfall\" ;
> PRECSNO:units = \"kg m-2 s-1\" ; PRECSNO:~FillValue~ = 1.e+15f ;
> PRECSNO:missing~value~ = 1.e+15f ; PRECSNO:fmissing~value~ = 1.e+15f ;
> PRECSNO:scale~factor~ = 1.f ; PRECSNO:add~offset~ = 0.f ;
> PRECSNO:standard~name~ = \"snowfall\" ; PRECSNO:vmax = 1.e+15f ;
> PRECSNO:vmin = -1.e+15f ; PRECSNO:valid~range~ = -1.e+15f, 1.e+15f ;
> float PRECTOT(time, lat, lon) ; PRECTOT:long~name~ =
> \"total~precipitation~\" ; PRECTOT:units = \"kg m-2 s-1\" ;
> PRECTOT:~FillValue~ = 1.e+15f ; PRECTOT:missing~value~ = 1.e+15f ;
> PRECTOT:fmissing~value~ = 1.e+15f ; PRECTOT:scale~factor~ = 1.f ;
> PRECTOT:add~offset~ = 0.f ; PRECTOT:standard~name~ =
> \"total~precipitation~\" ; PRECTOT:vmax = 1.e+15f ; PRECTOT:vmin =
> -1.e+15f ; PRECTOT:valid~range~ = -1.e+15f, 1.e+15f ; float
> PRECTOTCORR(time, lat, lon) ; PRECTOTCORR:long~name~ =
> \"total~precipitation~\" ; PRECTOTCORR:units = \"kg m-2 s-1\" ;
> PRECTOTCORR:~FillValue~ = 1.e+15f ; PRECTOTCORR:missing~value~ =
> 1.e+15f ; PRECTOTCORR:fmissing~value~ = 1.e+15f ;
> PRECTOTCORR:scale~factor~ = 1.f ; PRECTOTCORR:add~offset~ = 0.f ;
> PRECTOTCORR:standard~name~ = \"total~precipitation~\" ;
> PRECTOTCORR:vmax = 1.e+15f ; PRECTOTCORR:vmin = -1.e+15f ;
> PRECTOTCORR:valid~range~ = -1.e+15f, 1.e+15f ; float PREVTOT(time,
> lat, lon) ; PREVTOT:long~name~ =
> \"Total~columnre~-evap/subl~ofprecipitation~\" ; PREVTOT:units = \"kg
> m-2 s-1\" ; PREVTOT:~FillValue~ = 1.e+15f ; PREVTOT:missing~value~ =
> 1.e+15f ; PREVTOT:fmissing~value~ = 1.e+15f ; PREVTOT:scale~factor~ =
> 1.f ; PREVTOT:add~offset~ = 0.f ; PREVTOT:standard~name~ =
> \"Total~columnre~-evap/subl~ofprecipitation~\" ; PREVTOT:vmax =
> 1.e+15f ; PREVTOT:vmin = -1.e+15f ; PREVTOT:valid~range~ = -1.e+15f,
> 1.e+15f ; float QLML(time, lat, lon) ; QLML:long~name~ =
> \"surface~specifichumidity~\" ; QLML:units = \"1\" ; QLML:~FillValue~
> = 1.e+15f ; QLML:missing~value~ = 1.e+15f ; QLML:fmissing~value~ =
> 1.e+15f ; QLML:scale~factor~ = 1.f ; QLML:add~offset~ = 0.f ;
> QLML:standard~name~ = \"surface~specifichumidity~\" ; QLML:vmax =
> 1.e+15f ; QLML:vmin = -1.e+15f ; QLML:valid~range~ = -1.e+15f, 1.e+15f
> ; float QSH(time, lat, lon) ; QSH:long~name~ =
> \"effective~surfacespecifichumidity~\" ; QSH:units = \"kg kg-1\" ;
> QSH:~FillValue~ = 1.e+15f ; QSH:missing~value~ = 1.e+15f ;
> QSH:fmissing~value~ = 1.e+15f ; QSH:scale~factor~ = 1.f ;
> QSH:add~offset~ = 0.f ; QSH:standard~name~ =
> \"effective~surfacespecifichumidity~\" ; QSH:vmax = 1.e+15f ; QSH:vmin
> = -1.e+15f ; QSH:valid~range~ = -1.e+15f, 1.e+15f ; float QSTAR(time,
> lat, lon) ; QSTAR:long~name~ = \"surface~moisturescale~\" ;
> QSTAR:units = \"kg kg-1\" ; QSTAR:~FillValue~ = 1.e+15f ;
> QSTAR:missing~value~ = 1.e+15f ; QSTAR:fmissing~value~ = 1.e+15f ;
> QSTAR:scale~factor~ = 1.f ; QSTAR:add~offset~ = 0.f ;
> QSTAR:standard~name~ = \"surface~moisturescale~\" ; QSTAR:vmax =
> 1.e+15f ; QSTAR:vmin = -1.e+15f ; QSTAR:valid~range~ = -1.e+15f,
> 1.e+15f ; float RHOA(time, lat, lon) ; RHOA:long~name~ =
> \"air~densityatsurface~\" ; RHOA:units = \"kg m-3\" ; RHOA:~FillValue~
> = 1.e+15f ; RHOA:missing~value~ = 1.e+15f ; RHOA:fmissing~value~ =
> 1.e+15f ; RHOA:scale~factor~ = 1.f ; RHOA:add~offset~ = 0.f ;
> RHOA:standard~name~ = \"air~densityatsurface~\" ; RHOA:vmax = 1.e+15f
> ; RHOA:vmin = -1.e+15f ; RHOA:valid~range~ = -1.e+15f, 1.e+15f ; float
> RISFC(time, lat, lon) ; RISFC:long~name~ =
> \"surface~bulkrichardsonnumber~\" ; RISFC:units = \"1\" ;
> RISFC:~FillValue~ = 1.e+15f ; RISFC:missing~value~ = 1.e+15f ;
> RISFC:fmissing~value~ = 1.e+15f ; RISFC:scale~factor~ = 1.f ;
> RISFC:add~offset~ = 0.f ; RISFC:standard~name~ =
> \"surface~bulkrichardsonnumber~\" ; RISFC:vmax = 1.e+15f ; RISFC:vmin
> = -1.e+15f ; RISFC:valid~range~ = -1.e+15f, 1.e+15f ; float
> SPEED(time, lat, lon) ; SPEED:long~name~ = \"surface~windspeed~\" ;
> SPEED:units = \"m s-1\" ; SPEED:~FillValue~ = 1.e+15f ;
> SPEED:missing~value~ = 1.e+15f ; SPEED:fmissing~value~ = 1.e+15f ;
> SPEED:scale~factor~ = 1.f ; SPEED:add~offset~ = 0.f ;
> SPEED:standard~name~ = \"surface~windspeed~\" ; SPEED:vmax = 1.e+15f ;
> SPEED:vmin = -1.e+15f ; SPEED:valid~range~ = -1.e+15f, 1.e+15f ; float
> SPEEDMAX(time, lat, lon) ; SPEEDMAX:long~name~ =
> \"surface~windspeed~\" ; SPEEDMAX:units = \"m s-1\" ;
> SPEEDMAX:~FillValue~ = 1.e+15f ; SPEEDMAX:missing~value~ = 1.e+15f ;
> SPEEDMAX:fmissing~value~ = 1.e+15f ; SPEEDMAX:scale~factor~ = 1.f ;
> SPEEDMAX:add~offset~ = 0.f ; SPEEDMAX:standard~name~ =
> \"surface~windspeed~\" ; SPEEDMAX:vmax = 1.e+15f ; SPEEDMAX:vmin =
> -1.e+15f ; SPEEDMAX:valid~range~ = -1.e+15f, 1.e+15f ; float
> TAUGWX(time, lat, lon) ; TAUGWX:long~name~ =
> \"surface~eastwardgravitywavestress~\" ; TAUGWX:units = \"N m-2\" ;
> TAUGWX:~FillValue~ = 1.e+15f ; TAUGWX:missing~value~ = 1.e+15f ;
> TAUGWX:fmissing~value~ = 1.e+15f ; TAUGWX:scale~factor~ = 1.f ;
> TAUGWX:add~offset~ = 0.f ; TAUGWX:standard~name~ =
> \"surface~eastwardgravitywavestress~\" ; TAUGWX:vmax = 1.e+15f ;
> TAUGWX:vmin = -1.e+15f ; TAUGWX:valid~range~ = -1.e+15f, 1.e+15f ;
> float TAUGWY(time, lat, lon) ; TAUGWY:long~name~ =
> \"surface~northwardgravitywavestress~\" ; TAUGWY:units = \"N m-2\" ;
> TAUGWY:~FillValue~ = 1.e+15f ; TAUGWY:missing~value~ = 1.e+15f ;
> TAUGWY:fmissing~value~ = 1.e+15f ; TAUGWY:scale~factor~ = 1.f ;
> TAUGWY:add~offset~ = 0.f ; TAUGWY:standard~name~ =
> \"surface~northwardgravitywavestress~\" ; TAUGWY:vmax = 1.e+15f ;
> TAUGWY:vmin = -1.e+15f ; TAUGWY:valid~range~ = -1.e+15f, 1.e+15f ;
> float TAUX(time, lat, lon) ; TAUX:long~name~ =
> \"eastward~surfacestress~\" ; TAUX:units = \"N m-2\" ;
> TAUX:~FillValue~ = 1.e+15f ; TAUX:missing~value~ = 1.e+15f ;
> TAUX:fmissing~value~ = 1.e+15f ; TAUX:scale~factor~ = 1.f ;
> TAUX:add~offset~ = 0.f ; TAUX:standard~name~ =
> \"eastward~surfacestress~\" ; TAUX:vmax = 1.e+15f ; TAUX:vmin =
> -1.e+15f ; TAUX:valid~range~ = -1.e+15f, 1.e+15f ; float TAUY(time,
> lat, lon) ; TAUY:long~name~ = \"northward~surfacestress~\" ;
> TAUY:units = \"N m-2\" ; TAUY:~FillValue~ = 1.e+15f ;
> TAUY:missing~value~ = 1.e+15f ; TAUY:fmissing~value~ = 1.e+15f ;
> TAUY:scale~factor~ = 1.f ; TAUY:add~offset~ = 0.f ;
> TAUY:standard~name~ = \"northward~surfacestress~\" ; TAUY:vmax =
> 1.e+15f ; TAUY:vmin = -1.e+15f ; TAUY:valid~range~ = -1.e+15f, 1.e+15f
> ; float TCZPBL(time, lat, lon) ; TCZPBL:long~name~ =
> \"transcom~planetaryboundarylayerheight~\" ; TCZPBL:units = \"m\" ;
> TCZPBL:~FillValue~ = 1.e+15f ; TCZPBL:missing~value~ = 1.e+15f ;
> TCZPBL:fmissing~value~ = 1.e+15f ; TCZPBL:scale~factor~ = 1.f ;
> TCZPBL:add~offset~ = 0.f ; TCZPBL:standard~name~ =
> \"transcom~planetaryboundarylayerheight~\" ; TCZPBL:vmax = 1.e+15f ;
> TCZPBL:vmin = -1.e+15f ; TCZPBL:valid~range~ = -1.e+15f, 1.e+15f ;
> float TLML(time, lat, lon) ; TLML:long~name~ =
> \"surface~airtemperature~\" ; TLML:units = \"K\" ; TLML:~FillValue~ =
> 1.e+15f ; TLML:missing~value~ = 1.e+15f ; TLML:fmissing~value~ =
> 1.e+15f ; TLML:scale~factor~ = 1.f ; TLML:add~offset~ = 0.f ;
> TLML:standard~name~ = \"surface~airtemperature~\" ; TLML:vmax =
> 1.e+15f ; TLML:vmin = -1.e+15f ; TLML:valid~range~ = -1.e+15f, 1.e+15f
> ; float TSH(time, lat, lon) ; TSH:long~name~ =
> \"effective~surfaceskintemperature~\" ; TSH:units = \"K\" ;
> TSH:~FillValue~ = 1.e+15f ; TSH:missing~value~ = 1.e+15f ;
> TSH:fmissing~value~ = 1.e+15f ; TSH:scale~factor~ = 1.f ;
> TSH:add~offset~ = 0.f ; TSH:standard~name~ =
> \"effective~surfaceskintemperature~\" ; TSH:vmax = 1.e+15f ; TSH:vmin
> = -1.e+15f ; TSH:valid~range~ = -1.e+15f, 1.e+15f ; float TSTAR(time,
> lat, lon) ; TSTAR:long~name~ = \"surface~temperaturescale~\" ;
> TSTAR:units = \"K\" ; TSTAR:~FillValue~ = 1.e+15f ;
> TSTAR:missing~value~ = 1.e+15f ; TSTAR:fmissing~value~ = 1.e+15f ;
> TSTAR:scale~factor~ = 1.f ; TSTAR:add~offset~ = 0.f ;
> TSTAR:standard~name~ = \"surface~temperaturescale~\" ; TSTAR:vmax =
> 1.e+15f ; TSTAR:vmin = -1.e+15f ; TSTAR:valid~range~ = -1.e+15f,
> 1.e+15f ; float ULML(time, lat, lon) ; ULML:long~name~ =
> \"surface~eastwardwind~\" ; ULML:units = \"m s-1\" ; ULML:~FillValue~
> = 1.e+15f ; ULML:missing~value~ = 1.e+15f ; ULML:fmissing~value~ =
> 1.e+15f ; ULML:scale~factor~ = 1.f ; ULML:add~offset~ = 0.f ;
> ULML:standard~name~ = \"surface~eastwardwind~\" ; ULML:vmax = 1.e+15f
> ; ULML:vmin = -1.e+15f ; ULML:valid~range~ = -1.e+15f, 1.e+15f ; float
> USTAR(time, lat, lon) ; USTAR:long~name~ = \"surface~velocityscale~\"
> ; USTAR:units = \"m s-1\" ; USTAR:~FillValue~ = 1.e+15f ;
> USTAR:missing~value~ = 1.e+15f ; USTAR:fmissing~value~ = 1.e+15f ;
> USTAR:scale~factor~ = 1.f ; USTAR:add~offset~ = 0.f ;
> USTAR:standard~name~ = \"surface~velocityscale~\" ; USTAR:vmax =
> 1.e+15f ; USTAR:vmin = -1.e+15f ; USTAR:valid~range~ = -1.e+15f,
> 1.e+15f ; float VLML(time, lat, lon) ; VLML:long~name~ =
> \"surface~northwardwind~\" ; VLML:units = \"m s-1\" ; VLML:~FillValue~
> = 1.e+15f ; VLML:missing~value~ = 1.e+15f ; VLML:fmissing~value~ =
> 1.e+15f ; VLML:scale~factor~ = 1.f ; VLML:add~offset~ = 0.f ;
> VLML:standard~name~ = \"surface~northwardwind~\" ; VLML:vmax = 1.e+15f
> ; VLML:vmin = -1.e+15f ; VLML:valid~range~ = -1.e+15f, 1.e+15f ; float
> Z0H(time, lat, lon) ; Z0H:long~name~ = \"surface~roughnessforheat~\" ;
> Z0H:units = \"m\" ; Z0H:~FillValue~ = 1.e+15f ; Z0H:missing~value~ =
> 1.e+15f ; Z0H:fmissing~value~ = 1.e+15f ; Z0H:scale~factor~ = 1.f ;
> Z0H:add~offset~ = 0.f ; Z0H:standard~name~ =
> \"surface~roughnessforheat~\" ; Z0H:vmax = 1.e+15f ; Z0H:vmin =
> -1.e+15f ; Z0H:valid~range~ = -1.e+15f, 1.e+15f ; float Z0M(time, lat,
> lon) ; Z0M:long~name~ = \"surface~roughness~\" ; Z0M:units = \"m\" ;
> Z0M:~FillValue~ = 1.e+15f ; Z0M:missing~value~ = 1.e+15f ;
> Z0M:fmissing~value~ = 1.e+15f ; Z0M:scale~factor~ = 1.f ;
> Z0M:add~offset~ = 0.f ; Z0M:standard~name~ = \"surface~roughness~\" ;
> Z0M:vmax = 1.e+15f ; Z0M:vmin = -1.e+15f ; Z0M:valid~range~ =
> -1.e+15f, 1.e+15f ;
>
> // global attributes: :History = \"Original file generated: Mon Mar 23
> 03:29:05 2015 GMT\" ; :Comment = \"GMAO filename:
> d5124~m2jan00~.tavg1~2dflxNx~.20100101.nc4\" ; :Filename =
> \"MERRA2~300~.tavg1~2dflxNx~.20100101.nc4\" ; :Conventions = \"CF-1\"
> ; :Institution = \"NASA Global Modeling and Assimilation Office\" ;
> :References = \"<http://gmao.gsfc.nasa.gov>\" ; :Format =
> \"NetCDF-4/HDF-5\" ; :SpatialCoverage = \"global\" ; :VersionID =
> \"5.12.4\" ; :TemporalRange = \"1980-01-01 -\> 2016-12-31\" ;
> :identifier~productdoiauthority~ = \"<http://dx.doi.org/>\" ;
> :ShortName = \"M2T1NXFLX\" ; :GranuleID =
> \"MERRA2~300~.tavg1~2dflxNx~.20100101.nc4\" ; :ProductionDateTime =
> \"Original file generated: Mon Mar 23 03:29:05 2015 GMT\" ; :LongName
> = \"MERRA2 tavg1~2dflxNx~:
> 2d,1-Hourly,Time-Averaged,Single-Level,Assimilation,Surface Flux
> Diagnostics\" ; :Title = \"MERRA2 tavg1~2dflxNx~:
> 2d,1-Hourly,Time-Averaged,Single-Level,Assimilation,Surface Flux
> Diagnostics\" ; :SouthernmostLatitude = \"-90.0\" ;
> :NorthernmostLatitude = \"90.0\" ; :WesternmostLongitude = \"-180.0\"
> ; :EasternmostLongitude = \"179.375\" ; :LatitudeResolution = \"0.5\"
> ; :LongitudeResolution = \"0.625\" ; :DataResolution = \"0.5 x 0.625\"
> ; :Source = \"CVS tag: GEOSadas-5~124~\" ; :Contact =
> \"<http://gmao.gsfc.nasa.gov>\" ; :identifier~productdoi~ =
> \"10.5067/7MCPBJ41Y0K6\" ; :RangeBeginningDate = \"2010-01-01\" ;
> :RangeBeginningTime = \"00:00:00.000000\" ; :RangeEndingDate =
> \"2010-01-01\" ; :RangeEndingTime = \"23:59:59.000000\" ; }

### Ocean Surface Diagnostics (\'ocn\')

Not presently on any projects associated with gadi. I have inquired with
the one of the project leaders at UTas if there is any scope to \'host\'
(for the lack of a better term) this portion of MERRAv2 in gadi project
\'ua8\'

### Radiation Diagnostics (\'rad\')

I have downloaded years 2010-2019 files from gadi to my [local
repository](file:///volumes/ioa01/reanalysis/MERRAv2/rad). Here is the
header of the first file:

> netcdf MERRA2~300~.tavg1~2dradNx~.20100101 { dimensions: lon = 576 ;
> lat = 361 ; time = UNLIMITED ; // (24 currently) variables: double
> lon(lon) ; lon:long~name~ = \"longitude\" ; lon:units =
> \"degrees~east~\" ; lon:vmax = 1.e+15f ; lon:vmin = -1.e+15f ;
> lon:valid~range~ = -1.e+15f, 1.e+15f ; double lat(lat) ;
> lat:long~name~ = \"latitude\" ; lat:units = \"degrees~north~\" ;
> lat:vmax = 1.e+15f ; lat:vmin = -1.e+15f ; lat:valid~range~ =
> -1.e+15f, 1.e+15f ; int time(time) ; time:long~name~ = \"time\" ;
> time:units = \"minutes since 2010-01-01 00:30:00\" ;
> time:time~increment~ = 10000 ; time:begin~date~ = 20100101 ;
> time:begin~time~ = 3000 ; time:vmax = 1.e+15f ; time:vmin = -1.e+15f ;
> time:valid~range~ = -1.e+15f, 1.e+15f ; float ALBEDO(time, lat, lon) ;
> ALBEDO:long~name~ = \"surface~albedo~\" ; ALBEDO:units = \"1\" ;
> ALBEDO:~FillValue~ = 1.e+15f ; ALBEDO:missing~value~ = 1.e+15f ;
> ALBEDO:fmissing~value~ = 1.e+15f ; ALBEDO:scale~factor~ = 1.f ;
> ALBEDO:add~offset~ = 0.f ; ALBEDO:standard~name~ = \"surface~albedo~\"
> ; ALBEDO:vmax = 1.e+15f ; ALBEDO:vmin = -1.e+15f ; ALBEDO:valid~range~
> = -1.e+15f, 1.e+15f ; float ALBNIRDF(time, lat, lon) ;
> ALBNIRDF:long~name~ = \"surface~albedofornearinfrareddiffuse~\" ;
> ALBNIRDF:units = \"1\" ; ALBNIRDF:~FillValue~ = 1.e+15f ;
> ALBNIRDF:missing~value~ = 1.e+15f ; ALBNIRDF:fmissing~value~ = 1.e+15f
> ; ALBNIRDF:scale~factor~ = 1.f ; ALBNIRDF:add~offset~ = 0.f ;
> ALBNIRDF:standard~name~ = \"surface~albedofornearinfrareddiffuse~\" ;
> ALBNIRDF:vmax = 1.e+15f ; ALBNIRDF:vmin = -1.e+15f ;
> ALBNIRDF:valid~range~ = -1.e+15f, 1.e+15f ; float ALBNIRDR(time, lat,
> lon) ; ALBNIRDR:long~name~ = \"surface~albedofornearinfraredbeam~\" ;
> ALBNIRDR:units = \"1\" ; ALBNIRDR:~FillValue~ = 1.e+15f ;
> ALBNIRDR:missing~value~ = 1.e+15f ; ALBNIRDR:fmissing~value~ = 1.e+15f
> ; ALBNIRDR:scale~factor~ = 1.f ; ALBNIRDR:add~offset~ = 0.f ;
> ALBNIRDR:standard~name~ = \"surface~albedofornearinfraredbeam~\" ;
> ALBNIRDR:vmax = 1.e+15f ; ALBNIRDR:vmin = -1.e+15f ;
> ALBNIRDR:valid~range~ = -1.e+15f, 1.e+15f ; float ALBVISDF(time, lat,
> lon) ; ALBVISDF:long~name~ = \"surface~albedoforvisiblediffuse~\" ;
> ALBVISDF:units = \"1\" ; ALBVISDF:~FillValue~ = 1.e+15f ;
> ALBVISDF:missing~value~ = 1.e+15f ; ALBVISDF:fmissing~value~ = 1.e+15f
> ; ALBVISDF:scale~factor~ = 1.f ; ALBVISDF:add~offset~ = 0.f ;
> ALBVISDF:standard~name~ = \"surface~albedoforvisiblediffuse~\" ;
> ALBVISDF:vmax = 1.e+15f ; ALBVISDF:vmin = -1.e+15f ;
> ALBVISDF:valid~range~ = -1.e+15f, 1.e+15f ; float ALBVISDR(time, lat,
> lon) ; ALBVISDR:long~name~ = \"surface~albedoforvisiblebeam~\" ;
> ALBVISDR:units = \"1\" ; ALBVISDR:~FillValue~ = 1.e+15f ;
> ALBVISDR:missing~value~ = 1.e+15f ; ALBVISDR:fmissing~value~ = 1.e+15f
> ; ALBVISDR:scale~factor~ = 1.f ; ALBVISDR:add~offset~ = 0.f ;
> ALBVISDR:standard~name~ = \"surface~albedoforvisiblebeam~\" ;
> ALBVISDR:vmax = 1.e+15f ; ALBVISDR:vmin = -1.e+15f ;
> ALBVISDR:valid~range~ = -1.e+15f, 1.e+15f ; float CLDHGH(time, lat,
> lon) ; CLDHGH:long~name~ = \"cloud~areafractionforhighclouds~\" ;
> CLDHGH:units = \"1\" ; CLDHGH:~FillValue~ = 1.e+15f ;
> CLDHGH:missing~value~ = 1.e+15f ; CLDHGH:fmissing~value~ = 1.e+15f ;
> CLDHGH:scale~factor~ = 1.f ; CLDHGH:add~offset~ = 0.f ;
> CLDHGH:standard~name~ = \"cloud~areafractionforhighclouds~\" ;
> CLDHGH:vmax = 1.e+15f ; CLDHGH:vmin = -1.e+15f ; CLDHGH:valid~range~ =
> -1.e+15f, 1.e+15f ; float CLDLOW(time, lat, lon) ; CLDLOW:long~name~ =
> \"cloud~areafractionforlowclouds~\" ; CLDLOW:units = \"1\" ;
> CLDLOW:~FillValue~ = 1.e+15f ; CLDLOW:missing~value~ = 1.e+15f ;
> CLDLOW:fmissing~value~ = 1.e+15f ; CLDLOW:scale~factor~ = 1.f ;
> CLDLOW:add~offset~ = 0.f ; CLDLOW:standard~name~ =
> \"cloud~areafractionforlowclouds~\" ; CLDLOW:vmax = 1.e+15f ;
> CLDLOW:vmin = -1.e+15f ; CLDLOW:valid~range~ = -1.e+15f, 1.e+15f ;
> float CLDMID(time, lat, lon) ; CLDMID:long~name~ =
> \"cloud~areafractionformiddleclouds~\" ; CLDMID:units = \"1\" ;
> CLDMID:~FillValue~ = 1.e+15f ; CLDMID:missing~value~ = 1.e+15f ;
> CLDMID:fmissing~value~ = 1.e+15f ; CLDMID:scale~factor~ = 1.f ;
> CLDMID:add~offset~ = 0.f ; CLDMID:standard~name~ =
> \"cloud~areafractionformiddleclouds~\" ; CLDMID:vmax = 1.e+15f ;
> CLDMID:vmin = -1.e+15f ; CLDMID:valid~range~ = -1.e+15f, 1.e+15f ;
> float CLDTOT(time, lat, lon) ; CLDTOT:long~name~ =
> \"total~cloudareafraction~\" ; CLDTOT:units = \"1\" ;
> CLDTOT:~FillValue~ = 1.e+15f ; CLDTOT:missing~value~ = 1.e+15f ;
> CLDTOT:fmissing~value~ = 1.e+15f ; CLDTOT:scale~factor~ = 1.f ;
> CLDTOT:add~offset~ = 0.f ; CLDTOT:standard~name~ =
> \"total~cloudareafraction~\" ; CLDTOT:vmax = 1.e+15f ; CLDTOT:vmin =
> -1.e+15f ; CLDTOT:valid~range~ = -1.e+15f, 1.e+15f ; float EMIS(time,
> lat, lon) ; EMIS:long~name~ = \"surface~emissivity~\" ; EMIS:units =
> \"1\" ; EMIS:~FillValue~ = 1.e+15f ; EMIS:missing~value~ = 1.e+15f ;
> EMIS:fmissing~value~ = 1.e+15f ; EMIS:scale~factor~ = 1.f ;
> EMIS:add~offset~ = 0.f ; EMIS:standard~name~ = \"surface~emissivity~\"
> ; EMIS:vmax = 1.e+15f ; EMIS:vmin = -1.e+15f ; EMIS:valid~range~ =
> -1.e+15f, 1.e+15f ; float LWGAB(time, lat, lon) ; LWGAB:long~name~ =
> \"surface~absorbedlongwaveradiation~\" ; LWGAB:units = \"W m-2\" ;
> LWGAB:~FillValue~ = 1.e+15f ; LWGAB:missing~value~ = 1.e+15f ;
> LWGAB:fmissing~value~ = 1.e+15f ; LWGAB:scale~factor~ = 1.f ;
> LWGAB:add~offset~ = 0.f ; LWGAB:standard~name~ =
> \"surface~absorbedlongwaveradiation~\" ; LWGAB:vmax = 1.e+15f ;
> LWGAB:vmin = -1.e+15f ; LWGAB:valid~range~ = -1.e+15f, 1.e+15f ; float
> LWGABCLR(time, lat, lon) ; LWGABCLR:long~name~ =
> \"surface~absorbedlongwaveradiationassumingclearsky~\" ;
> LWGABCLR:units = \"W m-2\" ; LWGABCLR:~FillValue~ = 1.e+15f ;
> LWGABCLR:missing~value~ = 1.e+15f ; LWGABCLR:fmissing~value~ = 1.e+15f
> ; LWGABCLR:scale~factor~ = 1.f ; LWGABCLR:add~offset~ = 0.f ;
> LWGABCLR:standard~name~ =
> \"surface~absorbedlongwaveradiationassumingclearsky~\" ; LWGABCLR:vmax
> = 1.e+15f ; LWGABCLR:vmin = -1.e+15f ; LWGABCLR:valid~range~ =
> -1.e+15f, 1.e+15f ; float LWGABCLRCLN(time, lat, lon) ;
> LWGABCLRCLN:long~name~ =
> \"surface~absorbedlongwaveradiationassumingclearskyandnoaerosol~\" ;
> LWGABCLRCLN:units = \"W m-2\" ; LWGABCLRCLN:~FillValue~ = 1.e+15f ;
> LWGABCLRCLN:missing~value~ = 1.e+15f ; LWGABCLRCLN:fmissing~value~ =
> 1.e+15f ; LWGABCLRCLN:scale~factor~ = 1.f ; LWGABCLRCLN:add~offset~ =
> 0.f ; LWGABCLRCLN:standard~name~ =
> \"surface~absorbedlongwaveradiationassumingclearskyandnoaerosol~\" ;
> LWGABCLRCLN:vmax = 1.e+15f ; LWGABCLRCLN:vmin = -1.e+15f ;
> LWGABCLRCLN:valid~range~ = -1.e+15f, 1.e+15f ; float LWGEM(time, lat,
> lon) ; LWGEM:long~name~ = \"longwave~fluxemittedfromsurface~\" ;
> LWGEM:units = \"W m-2\" ; LWGEM:~FillValue~ = 1.e+15f ;
> LWGEM:missing~value~ = 1.e+15f ; LWGEM:fmissing~value~ = 1.e+15f ;
> LWGEM:scale~factor~ = 1.f ; LWGEM:add~offset~ = 0.f ;
> LWGEM:standard~name~ = \"longwave~fluxemittedfromsurface~\" ;
> LWGEM:vmax = 1.e+15f ; LWGEM:vmin = -1.e+15f ; LWGEM:valid~range~ =
> -1.e+15f, 1.e+15f ; float LWGNT(time, lat, lon) ; LWGNT:long~name~ =
> \"surface~netdownwardlongwaveflux~\" ; LWGNT:units = \"W m-2\" ;
> LWGNT:~FillValue~ = 1.e+15f ; LWGNT:missing~value~ = 1.e+15f ;
> LWGNT:fmissing~value~ = 1.e+15f ; LWGNT:scale~factor~ = 1.f ;
> LWGNT:add~offset~ = 0.f ; LWGNT:standard~name~ =
> \"surface~netdownwardlongwaveflux~\" ; LWGNT:vmax = 1.e+15f ;
> LWGNT:vmin = -1.e+15f ; LWGNT:valid~range~ = -1.e+15f, 1.e+15f ; float
> LWGNTCLR(time, lat, lon) ; LWGNTCLR:long~name~ =
> \"surface~netdownwardlongwavefluxassumingclearsky~\" ; LWGNTCLR:units
> = \"W m-2\" ; LWGNTCLR:~FillValue~ = 1.e+15f ; LWGNTCLR:missing~value~
> = 1.e+15f ; LWGNTCLR:fmissing~value~ = 1.e+15f ;
> LWGNTCLR:scale~factor~ = 1.f ; LWGNTCLR:add~offset~ = 0.f ;
> LWGNTCLR:standard~name~ =
> \"surface~netdownwardlongwavefluxassumingclearsky~\" ; LWGNTCLR:vmax =
> 1.e+15f ; LWGNTCLR:vmin = -1.e+15f ; LWGNTCLR:valid~range~ = -1.e+15f,
> 1.e+15f ; float LWGNTCLRCLN(time, lat, lon) ; LWGNTCLRCLN:long~name~ =
> \"surface~netdownwardlongwavefluxassumingclearskyandnoaerosol~\" ;
> LWGNTCLRCLN:units = \"W m-2\" ; LWGNTCLRCLN:~FillValue~ = 1.e+15f ;
> LWGNTCLRCLN:missing~value~ = 1.e+15f ; LWGNTCLRCLN:fmissing~value~ =
> 1.e+15f ; LWGNTCLRCLN:scale~factor~ = 1.f ; LWGNTCLRCLN:add~offset~ =
> 0.f ; LWGNTCLRCLN:standard~name~ =
> \"surface~netdownwardlongwavefluxassumingclearskyandnoaerosol~\" ;
> LWGNTCLRCLN:vmax = 1.e+15f ; LWGNTCLRCLN:vmin = -1.e+15f ;
> LWGNTCLRCLN:valid~range~ = -1.e+15f, 1.e+15f ; float LWTUP(time, lat,
> lon) ; LWTUP:long~name~ = \"upwelling~longwavefluxattoa~\" ;
> LWTUP:units = \"W m-2\" ; LWTUP:~FillValue~ = 1.e+15f ;
> LWTUP:missing~value~ = 1.e+15f ; LWTUP:fmissing~value~ = 1.e+15f ;
> LWTUP:scale~factor~ = 1.f ; LWTUP:add~offset~ = 0.f ;
> LWTUP:standard~name~ = \"upwelling~longwavefluxattoa~\" ; LWTUP:vmax =
> 1.e+15f ; LWTUP:vmin = -1.e+15f ; LWTUP:valid~range~ = -1.e+15f,
> 1.e+15f ; float LWTUPCLR(time, lat, lon) ; LWTUPCLR:long~name~ =
> \"upwelling~longwavefluxattoaassumingclearsky~\" ; LWTUPCLR:units =
> \"W m-2\" ; LWTUPCLR:~FillValue~ = 1.e+15f ; LWTUPCLR:missing~value~ =
> 1.e+15f ; LWTUPCLR:fmissing~value~ = 1.e+15f ; LWTUPCLR:scale~factor~
> = 1.f ; LWTUPCLR:add~offset~ = 0.f ; LWTUPCLR:standard~name~ =
> \"upwelling~longwavefluxattoaassumingclearsky~\" ; LWTUPCLR:vmax =
> 1.e+15f ; LWTUPCLR:vmin = -1.e+15f ; LWTUPCLR:valid~range~ = -1.e+15f,
> 1.e+15f ; float LWTUPCLRCLN(time, lat, lon) ; LWTUPCLRCLN:long~name~ =
> \"upwelling~longwavefluxattoaassumingclearskyandnoaerosol~\" ;
> LWTUPCLRCLN:units = \"W m-2\" ; LWTUPCLRCLN:~FillValue~ = 1.e+15f ;
> LWTUPCLRCLN:missing~value~ = 1.e+15f ; LWTUPCLRCLN:fmissing~value~ =
> 1.e+15f ; LWTUPCLRCLN:scale~factor~ = 1.f ; LWTUPCLRCLN:add~offset~ =
> 0.f ; LWTUPCLRCLN:standard~name~ =
> \"upwelling~longwavefluxattoaassumingclearskyandnoaerosol~\" ;
> LWTUPCLRCLN:vmax = 1.e+15f ; LWTUPCLRCLN:vmin = -1.e+15f ;
> LWTUPCLRCLN:valid~range~ = -1.e+15f, 1.e+15f ; float SWGDN(time, lat,
> lon) ; SWGDN:long~name~ = \"surface~incomingshortwaveflux~\" ;
> SWGDN:units = \"W m-2\" ; SWGDN:~FillValue~ = 1.e+15f ;
> SWGDN:missing~value~ = 1.e+15f ; SWGDN:fmissing~value~ = 1.e+15f ;
> SWGDN:scale~factor~ = 1.f ; SWGDN:add~offset~ = 0.f ;
> SWGDN:standard~name~ = \"surface~incomingshortwaveflux~\" ; SWGDN:vmax
> = 1.e+15f ; SWGDN:vmin = -1.e+15f ; SWGDN:valid~range~ = -1.e+15f,
> 1.e+15f ; float SWGDNCLR(time, lat, lon) ; SWGDNCLR:long~name~ =
> \"surface~incomingshortwavefluxassumingclearsky~\" ; SWGDNCLR:units =
> \"W m-2\" ; SWGDNCLR:~FillValue~ = 1.e+15f ; SWGDNCLR:missing~value~ =
> 1.e+15f ; SWGDNCLR:fmissing~value~ = 1.e+15f ; SWGDNCLR:scale~factor~
> = 1.f ; SWGDNCLR:add~offset~ = 0.f ; SWGDNCLR:standard~name~ =
> \"surface~incomingshortwavefluxassumingclearsky~\" ; SWGDNCLR:vmax =
> 1.e+15f ; SWGDNCLR:vmin = -1.e+15f ; SWGDNCLR:valid~range~ = -1.e+15f,
> 1.e+15f ; float SWGNT(time, lat, lon) ; SWGNT:long~name~ =
> \"surface~netdownwardshortwaveflux~\" ; SWGNT:units = \"W m-2\" ;
> SWGNT:~FillValue~ = 1.e+15f ; SWGNT:missing~value~ = 1.e+15f ;
> SWGNT:fmissing~value~ = 1.e+15f ; SWGNT:scale~factor~ = 1.f ;
> SWGNT:add~offset~ = 0.f ; SWGNT:standard~name~ =
> \"surface~netdownwardshortwaveflux~\" ; SWGNT:vmax = 1.e+15f ;
> SWGNT:vmin = -1.e+15f ; SWGNT:valid~range~ = -1.e+15f, 1.e+15f ; float
> SWGNTCLN(time, lat, lon) ; SWGNTCLN:long~name~ =
> \"surface~netdownwardshortwavefluxassumingnoaerosol~\" ;
> SWGNTCLN:units = \"W m-2\" ; SWGNTCLN:~FillValue~ = 1.e+15f ;
> SWGNTCLN:missing~value~ = 1.e+15f ; SWGNTCLN:fmissing~value~ = 1.e+15f
> ; SWGNTCLN:scale~factor~ = 1.f ; SWGNTCLN:add~offset~ = 0.f ;
> SWGNTCLN:standard~name~ =
> \"surface~netdownwardshortwavefluxassumingnoaerosol~\" ; SWGNTCLN:vmax
> = 1.e+15f ; SWGNTCLN:vmin = -1.e+15f ; SWGNTCLN:valid~range~ =
> -1.e+15f, 1.e+15f ; float SWGNTCLR(time, lat, lon) ;
> SWGNTCLR:long~name~ =
> \"surface~netdownwardshortwavefluxassumingclearsky~\" ; SWGNTCLR:units
> = \"W m-2\" ; SWGNTCLR:~FillValue~ = 1.e+15f ; SWGNTCLR:missing~value~
> = 1.e+15f ; SWGNTCLR:fmissing~value~ = 1.e+15f ;
> SWGNTCLR:scale~factor~ = 1.f ; SWGNTCLR:add~offset~ = 0.f ;
> SWGNTCLR:standard~name~ =
> \"surface~netdownwardshortwavefluxassumingclearsky~\" ; SWGNTCLR:vmax
> = 1.e+15f ; SWGNTCLR:vmin = -1.e+15f ; SWGNTCLR:valid~range~ =
> -1.e+15f, 1.e+15f ; float SWGNTCLRCLN(time, lat, lon) ;
> SWGNTCLRCLN:long~name~ =
> \"surface~netdownwardshortwavefluxassumingclearskyandnoaerosol~\" ;
> SWGNTCLRCLN:units = \"W m-2\" ; SWGNTCLRCLN:~FillValue~ = 1.e+15f ;
> SWGNTCLRCLN:missing~value~ = 1.e+15f ; SWGNTCLRCLN:fmissing~value~ =
> 1.e+15f ; SWGNTCLRCLN:scale~factor~ = 1.f ; SWGNTCLRCLN:add~offset~ =
> 0.f ; SWGNTCLRCLN:standard~name~ =
> \"surface~netdownwardshortwavefluxassumingclearskyandnoaerosol~\" ;
> SWGNTCLRCLN:vmax = 1.e+15f ; SWGNTCLRCLN:vmin = -1.e+15f ;
> SWGNTCLRCLN:valid~range~ = -1.e+15f, 1.e+15f ; float SWTDN(time, lat,
> lon) ; SWTDN:long~name~ = \"toa~incomingshortwaveflux~\" ; SWTDN:units
> = \"W m-2\" ; SWTDN:~FillValue~ = 1.e+15f ; SWTDN:missing~value~ =
> 1.e+15f ; SWTDN:fmissing~value~ = 1.e+15f ; SWTDN:scale~factor~ = 1.f
> ; SWTDN:add~offset~ = 0.f ; SWTDN:standard~name~ =
> \"toa~incomingshortwaveflux~\" ; SWTDN:vmax = 1.e+15f ; SWTDN:vmin =
> -1.e+15f ; SWTDN:valid~range~ = -1.e+15f, 1.e+15f ; float SWTNT(time,
> lat, lon) ; SWTNT:long~name~ = \"toa~netdownwardshortwaveflux~\" ;
> SWTNT:units = \"W m-2\" ; SWTNT:~FillValue~ = 1.e+15f ;
> SWTNT:missing~value~ = 1.e+15f ; SWTNT:fmissing~value~ = 1.e+15f ;
> SWTNT:scale~factor~ = 1.f ; SWTNT:add~offset~ = 0.f ;
> SWTNT:standard~name~ = \"toa~netdownwardshortwaveflux~\" ; SWTNT:vmax
> = 1.e+15f ; SWTNT:vmin = -1.e+15f ; SWTNT:valid~range~ = -1.e+15f,
> 1.e+15f ; float SWTNTCLN(time, lat, lon) ; SWTNTCLN:long~name~ =
> \"toa~netdownwardshortwavefluxassumingnoaerosol~\" ; SWTNTCLN:units =
> \"W m-2\" ; SWTNTCLN:~FillValue~ = 1.e+15f ; SWTNTCLN:missing~value~ =
> 1.e+15f ; SWTNTCLN:fmissing~value~ = 1.e+15f ; SWTNTCLN:scale~factor~
> = 1.f ; SWTNTCLN:add~offset~ = 0.f ; SWTNTCLN:standard~name~ =
> \"toa~netdownwardshortwavefluxassumingnoaerosol~\" ; SWTNTCLN:vmax =
> 1.e+15f ; SWTNTCLN:vmin = -1.e+15f ; SWTNTCLN:valid~range~ = -1.e+15f,
> 1.e+15f ; float SWTNTCLR(time, lat, lon) ; SWTNTCLR:long~name~ =
> \"toa~netdownwardshortwavefluxassumingclearsky~\" ; SWTNTCLR:units =
> \"W m-2\" ; SWTNTCLR:~FillValue~ = 1.e+15f ; SWTNTCLR:missing~value~ =
> 1.e+15f ; SWTNTCLR:fmissing~value~ = 1.e+15f ; SWTNTCLR:scale~factor~
> = 1.f ; SWTNTCLR:add~offset~ = 0.f ; SWTNTCLR:standard~name~ =
> \"toa~netdownwardshortwavefluxassumingclearsky~\" ; SWTNTCLR:vmax =
> 1.e+15f ; SWTNTCLR:vmin = -1.e+15f ; SWTNTCLR:valid~range~ = -1.e+15f,
> 1.e+15f ; float SWTNTCLRCLN(time, lat, lon) ; SWTNTCLRCLN:long~name~ =
> \"toa~netdownwardshortwavefluxassumingclearskyandnoaerosol~\" ;
> SWTNTCLRCLN:units = \"W m-2\" ; SWTNTCLRCLN:~FillValue~ = 1.e+15f ;
> SWTNTCLRCLN:missing~value~ = 1.e+15f ; SWTNTCLRCLN:fmissing~value~ =
> 1.e+15f ; SWTNTCLRCLN:scale~factor~ = 1.f ; SWTNTCLRCLN:add~offset~ =
> 0.f ; SWTNTCLRCLN:standard~name~ =
> \"toa~netdownwardshortwavefluxassumingclearskyandnoaerosol~\" ;
> SWTNTCLRCLN:vmax = 1.e+15f ; SWTNTCLRCLN:vmin = -1.e+15f ;
> SWTNTCLRCLN:valid~range~ = -1.e+15f, 1.e+15f ; float TAUHGH(time, lat,
> lon) ; TAUHGH:long~name~ =
> \"in~cloudopticalthicknessofhighclouds~(EXPORT)\" ; TAUHGH:units =
> \"1\" ; TAUHGH:~FillValue~ = 1.e+15f ; TAUHGH:missing~value~ = 1.e+15f
> ; TAUHGH:fmissing~value~ = 1.e+15f ; TAUHGH:scale~factor~ = 1.f ;
> TAUHGH:add~offset~ = 0.f ; TAUHGH:standard~name~ =
> \"in~cloudopticalthicknessofhighclouds~(EXPORT)\" ; TAUHGH:vmax =
> 1.e+15f ; TAUHGH:vmin = -1.e+15f ; TAUHGH:valid~range~ = -1.e+15f,
> 1.e+15f ; float TAULOW(time, lat, lon) ; TAULOW:long~name~ =
> \"in~cloudopticalthicknessoflowclouds~\" ; TAULOW:units = \"1\" ;
> TAULOW:~FillValue~ = 1.e+15f ; TAULOW:missing~value~ = 1.e+15f ;
> TAULOW:fmissing~value~ = 1.e+15f ; TAULOW:scale~factor~ = 1.f ;
> TAULOW:add~offset~ = 0.f ; TAULOW:standard~name~ =
> \"in~cloudopticalthicknessoflowclouds~\" ; TAULOW:vmax = 1.e+15f ;
> TAULOW:vmin = -1.e+15f ; TAULOW:valid~range~ = -1.e+15f, 1.e+15f ;
> float TAUMID(time, lat, lon) ; TAUMID:long~name~ =
> \"in~cloudopticalthicknessofmiddleclouds~\" ; TAUMID:units = \"1\" ;
> TAUMID:~FillValue~ = 1.e+15f ; TAUMID:missing~value~ = 1.e+15f ;
> TAUMID:fmissing~value~ = 1.e+15f ; TAUMID:scale~factor~ = 1.f ;
> TAUMID:add~offset~ = 0.f ; TAUMID:standard~name~ =
> \"in~cloudopticalthicknessofmiddleclouds~\" ; TAUMID:vmax = 1.e+15f ;
> TAUMID:vmin = -1.e+15f ; TAUMID:valid~range~ = -1.e+15f, 1.e+15f ;
> float TAUTOT(time, lat, lon) ; TAUTOT:long~name~ =
> \"in~cloudopticalthicknessofallclouds~\" ; TAUTOT:units = \"1\" ;
> TAUTOT:~FillValue~ = 1.e+15f ; TAUTOT:missing~value~ = 1.e+15f ;
> TAUTOT:fmissing~value~ = 1.e+15f ; TAUTOT:scale~factor~ = 1.f ;
> TAUTOT:add~offset~ = 0.f ; TAUTOT:standard~name~ =
> \"in~cloudopticalthicknessofallclouds~\" ; TAUTOT:vmax = 1.e+15f ;
> TAUTOT:vmin = -1.e+15f ; TAUTOT:valid~range~ = -1.e+15f, 1.e+15f ;
> float TS(time, lat, lon) ; TS:long~name~ =
> \"surface~skintemperature~\" ; TS:units = \"K\" ; TS:~FillValue~ =
> 1.e+15f ; TS:missing~value~ = 1.e+15f ; TS:fmissing~value~ = 1.e+15f ;
> TS:scale~factor~ = 1.f ; TS:add~offset~ = 0.f ; TS:standard~name~ =
> \"surface~skintemperature~\" ; TS:vmax = 1.e+15f ; TS:vmin = -1.e+15f
> ; TS:valid~range~ = -1.e+15f, 1.e+15f ;
>
> // global attributes: :History = \"Original file generated: Mon Mar 23
> 03:29:51 2015 GMT\" ; :Comment = \"GMAO filename:
> d5124~m2jan00~.tavg1~2dradNx~.20100101.nc4\" ; :Filename =
> \"MERRA2~300~.tavg1~2dradNx~.20100101.nc4\" ; :Conventions = \"CF-1\"
> ; :Institution = \"NASA Global Modeling and Assimilation Office\" ;
> :References = \"<http://gmao.gsfc.nasa.gov>\" ; :Format =
> \"NetCDF-4/HDF-5\" ; :SpatialCoverage = \"global\" ; :VersionID =
> \"5.12.4\" ; :TemporalRange = \"1980-01-01 -\> 2016-12-31\" ;
> :identifier~productdoiauthority~ = \"<http://dx.doi.org/>\" ;
> :ShortName = \"M2T1NXRAD\" ; :GranuleID =
> \"MERRA2~300~.tavg1~2dradNx~.20100101.nc4\" ; :ProductionDateTime =
> \"Original file generated: Mon Mar 23 03:29:51 2015 GMT\" ; :LongName
> = \"MERRA2 tavg1~2dradNx~:
> 2d,1-Hourly,Time-Averaged,Single-Level,Assimilation,Radiation
> Diagnostics\" ; :Title = \"MERRA2 tavg1~2dradNx~:
> 2d,1-Hourly,Time-Averaged,Single-Level,Assimilation,Radiation
> Diagnostics\" ; :SouthernmostLatitude = \"-90.0\" ;
> :NorthernmostLatitude = \"90.0\" ; :WesternmostLongitude = \"-180.0\"
> ; :EasternmostLongitude = \"179.375\" ; :LatitudeResolution = \"0.5\"
> ; :LongitudeResolution = \"0.625\" ; :DataResolution = \"0.5 x 0.625\"
> ; :Source = \"CVS tag: GEOSadas-5~124~\" ; :Contact =
> \"<http://gmao.gsfc.nasa.gov>\" ; :identifier~productdoi~ =
> \"10.5067/Q9QMY5PBNV1T\" ; :RangeBeginningDate = \"2010-01-01\" ;
> :RangeBeginningTime = \"00:00:00.000000\" ; :RangeEndingDate =
> \"2010-01-01\" ; :RangeEndingTime = \"23:59:59.000000\" ; }

### Single-level Diagnostics (\'slv\')

## SOSE

The following files were downloaded from SOSE [Iteration
122](http://sose.ucsd.edu/BSOSE6_iter122_solution.html) (2013-2017).

### Bottom pressure

> netcdf bsose~i1222013to20171dayBottomPres~ { dimensions: time = 1826 ;
> XC = 2160 ; YC = 588 ; variables: int64 iter(time) ; iter:long~name~ =
> \"model timestep number\" ; iter:standard~name~ = \"timestep\" ; int64
> time(time) ; time:long~name~ = \"Time\" ; time:standard~name~ =
> \"time\" ; time:axis = \"T\" ; time:units = \"seconds since
> 2012-12-01\" ; time:calendar = \"proleptic~gregorian~\" ; float XC(XC)
> ; XC:~FillValue~ = NaNf ; XC:coordinate = \"YC XC\" ; XC:units =
> \"degrees~east~\" ; XC:standard~name~ = \"longitude\" ; XC:long~name~
> = \"longitude\" ; XC:axis = \"X\" ; float YC(YC) ; YC:~FillValue~ =
> NaNf ; YC:coordinate = \"YC XC\" ; YC:units = \"degrees~north~\" ;
> YC:standard~name~ = \"latitude\" ; YC:long~name~ = \"latitude\" ;
> YC:axis = \"Y\" ; float Depth(YC, XC) ; Depth:~FillValue~ = NaNf ;
> Depth:coordinate = \"XC YC\" ; Depth:units = \"m\" ;
> Depth:standard~name~ = \"ocean~depth~\" ; Depth:long~name~ = \"ocean
> depth\" ; float rA(YC, XC) ; rA:~FillValue~ = NaNf ; rA:coordinate =
> \"YC XC\" ; rA:units = \"m2\" ; rA:standard~name~ = \"cell~area~\" ;
> rA:long~name~ = \"cell area\" ; float PHIBOT(time, YC, XC) ;
> PHIBOT:~FillValue~ = NaNf ; PHIBOT:units = \"m^2^/s^2^\" ;
> PHIBOT:long~name~ = \"Bottom Pressure Pot.(p/rho) Anomaly\" ;
> PHIBOT:standard~name~ = \"PHIBOT\" ; PHIBOT:coordinates = \"Depth rA
> iter\" ; }

### Fresh water from atmosphere and land

> netcdf bsose~i1222013to20171dayFWfrmAtmLnd~ { dimensions: time = 1826
> ; XC = 2160 ; YC = 588 ; variables: int64 iter(time) ; iter:long~name~
> = \"model timestep number\" ; iter:standard~name~ = \"timestep\" ;
> int64 time(time) ; time:long~name~ = \"Time\" ; time:standard~name~ =
> \"time\" ; time:axis = \"T\" ; time:units = \"seconds since
> 2012-12-01\" ; time:calendar = \"proleptic~gregorian~\" ; float XC(XC)
> ; XC:~FillValue~ = NaNf ; XC:coordinate = \"YC XC\" ; XC:units =
> \"degrees~east~\" ; XC:standard~name~ = \"longitude\" ; XC:long~name~
> = \"longitude\" ; XC:axis = \"X\" ; float YC(YC) ; YC:~FillValue~ =
> NaNf ; YC:coordinate = \"YC XC\" ; YC:units = \"degrees~north~\" ;
> YC:standard~name~ = \"latitude\" ; YC:long~name~ = \"latitude\" ;
> YC:axis = \"Y\" ; float Depth(YC, XC) ; Depth:~FillValue~ = NaNf ;
> Depth:coordinate = \"XC YC\" ; Depth:units = \"m\" ;
> Depth:standard~name~ = \"ocean~depth~\" ; Depth:long~name~ = \"ocean
> depth\" ; float rA(YC, XC) ; rA:~FillValue~ = NaNf ; rA:coordinate =
> \"YC XC\" ; rA:units = \"m2\" ; rA:standard~name~ = \"cell~area~\" ;
> rA:long~name~ = \"cell area\" ; float SIatmFW(time, YC, XC) ;
> SIatmFW:~FillValue~ = NaNf ; SIatmFW:units = \"kg/m^2^/s\" ;
> SIatmFW:long~name~ = \"Net freshwater flux from atmosphere & land
> (+=down)\" ; SIatmFW:standard~name~ = \"SIatmFW\" ;
> SIatmFW:coordinates = \"Depth rA iter\" ; }

### Mixed layer depth

> netcdf bsose~i1222013to20171dayMLD~ { dimensions: time = 1826 ; XC =
> 2160 ; YC = 588 ; variables: int64 iter(time) ; iter:long~name~ =
> \"model timestep number\" ; iter:standard~name~ = \"timestep\" ; int64
> time(time) ; time:long~name~ = \"Time\" ; time:standard~name~ =
> \"time\" ; time:axis = \"T\" ; time:units = \"seconds since
> 2012-12-01\" ; time:calendar = \"proleptic~gregorian~\" ; float XC(XC)
> ; XC:~FillValue~ = NaNf ; XC:coordinate = \"YC XC\" ; XC:units =
> \"degrees~east~\" ; XC:standard~name~ = \"longitude\" ; XC:long~name~
> = \"longitude\" ; XC:axis = \"X\" ; float YC(YC) ; YC:~FillValue~ =
> NaNf ; YC:coordinate = \"YC XC\" ; YC:units = \"degrees~north~\" ;
> YC:standard~name~ = \"latitude\" ; YC:long~name~ = \"latitude\" ;
> YC:axis = \"Y\" ; float Depth(YC, XC) ; Depth:~FillValue~ = NaNf ;
> Depth:coordinate = \"XC YC\" ; Depth:units = \"m\" ;
> Depth:standard~name~ = \"ocean~depth~\" ; Depth:long~name~ = \"ocean
> depth\" ; float rA(YC, XC) ; rA:~FillValue~ = NaNf ; rA:coordinate =
> \"YC XC\" ; rA:units = \"m2\" ; rA:standard~name~ = \"cell~area~\" ;
> rA:long~name~ = \"cell area\" ; float BLGMLD(time, YC, XC) ;
> BLGMLD:~FillValue~ = NaNf ; BLGMLD:units = \"m\" ; BLGMLD:long~name~ =
> \"Diagnosed mixed layer depth\" ; BLGMLD:standard~name~ = \"BLGMLD\" ;
> BLGMLD:coordinates = \"Depth rA iter\" ; }

### Sea surface height anomaly

> netcdf bsose~i1222013to20171daySSH~ { dimensions: time = 1826 ; XC =
> 2160 ; YC = 588 ; variables: int64 iter(time) ; iter:long~name~ =
> \"model timestep number\" ; iter:standard~name~ = \"timestep\" ; int64
> time(time) ; time:long~name~ = \"Time\" ; time:standard~name~ =
> \"time\" ; time:axis = \"T\" ; time:units = \"seconds since
> 2012-12-01\" ; time:calendar = \"proleptic~gregorian~\" ; float XC(XC)
> ; XC:~FillValue~ = NaNf ; XC:coordinate = \"YC XC\" ; XC:units =
> \"degrees~east~\" ; XC:standard~name~ = \"longitude\" ; XC:long~name~
> = \"longitude\" ; XC:axis = \"X\" ; float YC(YC) ; YC:~FillValue~ =
> NaNf ; YC:coordinate = \"YC XC\" ; YC:units = \"degrees~north~\" ;
> YC:standard~name~ = \"latitude\" ; YC:long~name~ = \"latitude\" ;
> YC:axis = \"Y\" ; float Depth(YC, XC) ; Depth:~FillValue~ = NaNf ;
> Depth:coordinate = \"XC YC\" ; Depth:units = \"m\" ;
> Depth:standard~name~ = \"ocean~depth~\" ; Depth:long~name~ = \"ocean
> depth\" ; float rA(YC, XC) ; rA:~FillValue~ = NaNf ; rA:coordinate =
> \"YC XC\" ; rA:units = \"m2\" ; rA:standard~name~ = \"cell~area~\" ;
> rA:long~name~ = \"cell area\" ; float ETAN(time, YC, XC) ;
> ETAN:~FillValue~ = NaNf ; ETAN:units = \"m\" ; ETAN:long~name~ =
> \"Surface Height Anomaly\" ; ETAN:standard~name~ = \"ETAN\" ;
> ETAN:coordinates = \"Depth rA iter\" ; }

### Sea ice area

> netcdf bsose~i1222013to20171daySeaIceArea~ { dimensions: time = 1826 ;
> XC = 2160 ; YC = 588 ; variables: int64 iter(time) ; iter:long~name~ =
> \"model timestep number\" ; iter:standard~name~ = \"timestep\" ; int64
> time(time) ; time:long~name~ = \"Time\" ; time:standard~name~ =
> \"time\" ; time:axis = \"T\" ; time:units = \"seconds since
> 2012-12-01\" ; time:calendar = \"proleptic~gregorian~\" ; float XC(XC)
> ; XC:~FillValue~ = NaNf ; XC:coordinate = \"YC XC\" ; XC:units =
> \"degrees~east~\" ; XC:standard~name~ = \"longitude\" ; XC:long~name~
> = \"longitude\" ; XC:axis = \"X\" ; float YC(YC) ; YC:~FillValue~ =
> NaNf ; YC:coordinate = \"YC XC\" ; YC:units = \"degrees~north~\" ;
> YC:standard~name~ = \"latitude\" ; YC:long~name~ = \"latitude\" ;
> YC:axis = \"Y\" ; float Depth(YC, XC) ; Depth:~FillValue~ = NaNf ;
> Depth:coordinate = \"XC YC\" ; Depth:units = \"m\" ;
> Depth:standard~name~ = \"ocean~depth~\" ; Depth:long~name~ = \"ocean
> depth\" ; float rA(YC, XC) ; rA:~FillValue~ = NaNf ; rA:coordinate =
> \"YC XC\" ; rA:units = \"m2\" ; rA:standard~name~ = \"cell~area~\" ;
> rA:long~name~ = \"cell area\" ; float SIarea(time, YC, XC) ;
> SIarea:~FillValue~ = NaNf ; SIarea:units = \"m^2^/m^2^\" ;
> SIarea:long~name~ = \"SEAICE fractional ice-covered area \[0 to 1\]\"
> ; SIarea:standard~name~ = \"SIarea\" ; SIarea:coordinates = \"Depth rA
> iter\" ; }

### Sea ice covered part of SIqnet

> netcdf bsose~i1222013to20171daySeaIceCvrdQnet~ { dimensions: time =
> 1826 ; XC = 2160 ; YC = 588 ; variables: int64 iter(time) ;
> iter:long~name~ = \"model timestep number\" ; iter:standard~name~ =
> \"timestep\" ; int64 time(time) ; time:long~name~ = \"Time\" ;
> time:standard~name~ = \"time\" ; time:axis = \"T\" ; time:units =
> \"seconds since 2012-12-01\" ; time:calendar =
> \"proleptic~gregorian~\" ; float XC(XC) ; XC:~FillValue~ = NaNf ;
> XC:coordinate = \"YC XC\" ; XC:units = \"degrees~east~\" ;
> XC:standard~name~ = \"longitude\" ; XC:long~name~ = \"longitude\" ;
> XC:axis = \"X\" ; float YC(YC) ; YC:~FillValue~ = NaNf ; YC:coordinate
> = \"YC XC\" ; YC:units = \"degrees~north~\" ; YC:standard~name~ =
> \"latitude\" ; YC:long~name~ = \"latitude\" ; YC:axis = \"Y\" ; float
> Depth(YC, XC) ; Depth:~FillValue~ = NaNf ; Depth:coordinate = \"XC
> YC\" ; Depth:units = \"m\" ; Depth:standard~name~ = \"ocean~depth~\" ;
> Depth:long~name~ = \"ocean depth\" ; float rA(YC, XC) ; rA:~FillValue~
> = NaNf ; rA:coordinate = \"YC XC\" ; rA:units = \"m2\" ;
> rA:standard~name~ = \"cell~area~\" ; rA:long~name~ = \"cell area\" ;
> float SIqneti(time, YC, XC) ; SIqneti:~FillValue~ = NaNf ;
> SIqneti:units = \"W/m^2^\" ; SIqneti:long~name~ = \"Ice Covered Part
> of SIqnet, turb+rad, \>0 decr theta\" ; SIqneti:standard~name~ =
> \"SIqneti\" ; SIqneti:coordinates = \"Depth rA iter\" ; }

### Ocean surface freshwater flux

> netcdf bsose~i1222013to20171daySeaIceEmPmR~ { dimensions: time = 1826
> ; XC = 2160 ; YC = 588 ; variables: int64 iter(time) ; iter:long~name~
> = \"model timestep number\" ; iter:standard~name~ = \"timestep\" ;
> int64 time(time) ; time:long~name~ = \"Time\" ; time:standard~name~ =
> \"time\" ; time:axis = \"T\" ; time:units = \"seconds since
> 2012-12-01\" ; time:calendar = \"proleptic~gregorian~\" ; float XC(XC)
> ; XC:~FillValue~ = NaNf ; XC:coordinate = \"YC XC\" ; XC:units =
> \"degrees~east~\" ; XC:standard~name~ = \"longitude\" ; XC:long~name~
> = \"longitude\" ; XC:axis = \"X\" ; float YC(YC) ; YC:~FillValue~ =
> NaNf ; YC:coordinate = \"YC XC\" ; YC:units = \"degrees~north~\" ;
> YC:standard~name~ = \"latitude\" ; YC:long~name~ = \"latitude\" ;
> YC:axis = \"Y\" ; float Depth(YC, XC) ; Depth:~FillValue~ = NaNf ;
> Depth:coordinate = \"XC YC\" ; Depth:units = \"m\" ;
> Depth:standard~name~ = \"ocean~depth~\" ; Depth:long~name~ = \"ocean
> depth\" ; float rA(YC, XC) ; rA:~FillValue~ = NaNf ; rA:coordinate =
> \"YC XC\" ; rA:units = \"m2\" ; rA:standard~name~ = \"cell~area~\" ;
> rA:long~name~ = \"cell area\" ; float SIempmr(time, YC, XC) ;
> SIempmr:~FillValue~ = NaNf ; SIempmr:units = \"kg/m^2^/s\" ;
> SIempmr:long~name~ = \"Ocean surface freshwater flux, \> 0 increases
> salt\" ; SIempmr:standard~name~ = \"SIempmr\" ; SIempmr:coordinates =
> \"Depth rA iter\" ; }

### Sea ice effective ice thickness

> netcdf bsose~i1222013to20171daySeaIceHeff~ { dimensions: time = 1826 ;
> XC = 2160 ; YC = 588 ; variables: int64 iter(time) ; iter:long~name~ =
> \"model timestep number\" ; iter:standard~name~ = \"timestep\" ; int64
> time(time) ; time:long~name~ = \"Time\" ; time:standard~name~ =
> \"time\" ; time:axis = \"T\" ; time:units = \"seconds since
> 2012-12-01\" ; time:calendar = \"proleptic~gregorian~\" ; float XC(XC)
> ; XC:~FillValue~ = NaNf ; XC:coordinate = \"YC XC\" ; XC:units =
> \"degrees~east~\" ; XC:standard~name~ = \"longitude\" ; XC:long~name~
> = \"longitude\" ; XC:axis = \"X\" ; float YC(YC) ; YC:~FillValue~ =
> NaNf ; YC:coordinate = \"YC XC\" ; YC:units = \"degrees~north~\" ;
> YC:standard~name~ = \"latitude\" ; YC:long~name~ = \"latitude\" ;
> YC:axis = \"Y\" ; float Depth(YC, XC) ; Depth:~FillValue~ = NaNf ;
> Depth:coordinate = \"XC YC\" ; Depth:units = \"m\" ;
> Depth:standard~name~ = \"ocean~depth~\" ; Depth:long~name~ = \"ocean
> depth\" ; float rA(YC, XC) ; rA:~FillValue~ = NaNf ; rA:coordinate =
> \"YC XC\" ; rA:units = \"m2\" ; rA:standard~name~ = \"cell~area~\" ;
> rA:long~name~ = \"cell area\" ; float SIheff(time, YC, XC) ;
> SIheff:~FillValue~ = NaNf ; SIheff:units = \"m\" ; SIheff:long~name~ =
> \"SEAICE effective ice thickness\" ; SIheff:standard~name~ =
> \"SIheff\" ; SIheff:coordinates = \"Depth rA iter\" ; }

### Sea ice surface temperature over sea ice

> netcdf bsose~i1222013to20171daySeaIceSurfT~ { dimensions: time = 1826
> ; XC = 2160 ; YC = 588 ; variables: int64 iter(time) ; iter:long~name~
> = \"model timestep number\" ; iter:standard~name~ = \"timestep\" ;
> int64 time(time) ; time:long~name~ = \"Time\" ; time:standard~name~ =
> \"time\" ; time:axis = \"T\" ; time:units = \"seconds since
> 2012-12-01\" ; time:calendar = \"proleptic~gregorian~\" ; float XC(XC)
> ; XC:~FillValue~ = NaNf ; XC:coordinate = \"YC XC\" ; XC:units =
> \"degrees~east~\" ; XC:standard~name~ = \"longitude\" ; XC:long~name~
> = \"longitude\" ; XC:axis = \"X\" ; float YC(YC) ; YC:~FillValue~ =
> NaNf ; YC:coordinate = \"YC XC\" ; YC:units = \"degrees~north~\" ;
> YC:standard~name~ = \"latitude\" ; YC:long~name~ = \"latitude\" ;
> YC:axis = \"Y\" ; float Depth(YC, XC) ; Depth:~FillValue~ = NaNf ;
> Depth:coordinate = \"XC YC\" ; Depth:units = \"m\" ;
> Depth:standard~name~ = \"ocean~depth~\" ; Depth:long~name~ = \"ocean
> depth\" ; float rA(YC, XC) ; rA:~FillValue~ = NaNf ; rA:coordinate =
> \"YC XC\" ; rA:units = \"m2\" ; rA:standard~name~ = \"cell~area~\" ;
> rA:long~name~ = \"cell area\" ; float SItices(time, YC, XC) ;
> SItices:~FillValue~ = NaNf ; SItices:units = \"K\" ;
> SItices:long~name~ = \"Surface Temperature over Sea-Ice (area
> weighted)\" ; SItices:standard~name~ = \"SItices\" ; SItices:mate =
> \"SIarea\" ; SItices:coordinates = \"Depth rA iter\" ; }

### Sea ice zonal transport effective thickness

> netcdf bsose~i1222013to20171daySeaIceUHeff~ { dimensions: time = 1826
> ; YC = 588 ; XG = 2160 ; variables: int64 iter(time) ; iter:long~name~
> = \"model timestep number\" ; iter:standard~name~ = \"timestep\" ;
> int64 time(time) ; time:long~name~ = \"Time\" ; time:standard~name~ =
> \"time\" ; time:axis = \"T\" ; time:units = \"seconds since
> 2012-12-01\" ; time:calendar = \"proleptic~gregorian~\" ; float YC(YC)
> ; YC:~FillValue~ = NaNf ; YC:coordinate = \"YC XC\" ; YC:units =
> \"degrees~north~\" ; YC:standard~name~ = \"latitude\" ; YC:long~name~
> = \"latitude\" ; YC:axis = \"Y\" ; float XG(XG) ; XG:~FillValue~ =
> NaNf ; XG:coordinate = \"YG XG\" ; XG:units = \"degrees~east~\" ;
> XG:standard~name~ = \"longitude~atflocation~\" ; XG:long~name~ =
> \"longitude\" ; XG:axis = \"X\" ; XG:c~gridaxisshift~ = -0.5 ; float
> dxC(YC, XG) ; dxC:~FillValue~ = NaNf ; dxC:coordinate = \"YC XG\" ;
> dxC:units = \"m\" ; dxC:standard~name~ = \"cell~xsizeatulocation~\" ;
> dxC:long~name~ = \"cell x size\" ; float rAw(YC, XG) ; rAw:~FillValue~
> = NaNf ; rAw:coordinate = \"YG XC\" ; rAw:units = \"m2\" ;
> rAw:standard~name~ = \"cell~areaatulocation~\" ; rAw:long~name~ =
> \"cell area\" ; float dyG(YC, XG) ; dyG:~FillValue~ = NaNf ;
> dyG:coordinate = \"YC XG\" ; dyG:units = \"m\" ; dyG:standard~name~ =
> \"cell~ysizeatulocation~\" ; dyG:long~name~ = \"cell y size\" ; float
> SIuheff(time, YC, XG) ; SIuheff:~FillValue~ = NaNf ; SIuheff:units =
> \"m^2^/s\" ; SIuheff:long~name~ = \"Zonal Transport of eff ice thickn
> (centered)\" ; SIuheff:standard~name~ = \"SIuheff\" ; SIuheff:mate =
> \"SIvheff\" ; SIuheff:coordinates = \"rAw dyG dxC iter\" ; }

### Sea ice zonal velocity

> netcdf bsose~i1222013to20171daySeaIceUvel~ { dimensions: time = 1826 ;
> YC = 588 ; XG = 2160 ; variables: int64 iter(time) ; iter:long~name~ =
> \"model timestep number\" ; iter:standard~name~ = \"timestep\" ; int64
> time(time) ; time:long~name~ = \"Time\" ; time:standard~name~ =
> \"time\" ; time:axis = \"T\" ; time:units = \"seconds since
> 2012-12-01\" ; time:calendar = \"proleptic~gregorian~\" ; float YC(YC)
> ; YC:~FillValue~ = NaNf ; YC:coordinate = \"YC XC\" ; YC:units =
> \"degrees~north~\" ; YC:standard~name~ = \"latitude\" ; YC:long~name~
> = \"latitude\" ; YC:axis = \"Y\" ; float XG(XG) ; XG:~FillValue~ =
> NaNf ; XG:coordinate = \"YG XG\" ; XG:units = \"degrees~east~\" ;
> XG:standard~name~ = \"longitude~atflocation~\" ; XG:long~name~ =
> \"longitude\" ; XG:axis = \"X\" ; XG:c~gridaxisshift~ = -0.5 ; float
> dxC(YC, XG) ; dxC:~FillValue~ = NaNf ; dxC:coordinate = \"YC XG\" ;
> dxC:units = \"m\" ; dxC:standard~name~ = \"cell~xsizeatulocation~\" ;
> dxC:long~name~ = \"cell x size\" ; float rAw(YC, XG) ; rAw:~FillValue~
> = NaNf ; rAw:coordinate = \"YG XC\" ; rAw:units = \"m2\" ;
> rAw:standard~name~ = \"cell~areaatulocation~\" ; rAw:long~name~ =
> \"cell area\" ; float dyG(YC, XG) ; dyG:~FillValue~ = NaNf ;
> dyG:coordinate = \"YC XG\" ; dyG:units = \"m\" ; dyG:standard~name~ =
> \"cell~ysizeatulocation~\" ; dyG:long~name~ = \"cell y size\" ; float
> SIuice(time, YC, XG) ; SIuice:~FillValue~ = NaNf ; SIuice:units =
> \"m/s\" ; SIuice:long~name~ = \"SEAICE zonal ice velocity, \>0 from
> West to East\" ; SIuice:standard~name~ = \"SIuice\" ; SIuice:mate =
> \"SIvice\" ; SIuice:coordinates = \"rAw dyG dxC iter\" ; }

### Sea ice meridional transport effective thickness

> netcdf bsose~i1222013to20171daySeaIceVHeff~ { dimensions: time = 1826
> ; XC = 2160 ; YG = 588 ; variables: int64 iter(time) ; iter:long~name~
> = \"model timestep number\" ; iter:standard~name~ = \"timestep\" ;
> int64 time(time) ; time:long~name~ = \"Time\" ; time:standard~name~ =
> \"time\" ; time:axis = \"T\" ; time:units = \"seconds since
> 2012-12-01\" ; time:calendar = \"proleptic~gregorian~\" ; float XC(XC)
> ; XC:~FillValue~ = NaNf ; XC:coordinate = \"YC XC\" ; XC:units =
> \"degrees~east~\" ; XC:standard~name~ = \"longitude\" ; XC:long~name~
> = \"longitude\" ; XC:axis = \"X\" ; float YG(YG) ; YG:~FillValue~ =
> NaNf ; YG:units = \"degrees~north~\" ; YG:long~name~ = \"latitude\" ;
> YG:standard~name~ = \"latitude~atflocation~\" ; YG:axis = \"Y\" ;
> YG:c~gridaxisshift~ = -0.5 ; float rAs(YG, XC) ; rAs:~FillValue~ =
> NaNf ; rAs:units = \"m2\" ; rAs:long~name~ = \"cell area\" ;
> rAs:standard~name~ = \"cell~areaatvlocation~\" ; float dxG(YG, XC) ;
> dxG:~FillValue~ = NaNf ; dxG:coordinate = \"YG XC\" ; dxG:units =
> \"m\" ; dxG:standard~name~ = \"cell~xsizeatvlocation~\" ;
> dxG:long~name~ = \"cell x size\" ; float dyC(YG, XC) ; dyC:~FillValue~
> = NaNf ; dyC:coordinate = \"YG XC\" ; dyC:units = \"m\" ;
> dyC:standard~name~ = \"cell~ysizeatvlocation~\" ; dyC:long~name~ =
> \"cell y size\" ; float SIvheff(time, YG, XC) ; SIvheff:~FillValue~ =
> NaNf ; SIvheff:units = \"m^2^/s\" ; SIvheff:long~name~ = \"Meridional
> Transport of eff ice thickn (centered)\" ; SIvheff:standard~name~ =
> \"SIvheff\" ; SIvheff:mate = \"SIuheff\" ; SIvheff:coordinates = \"dxG
> rAs iter dyC\" ; }

### Sea ice meridional velocity

> netcdf bsose~i1222013to20171daySeaIceVvel~ { dimensions: time = 1826 ;
> XC = 2160 ; YG = 588 ; variables: int64 iter(time) ; iter:long~name~ =
> \"model timestep number\" ; iter:standard~name~ = \"timestep\" ; int64
> time(time) ; time:long~name~ = \"Time\" ; time:standard~name~ =
> \"time\" ; time:axis = \"T\" ; time:units = \"seconds since
> 2012-12-01\" ; time:calendar = \"proleptic~gregorian~\" ; float XC(XC)
> ; XC:~FillValue~ = NaNf ; XC:coordinate = \"YC XC\" ; XC:units =
> \"degrees~east~\" ; XC:standard~name~ = \"longitude\" ; XC:long~name~
> = \"longitude\" ; XC:axis = \"X\" ; float YG(YG) ; YG:~FillValue~ =
> NaNf ; YG:units = \"degrees~north~\" ; YG:long~name~ = \"latitude\" ;
> YG:standard~name~ = \"latitude~atflocation~\" ; YG:axis = \"Y\" ;
> YG:c~gridaxisshift~ = -0.5 ; float rAs(YG, XC) ; rAs:~FillValue~ =
> NaNf ; rAs:units = \"m2\" ; rAs:long~name~ = \"cell area\" ;
> rAs:standard~name~ = \"cell~areaatvlocation~\" ; float dxG(YG, XC) ;
> dxG:~FillValue~ = NaNf ; dxG:coordinate = \"YG XC\" ; dxG:units =
> \"m\" ; dxG:standard~name~ = \"cell~xsizeatvlocation~\" ;
> dxG:long~name~ = \"cell x size\" ; float dyC(YG, XC) ; dyC:~FillValue~
> = NaNf ; dyC:coordinate = \"YG XC\" ; dyC:units = \"m\" ;
> dyC:standard~name~ = \"cell~ysizeatvlocation~\" ; dyC:long~name~ =
> \"cell y size\" ; float SIvice(time, YG, XC) ; SIvice:~FillValue~ =
> NaNf ; SIvice:units = \"m/s\" ; SIvice:long~name~ = \"SEAICE merid.
> ice velocity, \>0 from South to North\" ; SIvice:standard~name~ =
> \"SIvice\" ; SIvice:mate = \"SIuice\" ; SIvice:coordinates = \"dxG rAs
> iter dyC\" ; }

### Sea ice effective snow thickness

> netcdf bsose~i1222013to20171daySnowHeff~ { dimensions: time = 1826 ;
> XC = 2160 ; YC = 588 ; variables: int64 iter(time) ; iter:long~name~ =
> \"model timestep number\" ; iter:standard~name~ = \"timestep\" ; int64
> time(time) ; time:long~name~ = \"Time\" ; time:standard~name~ =
> \"time\" ; time:axis = \"T\" ; time:units = \"seconds since
> 2012-12-01\" ; time:calendar = \"proleptic~gregorian~\" ; float XC(XC)
> ; XC:~FillValue~ = NaNf ; XC:coordinate = \"YC XC\" ; XC:units =
> \"degrees~east~\" ; XC:standard~name~ = \"longitude\" ; XC:long~name~
> = \"longitude\" ; XC:axis = \"X\" ; float YC(YC) ; YC:~FillValue~ =
> NaNf ; YC:coordinate = \"YC XC\" ; YC:units = \"degrees~north~\" ;
> YC:standard~name~ = \"latitude\" ; YC:long~name~ = \"latitude\" ;
> YC:axis = \"Y\" ; float Depth(YC, XC) ; Depth:~FillValue~ = NaNf ;
> Depth:coordinate = \"XC YC\" ; Depth:units = \"m\" ;
> Depth:standard~name~ = \"ocean~depth~\" ; Depth:long~name~ = \"ocean
> depth\" ; float rA(YC, XC) ; rA:~FillValue~ = NaNf ; rA:coordinate =
> \"YC XC\" ; rA:units = \"m2\" ; rA:standard~name~ = \"cell~area~\" ;
> rA:long~name~ = \"cell area\" ; float SIhsnow(time, YC, XC) ;
> SIhsnow:~FillValue~ = NaNf ; SIhsnow:units = \"m\" ;
> SIhsnow:long~name~ = \"SEAICE effective snow thickness\" ;
> SIhsnow:standard~name~ = \"SIhsnow\" ; SIhsnow:coordinates = \"Depth
> rA iter\" ; }

### Sea ice snow precipitation over sea ice

> netcdf bsose~i1222013to20171daySnowPrecip~ { dimensions: time = 1826 ;
> XC = 2160 ; YC = 588 ; variables: int64 iter(time) ; iter:long~name~ =
> \"model timestep number\" ; iter:standard~name~ = \"timestep\" ; int64
> time(time) ; time:long~name~ = \"Time\" ; time:standard~name~ =
> \"time\" ; time:axis = \"T\" ; time:units = \"seconds since
> 2012-12-01\" ; time:calendar = \"proleptic~gregorian~\" ; float XC(XC)
> ; XC:~FillValue~ = NaNf ; XC:coordinate = \"YC XC\" ; XC:units =
> \"degrees~east~\" ; XC:standard~name~ = \"longitude\" ; XC:long~name~
> = \"longitude\" ; XC:axis = \"X\" ; float YC(YC) ; YC:~FillValue~ =
> NaNf ; YC:coordinate = \"YC XC\" ; YC:units = \"degrees~north~\" ;
> YC:standard~name~ = \"latitude\" ; YC:long~name~ = \"latitude\" ;
> YC:axis = \"Y\" ; float Depth(YC, XC) ; Depth:~FillValue~ = NaNf ;
> Depth:coordinate = \"XC YC\" ; Depth:units = \"m\" ;
> Depth:standard~name~ = \"ocean~depth~\" ; Depth:long~name~ = \"ocean
> depth\" ; float rA(YC, XC) ; rA:~FillValue~ = NaNf ; rA:coordinate =
> \"YC XC\" ; rA:units = \"m2\" ; rA:standard~name~ = \"cell~area~\" ;
> rA:long~name~ = \"cell area\" ; float SIsnPrcp(time, YC, XC) ;
> SIsnPrcp:~FillValue~ = NaNf ; SIsnPrcp:units = \"kg/m^2^/s\" ;
> SIsnPrcp:long~name~ = \"Snow precip. (+=dw) over Sea-Ice (area
> weighted)\" ; SIsnPrcp:standard~name~ = \"SIsnPrcp\" ;
> SIsnPrcp:coordinates = \"Depth rA iter\" ; }

### Zonal wind stress at surface

> netcdf bsose~i1222013to20171dayoceTAUX~ { dimensions: time = 1826 ; YC
> = 588 ; XG = 2160 ; variables: int64 iter(time) ; iter:long~name~ =
> \"model timestep number\" ; iter:standard~name~ = \"timestep\" ; int64
> time(time) ; time:long~name~ = \"Time\" ; time:standard~name~ =
> \"time\" ; time:axis = \"T\" ; time:units = \"seconds since
> 2012-12-01\" ; time:calendar = \"proleptic~gregorian~\" ; float YC(YC)
> ; YC:~FillValue~ = NaNf ; YC:coordinate = \"YC XC\" ; YC:units =
> \"degrees~north~\" ; YC:standard~name~ = \"latitude\" ; YC:long~name~
> = \"latitude\" ; YC:axis = \"Y\" ; float XG(XG) ; XG:~FillValue~ =
> NaNf ; XG:coordinate = \"YG XG\" ; XG:units = \"degrees~east~\" ;
> XG:standard~name~ = \"longitude~atflocation~\" ; XG:long~name~ =
> \"longitude\" ; XG:axis = \"X\" ; XG:c~gridaxisshift~ = -0.5 ; float
> dxC(YC, XG) ; dxC:~FillValue~ = NaNf ; dxC:coordinate = \"YC XG\" ;
> dxC:units = \"m\" ; dxC:standard~name~ = \"cell~xsizeatulocation~\" ;
> dxC:long~name~ = \"cell x size\" ; float rAw(YC, XG) ; rAw:~FillValue~
> = NaNf ; rAw:coordinate = \"YG XC\" ; rAw:units = \"m2\" ;
> rAw:standard~name~ = \"cell~areaatulocation~\" ; rAw:long~name~ =
> \"cell area\" ; float dyG(YC, XG) ; dyG:~FillValue~ = NaNf ;
> dyG:coordinate = \"YC XG\" ; dyG:units = \"m\" ; dyG:standard~name~ =
> \"cell~ysizeatulocation~\" ; dyG:long~name~ = \"cell y size\" ; float
> oceTAUX(time, YC, XG) ; oceTAUX:~FillValue~ = NaNf ; oceTAUX:units =
> \"N/m^2^\" ; oceTAUX:long~name~ = \"zonal surface wind stress, \>0
> increases uVel\" ; oceTAUX:standard~name~ = \"oceTAUX\" ; oceTAUX:mate
> = \"oceTAUY\" ; oceTAUX:coordinates = \"rAw dyG dxC iter\" ; }

### Meridional wind stress at surface

> netcdf bsose~i1222013to20171dayoceTAUY~ { dimensions: time = 1826 ; XC
> = 2160 ; YG = 588 ; variables: int64 iter(time) ; iter:long~name~ =
> \"model timestep number\" ; iter:standard~name~ = \"timestep\" ; int64
> time(time) ; time:long~name~ = \"Time\" ; time:standard~name~ =
> \"time\" ; time:axis = \"T\" ; time:units = \"seconds since
> 2012-12-01\" ; time:calendar = \"proleptic~gregorian~\" ; float XC(XC)
> ; XC:~FillValue~ = NaNf ; XC:coordinate = \"YC XC\" ; XC:units =
> \"degrees~east~\" ; XC:standard~name~ = \"longitude\" ; XC:long~name~
> = \"longitude\" ; XC:axis = \"X\" ; float YG(YG) ; YG:~FillValue~ =
> NaNf ; YG:units = \"degrees~north~\" ; YG:long~name~ = \"latitude\" ;
> YG:standard~name~ = \"latitude~atflocation~\" ; YG:axis = \"Y\" ;
> YG:c~gridaxisshift~ = -0.5 ; float rAs(YG, XC) ; rAs:~FillValue~ =
> NaNf ; rAs:units = \"m2\" ; rAs:long~name~ = \"cell area\" ;
> rAs:standard~name~ = \"cell~areaatvlocation~\" ; float dxG(YG, XC) ;
> dxG:~FillValue~ = NaNf ; dxG:coordinate = \"YG XC\" ; dxG:units =
> \"m\" ; dxG:standard~name~ = \"cell~xsizeatvlocation~\" ;
> dxG:long~name~ = \"cell x size\" ; float dyC(YG, XC) ; dyC:~FillValue~
> = NaNf ; dyC:coordinate = \"YG XC\" ; dyC:units = \"m\" ;
> dyC:standard~name~ = \"cell~ysizeatvlocation~\" ; dyC:long~name~ =
> \"cell y size\" ; float oceTAUY(time, YG, XC) ; oceTAUY:~FillValue~ =
> NaNf ; oceTAUY:units = \"N/m^2^\" ; oceTAUY:long~name~ = \"meridional
> surf. wind stress, \>0 increases vVel\" ; oceTAUY:standard~name~ =
> \"oceTAUY\" ; oceTAUY:mate = \"oceTAUX\" ; oceTAUY:coordinates = \"dxG
> rAs iter dyC\" ; }

### Total salt flux

> netcdf bsose~i1222013to20171daysurfSflx~ { dimensions: time = 1826 ;
> XC = 2160 ; YC = 588 ; variables: int64 iter(time) ; iter:long~name~ =
> \"model timestep number\" ; iter:standard~name~ = \"timestep\" ; int64
> time(time) ; time:long~name~ = \"Time\" ; time:standard~name~ =
> \"time\" ; time:axis = \"T\" ; time:units = \"seconds since
> 2012-12-01\" ; time:calendar = \"proleptic~gregorian~\" ; float XC(XC)
> ; XC:~FillValue~ = NaNf ; XC:coordinate = \"YC XC\" ; XC:units =
> \"degrees~east~\" ; XC:standard~name~ = \"longitude\" ; XC:long~name~
> = \"longitude\" ; XC:axis = \"X\" ; float YC(YC) ; YC:~FillValue~ =
> NaNf ; YC:coordinate = \"YC XC\" ; YC:units = \"degrees~north~\" ;
> YC:standard~name~ = \"latitude\" ; YC:long~name~ = \"latitude\" ;
> YC:axis = \"Y\" ; float Depth(YC, XC) ; Depth:~FillValue~ = NaNf ;
> Depth:coordinate = \"XC YC\" ; Depth:units = \"m\" ;
> Depth:standard~name~ = \"ocean~depth~\" ; Depth:long~name~ = \"ocean
> depth\" ; float rA(YC, XC) ; rA:~FillValue~ = NaNf ; rA:coordinate =
> \"YC XC\" ; rA:units = \"m2\" ; rA:standard~name~ = \"cell~area~\" ;
> rA:long~name~ = \"cell area\" ; float SFLUX(time, YC, XC) ;
> SFLUX:~FillValue~ = NaNf ; SFLUX:units = \"g/m^2^/s\" ;
> SFLUX:long~name~ = \"total salt flux (match salt-content variations),
> \>0 increases salt\" ; SFLUX:standard~name~ = \"SFLUX\" ;
> SFLUX:coordinates = \"Depth rA iter\" ; }

### Total heat flux

> netcdf bsose~i1222013to20171daysurfTflx~ { dimensions: time = 1826 ;
> XC = 2160 ; YC = 588 ; variables: int64 iter(time) ; iter:long~name~ =
> \"model timestep number\" ; iter:standard~name~ = \"timestep\" ; int64
> time(time) ; time:long~name~ = \"Time\" ; time:standard~name~ =
> \"time\" ; time:axis = \"T\" ; time:units = \"seconds since
> 2012-12-01\" ; time:calendar = \"proleptic~gregorian~\" ; float XC(XC)
> ; XC:~FillValue~ = NaNf ; XC:coordinate = \"YC XC\" ; XC:units =
> \"degrees~east~\" ; XC:standard~name~ = \"longitude\" ; XC:long~name~
> = \"longitude\" ; XC:axis = \"X\" ; float YC(YC) ; YC:~FillValue~ =
> NaNf ; YC:coordinate = \"YC XC\" ; YC:units = \"degrees~north~\" ;
> YC:standard~name~ = \"latitude\" ; YC:long~name~ = \"latitude\" ;
> YC:axis = \"Y\" ; float Depth(YC, XC) ; Depth:~FillValue~ = NaNf ;
> Depth:coordinate = \"XC YC\" ; Depth:units = \"m\" ;
> Depth:standard~name~ = \"ocean~depth~\" ; Depth:long~name~ = \"ocean
> depth\" ; float rA(YC, XC) ; rA:~FillValue~ = NaNf ; rA:coordinate =
> \"YC XC\" ; rA:units = \"m2\" ; rA:standard~name~ = \"cell~area~\" ;
> rA:long~name~ = \"cell area\" ; float TFLUX(time, YC, XC) ;
> TFLUX:~FillValue~ = NaNf ; TFLUX:units = \"W/m^2^\" ; TFLUX:long~name~
> = \"total heat flux (match heat-content variations), \>0 increases
> theta\" ; TFLUX:standard~name~ = \"TFLUX\" ; TFLUX:coordinates =
> \"Depth rA iter\" ; }

### Atmospheric surface pressure

> netcdf bsose~i1222013to2017dailyEXFpress~ { dimensions: time = 1826 ;
> XC = 2160 ; YC = 588 ; variables: int64 iter(time) ; iter:long~name~ =
> \"model timestep number\" ; iter:standard~name~ = \"timestep\" ; int64
> time(time) ; time:long~name~ = \"Time\" ; time:standard~name~ =
> \"time\" ; time:axis = \"T\" ; time:units = \"seconds since
> 2012-12-01\" ; time:calendar = \"proleptic~gregorian~\" ; float XC(XC)
> ; XC:~FillValue~ = NaNf ; XC:coordinate = \"YC XC\" ; XC:units =
> \"degrees~east~\" ; XC:standard~name~ = \"longitude\" ; XC:long~name~
> = \"longitude\" ; XC:axis = \"X\" ; float YC(YC) ; YC:~FillValue~ =
> NaNf ; YC:coordinate = \"YC XC\" ; YC:units = \"degrees~north~\" ;
> YC:standard~name~ = \"latitude\" ; YC:long~name~ = \"latitude\" ;
> YC:axis = \"Y\" ; float Depth(YC, XC) ; Depth:~FillValue~ = NaNf ;
> Depth:coordinate = \"XC YC\" ; Depth:units = \"m\" ;
> Depth:standard~name~ = \"ocean~depth~\" ; Depth:long~name~ = \"ocean
> depth\" ; float rA(YC, XC) ; rA:~FillValue~ = NaNf ; rA:coordinate =
> \"YC XC\" ; rA:units = \"m2\" ; rA:standard~name~ = \"cell~area~\" ;
> rA:long~name~ = \"cell area\" ; float EXFpress(time, YC, XC) ;
> EXFpress:~FillValue~ = NaNf ; EXFpress:units = \"N/m^2^\" ;
> EXFpress:long~name~ = \"atmospheric pressure field\" ;
> EXFpress:standard~name~ = \"EXFpress\" ; EXFpress:coordinates =
> \"Depth rA iter\" ; }

## BRAN

[CSIRO maintained
dataset](https://geonetwork.nci.org.au/geonetwork/srv/eng/catalog.search#/metadata/f9372_7752_2015_3718),
but available on gadi via the [BRAN
data](https://dapds00.nci.org.au/thredds/catalog/gb6/BRAN/BRAN2020/catalog.html)
repository.

## CAWCR

The [CICE GitHub online setup
guide](https://github.com/CICE-Consortium/CICE/wiki/CICE-Input-Data#descriptions-of-tar-files)
explains that ocean surface gravity wave information is necessary for
floe size distribution (FSD). Explanation is given that illustrates from
just one day with the only variable present from [WaveWatch
III,](https://github.com/NOAA-EMC/WW3) `efreq` ($m^2/s$) the importance
of including wave information in sea ice modelling. The [CICE
Documentation](https://cice-consortium-cice.readthedocs.io/en/main/#)
states the basic problem in sea ice modelling, the one CICE aims to
model, is ice thickness distribution (ITD). However, CICE can optionally
attempt to model FSD. It does this via a probability distribution
function that charcacterises the variability of the horizontal floe size
through vertical and lateral growth, new freezing, breaking waves and
welding of floes. This work is based on
\\[@horvatPrognosticModelSeaice2015]

-   Available to download from a
    [THREDDS](https://data-cbr.csiro.au/thredds/catalog/catch_all/CMAR_CAWCR-Wave_archive/CAWCR_Wave_Hindcast_aggregate/gridded/catalog.html)

# Pre-processing Forcing and Initial Conditions

## Grids

There are two approaches to providing a stand-alone grid to CICE6: using
someone elses or creating one from scratch. Methods accounting both for
my research with CICE6 stand-alone are given below.

### Pilfering ACCESS-OM2

[ACCESS-OM2](http://cosima.org.au/index.php/models/access-om2/) has
$1$-degree, $1/4$-degree, and $1/10$-degree grid files that are
Arakawa-B grid types. A good description of the varying types of grids
is given in [@nGridsNumericalWeather2013]. [CICE6 allows for both
Arakawa-B and
Arakawa-C](https://cice-consortium-cice.readthedocs.io/en/main/user_guide/ug_implementation.html#grid-boundary-conditions-and-masks)
grids.

Regridding of [*ERA5*]{.spurious-link target="ERA5"} and
[*BRAN*]{.spurious-link target="BRAN"} datasets is required to map the
essential forcing variables onto the same GRID that will be used for
CICE6 studies herein.
[afim.py](file:///Users/dpath2o/PHD/src/python/afim.py) was written for
this purpose. Currently it uses [*CDO*]{.spurious-link
target="file/../admin/computing.org:Climate Data Operator (CDO)"} but
there is a thought here to make it optionally use [NetCDF Operators
(NCO)](../admin/computing.org).

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
