# External Systems
## [GitHub AFIM](https://github.com/dpath2o/AFIM)
## Kunanyi
## gadi 
* [Account Created on 30 Mar 22](https://opus.nci.org.au/login.action?os_destination=%2Fdisplay%2FHelp%2F0.%2BWelcome%2Bto%2BGadi)
* [helpful information on GADI usage](https://opus.nci.org.au/display/Help/0.+Welcome+to+Gadi#id-0.WelcometoGadi-FileTransferto/fromGadi)
* [Queues and Limits](https://opus.nci.org.au/display/Help/Queue+Limits)
### Project and Group Memberships
1.  [gb6 (Bluelink)](https://dx.doi.org/10.25914/6009627c7af03)
2.  [hh5 (Climate-LIEF)](https://orcid.org/0000-0002-0164-5105)
3.  [ik11 (COSIMA)](https://orcid.org/0000-0001-5898-7635)
4.  [jk72 (Antarctic systems and future change)](https://orcid.org/0000-0003-1404-4103)
5.  [qv56 (Reference datasets for climate model analysis/forcing)](https://geonetwork.nci.org.au/geonetwork/srv/eng/catalog.search#/metadata/f4756_8918_5098_5424)
6.  [rr7 (Reanalysis reference data)]()
7.  [rt52 (ERA5 Replicated Data: Single and pressure-levels data)](https://geonetwork.nci.org.au/geonetwork/srv/eng/catalog.search#/metadata/f2836_1346_3774_9763)
8.  [ua8 (ARC Centre of Excellence for Climate System Science data)](https://orcid.org/0000-0002-0164-5105)
## CloudStor
## Remote Desktop 
### X2Go [X2Go](https://wiki.x2go.org/doku.php/doc:newtox2go)
## [Tasmanian Partnership for Advancing Computing (TPAC)](https://www.tpac.org.au)
## [National eResearch Colloboration Tools And Resources (NECTAR)](https://nectar.org.au)
### [Nectar Cloud Starter Training](https://www.tpac.org.au/training/)

# High Performance Computing Software
## NetCDF
### [NetCDF Operator Library (NCO)](https://nco.sourceforge.net)
#### Helpful NCO commands:
* Concatenating along the temporal dimension: [link](https://disc.gsfc.nasa.gov/information/howto?title=How%20to%20Concatenate%20the%20Time%20Dimension%20of%20netCDF%20Files%20with%20NCO)
* Change variable name: `ncrename -h -O -v old_variable_name,new_variable_name filename.nc`
* Delete variable “lev” from in.nc: `ncks -C -O -x -v lev in.nc out.nc`
* Delete dimension “lev” from in.nc: `ncwa -a lev in.nc out.nc`
* delete a global attribute: `ncatted -a global_attr_name,d,, infile.nc outfile.nc`
* convert a variable type: `ncap2 -s 'time=float(time)'`
* add a variable over a mapped dimension: `ncap2 -s 'time[$time]=array(54760,30,$time)' infile.nc outfile.nc`
* add an attribute one at a time: `ncatted -a attribute,variable,a,c,"Atrribute Value" infile.nc`
### [Climate Data Operators (CDO)](https://code.mpimet.mpg.de/projects/cdo/)
## Python
### [COSIMA Cookbook](https://cosima-recipes.readthedocs.io/en/latest/index.htmlx)
## Fortran
## General Software Applications
### [Zotero](https://www.zotero.org): reference management

# Hardware
## Primary Laptop
-   Macbook Pro (14-inch,2021)
-   Memory :: 32GB
-   Storage :: 1TB
## Primary Notepad
-   iPad Pro 11-inch