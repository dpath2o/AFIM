import os
import sys
import afim
import cmocean
import numpy  as np
import xarray as xr
from datetime import datetime

IC_AOM2_OCN      = xr.open_dataset("/g/data/cj50/access-om2/raw-output/access-om2-025/025deg_jra55v13_ryf9091_gmredi6/output052/ocean/ocean_month.nc")
D_graph_base     = "/g/data/jk72/da1339/GRAPHICAL/initial_conditions/2004_12/"
IC_AOM2_OCN_vars = {"sea_level"             : {"long_name" : "effective sea level (eta_t + patm/(rho0*g)) on T cells",
                                               "lon_name"  : "xt_ocean",
                                               "lat_name"  : "yt_ocean",
                                               "units"     : "m",
                                               "vmin"      : -3,
                                               "vmax"      : 0,
                                               "cmap"      : cmocean.cm.amp}, 
                    "eta_t"                 : {"long_name" : "surface height on T cells",
                                               "lon_name"  : "xt_ocean",
                                               "lat_name"  : "yt_ocean",
                                               "units"     : "m",
                                               "vmin"      : -3,
                                               "vmax"      : 0,
                                               "cmap"      : cmocean.cm.amp}, 
                    "sea_levelsq"           : {"long_name" : "square of effective sea level (eta_t + patm/(rho0*g)) on T cells",
                                               "lon_name"  : "xt_ocean",
                                               "lat_name"  : "yt_ocean",
                                               "units"     : "m^2",
                                               "vmin"      : 0,
                                               "vmax"      : 6,
                                               "cmap"      : cmocean.cm.amp}, 
                    "mld"                   : {"long_name" : "mixed layer depth determined by density criteria",
                                               "lon_name"  : "xt_ocean",
                                               "lat_name"  : "yt_ocean",
                                               "units"     : "m",
                                               "vmin"      : 0,
                                               "vmax"      : 200,
                                               "cmap"      : cmocean.cm.deep}, 
                    "surface_temp"          : {"long_name" : "Conservative temperature",
                                               "lon_name"  : "xt_ocean",
                                               "lat_name"  : "yt_ocean",
                                               "units"     : "C",
                                               "vmin"      : -3,
                                               "vmax"      : 10,
                                               "cmap"      : cmocean.cm.thermal}, 
                    "surface_salt"          : {"long_name" : "Practical Salinity",
                                               "lon_name"  : "xt_ocean",
                                               "lat_name"  : "yt_ocean",
                                               "units"     : "psu",
                                               "vmin"      : 34,
                                               "vmax"      : 36,
                                               "cmap"      : cmocean.cm.haline}, 
                    "evap"                  : {"long_name" : "mass flux from evaporation/condensation (>0 enters ocean)",
                                               "lon_name"  : "xt_ocean",
                                               "lat_name"  : "yt_ocean",
                                               "units"     : "(kg/m^3)*(m/sec)",
                                               "vmin"      : -0.0002,
                                               "vmax"      : 0,
                                               "cmap"      : cmocean.cm.amp_r}, 
                    "melt"                  : {"long_name" : "water flux transferred with sea ice form/melt (>0 enters ocean)",
                                               "lon_name"  : "xt_ocean",
                                               "lat_name"  : "yt_ocean",
                                               "units"     : "(kg/m^3)*(m/sec)",
                                               "vmin"      : 0,
                                               "vmax"      : 0.001,
                                               "cmap"      : cmocean.cm.amp}, 
                    "sfc_salt_flux_restore" : {"long_name" : "sfc_salt_flux_restore: flux from restoring term",
                                               "lon_name"  : "xt_ocean",
                                               "lat_name"  : "yt_ocean",
                                               "units"     : "kg/(m^2*sec)",
                                               "vmin"      : -6e-7,
                                               "vmax"      : 6e-7,
                                               "cmap"      : cmocean.cm.balance},
                    "sfc_salt_flux_ice"     : {"long_name" : "sfc_salt_flux_ice",
                                               "lon_name"  : "xt_ocean",
                                               "lat_name"  : "yt_ocean",
                                               "units"     : "kg/(m^2*sec)",
                                               "vmin"      : -8e-6,
                                               "vmax"      : 8e-6,
                                              "cmap"      : cmocean.cm.balance},
                    "sfc_salt_flux_coupler" : {"long_name" : "sfc_salt_flux_coupler: flux from the coupler",
                                               "lon_name"  : "xt_ocean",
                                               "lat_name"  : "yt_ocean",
                                               "units"     : "kg/(m^2*sec)",
                                               "vmin"      : -8e-6,
                                               "vmax"      : 8e-6,
                                              "cmap"      : cmocean.cm.balance}, 
                    "net_sfc_heating"       : {"long_name" : "surface ocean heat flux coming through coupler and mass transfer",
                                               "lon_name"  : "xt_ocean",
                                               "lat_name"  : "yt_ocean",
                                               "units"     : "W/m^2",
                                               "vmin"      : 10,
                                               "vmax"      : 1000,
                                               "cmap"      : cmocean.cm.amp},
                    "tau_x"                 : {"long_name" : "i-directed wind stress forcing u-velocity",
                                               "lon_name"  : "xu_ocean",
                                               "lat_name"  : "yu_ocean",
                                               "units"     : "N/m^2",
                                               "vmin"      : -0.4,
                                               "vmax"      :  0.4,
                                               "cmap"      : cmocean.cm.balance},
                    "tau_y"                 : {"long_name" : "j-directed wind stress forcing v-velocity",
                                               "lon_name"  : "xu_ocean",
                                               "lat_name"  : "yu_ocean",
                                               "units"     : "N/m^2",
                                               "vmin"      : -0.4,
                                               "vmax"      :  0.4,
                                               "cmap"      : cmocean.cm.balance},
                    "bmf_u"                 : {"long_name" : "Bottom u-stress via bottom drag",
                                               "lon_name"  : "xu_ocean",
                                               "lat_name"  : "yu_ocean",
                                               "units"     : "N/m^2",
                                               "vmin"      : -0.5,
                                               "vmax"      :  0.5,
                                               "cmap"      : cmocean.cm.balance},
                    "bmf_v"                 : {"long_name" : "Bottom v-stress via bottom drag",
                                               "lon_name"  : "xu_ocean",
                                               "lat_name"  : "yu_ocean",
                                               "units"     : "N/m^2",
                                               "vmin"      : -0.5,
                                               "vmax"      :  0.5,
                                               "cmap"      : cmocean.cm.balance},
                    "tx_trans_int_z"        : {"long_name" : "T-cell i-mass transport vertically summed",
                                               "lon_name"  : "xu_ocean",
                                               "lat_name"  : "yt_ocean",
                                               "units"     : "kg/s",
                                               "vmin"      : -3e10,
                                               "vmax"      :  3e10,
                                               "cmap"      : cmocean.cm.balance},
                    "ty_trans_int_z"        : {"long_name" : "T-cell j-mass transport vertically summed",
                                               "lon_name"  : "xt_ocean",
                                               "lat_name"  : "yu_ocean",
                                               "units"     : "kg/s",
                                               "vmin"      : -3e10,
                                               "vmax"      :  3e10,
                                               "cmap"      : cmocean.cm.balance},
                    "aredi"                 : {"long_name" : "neutral diffusivity at k=1",
                                               "lon_name"  : "xt_ocean",
                                               "lat_name"  : "yt_ocean",
                                               "units"     : "m^2/sec",
                                               "vmin"      : 0,
                                               "vmax"      : 200,
                                               "cmap"      : cmocean.cm.amp},
                    "agm"                   : {"long_name" : "GM diffusivity at surface",
                                               "lon_name"  : "xt_ocean",
                                               "lat_name"  : "yt_ocean",
                                               "units"     : "m^2/sec",
                                               "vmin"      : 0,
                                               "vmax"      : 200,
                                               "cmap"      : cmocean.cm.amp},
                    "frazil_3d"             : {"long_name" : "ocn frazil heat flux over time step",
                                               "lon_name"  : "xt_ocean",
                                               "lat_name"  : "yt_ocean",
                                               "units"     : "W/m^2",
                                               "vmin"      : -200,
                                               "vmax"      :  200,
                                               "cmap"      : cmocean.cm.balance},
                    "swflx"                 : {"long_name" : "shortwave flux into ocean (>0 heats ocean)",
                                               "lon_name"  : "xt_ocean",
                                               "lat_name"  : "yt_ocean",
                                               "units"     : "W/m^2",
                                               "vmin"      : 0,
                                               "vmax"      : 350,
                                               "cmap"      : cmocean.cm.amp},
                    "sw_heat"               : {"long_name" : "penetrative shortwave heating",
                                               "lon_name"  : "xt_ocean",
                                               "lat_name"  : "yt_ocean",
                                               "units"     : "W/m^2",
                                               "vmin"      : -150,
                                               "vmax"      : 0,
                                               "cmap"      : cmocean.cm.amp_r}}
for var_name, attributes in IC_AOM2_OCN_vars.items():
    if var_name=="frazil_3d" or var_name=="sw_heat": 
        da = IC_AOM2_OCN[var_name].isel(time=11,st_ocean=0)
    elif var_name=="surface_temp":
        da = da = IC_AOM2_OCN[var_name].isel(time=11)-273.15
    else:
        da = IC_AOM2_OCN[var_name].isel(time=11)
    t_str    = da.time.values.astype('datetime64[M]').astype(str)
    tit_str  = f"{attributes['long_name']}"
    F_name   = f"{t_str}_{var_name}_ACCESS-OM2-025_SH.png"
    cbar_lab = f"{attributes['units']}"
    if not os.path.exists(os.path.join(D_graph_base,F_name)):
        print(f"attempting to plot {var_name}")
        afim.plot_cartopy_pcolormesh(da=da, lon_name=attributes['lon_name'],
                                     lat_name=attributes['lat_name'],t_str=t_str,cmap=attributes['cmap'], 
                                     vmin=attributes['vmin'], vmax=attributes['vmax'],cbar_label=cbar_lab,
                                     title=tit_str,D_save=D_graph_base,F_save=F_name)

IC_AOM2_I2O      = xr.open_dataset("/g/data/cj50/access-om2/raw-output/access-om2-025/025deg_jra55v13_ryf9091_gmredi6/output052/ocean/o2i.nc", decode_times=False)
D_graph_base     = "/g/data/jk72/da1339/GRAPHICAL/initial_conditions/2004_12/"
IC_AOM2_I2O_vars = {"sst_i"    : {"long_name" : "Sea Surface Temperature",
                                  "lon_name"  : "nx",
                                   "lat_name"  : "ny",
                                   "units"     : "C",
                                   "vmin"      : -3,
                                   "vmax"      : 10,
                                   "cmap"      : cmocean.cm.thermal}, 
                    "sss_i"    : {"long_name" : "Sea Surface Salinity",
                                  "lon_name"  : "nx",
                                   "lat_name"  : "ny",
                                   "units"     : "psu",
                                   "vmin"      : 34,
                                   "vmax"      : 36,
                                   "cmap"      : cmocean.cm.haline}, 
                    "ssu_i"    : {"long_name" : "Sea Surface Zonal Flow",
                                  "lon_name"  : "nx",
                                   "lat_name"  : "ny",
                                   "units"     : "m/s",
                                   "vmin"      : -1.5,
                                   "vmax"      : 1.5,
                                   "cmap"      : cmocean.cm.balance}, 
                    "ssv_i"    : {"long_name" : "Sea Surface Meridional Flow",
                                  "lon_name"  : "nx",
                                   "lat_name"  : "ny",
                                   "units"     : "m/s",
                                   "vmin"      : -1.5,
                                   "vmax"      : 1.5,
                                   "cmap"      : cmocean.cm.balance}, 
                    "sslx_i"   : {"long_name" : "Sea Surface Zonal Tilt",
                                  "lon_name"  : "nx",
                                   "lat_name"  : "ny",
                                   "units"     : "m",
                                   "vmin"      : -8e-6,
                                   "vmax"      :  8e-6,
                                   "cmap"      : cmocean.cm.balance}, 
                    "ssly_i"   : {"long_name" : "Sea Surface Meridional Tilt",
                                  "lon_name"  : "nx",
                                   "lat_name"  : "ny",
                                   "units"     : "m",
                                   "vmin"      : -8e-6,
                                   "vmax"      :  8e-6,
                                   "cmap"      : cmocean.cm.balance}, 
                    "pfmice_i" : {"long_name" : "Freeze Melt",
                                  "lon_name"  : "nx",
                                   "lat_name"  : "ny",
                                   "units"     : "[UNKNOWN]",
                                   "vmin"      : -1.5e6,
                                   "vmax"      : -5.0e4,
                                   "cmap"      : cmocean.cm.amp_r}}
for var_name, attributes in IC_AOM2_I2O_vars.items():
    if var_name=="sst_i":
        da = IC_AOM2_I2O[var_name].isel(time=0)-273.15
    else:
        da = IC_AOM2_I2O[var_name].isel(time=0)
    tit_str  = f"{attributes['long_name']}"
    F_name   = f"2004-12_{var_name}_ACCESS-OM2-025_I2O_SH.png"
    cbar_lab = f"{attributes['units']}"
    if not os.path.exists(os.path.join(D_graph_base,F_name)):
        print(f"attempting to plot {var_name}")
        afim.plot_cartopy_pcolormesh(da=da, lon_name=attributes['lon_name'],
                                     lat_name=attributes['lat_name'],t_str="2004-12",cmap=attributes['cmap'], 
                                     vmin=attributes['vmin'], vmax=attributes['vmax'],cbar_label=cbar_lab,
                                     title=tit_str,D_save=D_graph_base,F_save=F_name)