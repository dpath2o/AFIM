# #!/usr/bin/env python

#!/usr/bin/env python
import sys, os, time
import numpy             as np
import pandas            as pd
import xarray            as xr
import xesmf             as xe
import gsw

# =======================
# User controls / paths
# =======================
compute_davg = False   # True: do your Southern Ocean depth-avg workflow
ow_nc        = False
yr0, yrN     = 1993, 2025

# Optional: rotate (east,north) → model (i,j) at E/N points
ROTATE_TO_GRID = False   # <— set True if you want i/j-aligned currents saved

# ORAS location
D_ORAS_base = os.path.join("/","g","data","gv90","da1339","afim_input","ORAS","daily")
D_ORAS_org  = os.path.join(D_ORAS_base,"org")
if compute_davg:
    D_ORAS_tmp = os.path.join(D_ORAS_base,f"d5","tmp")     # matches your prior dmax=5 branch
    D_ORAS_fin = os.path.join(D_ORAS_base,f"d5","Cgrd")
else:
    D_ORAS_tmp = os.path.join(D_ORAS_base,"sfc","tmp")
    D_ORAS_fin = os.path.join(D_ORAS_base,"sfc","Cgrd")
os.makedirs(D_ORAS_tmp, exist_ok=True)
os.makedirs(D_ORAS_fin, exist_ok=True)

# Variable names in ORAS files
pt_name     = "thetao_oras" # potential temperature
s_name      = "so_oras"     # salinity
d_name      = "depth"       # depth coordinate and dimension name
u_name      = "uo_oras"     # eastward current
v_name      = "vo_oras"     # northward current
t_dim_name  = "time"
x_dim_name  = "ni"
y_dim_name  = "nj"
z_dim_name  = "nk"

# Output dims
dim4_out    = [t_dim_name, z_dim_name, y_dim_name, x_dim_name]
dim3_out    = [t_dim_name, y_dim_name, x_dim_name]
dim2_out    = [y_dim_name, x_dim_name]

# Regridding
reG_meth    = "bilinear"
extrap_meth = "inverse_dist"

# Weights (separate for T/E/N)
W_BASE = os.path.join("/","g","data","gv90","da1339","grids","weights")
os.makedirs(W_BASE, exist_ok=True)
W_T = os.path.join(W_BASE, "map_ORAS_to_AOM3-025_T_xesmf-bil.nc")
W_E = os.path.join(W_BASE, "map_ORAS_to_AOM3-025_E_xesmf-bil.nc")
W_N = os.path.join(W_BASE, "map_ORAS_to_AOM3-025_N_xesmf-bil.nc")

# Land mask & grid (your new radians C-grid)
P_kmt = os.path.join("/","g","data","gv90","da1339","grids","ACCESS-OM3-025_kmt.nc")
P_G   = os.path.join("/","g","data","gv90","da1339","grids","ACCESS-OM3-025_Cgrid.nc")

# Depth-avg controls
dmax = 5.0     # top 5 m for velocity
D100 = 100.0   # top 100 m for CT/SSS


# =======================
# Helpers
# =======================
R_EARTH = 6378137.0

def _wrap_dlon_deg(dlon):
    return (dlon + 180.0) % 360.0 - 180.0

def _bearing_east_angle_deg(latA, lonA, latB, lonB):
    """Angle (radians) of A->B vs true east (positive CCW toward north). Inputs in degrees."""
    latm = np.deg2rad(0.5*(latA + latB))
    dlon = np.deg2rad(_wrap_dlon_deg(lonB - lonA))
    dlat = np.deg2rad(latB - latA)
    dx = R_EARTH * np.cos(latm) * dlon
    dy = R_EARTH * dlat
    return np.arctan2(dy, dx)

def guess_lat_lon(ds):
    """Return (lat_name, lon_name) from a dataset."""
    cands_lat = ["lat", "latitude", "nav_lat", "yt_ocean"]
    cands_lon = ["lon", "longitude", "nav_lon", "xt_ocean"]
    lat = next((n for n in cands_lat if n in ds.variables), None)
    lon = next((n for n in cands_lon if n in ds.variables), None)
    if lat is None or lon is None:
        # try coords
        lat = lat or next((n for n in cands_lat if n in ds.coords), None)
        lon = lon or next((n for n in cands_lon if n in ds.coords), None)
    if lat is None or lon is None:
        raise ValueError("Could not find lat/lon in ORAS dataset.")
    return lat, lon

def sample_oras_file_for_year(y):
    """Pick the first half or full-year file for a given year to read lat/lon."""
    if y in [1994, 2000]:
        return os.path.join(D_ORAS_org, f"{y}0101_{y}1231_ORAS_org.nc")
    else:
        return os.path.join(D_ORAS_org, f"{y}0101_{y}0630_ORAS_org.nc")


# =======================
# Load grid / masks
# =======================
kmt      = xr.open_dataset(P_kmt).kmt
kmt_mask = xr.where(kmt == 0, True, False).compute()   # (nj, ni)

G_CICE   = xr.open_dataset(P_G)

# Radians → degrees for xESMF
TLAT2D = np.rad2deg(G_CICE.tlat.values)
TLON2D = np.rad2deg(G_CICE.tlon.values)
ELAT2D = np.rad2deg(G_CICE.elat.values)
ELON2D = np.rad2deg(G_CICE.elon.values)
NLAT2D = np.rad2deg(G_CICE.nlat.values)
NLON2D = np.rad2deg(G_CICE.nlon.values)

# xESMF target grids (2D curvilinear)
tgt_T = xr.Dataset({"lat": (("nj","ni"), TLAT2D),
                    "lon": (("nj","ni"), TLON2D)})
tgt_E = xr.Dataset({"lat": (("nj","ni"), ELAT2D),
                    "lon": (("nj","ni"), ELON2D)})
tgt_N = xr.Dataset({"lat": (("nj","ni"), NLAT2D),
                    "lon": (("nj","ni"), NLON2D)})

# Angles (radians) for optional rotation
# U angle at E points (from your grid; name may be 'angle' or 'uangle')
if "angle" in G_CICE:
    ANGLE_U = G_CICE["angle"].values
elif "uangle" in G_CICE:
    ANGLE_U = G_CICE["uangle"].values
else:
    raise ValueError("Angle at E points not found in grid file (expected 'angle' or 'uangle').")

# Angle at N points: compute from N locations along +i (east) direction
ANG_N = _bearing_east_angle_deg(NLAT2D, NLON2D,
                                NLAT2D, np.roll(NLON2D, -1, axis=1))   # radians

# Also hold T-point coords (degrees) for output convenience
TLAT = TLAT2D
TLON = TLON2D
ELAT = ELAT2D
ELON = ELON2D
NLAT = NLAT2D
NLON = NLON2D


# =======================
# Build regridders once
# =======================
# Use a sample ORAS file from yr0 to get src grid
P_ORAS_sample = sample_oras_file_for_year(yr0)
ORAS_sample   = xr.open_dataset(P_ORAS_sample, engine="netcdf4")
lat_name, lon_name = guess_lat_lon(ORAS_sample)
src_grid = xr.Dataset({"lat": ORAS_sample[lat_name], "lon": ORAS_sample[lon_name]})

reG_T = xe.Regridder(src_grid, tgt_T, method=reG_meth, periodic=True,
                     extrap_method=extrap_meth,
                     filename=W_T, reuse_weights=os.path.exists(W_T))
reG_E = xe.Regridder(src_grid, tgt_E, method=reG_meth, periodic=True,
                     extrap_method=extrap_meth,
                     filename=W_E, reuse_weights=os.path.exists(W_E))
reG_N = xe.Regridder(src_grid, tgt_N, method=reG_meth, periodic=True,
                     extrap_method=extrap_meth,
                     filename=W_N, reuse_weights=os.path.exists(W_N))

print("Regridders ready:", W_T, W_E, W_N)


# =======================
# Year loop
# =======================
for yr in range(yr0, yrN):
    year_t0 = time.time()

    # Determine file segments for the year
    if yr in [1994, 2000]:
        F_ORAS_orgs = [ f"{yr}0101_{yr}1231_ORAS_org.nc" ]
        F_ORAS_tmps = [ f"{yr}0101_{yr}1231_ORAS_tmp.nc" ]
    else:
        F_ORAS_orgs = [f"{yr}0101_{yr}0630_ORAS_org.nc", f"{yr}0701_{yr}1231_ORAS_org.nc"]
        F_ORAS_tmps = [f"{yr}0101_{yr}0630_ORAS_tmp.nc", f"{yr}0701_{yr}1231_ORAS_tmp.nc"]

    for i, F_ORAS_org in enumerate(F_ORAS_orgs):
        P_ORAS_org    = os.path.join(D_ORAS_org, F_ORAS_org)
        P_ORAS_tmp    = os.path.join(D_ORAS_tmp, F_ORAS_tmps[i])
        P_ORAS_fin    = os.path.join(D_ORAS_fin, f"ORAS_{yr:04d}.nc")
        P_davg_wldcrd = os.path.join(D_ORAS_tmp, f"{yr:04d}*_tmp.nc")

        if os.path.exists(P_ORAS_fin) and not ow_nc:
            print(f"file already exists and over-writing disabled; skipping {P_ORAS_fin}")
            continue
        else:
            print(f"working on creating {P_ORAS_fin}")

        if os.path.exists(P_ORAS_tmp) and not ow_nc:
            print(f"file already exists and over-writing disabled; skipping {P_ORAS_tmp}")
            continue
        else:
            print(f"working on creating {P_ORAS_tmp}")

        t0 = time.time()
        ORAS_mos = xr.open_dataset(P_ORAS_org, engine="netcdf4")
        PT       = ORAS_mos[pt_name]   # (time, depth, lat, lon)
        S        = ORAS_mos[s_name]
        d        = ORAS_mos[d_name]
        U        = ORAS_mos[u_name]
        V        = ORAS_mos[v_name]
        print(f"loading data took {time.time()-t0:.1f} s")

        # Regrid
        print("regridding PT (→ T)"); t1 = time.time()
        PT_reG = reG_T(PT)
        print(f"\t{time.time()-t1:.1f} s")

        print("regridding S  (→ T)"); t1 = time.time()
        S_reG  = reG_T(S)
        print(f"\t{time.time()-t1:.1f} s")

        print("regridding U eastward (→ E)"); t1 = time.time()
        U_E    = reG_E(U)
        print(f"\t{time.time()-t1:.1f} s")

        print("regridding V northward (→ N)"); t1 = time.time()
        V_N    = reG_N(V)
        print(f"\t{time.time()-t1:.1f} s")

        # Optional rotation to grid axes
        if ROTATE_TO_GRID:
            print("ROTATE_TO_GRID=True → rotating (east,north) to (i,j)")
            # Need the companion component at the same point
            V_E = reG_E(V)   # V at E locations
            U_N = reG_N(U)   # U at N locations

            # u_i at E:  u_i =  u_east*cos(theta_i) + v_north*sin(theta_i)
            U_grid = U_E * np.cos(ANGLE_U) + V_E * np.sin(ANGLE_U)

            # v_j at N:  v_j = -u_east*sin(theta_i_at_N) + v_north*cos(theta_i_at_N)
            V_grid = -U_N * np.sin(ANG_N) + V_N * np.cos(ANG_N)
        else:
            # keep geographic components on their native E/N points
            U_grid = U_E
            V_grid = V_N

        # Convert PT→CT (uses regridded S)
        print("convert potential temperature to Conservative Temperature (CT)")
        t1 = time.time()
        CT_reG = gsw.conversions.CT_from_pt(S_reG, PT_reG)
        print(f"\t{time.time()-t1:.1f} s")

        # ===== Depth treatment =====
        if compute_davg:
            print("Depth-averaging branch")
            t1 = time.time()

            # Southern Ocean slice (your indices); adjust if needed
            so_j_stop = 233  # keep as in your script
            CT_lt_slc = CT_reG.isel(nj=slice(0, so_j_stop))
            S_lt_slc  = S_reG .isel(nj=slice(0, so_j_stop))
            U_lt_slc  = U_grid.isel(nj=slice(0, so_j_stop))
            V_lt_slc  = V_grid.isel(nj=slice(0, so_j_stop))

            # Depth weights/top ranges
            d_vals  = d.values
            d_diff  = np.diff(d_vals, append=d_vals[-1])
            d_wghts = xr.DataArray(d_diff, dims=[d_name], coords={d_name: d_vals})

            # Top 100 m for CT/SSS
            d100_mask  = d_vals <= D100
            w_d100     = d_wghts.sel({d_name: d_vals[d100_mask]})
            ct_d100    = CT_lt_slc.sel({d_name: d_vals[d100_mask]})
            s_d100     = S_lt_slc .sel({d_name: d_vals[d100_mask]})
            sst_lt_slc = (ct_d100 * w_d100).sum(dim=d_name) / w_d100.sum()
            sss_lt_slc = (s_d100  * w_d100).sum(dim=d_name) / w_d100.sum()

            # Top 5 m for velocity
            d5_mask = d_vals <= dmax
            w_d5    = d_wghts.sel({d_name: d_vals[d5_mask]})
            u_d5    = U_lt_slc.sel({d_name: d_vals[d5_mask]})
            v_d5    = V_lt_slc.sel({d_name: d_vals[d5_mask]})
            u_lt_slc = (u_d5 * w_d5).sum(dim=d_name) / w_d5.sum()
            v_lt_slc = (v_d5 * w_d5).sum(dim=d_name) / w_d5.sum()

            # Build region datasets with correct coords
            TLON_lt   = TLON[:so_j_stop, :]
            TLAT_lt   = TLAT[:so_j_stop, :]
            sst_lt_ds = xr.DataArray(sst_lt_slc, dims=dim3_out,
                                     coords={"TLON": (dim2_out, TLON_lt),
                                             "TLAT": (dim2_out, TLAT_lt)})

            sss_lt_ds = xr.DataArray(sss_lt_slc, dims=dim3_out,
                                     coords={"TLON": (dim2_out, TLON_lt),
                                             "TLAT": (dim2_out, TLAT_lt)})

            u_lt_ds   = xr.DataArray(u_lt_slc,   dims=dim3_out,
                                     coords={"ELON": (dim2_out, ELON[:so_j_stop,:]),
                                             "ELAT": (dim2_out, ELAT[:so_j_stop,:])})

            v_lt_ds   = xr.DataArray(v_lt_slc,   dims=dim3_out,
                                     coords={"NLON": (dim2_out, NLON[:so_j_stop,:]),
                                             "NLAT": (dim2_out, NLAT[:so_j_stop,:])})

            # Non-depth-avg part (rest of globe, surface only)
            sst_no_davg = xr.DataArray(CT_reG.isel(depth=0, nj=slice(so_j_stop, None)),
                                       dims=dim3_out,
                                       coords={"TLON": (dim2_out, TLON[so_j_stop:,:]),
                                               "TLAT": (dim2_out, TLAT[so_j_stop:,:])})

            sss_no_davg = xr.DataArray(S_reG.isel(depth=0, nj=slice(so_j_stop, None)),
                                       dims=dim3_out,
                                       coords={"TLON": (dim2_out, TLON[so_j_stop:,:]),
                                               "TLAT": (dim2_out, TLAT[so_j_stop:,:])})

            u_no_davg   = xr.DataArray(U_grid.isel(depth=0, nj=slice(so_j_stop, None)),
                                       dims=dim3_out,
                                       coords={"ELON": (dim2_out, ELON[so_j_stop:,:]),
                                               "ELAT": (dim2_out, ELAT[so_j_stop:,:])})

            v_no_davg   = xr.DataArray(V_grid.isel(depth=0, nj=slice(so_j_stop, None)),
                                       dims=dim3_out,
                                       coords={"NLON": (dim2_out, NLON[so_j_stop:,:]),
                                               "NLAT": (dim2_out, NLAT[so_j_stop:,:])})

            # Concatenate back to full globe
            sst_global = xr.concat([sst_lt_ds, sst_no_davg], dim="nj")
            sss_global = xr.concat([sss_lt_ds, sss_no_davg], dim="nj")
            u_global   = xr.concat([u_lt_ds ,  u_no_davg  ], dim="nj")
            v_global   = xr.concat([v_lt_ds ,  v_no_davg  ], dim="nj")

            print(f"\tdepth-avg branch took {time.time()-t1:.1f} s")
        else:
            print("Surface-only branch (depth=0)")
            t1 = time.time()
            sst_global = xr.DataArray(CT_reG.isel(depth=0).values,
                                      dims=dim3_out,
                                      coords={"TLON": (dim2_out, TLON),
                                              "TLAT": (dim2_out, TLAT)})
            sss_global = xr.DataArray(S_reG.isel(depth=0).values,
                                      dims=dim3_out,
                                      coords={"TLON": (dim2_out, TLON),
                                              "TLAT": (dim2_out, TLAT)})
            u_global   = xr.DataArray(U_grid.isel(depth=0).values,
                                      dims=dim3_out,
                                      coords={"ELON": (dim2_out, ELON),
                                              "ELAT": (dim2_out, ELAT)})
            v_global   = xr.DataArray(V_grid.isel(depth=0).values,
                                      dims=dim3_out,
                                      coords={"NLON": (dim2_out, NLON),
                                              "NLAT": (dim2_out, NLAT)})
            print(f"\t{time.time()-t1:.1f} s")

        # Insurance landmask (T mask); okay to use for u/v too
        t1 = time.time()
        sst_mask = sst_global.where(~kmt_mask)
        sss_mask = sss_global.where(~kmt_mask)
        u_mask   = u_global.where(~kmt_mask)
        v_mask   = v_global.where(~kmt_mask)
        print("applied kmt mask in {:.1f} s".format(time.time()-t1))

        # Intermediate tmp dataset for this segment
        t1 = time.time()
        ORAS_tmp = xr.Dataset(
            {
                "sst": (dim3_out, sst_mask.values,
                        {"units": "degrees_C",
                         "long_name": ("Sea Surface Temperature (CT over top 100 m)"
                                       if compute_davg else "Sea Surface Temperature (CT, surface)"),
                         "coordinates": "TLAT TLON"}),
                "sss": (dim3_out, sss_mask.values,
                        {"units": "psu",
                         "long_name": ("Sea Surface Salinity (top 100 m)"
                                       if compute_davg else "Sea Surface Salinity (surface)"),
                         "coordinates": "TLAT TLON"}),
                "u":   (dim3_out, u_mask.values,
                        {"units": "m s-1",
                         "long_name": ("Current along i on E points (top 5 m)"
                                       if ROTATE_TO_GRID else "Eastward current on E points (top 5 m)"
                                       if compute_davg else
                                       ("Current along i on E points (surface)"
                                        if ROTATE_TO_GRID else "Eastward current on E points (surface)")),
                         "coordinates": "ELAT ELON"}),
                "v":   (dim3_out, v_mask.values,
                        {"units": "m s-1",
                         "long_name": ("Current along j on N points (top 5 m)"
                                       if ROTATE_TO_GRID else "Northward current on N points (top 5 m)"
                                       if compute_davg else
                                       ("Current along j on N points (surface)"
                                        if ROTATE_TO_GRID else "Northward current on N points (surface)")),
                         "coordinates": "NLAT NLON"}),
            },
            coords={
                "time": (t_dim_name, ORAS_mos.time.values),
                "TLAT": (dim2_out, TLAT, {"units":"degrees_north","long_name":"T-point latitude"}),
                "TLON": (dim2_out, TLON, {"units":"degrees_east" ,"long_name":"T-point longitude"}),
                "ELAT": (dim2_out, ELAT, {"units":"degrees_north","long_name":"E-point latitude"}),
                "ELON": (dim2_out, ELON, {"units":"degrees_east" ,"long_name":"E-point longitude"}),
                "NLAT": (dim2_out, NLAT, {"units":"degrees_north","long_name":"N-point latitude"}),
                "NLON": (dim2_out, NLON, {"units":"degrees_east" ,"long_name":"N-point longitude"}),
            },
            attrs={
                "title"   : "Processed Ocean Dataset (ORAS regridded to 0.25° tripolar C-grid)",
                "summary" : ("Regridded depth-averages of key variables for CICE-SA (Southern Ocean only merged with surface elsewhere)."
                             if compute_davg else
                             "Regridded surface fields for CICE-SA."),
                "author"  : "DP@H2O; daniel.atwater@utas.edu.au",
                "source"  : "ORAS (CMEMS); processed with TEOS-10 and xESMF; mask from CICE kmt.",
                "Conventions": "CF-1.7"
            }
        )
        ORAS_tmp.encoding["unlimited_dims"] = {"time"}
        print("built ORAS_tmp in {:.1f} s".format(time.time()-t1))

        # Write tmp file
        if os.path.exists(P_ORAS_tmp):
            os.remove(P_ORAS_tmp)
        t1 = time.time()
        ORAS_tmp.to_netcdf(P_ORAS_tmp)
        print(f"wrote tmp: {P_ORAS_tmp} in {time.time()-t1:.1f} s")

        # If this is the 2nd half or a full-year file, build the yearly file
        if i == 1 or yr in [1994, 2000]:
            print(f"open wildcard files: {P_davg_wldcrd}")
            ORAS_cat = xr.open_mfdataset(P_davg_wldcrd, combine="by_coords")

            t1 = time.time()
            ORAS = xr.Dataset(
                {
                    "sst": (dim3_out, ORAS_cat.sst.values,
                            {"units":"degrees_C","long_name":ORAS_tmp.sst.long_name,
                             "coordinates":"TLAT TLON"}),
                    "sss": (dim3_out, ORAS_cat.sss.values,
                            {"units":"psu","long_name":ORAS_tmp.sss.long_name,
                             "coordinates":"TLAT TLON"}),
                    "u":   (dim3_out, ORAS_cat.u.values,
                            {"units":"m s-1","long_name":ORAS_tmp.u.long_name,
                             "coordinates":"ELAT ELON"}),
                    "v":   (dim3_out, ORAS_cat.v.values,
                            {"units":"m s-1","long_name":ORAS_tmp.v.long_name,
                             "coordinates":"NLAT NLON"}),
                },
                coords={
                    "time": (t_dim_name, ORAS_cat.time.values),
                    "TLAT": (dim2_out, TLAT, {"units":"degrees_north","long_name":"T-point latitude"}),
                    "TLON": (dim2_out, TLON, {"units":"degrees_east" ,"long_name":"T-point longitude"}),
                    "ELAT": (dim2_out, ELAT, {"units":"degrees_north","long_name":"E-point latitude"}),
                    "ELON": (dim2_out, ELON, {"units":"degrees_east" ,"long_name":"E-point longitude"}),
                    "NLAT": (dim2_out, NLAT, {"units":"degrees_north","long_name":"N-point latitude"}),
                    "NLON": (dim2_out, NLON, {"units":"degrees_east" ,"long_name":"N-point longitude"}),
                },
                attrs={
                    "title"   : "ORAS for CICE-SA (0.25° tripolar C-grid)",
                    "summary" : ORAS_tmp.attrs["summary"],
                    "author"  : "DP@H2O; daniel.atwater@utas.edu.au",
                    "source"  : "ORAS dataset processed with TEOS-10 and xESMF.",
                    "Conventions": "CF-1.7"
                }
            )
            ORAS.encoding["unlimited_dims"] = {"time"}
            print(f"built yearly dataset in {time.time()-t1:.1f} s")

            if os.path.exists(P_ORAS_fin):
                os.remove(P_ORAS_fin)
            t1 = time.time()
            ORAS.to_netcdf(P_ORAS_fin)
            print(f"wrote yearly file: {P_ORAS_fin} in {time.time()-t1:.1f} s")

    print(f"year {yr} total time: {time.time()-year_t0:.1f} s\n")


# import sys, os, time
# import gsw
# import copernicusmarine
# import xesmf             as xe
# import numpy             as np
# import pandas            as pd
# import xarray            as xr

# # Initialise
# compute_davg= False
# ow_nc       = False
# yr0         = 1993
# yrN         = 2025
# dmin        = 0.505760014057159
# dmax        = 5
# D_ORAS_base = os.path.join("/","g","data","gv90","da1339","afim_input","ORAS","daily")
# D_ORAS_org  = os.path.join(D_ORAS_base,"org")
# if compute_davg:
#     D_ORAS_tmp = os.path.join(D_ORAS_base,f"d{int(dmax)}","tmp")
#     D_ORAS_fin = os.path.join(D_ORAS_base,f"d{int(dmax)}","Cgrd")
# else:
#     D_ORAS_tmp = os.path.join(D_ORAS_base,"sfc","tmp")
#     D_ORAS_fin = os.path.join(D_ORAS_base,"sfc","Cgrd")
# pt_name     = "thetao_oras" # potential temperature
# s_name      = "so_oras"     # salinity
# d_name      = "depth"       # depth coordinate and dimension name
# u_name      = "uo_oras"     # eastward current
# v_name      = "vo_oras"     # northward current
# t_dim_name  = "time"
# x_dim_name  = "ni"
# y_dim_name  = "nj"
# z_dim_name  = "nk"
# dim4_out    = [t_dim_name, z_dim_name, y_dim_name, x_dim_name]
# dim3_out    = [t_dim_name, y_dim_name, x_dim_name]
# dim2_out    = [y_dim_name, x_dim_name]
# var_list    = ["thetao_oras", "so_oras", "uo_oras", "vo_oras"]
# reG_meth    = "bilinear"
# extrap_meth = "inverse_dist"
# P_wghts     = os.path.join("/","g","data","gv90","da1339","grids",'weights',"map_ORAS_to_AOM3_0p25_xesmf-bil_20251009.nc")
# P_kmt       = os.path.join("/","g","data","gv90","da1339","grids","ACCESS-OM3-025_kmt.nc")
# P_G         = os.path.join("/","g","data","gv90","da1339","grids","ACCESS-OM3-025_Cgrid.nc")
# kmt         = xr.open_dataset(P_kmt).kmt
# kmt_mask    = xr.where(kmt == 0, True, False).compute()#.rename({'ny':'nj','nx':'ni'}).compute()
# G_CICE      = xr.open_dataset(P_G)
# TLAT        = G_CICE.tlat.values * 180 / np.pi
# TLON        = G_CICE.tlon.values * 180 / np.pi
# ds_tcice    = xr.Dataset({'mask' : (dim2_out, kmt.values)},
#                           coords = dict(lon = (["ni"], TLON[0,:]),
#                                         lat = (["nj"], TLAT[:,0])))
# # loop over years
# for yr in range(yr0,yrN):
#     t0 = time.time()
#     if yr in [1994, 2000]:
#         F_ORAS_orgs = [ f"{yr}0101_{yr}1231_ORAS_org.nc" ]
#         F_ORAS_tmps = [ f"{yr}0101_{yr}1231_ORAS_tmp.nc" ]
#     else:
#         F_ORAS_orgs = [f"{yr}0101_{yr}0630_ORAS_org.nc",f"{yr}0701_{yr}1231_ORAS_org.nc"]
#         F_ORAS_tmps = [f"{yr}0101_{yr}0630_ORAS_tmp.nc",f"{yr}0701_{yr}1231_ORAS_tmp.nc"]
#     for i,F_ORAS_org in enumerate(F_ORAS_orgs):
#         P_ORAS_org    = os.path.join(D_ORAS_org,F_ORAS_org)
#         P_ORAS_tmp    = os.path.join(D_ORAS_tmp,F_ORAS_tmps[i])
#         F_ORAS_fin    = f"ORAS_{yr}.nc"
#         P_ORAS_fin    = os.path.join(D_ORAS_fin,F_ORAS_fin)
#         P_davg_wldcrd = os.path.join(D_ORAS_tmp,f"{yr:04d}*_tmp.nc")
#         P_ORAS_fin    = os.path.join(D_ORAS_fin,f"ORAS_{yr:04d}.nc")
#         if os.path.exists(P_ORAS_fin) and not ow_nc:
#             print(f"file already exists and over-writing enabled; skipping {P_ORAS_fin}")
#             continue
#         else:
#             print(f"working on creating {P_ORAS_fin}")
#         if os.path.exists(P_ORAS_tmp) and not ow_nc:
#             print(f"file already exists and over-writing enabled; skipping {P_ORAS_tmp}")
#             continue
#         else:
#             print(f"working on creating {P_ORAS_tmp}")
#         t0        = time.time()
#         ORAS_mos  = xr.open_dataset(P_ORAS_org, engine="netcdf4")
#         PT        = ORAS_mos[pt_name]
#         S         = ORAS_mos[s_name]
#         d         = ORAS_mos[d_name]
#         U         = ORAS_mos[u_name]
#         V         = ORAS_mos[v_name]
#         print(f"loading data took {time.time()-t0:.1f} seconds")
#         ds_t         = PT.isel(time=0,depth=0)
#         ds_t["mask"] = xr.where(~np.isnan(ds_t), 1, 0)
#         print("loading weights:")
#         t1       = time.time()
#         if os.path.exists(P_wghts):
#             reG_tbil = xe.Regridder(ds_t, ds_tcice, method=reG_meth, periodic=True, extrap_method=extrap_meth, weights=P_wghts)
#         else:
#             reG_tbil = xe.Regridder(ds_t, ds_tcice, method=reG_meth, periodic=True, extrap_method=extrap_meth, filename=P_wghts)
#         print(f"\ttime taken: {time.time()-t1:.1f} seconds")
#         print("regridding PT")
#         t1     = time.time()
#         PT_reG = reG_tbil(PT)
#         print(f"\ttime taken: {time.time()-t1:.1f} seconds")
#         print("regridding S")
#         t1     = time.time()
#         S_reG  = reG_tbil(S)
#         print(f"\ttime taken: {time.time()-t1:.1f} seconds")
#         print("regridding U")
#         t1     = time.time()
#         U_reG  = reG_tbil(U)
#         print(f"\ttime taken: {time.time()-t1:.1f} seconds")
#         print("regridding V")
#         t1     = time.time()
#         V_reG  = reG_tbil(V)
#         print(f"\ttime taken: {time.time()-t1:.1f} seconds")
#         print("convert potential temperature into conservative temperature")
#         t1 = time.time()
#         CT = gsw.conversions.CT_from_pt( S_reG , PT_reG )
#         print(f"\ttime taken: {time.time()-t1:.1f} seconds")
#         if compute_davg:
#             print("slicing-out Southern Ocean")
#             t1          = time.time()
#             CT_lt_slc   = CT.isel(nj=slice(0,233))
#             S_lt_slc    = S_reG.isel(nj=slice(0,233))
#             U_lt_slc    = U_reG.isel(nj=slice(0,233))
#             V_lt_slc    = V_reG.isel(nj=slice(0,233))
#             TLON_lt_slc = G_CICE.isel(ny=slice(0,233)).tlon.values * 180 / np.pi 
#             TLAT_lt_slc = G_CICE.isel(ny=slice(0,233)).tlat.values * 180 / np.pi 
#             print(f"\ttime taken: {time.time()-t1:.1f} seconds")
#             print("compute depth-averages for Southern Ocean")
#             t1         = time.time()
#             d_diff     = np.diff(d, append=d[-1])  # Depth layer thickness
#             d_wghts    = xr.DataArray(d_diff, dims=[d_name], coords={d_name: d})
#             d100_mask  = d<=dmax
#             wghts_d100 = d_wghts.sel(depth=d[d100_mask])
#             ct_d100    = CT_lt_slc.sel(depth=d[d100_mask])
#             s_d100     = S_lt_slc.sel(depth=d[d100_mask])
#             sst_lt_slc = (ct_d100 * wghts_d100).sum(dim=d_name) / wghts_d100.sum()
#             sss_lt_slc = (s_d100  * wghts_d100).sum(dim=d_name) / wghts_d100.sum()
#             d_diff     = np.diff(d[:4],append=d[4])
#             d_wghts    = xr.DataArray(d_diff, dims=[d_name], coords={d_name: d[:4]})
#             wghts_d5   = d_wghts.sel(depth=d[:4])
#             u_d5       = U_lt_slc.sel(depth=d[:4])
#             v_d5       = V_lt_slc.sel(depth=d[:4])
#             u_lt_slc   = (u_d5 * wghts_d5).sum(dim=d_name) / wghts_d5.sum()
#             v_lt_slc   = (v_d5 * wghts_d5).sum(dim=d_name) / wghts_d5.sum()
#             print(f"\ttime taken: {time.time()-t1:.1f} seconds")
#             print("Create datasets for the depth-averaged region (nj=0:233)")
#             t1            = time.time()
#             sst_lt_slc_ds = xr.DataArray(sst_lt_slc, 
#                                          dims=dim3_out,
#                                          coords={"TLON": (dim2_out,TLON_lt_slc),
#                                                  "TLAT": (dim2_out,TLAT_lt_slc)})
#             sss_lt_slc_ds = xr.DataArray(sss_lt_slc,
#                                          dims=dim3_out, 
#                                          coords={"TLON": (dim2_out,TLON_lt_slc),
#                                                  "TLAT": (dim2_out,TLAT_lt_slc)})
#             u_lt_slc_ds   = xr.DataArray(u_lt_slc,
#                                          dims=dim3_out,
#                                          coords={"TLON": (dim2_out,TLON_lt_slc),
#                                                  "TLAT": (dim2_out,TLAT_lt_slc)})
#             v_lt_slc_ds   = xr.DataArray(v_lt_slc,
#                                          dims=dim3_out,
#                                          coords={"TLON": (dim2_out,TLON_lt_slc),
#                                                  "TLAT": (dim2_out,TLAT_lt_slc)})
#             print(f"\ttime taken: {time.time()-t1:.1f} seconds")
#             print("Slice non-depth-averaged data (nj=234:-1, depth=0)")
#             t1           = time.time()
#             TLON_no_davg = G_CICE.isel(ny=slice(233,1080)).tlon.values * 180 / np.pi 
#             TLAT_no_davg = G_CICE.isel(ny=slice(233,1080)).tlat.values * 180 / np.pi 
#             sst_no_davg  = xr.DataArray(CT.isel(depth=0, nj=slice(233,1080)).drop_vars(['depth','lon','lat']),
#                                         dims=dim3_out, 
#                                         coords={"TLON": (dim2_out,TLON_no_davg),
#                                                 "TLAT": (dim2_out,TLAT_no_davg)})
#             sss_no_davg  = xr.DataArray(S_reG.isel(depth=0, nj=slice(233,1080)).drop_vars(['depth','lon','lat']),
#                                         dims=dim3_out, 
#                                         coords={"TLON": (dim2_out,TLON_no_davg),
#                                                 "TLAT": (dim2_out,TLAT_no_davg)})
#             u_no_davg    = xr.DataArray(U_reG.isel(depth=0, nj=slice(233,1080)).drop_vars(['depth','lon','lat']),
#                                         dims=dim3_out, 
#                                         coords={"TLON": (dim2_out,TLON_no_davg),
#                                                 "TLAT": (dim2_out,TLAT_no_davg)})
#             v_no_davg    = xr.DataArray(V_reG.isel(depth=0, nj=slice(233,1080)).drop_vars(['depth','lon','lat']),
#                                         dims=dim3_out, 
#                                         coords={"TLON": (dim2_out,TLON_no_davg),
#                                                 "TLAT": (dim2_out,TLAT_no_davg)})
#             print(f"\ttime taken: {time.time()-t1:.1f} seconds")
#             print("Concatenate depth-averaged and non-depth-averaged data")
#             t1         = time.time()
#             sst_global = xr.concat([sst_lt_slc_ds, sst_no_davg], dim="nj")
#             sss_global = xr.concat([sss_lt_slc_ds, sss_no_davg], dim="nj")
#             u_global   = xr.concat([u_lt_slc_ds  , u_no_davg  ], dim="nj")
#             v_global   = xr.concat([v_lt_slc_ds  , v_no_davg  ], dim="nj")
#             print(f"\ttime taken: {time.time()-t1:.1f} seconds")
#         else:
#             print("NO DEPTH AVERAGING")
#             print("creating global data arrays of just surface values (i.e. isel( depth=0 ))")
#             t1          = time.time()
#             sst_global  = xr.DataArray(CT.isel(depth=0).values,
#                                         dims=dim3_out, 
#                                         coords={"TLON": (dim2_out,TLON),
#                                                 "TLAT": (dim2_out,TLAT)})
#             sss_global  = xr.DataArray(S_reG.isel(depth=0).values,
#                                         dims=dim3_out, 
#                                         coords={"TLON": (dim2_out,TLON),
#                                                 "TLAT": (dim2_out,TLAT)})
#             u_global    = xr.DataArray(U_reG.isel(depth=0).values,
#                                         dims=dim3_out, 
#                                         coords={"TLON": (dim2_out,TLON),
#                                                 "TLAT": (dim2_out,TLAT)})
#             v_global    = xr.DataArray(V_reG.isel(depth=0).values,
#                                         dims=dim3_out, 
#                                         coords={"TLON": (dim2_out,TLON),
#                                                 "TLAT": (dim2_out,TLAT)})
#             print(f"\ttime taken: {time.time()-t1:.1f} seconds")
#         print("insurance landmask")
#         t1       = time.time()
#         sst_mask = sst_global.where(~kmt_mask)
#         sss_mask = sss_global.where(~kmt_mask)
#         u_mask   = u_global.where(~kmt_mask)
#         v_mask   = v_global.where(~kmt_mask)
#         print(f"\ttime taken: {time.time()-t1:.1f} seconds")
#         print("creating intermediate ('tmp') dataset")
#         t1 = time.time()
#         ORAS_tmp = xr.Dataset({"sst": (dim3_out, sst_mask.values , {"units": "degrees_C", "long_name": "Sea Surface Temperature (CT over top 100m)"}),
#                                "sss": (dim3_out, sss_mask.values , {"units": "psu"      , "long_name": "Sea Surface Salinity (over top 100m)"}      ),
#                                "u"  : (dim3_out, u_mask.values   , {"units": "m/s"      , "long_name": "Surface Eastward Velocity (over top 5m)"}   ),
#                                "v"  : (dim3_out, v_mask.values   , {"units": "m/s"      , "long_name": "Surface Northward Velocity (over top 5m)"}  )},
#                                coords={"time" : (t_dim_name, ORAS_mos.time.values),
#                                        "TLAT" : (dim2_out  , TLAT , {"units" : "degrees_north", "long_name": "grid cell centre latitudes"}),
#                                        "TLON" : (dim2_out  , TLON , {"units" : "degrees_east" , "long_name": "grid cell centre longitudes"})},
#                                attrs={"title"   : "Processed Ocean Dataset",
#                                       "summary" : "re-gridded depth-averages of key variables for CICE-SA.",
#                                       "author"  : "DP@H2O; daniel.atwater@utas.edu.au",
#                                       "source"  : "Original CMEMS dataset, processed with TEOS-10 and weighted averaging.",
#                                       "Conventions": "CF-1.7"})
#         ORAS_tmp.encoding["unlimited_dims"] = {"time"}
#         print(f"\ttime taken: {time.time()-t1:.1f} seconds")
#         print(f"writing file to disk: {P_ORAS_tmp}")
#         if os.path.exists(P_ORAS_tmp):
#             os.remove(P_ORAS_tmp)
#         t1 = time.time()
#         ORAS_tmp.to_netcdf(P_ORAS_tmp)
#         print(f"\ttime taken: {time.time()-t1:.1f} seconds")
#         if i==1 or yr==1994 or yr==2000:
#             print(f"open wildcard files: {P_davg_wldcrd}")
#             ORAS_cat = xr.open_mfdataset(P_davg_wldcrd)
#             print("creating yearly dataset")
#             t1   = time.time() 
#             ORAS = xr.Dataset({"sst": (dim3_out, ORAS_cat.sst.values , {"units": "degrees_C", "long_name": "Sea Surface Temperature (CT over top 100m)"}),
#                                "sss": (dim3_out, ORAS_cat.sss.values , {"units": "psu"      , "long_name": "Sea Surface Salinity (over top 100m)"}      ),
#                                "u"  : (dim3_out, ORAS_cat.u.values   , {"units": "m/s"      , "long_name": "Surface Eastward Velocity (over top 5m)"}   ),
#                                "v"  : (dim3_out, ORAS_cat.v.values   , {"units": "m/s"      , "long_name": "Surface Northward Velocity (over top 5m)"}  )},
#                               coords={"time" : (t_dim_name, ORAS_cat.time.values)    ,
#                                       "TLAT" : (dim2_out  , TLAT , {"units" : "degrees_north", "long_name": "grid cell centre latitudes"}),
#                                       "TLON" : (dim2_out  , TLON , {"units" : "degrees_east" , "long_name": "grid cell centre longitudes"})},
#                               attrs={"title"   : "ORAS for CICE-SA",
#                                      "summary" : "Re-gridded and depth-averaged ORAS key variables for CICE-SA; masked with CICE-SA kmt file",
#                                      "author"  : "DP@H2O; daniel.atwater@utas.edu.au",
#                                      "source"  : "ORAS dataset (https://os.copernicus.org/articles/15/779/2019/), processed with TEOS-10 and weighted averaging.",
#                                      "Conventions": "CF-1.7"})
#             ORAS.encoding["unlimited_dims"] = {"time"}
#             print(f"\ttime taken: {time.time()-t1:.1f} seconds")
#             print(f"writing file to disk: {P_ORAS_fin}")
#             if os.path.exists(P_ORAS_fin):
#                 os.remove(P_ORAS_fin)
#             t1 = time.time()
#             ORAS.to_netcdf(P_ORAS_fin)
#             print(f"\ttime taken: {time.time()-t1:.1f} seconds")
#     print(f"year took total time of {time.time()-t0:.1f} seconds\n\n")
