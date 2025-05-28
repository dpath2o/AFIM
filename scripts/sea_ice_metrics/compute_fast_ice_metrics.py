import sys, os, argparse, pygmt
sys.path.insert(0, '/home/581/da1339/AFIM/src/AFIM/src')
from sea_ice_processor import SeaIceProcessor
from sea_ice_plotter   import SeaIcePlotter
from pathlib           import Path
import xarray          as xr
import numpy           as np

def plot_cice_field(x1, y1, z1, title, region, 
                    cmap="cmocean/haline", cmap_series=[0.9, 1],
                    GI_coords=None, GI_color="red", 
                    var_style="s0.5c", GI_style="c0.3c", 
                    projection="S75/-90/30c", 
                    cbar_label="sea ice concentration", units="1/100"):
    fig = pygmt.Figure()
    pygmt.config(FORMAT_GEO_MAP="ddd.x", FONT_TITLE="16p,Courier-Bold", FONT_ANNOT_PRIMARY="14p,Helvetica")
    pygmt.makecpt(cmap=cmap, series=cmap_series)
    fig.basemap(projection=projection, region=region, frame=["af", f"+t{title}"])
    fig.plot(x=x1, y=y1, fill=z1, cmap=True, style=var_style)
    if GI_coords:
        fig.plot(x=GI_coords[0], y=GI_coords[1], fill=GI_color, style=GI_style)
    fig.coast(land='gray')
    fig.colorbar(position="JBC+o0c/1.5c", frame=[f"x+l{cbar_label}", f"y+l{units}"])
    return fig

def plot_persistence_map(DA, sim_name, itype, ispd_str, ktens, elps, GI_thin, SI_plot, dt_range_str, overwrite_png):
    da   = DA.values.ravel()
    lon  = DA["TLON"].values.ravel()
    lat  = DA["TLAT"].values.ravel()
    mask = ~np.isnan(da)
    lat_plt, lon_plt, da_plt = lat[mask], lon[mask], da[mask]
    for reg_name, reg_cfg in SI_plot.reg_dict.items():
        D_plt = Path(SI_plot.D_graph, sim_name, reg_name, "FIP", f"ispd_thresh_{ispd_str}", itype)
        D_plt.mkdir(parents=True, exist_ok=True)
        P_plt = Path(D_plt,f"FIP_{dt_range_str}.png")
        if P_plt.exists() and not overwrite_png:
            print(f"figure {P_plt} exists and not overwriting")
            continue
        else:
            MC  = SI_plot.get_meridian_center_from_geographic_extent(reg_cfg['plt_ext'])
            fig = plot_cice_field(lon_plt, lat_plt, da_plt,
                                  f"{sim_name} {itype} ispd_thresh={ispd_str}: ktens={ktens}, elps={elps}, GI-thin={GI_thin:.2f}",
                                  region=reg_cfg['plt_ext'],
                                  cmap="/home/581/da1339/graphical/CPTs/AF2020_YlGnBu.cpt",
                                  cmap_series=[0, 1],
                                  GI_coords=(SI_plot.GI_proc.G_t['GI_lon'], SI_plot.GI_proc.G_t['GI_lat']),
                                  GI_color="red", var_style="s0.3c", GI_style="c0.1c", 
                                  projection=f"S{MC}/-90/30c", cbar_label="fast ice persistence")
            print(f"📸 Saving figure: {P_plt}")
            fig.savefig(P_plt)

def compute_and_save_metrics(SI_proc, DS_FI, P_METS, obs_clim=None):
    METS = {}
    # 3D + 1D metrics
    FIA = SI_proc.compute_ice_area(DS_FI['aice'], DS_FI['tarea']).compute()
    FIP = SI_proc.compute_variable_aggregate(DS_FI['aice']).compute()
    METS["FIA"] = FIA
    METS["FIP"] = FIP
    # Scalar / 1D metrics
    summary = {}
    summary["onset_doy"]    = SI_proc.compute_fia_onset_doy(FIA)
    summary["growth_rate"]  = SI_proc.compute_fia_growth_rate(FIA)
    summary["max_growth"]   = SI_proc.compute_fia_max_growth(FIA)
    summary["doy_max"]      = SI_proc.compute_doy_max(FIA)
    summary["duration"]     = SI_proc.compute_fast_ice_duration(FIA)
    fip_mean, fip_std       = SI_proc.compute_fip_spatial_stats(FIP)
    summary["FIP_mean"]     = fip_mean
    summary["FIP_std"]      = fip_std
    mean_dist, max_dist     = SI_proc.compute_fast_ice_distance_extent(DS_FI['FI_mask'])
    summary["mean_FI_dist"] = mean_dist
    summary["max_FI_dist"]  = max_dist

    if obs_clim is not None:
        model_doy = FIA["time"].dt.dayofyear.values
        obs_vals = np.interp(model_doy, obs_clim.coords["doy"].values, obs_clim.values)
        summary["rmse_to_obs"] = SI_proc.compute_fia_rmse(FIA, xr.DataArray(obs_vals, coords=[("time", FIA["time"])]))
    # Combine all into a dataset
    DS_METS = xr.Dataset(summary)
    DS_METS["FIA"] = FIA
    DS_METS["FIP"] = FIP
    DS_METS.to_zarr(P_METS, mode="w", consolidated=True)
    print(f"📊 Metrics written to {P_METS}")
    return DS_METS


def save_metrics_csv(metrics_dict, sim_name, i_type, ispd_str, D_out):
    from pandas import DataFrame
    df = DataFrame([metrics_dict])
    df["sim_name"]     = sim_name
    df["ice_type"]     = i_type
    df["ispd_thresh"]  = ispd_str
    df.to_csv(Path(D_out, f"metrics_summary_{i_type}.csv"), index=False)


def main(sim_name, ispd_thresh, ice_type, compute_boolean, smooth_FIA_days, overwrite_zarr, overwrite_png):
    print(f"ice_type passed is: {ice_type}")
    ispd_str = f"{ispd_thresh:.1e}".replace("e-0", "e-")
    SI_proc  = SeaIceProcessor(sim_name = sim_name)
    SI_plot  = SeaIcePlotter(sim_name  = sim_name,
                             ice_type  ='FI', 
                             plot_type = 'timeseries',
                             var_name  = 'FIA',
                             show_fig  = False,
                             save_fig  = True,
                             overwrite = True)
    dt0_str, dtN_str = SI_proc.sim_config['dt0_str'], SI_proc.sim_config['dtN_str']
    dt_range_str     = f"{dt0_str[:4]}-{dtN_str[:4]}"
    cfg              = SI_proc.sim_config
    ktens, elps, GI_thin = cfg.get("Ktens", "?"), cfg.get("e_f", "?"), 1 - cfg.get("GI_thin_fact", 0)
    ice_types = [ice_type] if isinstance(ice_type, str) else list(ice_type)
    D_out = Path(SI_proc.config['D_dict']['AFIM_out'], sim_name, "zarr", f"ispd_thresh_{ispd_str}", "ice_metrics")
    D_out.mkdir(parents=True, exist_ok=True)
    FIA_comp, FIP_comp = {}, {}
    obs_clim = None
    try:
        af2020_df = SI_proc.load_AF2020_FI_area_CSV(doy_start=1)
        obs_clim  = SI_proc.interpolate_obs_fia(af2020_df)
    except Exception as e:
        print(f"⚠️ Could not load observational climatology: {e}")
    for itype in ice_type:
        for roll in [False,True]:
            i_type = f"{itype}_roll" if roll else itype
            P_METS = Path(D_out, f"metrics_{i_type}.zarr")
            if P_METS.exists() and not overwrite_zarr:
                print(f"{P_METS} exists and not overwriting--loading")
                METS = xr.open_zarr(P_METS)
            else:
                print(f"{P_METS} does NOT exists and/or overwriting--computing")
                DS_FI, CICE_SO = SI_proc.load_processed_cice(dt0_str     = dt0_str,
                                                             dtN_str     = dtN_str,
                                                             ispd_thresh = ispd_thresh,
                                                             ice_type    = itype,
                                                             zarr_CICE   = True,
                                                             rolling     = roll)
                if DS_FI is None: continue
                METS = compute_and_save_metrics(SI_proc, DS_FI, P_METS, obs_clim=obs_clim)
                save_metrics_csv(METS, sim_name=sim_name, i_type=i_type, ispd_str=ispd_str, D_out=D_out)
            FIA_comp[i_type] = METS['FIA']
            FIP_comp[i_type] = METS['FIP']
            plot_persistence_map(METS['FIP'], sim_name, i_type, ispd_str, ktens, elps, GI_thin, SI_plot, dt_range_str, overwrite_png)
            if compute_boolean and not roll:
                P_METS, METS = None, None
                i_type = f"{itype}_bool"
                P_METS = Path(D_out,f"metrics_{i_type}.zarr")
                if P_METS.exists() and not overwrite_zarr:
                    print(f"{P_METS} exists and not overwriting--loading")
                    METS = xr.open_zarr(P_METS)
                else:
                    print(f"{P_METS} does NOT exists and/or overwriting--computing")
                    bool_mask = SI_proc.boolean_fast_ice(DS_FI['FI_mask'], dim="time", window=7, min_count=6)
                    DS_bool   = CICE_SO.where(bool_mask)
                    METS      = compute_and_save_metrics(SI_proc, DS_bool, P_METS, obs_clim=obs_clim)
                    save_metrics_csv(METS, sim_name=sim_name, i_type=i_type, ispd_str=ispd_str, D_out=D_out)
                FIA_comp[i_type] = METS['FIA']
                FIP_comp[i_type] = METS['FIP']
                plot_persistence_map(METS['FIP'], sim_name, i_type, ispd_str, ktens, elps, GI_thin, SI_plot, dt_range_str, overwrite_png)
    tit_str = f"{sim_name} ispd_thresh={ispd_str}: ktens={ktens}, elps={elps}, GI-thin={GI_thin:.2f}"
    if "duration" in METS:
        tit_str += f", dur={METS['duration']}"
    if smooth_FIA_days>0:
        P_png = Path(SI_plot.D_graph, "timeseries", f"FIA_{sim_name}_{ispd_str}_smoothed_{dt_range_str}.png")
    else:
        P_png = Path(SI_plot.D_graph, "timeseries", f"FIA_{sim_name}_{ispd_str}_{dt_range_str}.png")
    if P_png.exists and not overwrite_png:
        print(f"{P_png} exists and not overwriting")
    else:
        print(f"{P_png} does NOT exists and/or overwriting--creating")
        SI_plot.plot_ice_area(FIA_comp, tit_str=tit_str, P_png=P_png, roll_days=smooth_FIA_days)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute FIA and FIP metrics, apply boolean mask, and plot spatial + temporal outputs.")
    parser.add_argument("--sim_name"        , type=str, required=True)
    parser.add_argument("--ispd_thresh"     , type=float, required=True)
    parser.add_argument("--ice_type"        , nargs="+", default=["FI_B", "FI_Ta", "FI_Tx", "FI_BT"])
    parser.add_argument("--compute_boolean" , action="store_true")
    parser.add_argument("--overwrite_zarr"  , action="store_true")
    parser.add_argument("--overwrite_png"   , action="store_true")
    parser.add_argument("--smooth_FIA_days" , type=int, default=0)
    parser.set_defaults(compute_boolean=False, overwrite_zarr=False, overwrite_png=False)
    args     = parser.parse_args()
    ice_type = args.ice_type
    # Flatten any space-separated or comma-separated strings into a proper list
    if isinstance(ice_type, list):
        if len(ice_type) == 1:
            ice_type = ice_type[0].replace(",", " ").split()
    main(args.sim_name,
         args.ispd_thresh,
         ice_type,
         args.compute_boolean,
         args.smooth_FIA_days,
         args.overwrite_zarr,
         args.overwrite_png)
