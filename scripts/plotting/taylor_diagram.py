#!/usr/bin/env python3
import os
import numpy             as np
import pandas            as pd
import xarray            as xr
import matplotlib.pyplot as plt

def align_and_subset(obs, cice, aom2):
    for da in [obs, cice, aom2]:
        da["time"] = pd.to_datetime(da.time.values).normalize()
    t_common = np.intersect1d(np.intersect1d(obs.time.values, cice.time.values), aom2.time.values)
    return obs.sel(time=t_common), cice.sel(time=t_common), aom2.sel(time=t_common)

def compute_taylor_stats_streaming(da1, da2, chunk_size=30):
    all_model, all_obs = [], []
    for i in range(0, da1.sizes['time'], chunk_size):
        a = da1.isel(time=slice(i, i+chunk_size)).compute().values.flatten()
        b = da2.isel(time=slice(i, i+chunk_size)).compute().values.flatten()
        valid = np.isfinite(a) & np.isfinite(b)
        if valid.sum() > 0:
            all_model.append(a[valid])
            all_obs.append(b[valid])
    model, obs = np.concatenate(all_model), np.concatenate(all_obs)
    corr = np.corrcoef(model, obs)[0, 1]
    std_model, std_obs = np.std(model), np.std(obs)
    rmsd = np.sqrt(np.mean((model - obs - (np.mean(model) - np.mean(obs)))**2))
    return {"corr": corr, "std_ratio": std_model / std_obs, "rmsd_ratio": rmsd / std_obs}

def plot_taylor(stats_dict, out_path):
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111, polar=True)
    ax.set_theta_zero_location('E')
    ax.set_theta_direction(-1)

    for rms in [0.2, 0.5, 1.0, 1.5]:
        rs = np.linspace(0.5, 2.0, 500)
        theta = np.arccos(np.clip(1 - (rms**2 - 1)/(2 * rs), -1, 1))
        ax.plot(theta, rs, '--', color='gray', lw=0.6)
    ax.plot([0], [1], 'ko', label='Reference')

    for label, stat in stats_dict.items():
        angle = np.arccos(stat["corr"])
        r = stat["std_ratio"]
        ax.plot(angle, r, 'o', label=label)

    ax.set_rmax(2)
    ax.set_rticks([0.5, 1.0, 1.5, 2.0])
    ax.set_rlabel_position(135)
    ax.set_title("Taylor Diagram: Sea Ice Speed Comparison", fontsize=12)
    ax.legend(loc='upper right', bbox_to_anchor=(1.45, 1.1))
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()

def main():
    mod_path = '/home/581/da1339/AFIM/src/AFIM/src'
    sys.path.insert(0, mod_path)
    from sea_ice_toolbox import SeaIceToolbox
    dt0_str = "1993-01-01"
    dtN_str = "1999-12-31"
    out_dir = "/g/data/gv90/da1339/GRAPHICAL/AFIM/taylor_diagram"
    sim_name = "gi-nil"
    os.makedirs(out_dir, exist_ok=True)

    # Load datasets
    SI_cice = SeaIceToolbox(sim_name, dt0_str, dtN_str, True, True, True, True)
    _, CICE = SI_cice.load_processed_cice(zarr_CICE=True)
    CICE = CICE.isel(nj=SI_cice.hemisphere_dict['nj_slice'])

    SI_aom2 = SeaIceToolbox("AOM2-ERA5", dt0_str, dtN_str, True, True, True, True)
    _, AOM2 = SI_aom2.load_processed_cice(zarr_CICE=True)
    AOM2 = AOM2.isel(nj=SI_aom2.hemisphere_dict['nj_slice'])

    OSI_SAF = xr.open_mfdataset("/home/581/da1339/seaice/OSI_SAF/ispd_reG_SH*")

    da_obs = OSI_SAF.isel(nj=slice(0, 540)).ice_speed * 1000
    da_cice = CICE.isel(nj=slice(0, 540)).ispd_BT
    da_aom2 = AOM2.isel(nj=slice(0, 540)).ispd_BT

    da_obs, da_cice, da_aom2 = align_and_subset(da_obs, da_cice, da_aom2)

    stats = {
        "CICE vs OSI-SAF": compute_taylor_stats_streaming(da_cice, da_obs),
        "AOM2 vs OSI-SAF": compute_taylor_stats_streaming(da_aom2, da_obs),
        "CICE vs AOM2":    compute_taylor_stats_streaming(da_cice, da_aom2)
    }

    plot_taylor(stats, os.path.join(out_dir, f"{sim_name}_taylor_diagram_seaice_speed.png"))

if __name__ == "__main__":
    from multiprocessing import freeze_support
    freeze_support()
    main()
