# create_animations.py
import sys, os
import numpy as np
import pandas as pd
from pathlib import Path

def main():
    mod_path = '/home/581/da1339/AFIM/src/AFIM/src'
    sys.path.insert(0, mod_path)
    from sea_ice_toolbox import SeaIceToolbox
    sims       = ["elps-min", "gi-nil", "AOM2-ERA5"]
    regions    = ["north", "south"]
    var_names  = ["aice", "hi"]
    video_paths = {sim: {region: {} for region in regions} for sim in sims}
    dt0_str     = "1999-08-01"
    dtN_str     = "1999-10-31"
    # STEP 1: Generate PNG frames
    for sim_name in sims:
        for var_name in var_names:
            SI_tools = SeaIceToolbox(sim_name             = sim_name,
                                     dt0_str              = dt0_str,
                                     dtN_str              = dtN_str,
                                     P_log                = f"/g/data/gv90/da1339/logs/SeaIceToolBox_{sim_name}_{var_name}_plotting.log",
                                     save_new_figs        = True,
                                     show_figs            = False,
                                     overwrite_saved_figs = True)
            CICE_all = SI_tools.load_iceh_zarr()
            CICE_slice = CICE_all.sel(time=slice(dt0_str,dtN_str))
            for i,dt in enumerate(CICE_slice['time'].values):
                ds_slice = CICE_slice.isel(time=i,nj=SI_tools.hemisphere_dict['nj_slice'])
                plt_da = ds_slice[var_name]
                dt = pd.Timestamp(plt_da.time.values)
                dt_str = f"{dt.year}-{dt.month:02d}-{dt.day:02d}"
                SI_tools.pygmt_map_plot_one_var(plt_da, var_name,
                                                plot_regions  = None,
                                                time_stamp    = dt_str,
                                                tit_str       = dt_str,
                                                plot_GI       = bool(SI_tools.use_gi),
                                                var_sq_size   = 0.2,
                                                GI_fill_color = 'red',
                                                GI_sq_size    = 0.01,
                                                extend_cbar   = True,
                                                show_fig      = False)
            SI_tools.client.close()

    # STEP 2: Convert PNGs to MP4s
    for sim_name in sims:
        for region in regions:
            for var_name in var_names:
                SI_tools = SeaIceToolbox(sim_name=sim_name)
                D_png = Path(SI_tools.D_graph, sim_name, region, var_name)
                D_ani = Path(SI_tools.D_graph, "animations", sim_name, var_name)
                P_mp4 = Path.home() / "AFIM" / "src" / "AFIM" / "docs" / "figures" / f"{sim_name}_{var_name}_{region}_1999-winter.mp4"
                D_ani.mkdir(parents=True, exist_ok=True)
                frames = sorted([f for f in os.listdir(D_png) if f.endswith(".png")])
                os.system(f"rm {SI_tools.D_tmp}/frame_*.png")
                for i, f in enumerate(frames):
                    src = D_png / f
                    dst = Path(SI_tools.D_tmp) / f"frame_{i:04d}.png"
                    if not dst.exists():
                        os.symlink(src, dst)
                os.system(f"ffmpeg -y -r 2 -i {SI_tools.D_tmp}/frame_%04d.png -vf \"scale=iw:ih+mod(2-ih\\,2)\" -c:v libx264 -pix_fmt yuv420p {P_mp4}")
                video_paths[sim_name][region][var_name] = P_mp4

if __name__ == '__main__':
    from multiprocessing import freeze_support
    freeze_support()  # safe on all platforms, even if not frozen
    main()

