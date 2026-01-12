#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compute Liu et al. (2022) coastal-drag form factors F2 from a high-res coastline shapefile.

Implements (per T-cell):
  F2x = sum(|l cos(theta)|) / dx
  F2y = sum(|l sin(theta)|) / dy
where theta is segment angle relative to model x-axis (i-direction).

Outputs NetCDF with F2x_T, F2y_T (+ optional F2u/F2v via averaging).
"""
import argparse, sys
from pathlib import Path
mod_path = '/home/581/da1339/AFIM/src/AFIM/src'
sys.path.insert(0, mod_path)
# Your local project import pattern may differ on Gadi; keep consistent with how you run other AFIM scripts
from sea_ice_toolbox import SeaIceToolbox, SeaIceToolboxManager

def main():
    p = argparse.ArgumentParser(description="Build Liu et al. (2022) F2 form factors from a high-resolution coastline shapefile.")
    p.add_argument("--sim_name", type=str, default="elps-min",
                   help="Simulation name used by SeaIceToolbox config initialisation (grid/coast are overridden).")
    p.add_argument("--P_json", type=str, default=None,
                   help="Path to sea_ice_config.json (defaults to SeaIceToolbox internal default if omitted).")
    p.add_argument("--P_CICE_grid", type=str, required=True,
                   help="Full path to the CICE C-grid NetCDF (overrides CICE_dict['P_G']).")
    p.add_argument("--P_high_res_coast", type=str, default=None,
                   help="Path to high-res coastline shapefile (defaults to CICE_dict['P_high_res_coast']).")
    p.add_argument("--P_out", type=str, default=None,
                   help="Output NetCDF path (defaults to CICE_dict['P_F2_coast']).")
    p.add_argument("--proj_crs", type=str, default="EPSG:3031",
                   help="Projection for KDTree mapping of segments to cells (default: EPSG:3031).")
    p.add_argument("--chunk_segments", type=int, default=2_000_000,
                   help="Segment chunk size for vectorised processing.")
    p.add_argument("--coast_write_stride", type=int, default=25,
                   help="Stride for writing coastline lon/lat to NetCDF (provenance only).")
    p.add_argument("--lat_subset_max", type=float, default=-30.0,
                   help="Only build KDTree using grid cells with tlat <= this value (deg).")
    p.add_argument("--P_log", type=str, default=None,
                   help="Optional log file path for SeaIceToolbox.")
    args = p.parse_args()
    # Initialise toolbox
    SI_mgr = SeaIceToolboxManager(P_log = args.P_log)
    SI     = SI_mgr.get_toolbox(sim_name    = args.sim_name,
                                P_json      = args.P_json,
                                P_CICE_grid = args.P_CICE_grid)
    # Optional overrides from CLI
    if args.P_high_res_coast is not None:
        SI.CICE_dict["P_high_res_coast"] = args.P_high_res_coast
    if args.P_out is not None:
        SI.CICE_dict["P_F2_coast"] = args.P_out
    SI.build_F2_form_factors_from_high_res_coast(proj_crs           = args.proj_crs,
                                                 chunk_segments     = args.chunk_segments,
                                                 coast_write_stride = args.coast_write_stride,
                                                 lat_subset_max     = args.lat_subset_max)

if __name__ == "__main__":
    main()

