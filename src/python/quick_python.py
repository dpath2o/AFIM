import xarray as xr
import numpy as np
import dask
import sys

def get_size(dataset, units='MB'):
    """Return the size of an object in specified units (default is MB)"""
    size_bytes = dataset.nbytes
    
    if units == 'KB':
        return f"{size_bytes / 1024:.2f} KB"
    elif units == 'MB':
        return f"{size_bytes / (1024**2):.2f} MB"
    elif units == 'GB':
        return f"{size_bytes / (1024**3):.2f} GB"
    else:
        return f"{size_bytes} bytes"

# Load the source data
src_file = "/g/data/ik11/inputs/access-om2/input_20200530/cice_025deg/grid.nc"
print(f"loading grid file: {src_file}")
ds = xr.open_dataset(src_file, chunks={'ny': 500, 'nx': 500})  # you might need to adjust the chunk sizes.
print(get_size(ds))
# Calculate grid_size from dimensions
ny, nx = ds.dims['ny'], ds.dims['nx']
grid_size = nx * ny

# Stack ny and nx dimensions into a new multi-index 'grid_size'
print("stacking")
ds_stacked = ds.stack(grid_size=('ny', 'nx'))

# Convert radians to degrees for lat/lon and assign to stacked dataset
print("converting radians to degrees")
ds_stacked['grid_center_lat'] = ds_stacked['tlat'] * 180.0 / np.pi
ds_stacked['grid_center_lon'] = ds_stacked['tlon'] * 180.0 / np.pi

# Create placeholder arrays for corner lat/lons in stacked dataset
print("placeholders")
ds_stacked['grid_corner_lat'] = ds_stacked['grid_center_lat'].expand_dims(grid_corners=4)
ds_stacked['grid_corner_lon'] = ds_stacked['grid_center_lon'].expand_dims(grid_corners=4)

# Add new dimensions for grid_corners and grid_rank to the original dataset
print("coordinates")
ds = ds.assign_coords(grid_corners=np.arange(4))
ds = ds.assign_coords(grid_rank=np.arange(2))

# Create grid_dims variable in the original dataset
print("grid dimensions")
ds['grid_dims'] = xr.DataArray([ny, nx], dims=['grid_rank'], coords={'grid_rank': [0, 1]})
ds['grid_dims'].attrs['long_name'] = "grid dimensions"

# Merge the original dataset and the stacked dataset
print("merging")
print("ds new size:",get_size(ds))
print("ds stack size:",get_size(ds_stacked))
merged_ds = xr.merge([ds, ds_stacked])

# Save to the destination file
print("writing out")
dst_file = "/home/581/da1339/DATA/grids/ERA5_grid_scrip.20231027.nc"
merged_ds.to_netcdf(dst_file)
print(f"Transformation complete. Output saved to {dst_file}")
