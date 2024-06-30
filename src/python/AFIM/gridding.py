import xarray as xr
#import xesmf  as xe

############################################################################
def xesmf_regrid_dataset(DS_src, DS_dst, F_DS_na,
                         method='bilinear',
                         periodic=True,
                         reuse_weights=True,
                         F_weights='',
                         multi_file=False):
    '''
    Regrid a given dataset from source grid to destination grid using xESMF.
    
    Parameters:
    DS_src (xarray.Dataset) : Source dataset
    DS_dst (xarray.Dataset) : Destination dataset
    F_DS_na (str)           : File name or path of the dataset to be regridded
    method (str)            : Regridding method; default is 'bilinear'
    periodic (bool)         : Whether the grid is periodic; default is True
    reuse_weights (bool)    : Whether to reuse weights; default is True
    F_weights (str)         : File name for saving/loading regridding weights; default is ''
    multi_file (bool)       : Whether to open multiple files; default is False
    
    Returns:
    xarray.Dataset : Regridded dataset
    '''
    print("regridding file: ",F_DS_na)
    rg = xe.Regridder(DS_src, DS_dst, 
                      method=method, 
                      periodic=periodic, 
                      filename=F_weights, 
                      reuse_weights=reuse_weights)
    if multi_file:
        DS_na = xr.open_mfdataset(F_DS_na, parallel=True, chunks={'time':1})
    else:
        DS_na = xr.open_dataset(F_DS_na)
    DS_rg = rg(DS_na)
    return DS_rg
