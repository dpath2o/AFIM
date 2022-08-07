import xarray as xr
import rioxarray
import pygmt
import os
F1 = os.path.join('/','Users','dpath2o','DATA','Satellite','Sentinel','CapeDarnley','EW','2015',
                  'S1A_EW_GRDM_1SSH_20150321T152013_20150321T152101_005131_006765_4867.SAFE',
                  'measurement',
                  's1a-ew-grd-hh-20150321t152013-20150321t152101-005131-006765-001.tiff')
ras1           = rasterio.open(F1)
gcps1,ras_crs1 = ras1.get_gcps()
da1            = xr.open_rasterio(F1)
ds1            = da1.to_dataset('band')
ds1            = ds1.rename({1:'sar'})
destination1   = ds1.sar.values
dest1,trans1   = rasterio.warp.reproject(ras1.read(1),
                                         destination1,
                                         gcps=gcps1,
                                         src_crs=ras_crs1,
                                         dst_crs=ras_crs1)
dn1            = xr.DataArray(dest1, dims=['y','x'])
dn1            = dn1.rio.write_crs(ras_crs1)
dn1            = dn1.rio.write_transform(trans1)
dn1            = dn1.assign_coords( x=(trans1*(dn1.x,dn1.y[0]))[0],
                                    y=(trans1*(dn1.x[0],dn1.y))[1] )
dn1            = dn1.astype(float)
region         = [dn1.x.min().values,dn1.x.max().values,dn1.y.min().values,dn1.y.max().values]
fig            = pygmt.Figure()
fig.basemap( region=region , projection="M10i" , frame=["af",'+t"S1A SAR GRDM 1SSH 2015-03-21 1520Z"'] )
pygmt.makecpt( cmap="gray" , series=(0,dn1.max().values,1) )
fig.grdimage( dn1 , cmap=True , transparency=20 )
fig.coast(shorelines=".5p,blue")
fig.savefig(os.path.join('/','Users','dpath2o','PHD','GRAPHICAL','grounded_icebergs','unmerged_tiffs',
                         'S1A_EW_GRDM_1SSH_20150321T152013_20150321T152101_005131_006765_4867_EC.png'))
F2 = os.path.join('/','Users','dpath2o','DATA','Satellite','Sentinel','CapeDarnley','EW','2015',
                  'S1A_EW_GRDM_1SSH_20150402T152013_20150402T152101_005306_006B78_0527.SAFE',
                  'measurement',
                  's1a-ew-grd-hh-20150402t152013-20150402t152101-005306-006b78-001.tiff')
ras2           = rasterio.open(F2)
gcps2,ras_crs2 = ras2.get_gcps()
da2            = xr.open_rasterio(F2)
ds2            = da2.to_dataset('band')
ds2            = ds2.rename({1:'sar'})
destination2   = ds2.sar.values
dest2,trans2   = rasterio.warp.reproject(ras2.read(1),
                                         destination2,
                                         gcps=gcps2,
                                         src_crs=ras_crs2,
                                         dst_crs=ras_crs2)
dn2            = xr.DataArray(dest2, dims=['y','x'])
dn2            = dn2.rio.write_crs(ras_crs2)
dn2            = dn2.rio.write_transform(trans2)
dn2            = dn2.assign_coords( x=(trans2*(dn2.x,dn2.y[0]))[0],
                                    y=(trans2*(dn2.x[0],dn2.y))[1] )
dn2            = dn2.astype(float)
region         = [dn2.x.min().values,dn2.x.max().values,dn2.y.min().values,dn2.y.max().values]
fig            = pygmt.Figure()
fig.basemap( region=region , projection="M10i" , frame=["af",'+t"S1A SAR GRDM 1SSH 2015-04-02 1520Z"'] )
pygmt.makecpt( cmap="gray" , series=(0,dn2.max().values,1) )
fig.grdimage( dn2 , cmap=True , transparency=20 )
fig.coast(shorelines=".5p,blue")
fig.savefig(os.path.join('/','Users','dpath2o','PHD','GRAPHICAL','grounded_icebergs','unmerged_tiffs',
                         'S1A_EW_GRDM_1SSH_20150402T152013_20150402T152101_005306_006B78_0527.png'))

tmp2 = dn2.isel(x=slice(0,10360),y=slice(0,8145))
tmp1 = dn1.isel(x=slice(0,10360),y=slice(0,8145))
dd   = xr.DataArray(tmp2.values-tmp1.values, dims=['y','x'])
fig  = pygmt.Figure()
fig.basemap( region=region , projection="M10i" , frame=["af",'+t"S1A Difference (2105-04-02T1520Z - 2015-03-21T1520Z)"'] )
pygmt.makecpt( cmap="gray" , series=(0,dd.max().values,1) )
fig.grdimage( dd , cmap=True , transparency=20 )
fig.coast(shorelines=".5p,blue")
fig.savefig(os.path.join('/','Users','dpath2o','PHD','GRAPHICAL','grounded_icebergs',
                         'CapeDarnley_S1A_EW_GRDM_1SSH_20150402_diff_20150321.png'))
