import os
import afim

cice_prep = afim.cice_prep('/Users/dpath2o/PHD/src/python/afim.json')
# do ERA5 regridding and derivations for forcing
#cice_prep.regrid_wrapper(ds_name='ERA5')
# do BRAN regridding and derivations for forcing
cice_prep.regrid_wrapper(ds_name='BRAN')


