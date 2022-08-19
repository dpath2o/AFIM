

import os
import afim

cice_prep = afim.cice_prep('/Users/dpath2o/PHD/src/python/afim.json')
cice_prep.regrid_wrapper(ds_name='ERA5')


