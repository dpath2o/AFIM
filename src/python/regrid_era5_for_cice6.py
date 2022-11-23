import os
import afim
cice_prep = afim.cice_prep(os.path.join(os.path.expanduser('~'), 'src','python', 'afim_on_gadi.json'))
ice_prep.regrid_era5_for_cice6()

