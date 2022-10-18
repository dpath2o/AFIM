import os
import afim
cice_prep = afim.cice_prep(os.path.join(os.path.expanduser('~'), 'src', 'python', 'afim_on_gadi.json'))
cice_prep.esmf_generate_weights(ds_name='BRAN',interp_type='patch',regrid_options='-p none -i')
