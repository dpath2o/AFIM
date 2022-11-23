import os
import afim
cice_prep = afim.cice_prep( os.path.join(os.path.expanduser('~') , 'src' , 'python' , 'afim_on_gadi.json') )
cice_prep.acom2_coallate_for_cice_forcing()
