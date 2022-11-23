import os
import afim
#from dask.distributed import Client
#client = Client(n_workers=24, threads_per_worker=2, memory_limit='765GB')
cice_prep = afim.cice_prep(os.path.join(os.path.expanduser('~'), 'src','python', 'afim_on_gadi.json'))
cice_prep.regrid_bran_for_cice6()
