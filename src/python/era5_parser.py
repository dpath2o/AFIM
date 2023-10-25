import sys
import os
import subprocess

def process_file(filepath, variable):
    fname = os.path.basename(filepath)
    Dout = f"/g/data/jk72/da1339/afim_input/ERA5/sfc_vars_sh/{variable}"
    if not os.path.exists(Dout): os.makedirs(Dout)
    dest_file = f"{Dout}/{fname}"
    if not os.path.exists(dest_file):
        subprocess.run(["ncks", "-O", "-d", "latitude,360,720,1", filepath, dest_file]) #"-d", "longitude,1120,1399,1", 

if __name__ == "__main__":
    variable = sys.argv[1]
    files = sys.argv[2:]
    for file in files:
        process_file(file, variable)
