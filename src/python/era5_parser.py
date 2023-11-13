import sys
import os
import subprocess
import logging

def process_file(filepath, variable):
    fname = os.path.basename(filepath)
    Dout = f"/g/data/jk72/da1339/afim_input/ERA5/sfc_vars_sh/{variable}"
    if not os.path.exists(Dout):
        os.makedirs(Dout)
    dest_file = f"{Dout}/{fname}"
    if not os.path.exists(dest_file):
        logging.info(f"Started processing {dest_file}")
        subprocess.run(["ncks", "-O", "-d", "latitude,360,720,1", filepath, dest_file])
        logging.info(f"Completed processing {dest_file}")

if __name__ == "__main__":
    logging.basicConfig(filename='/home/581/da1339/logs/era5_parser.log', level=logging.INFO, format='%(asctime)s - %(message)s')
    variable = sys.argv[1]
    files = sys.argv[2:]
    for file in files:
        process_file(file, variable)
