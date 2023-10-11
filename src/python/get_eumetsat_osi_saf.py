import requests
import netCDF4 as nc

# Define the OpenDAP URL
url = "https://thredds.met.no/thredds/dodsC/osisaf/met.no/ice/drift_lr_sh_agg.nc"

# Specify the variables you want to download
variables = [
    "Polar_Stereographic_Grid",
    "xc[0:1:124]",
    "yc[0:1:130]",
    "lat[0:1:0][0:1:0]",
    "lon[0:1:0][0:1:0]",
    "time[0:1:3837]",
    "time_bnds[0:1:3837][0:1:1]",
    "dt0[0:1:0][0:1:0][0:1:0]",
    "lon1[0:1:0][0:1:0][0:1:0]",
    "lat1[0:1:0][0:1:0][0:1:0]",
    "dt1[0:1:0][0:1:0][0:1:0]",
    "dX[0:1:0][0:1:0][0:1:0]",
    "dY[0:1:0][0:1:0][0:1:0]",
    "status_flag[0:1:0][0:1:0][0:1:0]",
    "uncert_dX_and_dY[0:1:0][0:1:0][0:1:0]",
]

# Create the OpenDAP URL with variables
dap_url = f"{url}{'?'.join(variables)}"

# Download the data
response = requests.get(dap_url)

# Save the data to a local NetCDF file
with open("drift_lr_sh_agg.nc", "wb") as f:
    f.write(response.content)

print("Data downloaded successfully.")

