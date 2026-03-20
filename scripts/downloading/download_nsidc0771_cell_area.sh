#!/bin/bash
# download_nsidc0771_cell_area.sh
set -euo pipefail

D_AREA="/g/data/gv90/da1339/SeaIce/NSIDC/NSIDC0771"
U_AREA="https://daacdata.apps.nsidc.org/pub/DATASETS/nsidc0771_polarstereo_anc_grid_info/"

mkdir -p "${D_AREA}"
cd "${D_AREA}"

# ~/.netrc must contain:
# machine urs.earthdata.nasa.gov login <uid> password <password>
# chmod 600 ~/.netrc

wget \
    --load-cookies "${HOME}/.urs_cookies" \
    --save-cookies "${HOME}/.urs_cookies" \
    --keep-session-cookies \
    --no-check-certificate \
    --auth-no-challenge=on \
    -r -l 3 -np -nd \
    -A "NSIDC0771_CellArea_PS_N25km_v1.0.nc,NSIDC0771_CellArea_PS_S25km_v1.0.nc" \
    "${U_AREA}"
