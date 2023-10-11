#!/bin/bash

for region in mawson olav sabrina; do
    for month in 01 02 03 04 05 06 07 08 09 10 11 12; do
        if [ -d "$region/$month" ]; then
        	montage \
		${region}/${month}/*baseline*.png \
            	${region}/${month}/*bathy_*.png \
            	${region}/${month}/*FI-noGI*.png \
            	${region}/${month}/*GI-10m*.png \
            	${region}/${month}/*GI-neg0p2m*.png \
            	-quality 100 -depth 16 -tile 2x -geometry 700x500! \
            	${region}/cice6_${month}_0veldays_${region}.png
        fi
    done
done

