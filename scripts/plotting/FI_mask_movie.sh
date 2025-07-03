#!/bin/bash
gmt begin
    datetime=$(date --date=$(gmt math -T2008-08-01T00:00:00/2008-08-31T00:00:00/1d -o0 -qo${MOVIE_FRAME} T =) "+%Y%m%d")
    #gmt set FONT_ANNOT=white
    #gmt set MAP_FRAME_TYPE=inside
    #gmt set PS_PAGE_COLOR=black
    gmt grdimage @earth_day_30s -R160/180/-80/-70 -JM10c -Baf -BWSen -t30
    #gmt coast -R160/180/-80/-70 -Df -Gwhite -JM10c -Baf -BWSen
    gmt plot ~/AFIM_gv90/tmp/FI_mask_txt/FI_mask_${datetime}.txt -Sc0.15c -Gred
    #gmt grdimage ~/AFIM_gv90/tmp/FI_mask_frames_regular/FI_mask_${datetime}.nc -Ccmocean/amp -t0 -B
    #gmt grdimage ~/AFIM_gv90/tmp/FI_mask_frames/FI_mask_${datetime}.nc -Creds -t0 -B
    #gmt colorbar -Ccmocean/amp -DJBC+10c/0.6c+w-6c/0.5c -Baf+l"fast ice mask binary-days" -F+p+gwhite+s+r
    echo $datetime | gmt text -F+cTC -D0c/-1.5c
gmt end
