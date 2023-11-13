begin
;---Input file
    srcFileName = "/g/data/jk72/da1339/grids/ERA5_grid.nc"

;---Output (and input) files
    dstGridName = "/g/data/jk72/da1339/grids/ERA5_grid_SCRIP.20231107.nc"

;----------------------------------------------------------------------
; Convert original source NCEP grid to a SCRIP convention file.
;----------------------------------------------------------------------
    src_file = addfile(srcFileName,"r")
    lon     = src_file->lon(0,:)
    lat     = src_file->lat(:,0)
    print(lon)
    print(lat)
    Opt                = True
    Opt@ForceOverwrite = True
    Opt@PrintTimings   = True
    Opt@Title          = "ERA5 Grid"

    rectilinear_to_SCRIP(dstGridName,lat,lon,Opt)
    
end