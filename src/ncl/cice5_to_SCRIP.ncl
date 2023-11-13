begin
    start_time = get_cpu_time()
;---Input file
    srcFileName = "/g/data/jk72/da1339/grids/cice5_025_grid.nc"

;---Output (and input) files
    dstGridName = "/g/data/jk72/da1339/grids/CICE5_grid_SCRIP.20231107.nc"

;----------------------------------------------------------------------
; Convert original source NCEP grid to a SCRIP convention file.
;----------------------------------------------------------------------
    src_file = addfile(srcFileName,"r")
    lon      = src_file->tlon(:,:)
    lat      = src_file->tlat(:,:)
    printVarSummary(lat)
    printVarSummary(lon)
    lon2      = lon*(get_r2d(1))
    lat2      = lat*(get_r2d(1))
    printVarSummary(lat2)
    printVarSummary(lon2)
    Opt                = True
    Opt@ForceOverwrite = True
    Opt@PrintTimings   = True
    Opt@Title          = "CICE5 Grid"

    curvilinear_to_SCRIP(dstGridName,lat,lon,Opt)
    
;-- get the resulting CPU time
    end_time = get_cpu_time()
    cpu_time = end_time - start_time
    print("Elapsed time: "+ cpu_time + "s")
    
end