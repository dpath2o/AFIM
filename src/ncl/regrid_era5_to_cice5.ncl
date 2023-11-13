begin

;--- file names
  F_G_ERA5        = "/g/data/jk72/da1339/grids/ERA5_grid.nc"
  F_G_ERA5_SCRIP  = "/g/data/jk72/da1339/grids/ERA5_SCRIP_native.nc"
  F_G_CICE5       = "/g/data/jk72/da1339/grids/CICE5_0p25_grid_with_t_deg.nc"
  F_G_CICE5_SCRIP = "/g/data/jk72/da1339/grids/CICE5_SCRIP_0p25.nc"
  F_t2m           = "/g/data/rt52/era5/single-levels/reanalysis/2t/2005/2t_era5_oper_sfc_20050101-20050131.nc"
  F_wgt           = "/g/data/jk72/da1339/grids/weights/map_ERA5_to_CICE5_0p25_patch_extrap_neareststod.nc"

;--- read files and define destination variable names
  SRC              = addfile(F_t2m,"r")
  t2m              = SRC->t2m(0,:,:)
  time_units       = SRC->time@units
  printVarSummary(t2m)
  src_lat          = SRC->latitude
  src_lon          = SRC->longitude
  DST              = addfile(F_G_CICE5,"r")
  dst_lat          = DST->lat
  dst_lon          = DST->lon
  dst_dims         = dimsizes(dst_lat)
  dst_nlat         = dst_dims(0)
  dst_nlon         = dst_dims(1)

;--- set ESMF resources for rectilinear to SCRIP creation
  ; Opt                = True
  ; Opt@GridMask       = where(.not.ismissing(t2m),1,0)
  ; Opt@ForceOverwrite = True
  ; Opt@PrintTimings   = True
  ; Opt@Title          = "ERA5 Rectilinear SCRIP"
  ; rectilinear_to_SCRIP(F_G_ERA5_SCRIP,src_lat,src_lon,Opt)
  ; delete([Opt])
;--- set ESMF resources for curvilinear to SCRIP creation
  ; Opt                = True
  ; Opt@ForceOverwrite = True
  ; Opt@PrintTimings   = True
  ; Opt@Title          = "CICE5 Curvilinear SCRIP"
  ; curvilinear_to_SCRIP(F_G_CICE5_SCRIP,dst_lat,dst_lon,Opt)
  ; delete([Opt])

;--- set ESMF resources for curvilinear to SCRIP creation
  Opt                = True
  Opt@InterpMethod   = "neareststod"
  Opt@SrcFileName    = F_G_ERA5_SCRIP
  Opt@DstFileName    = F_G_CICE5_SCRIP
  Opt@WgtFileName    = F_wgt
  Opt@ForceOverwrite = True
  Opt@DstMask2D      = where(ismissing(t2m),0,1)
  Opt@DstGridType    = "curvilinear"
  Opt@DstTitle       = "CICE5"
  Opt@DstGridLon     = dst_lon
  Opt@DstGridLat     = dst_lat
  ;ESMF_regrid_gen_weights(F_G_ERA5_SCRIP,F_G_CICE5_SCRIP,F_wgt,Opt)

;--- regrid with weights
  ;; t2m_reG = ESMF_regrid_with_weights(t2m,F_wgt,Opt)
  ;; printVarSummary(t2m_reG)

;--- RE-DEFINE t2m variable for looping over
  delete(t2m)
  t2m     = SRC->t2m(:, :, :)  ; Access all time steps
  ;; lat     = SRC->latitude
  ;; lon     = SRC->longitude
colors = read_colormap_file("BlWhRe")
nscan = 2000             ; usually *much* fewer
eps   = 0.001            ; variable depended
gtype = False            ; "gendat" does not generate cyclic fields
guess = 0                ; use zonal means
relc  = 0.6              ; standard relaxation coef
opt   = 0
; Loop over each time step
  do i = 0, dimsizes(t2m(:, 0, 0)) - 1
    t2m_i       = t2m(i, :, :)   ; Extract data for a specific time step
    t2m_i      := t2m_i - 273.15
    t2m_i@units = "degC"
    t2m_i!0     = "lat"
    t2m_i!1     = "lon"
    t2m_i&lat   = src_lat
    t2m_i&lon   = src_lon
    print("RE-gridding for time index {i}")
    t2m_reG = ESMF_regrid_with_weights(t2m_i,F_wgt,Opt)
    dst_lat = t2m_reG&lat
    dst_lon = t2m_reg&lon
    print("regridded:")
    printVarSummary(t2m_reG)
;******************************************************************
; set the poisson_grid_fill arguments
;******************************************************************
    zc = t2m_reG          ; save for demo plot
    poisson_grid_fill(zc, gtype, guess, nscan, eps, relc, opt)
;----------------------------------------------------------------------
; Plot the original and regridded data on a panel plot
;----------------------------------------------------------------------
    t_           = SRC->time(i)
    t_@units     = time_units
    time_date    = cd_calendar(t_, 0)
    wks          = gsn_open_wks("png", "t2m_reG_testing_"+time_date(0,0)+"_"+time_date(0,1)+"_"+time_date(0,2)+"_"+time_date(0,3)+"00")
    ;; gsn_define_colormap(wks,"gui_default")      ; choose colormap
    plot         = new(6,graphic)
  
;--- Plot 1
    res                      = True
    res@gsnMaximize          = True
    res@gsnDraw              = False
    res@gsnFrame             = False
    res@gsnPolar             = "NH"
    res@mpMinLatF            = 50
    res@mpOutlineOn          = True
    res@mpGridAndLimbOn      = False
    res@cnFillOn             = True
    res@cnLinesOn            = False
    res@cnLineLabelsOn       = False
    res@cnFillPalette        = colors(8:,:)   ; set color map
    res@cnLevelSelectionMode =  "ManualLevels"   ; manual levels
    res@cnMinLevelValF       = -50                 ; min level
    res@cnMaxLevelValF       = 20                ; max level
    res@cnLevelSpacingF      =  2                ; interval
    res@lbLabelBarOn         = False              ; Labelbar will be in panel
    res@tiMainString         = "Original t2m"
    cres@gsnCoordsLat          = src_lat(::4,::4)
    cres@gsnCoordsLon          = src_lon(::4,::4)
    cres@gsnCoordsAsLines      = True
    cres@gsLineThicknessF      = 1
    cres@gsnCoordsMissingColor = "red"
    cres@gsMarkerSizeF         = 1
    plot(0) = gsn_csm_contour_map_polar(wks, t2m_i, res)
    gsn_coordinates(wks,plot(0),t2m_i,cres)
    delete([/res@tiMainString,cres@gsnCoordsLat,cres@gsnCoordsLon/])

;--- Plot 2
    res@tiMainString = "Regridded t2m"
    cres@gsnCoordsLat          = dst_lat(::4,::4)
    cres@gsnCoordsLon          = dst_lon(::4,::4)
    plot(1) = gsn_csm_contour_map_polar(wks, t2m_reG, res)
    gsn_coordinates(wks,plot(1),t2m_reG,cres)
    delete([/res@tiMainString,cres@gsnCoordsLat,cres@gsnCoordsLon/])

;--- Plot 3
    res@tiMainString = "Regridded t2m -- FILLED"
    plot(2) = gsn_csm_contour_map_polar(wks, zc, res)
    delete([/res@tiMainString,res@gsnPolar,res@mpMinLatF/])

;--- Plot 4
    res@gsnPolar     = "SH"
    res@mpMinLatF    = 50
    res@tiMainString = "Original t2m"
    plot(3) = gsn_csm_contour_map_polar(wks, t2m_i, res)
    delete([/res@tiMainString/])

;--- Plot 5
    res@tiMainString = "Regridded t2m"
    plot(4) = gsn_csm_contour_map_polar(wks, t2m_reG, res)
    delete([/res@tiMainString/])

;--- Plot 5
    res@tiMainString = "Regridded t2m -- FILLED"
    plot(4) = gsn_csm_contour_map_polar(wks, zc, res)

; Panel the plots
    pres                    = True
    ; pres@gsnMaximize       = True
    ; pres@gsnPanelLabelBar  = True
    pres@gsnPanelMainSTring = time_date(0,0)+"_"+time_date(0,1)+"_"+time_date(0,2)+"_"+time_date(0,3)+"00"
    gsn_panel(wks,plot,(/2,3/), pres)
    delete([/t2m_i,t2m_reG,zc,t_,wks,res,pres/])
  end do

end
