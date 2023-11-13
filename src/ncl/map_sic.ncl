begin
  D_data = "/scratch/jk72/da1339/cice-dirs/runs/afim_era5/history/daily/"
  F_spat = "iceh.*-*-*.nc"
  F_list = systemfunc("ls "+D_data+F_spat)
  do i = 0, dimsizes(F_list) - 1
    F_ = F_list(i)
    print("Attempting to map 'aice' from: "+F_)
    in        = addfile(F_, "r")
    t_        = in->time
    t_       := t_ - 1
    t_@units  = in->time@units
    time_date = cd_calendar(t_, 0)
    yr        = time_date(0, 0)
    month     = time_date(0, 1)
    mo        = sprintf("%02.0f", month)
    day       = time_date(0, 2)
    dy        = sprintf("%02.0f", day)
    hour      = time_date(0, 3)
    hr        = sprintf("%02.0f", hour)
    d_str     = yr+"_"+mo+"_"+dy
    F_out     = "/g/data/jk72/da1339/GRAPHICAL/CICE/afim_era5/reG_patch/aice_"+d_str
    F_out_str = F_out+".png"
    if (fileexists(F_out_str)) then
      print("File exists: "+F_out_str+" ... skipping")
      continue
    end if
    print("creating: "+F_out_str)
    var         = in->aice(0,:,:)
    TLAT        = in->TLAT
    TLON        = in->TLON
    var         = var * 100.0
    var@lat2d   = TLAT
    var@lon2d   = TLON
    print("min/max values:")
    printMinMax(var, False)
    print("opening the figure")
    wks  = gsn_open_wks("png", F_out) 
    plot = new(2,graphic)
    print("defining a colormap")
    gsn_define_colormap(wks, "BlueDarkRed18")  
    print("defining northern hemisphere resources")
    mpres           = True
    mpres@gsnDraw   = False
    mpres@gsnFrame  = False
    mpres@gsnPolar  = "NH"
    mpres@mpMinLatF = 45.
    mpres@mpMaxLatF = 90.
    mpres@mpFillOn  = True
    print("plotting northern hemisphere basemap")
    map1 = gsn_csm_map_polar(wks,mpres)
    ; filled contours
    res                      = True            ; plot mods desire
    res@gsnDraw              = False           ; draw plot
    res@gsnFrame             = False           ; advance frame
    res@gsnLeftString        = "AICE (Line: 15%)"
    res@gsnCenterString      = ""
    res@gsnRightString       = ""
    res@gsnSpreadColors      = True            ; spread out color table
    res@gsnSpreadColorStart  = 2
    res@gsnSpreadColorEnd    = 11
    res@txFontHeightF        = 0.02  ; fontsize of the subtitles
    res@lbLabelBarOn         = False
    res@cnFillOn             = True            ; turn on color fill
    res@cnLinesOn            = False           ; turn off contour lines
    res@cnInfoLabelOn        = False           ; turn off info label
    res@cnLineLabelsOn       = False           ; turn off contour labels
    res@cnLevelSelectionMode = "ManualLevels"  ; set manual contour levels
    res@cnMinLevelValF       = 10.             ; set min contour level
    res@cnMaxLevelValF       = 90.             ; set max contour level
    res@cnLevelSpacingF      = 10.             ; set contour spacing
    print("plotting northern hemisphere filled contours")
    plotA = gsn_csm_contour(wks,var,res)
    ; 15% concentration contour
    res2                      = True   ; plot mods desire
    res2@gsnDraw              = False  ; draw plot
    res2@gsnFrame             = False  ; advance frame
    res2@cnFillOn             = False  ; turn on color fill
    res2@cnLinesOn            = True   ; turn on contour lines
    res2@cnLevelSelectionMode = "ExplicitLevels"  ; set explicit contour levels
    res2@cnLineColor          = "red"
    res2@cnLineThicknessF     = 3
    res2@cnLevels             = (/15./)
    print("plotting northern hemisphere contour lines")
    plotB = gsn_csm_contour(wks,var,res2)
    print("overlaying the northern hemisphere plots onto the basemap")
    overlay(map1,plotA)
    overlay(map1,plotB)
    plot(0) = map1
    delete(mpres)
    print("defining southern hemisphere resources")
    mpres           = True
    mpres@gsnDraw   = False
    mpres@gsnFrame  = False
    mpres@gsnPolar  = "SH"
    mpres@mpMinLatF = -90.
    mpres@mpMaxLatF = -55.
    mpres@mpFillOn  = True
    print("plotting southern hemisphere basemap")
    map2 = gsn_csm_map_polar(wks,mpres)
    print("plotting southern hemisphere filled contours")
    plotC = gsn_csm_contour(wks,var,res)
    print("plotting southern hemisphere contour lines")
    plotD = gsn_csm_contour(wks,var,res2)
    print("overlaying the southern hemisphere contour lines")
    overlay(map2,plotC)
    overlay(map2,plotD)
    plot(1) = map2
    print("creating the panel resources")
    pres                    = True
    pres@gsnMaximize        = True
    pres@gsnPanelLabelBar   = True
    pres@gsnPanelMainString = "CICE6-SA-ERA5 " + d_str 
    pres@lbTitleOn          = True
    pres@lbTitleString      = "(%)"
    pres@lbTitleFontHeightF = 0.018
    pres@lbTitlePosition    = "Bottom"
    pres@lbOrientation      = "horizontal"                    
    pres@lbBoxEndCapStyle   = "TriangleBothEnds"              
    pres@lbLabelStride      = 1
    print("plotting two hemispheres on one figure using 'gsn_panel'")
    gsn_panel(wks,plot,(/1,2/), pres)
    frame(wks)
    delete([/mpres,res,res2/])
    print("")
  end do
end
