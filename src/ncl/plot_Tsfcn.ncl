begin
  
    D_graph    = "/g/data/jk72/da1339/GRAPHICAL/CICE/IC/"
    F_ic       = "/scratch/jk72/da1339/cice-dirs/runs/afim_era5/restart/iced.2005-01-01-00000.AOM2_and_NCAR_with_o2i.nc"
    src_file   = addfile(F_ic, "r")
    Tsfc       = src_file->Tsfc(:, :)
    Tsfc@lat2d = src_file->TLAT
    Tsfc@lon2d = src_file->TLON
;--- FIG1: 
    wks                      = gsn_open_wks("png", D_graph+"Tsfcn_NCAR_nc0")
    plot                     = new(2,graphic)
    res                      = True            ; plot mods desire
    res@gsnDraw              = False           ; draw plot
    res@gsnFrame             = False           ; advance frame
    res@gsnPolar             = "NH"
    res@mpMinLatF            = 45.
    res@mpMaxLatF            = 90.
    res@mpFillOn             = True
    res@gsnLeftString        = "nc0"
    res@gsnCenterString      = ""
    res@gsnRightString       = ""
    res@gsnSpreadColors      = True            ; spread out color table
    res@gsnSpreadColorStart  = 2
    res@gsnSpreadColorEnd    = 11
    res@txFontHeightF        = 0.02            ; fontsize of the subtitles
    res@lbLabelBarOn         = False
    res@cnFillOn             = True            ; turn on color fill
    res@cnFillPalette       = "gui_default"
    res@cnLinesOn            = False           ; turn off contour lines
    res@cnInfoLabelOn        = False           ; turn off info label
    res@cnLineLabelsOn       = False           ; turn off contour labels
    res@cnLevelSelectionMode = "ManualLevels"  ; set manual contour levels
    res@cnMinLevelValF       = -2              ; set min contour level
    res@cnMaxLevelValF       = 1               ; set max contour level
    res@cnLevelSpacingF      = .1             ; set contour spacing
    plot(0)                  = gsn_csm_contour_map_polar(wks,Tsfc,res)
    delete([/res@gsnPolar,res@mpMinLatF,res@mpMaxLatF/])
    res@gsnPolar             = "SH"
    res@mpMinLatF            = -90.
    res@mpMaxLatF            = -55.
    plot(1)                  = gsn_csm_contour_map_polar(wks,Tsfc,res)
    pres                     = True
    pres@gsnMaximize         = True
    pres@gsnPanelLabelBar    = True
    pres@gsnPanelMainString  = "snow/ice surface temperature"
    pres@lbTitleOn           = True
    pres@lbTitleString       = "C"
    pres@lbTitleFontHeightF  = 0.018
    pres@lbTitlePosition     = "Bottom"
    pres@lbOrientation       = "horizontal"                    
    pres@lbBoxEndCapStyle    = "TriangleBothEnds"
    gsn_panel(wks,plot,(/1,2/), pres)
    delete([/wks,plot,res,pres,Tsfc/])

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    F_ic      = "/scratch/jk72/da1339/cice-dirs/runs/afim_era5/history/iceh_ic.2005-01-01-00600.nc"
    src_file  = addfile(F_ic, "r")
    Tsfc = src_file->Tsfc(0, :, :)
    Tsfc@lat2d = src_file->TLAT
    Tsfc@lon2d = src_file->TLON
;--- FIG1: NC(0)
    wks                      = gsn_open_wks("png", D_graph+"Tsfcn_CICE6_dt1_nc0")
    plot                     = new(2,graphic)
    res                      = True            ; plot mods desire
    res@gsnDraw              = False           ; draw plot
    res@gsnFrame             = False           ; advance frame
    res@gsnPolar             = "NH"
    res@mpMinLatF            = 45.
    res@mpMaxLatF            = 90.
    res@mpFillOn             = True
    res@gsnLeftString        = "nc0"
    res@gsnCenterString      = ""
    res@gsnRightString       = ""
    res@gsnSpreadColors      = True            ; spread out color table
    res@gsnSpreadColorStart  = 2
    res@gsnSpreadColorEnd    = 11
    res@txFontHeightF        = 0.02            ; fontsize of the subtitles
    res@lbLabelBarOn         = False
    res@cnFillOn             = True            ; turn on color fill
    res@cnFillPalette       = "gui_default"
    res@cnLinesOn            = False           ; turn off contour lines
    res@cnInfoLabelOn        = False           ; turn off info label
    res@cnLineLabelsOn       = False           ; turn off contour labels
    res@cnLevelSelectionMode = "ManualLevels"  ; set manual contour levels
    res@cnMinLevelValF       = -2             ; set min contour level
    res@cnMaxLevelValF       = 1             ; set max contour level
    res@cnLevelSpacingF      = .1             ; set contour spacing
    plot(0)                  = gsn_csm_contour_map_polar(wks,Tsfc,res)
    delete([/res@gsnPolar,res@mpMinLatF,res@mpMaxLatF/])
    res@gsnPolar             = "SH"
    res@mpMinLatF            = -90.
    res@mpMaxLatF            = -55.
    plot(1)                  = gsn_csm_contour_map_polar(wks,Tsfc,res)
    pres                     = True
    pres@gsnMaximize         = True
    pres@gsnPanelLabelBar    = True
    pres@gsnPanelMainString  = "snow/ice surface temperature"
    pres@lbTitleOn           = True
    pres@lbTitleString       = "C"
    pres@lbTitleFontHeightF  = 0.018
    pres@lbTitlePosition     = "Bottom"
    pres@lbOrientation       = "horizontal"                    
    pres@lbBoxEndCapStyle    = "TriangleBothEnds"
    gsn_panel(wks,plot,(/1,2/), pres)
    delete([/wks,plot,res,pres/])

end