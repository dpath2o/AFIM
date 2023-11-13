begin

;--- file names
  F_G_ERA5 = "/g/data/jk72/da1339/grids/ERA5_grid.nc"
  F_orig   = "/g/data/jk72/da1339/afim_input/ERA5/t2m_regrid_comparison/t2m_original.nc"
  F_nco    = "/g/data/jk72/da1339/afim_input/ERA5/t2m_regrid_comparison/airtmp_nco.nc"
  F_xesmf  = "/g/data/jk72/da1339/afim_input/ERA5/t2m_regrid_comparison/airtmp_xesmf.nc"
 
;--- read files and define destination variable names
  SRC        = addfile(F_orig,"r")
  NCO        = addfile(F_nco,"r")
  XESMF      = addfile(F_xesmf,"r")
  G_ERA5     = addfile(F_G_ERA5,"r")
  time_units = SRC->time@units
  t2m_orig   = SRC->t2m(:,:,:)  ; Access all time steps
  t2m_nco    = NCO->airtmp(:,:,:)
  t2m_xesmf  = XESMF->airtmp(:,:,:)
  lat_o      = G_ERA5->lat
  lon_o      = G_ERA5->lon
  lat_n      = NCO->lat
  lon_n      = NCO->lon
  lat_x      = XESMF->lat
  lon_x      = XESMF->lon

; Loop over each time step
  do i = 0, dimsizes(t2m_orig(:, 0, 0)) - 1
    t2m_o       = t2m_orig(i,:,:)   ; Extract data for a specific time step
    t2m_n       = t2m_nco(i,:,:)
    t2m_x       = t2m_xesmf(i,:,:)
    t2m_o      := t2m_o - 273.15
    t2m_n      := t2m_n - 273.15
    t2m_x      := t2m_x - 273.15
    t2m_o@units = "degC"
    t2m_n@units = "degC"
    t2m_x@units = "degC"
    t2m_o@lat2d = lat_o
    t2m_o@lon2d = lon_o
    t2m_n@lat2d = lat_n
    t2m_n@lon2d = lon_n
    t2m_x@lat2d = lat_x
    t2m_x@lon2d = lon_x
    t_          = SRC->time(i)
    t_@units    = time_units
    time_date   = cd_calendar(t_, 0)
    yr          = time_date(0, 0)
    month       = time_date(0, 1)
    mo          = sprintf("%02.0f", month)
    day         = time_date(0, 2)
    dy          = sprintf("%02.0f", day)
    hour        = time_date(0, 3)
    hr          = sprintf("%02.0f", hour)
    new_dt_label= yr+"_"+mo+"_"+dy+"_"+hr+"00"
    print(new_dt_label)
    wks         = gsn_open_wks("png", "/g/data/jk72/da1339/GRAPHICAL/ERA5/regrid_comparison/t2m_reG_"+new_dt_label)
    plot        = new(6,graphic)
  
;--- Plot 1
    res                      = True
    res@gsnMaximize          = True
    res@gsnDraw              = False
    res@gsnFrame             = False
    res@gsnPolar             = "NH"
    res@mpMinLatF            = 50
    res@mpOutlineOn          = True
    res@cnFillOn             = True
    res@cnLineLabelsOn       = False
    res@tiMainString         = "Original t2m"
    plot(0) = gsn_csm_contour_map_polar(wks, t2m_o, res)
    delete([/res@tiMainString/])

;--- Plot 2
    res@tiMainString = "Regridded t2m -- NCO method"
    plot(1) = gsn_csm_contour_map_polar(wks, t2m_n, res)
    delete([/res@tiMainString/])

;--- Plot 3
    res@tiMainString = "Regridded t2m -- XESMF method"
    plot(2) = gsn_csm_contour_map_polar(wks, t2m_x, res)
    delete([/res/])

;--- Plot 4
    res                      = True
    res@gsnMaximize          = True
    res@gsnDraw              = False
    res@gsnFrame             = False
    res@gsnPolar             = "SH"
    res@mpMaxLatF            = 50
    res@mpOutlineOn          = True
    res@cnFillOn             = True
    res@cnLineLabelsOn       = False
    res@tiMainString         = "Original t2m"
    plot(3) = gsn_csm_contour_map_polar(wks, t2m_o, res)
    delete([/res@tiMainString/])

;--- Plot 5
    res@tiMainString = "Regridded t2m -- NCO method"
    plot(4) = gsn_csm_contour_map_polar(wks, t2m_n, res)
    delete([/res@tiMainString/])

;--- Plot 5
    res@tiMainString = "Regridded t2m -- XESMF method"
    plot(5) = gsn_csm_contour_map_polar(wks, t2m_x, res)

; Panel the plots
    pres                    = True
    ; pres@gsnMaximize       = True
    pres@gsnPanelLabelBar   = False
    pres@gsnPanelMainSTring = new_dt_label
    gsn_panel(wks,plot,(/2,3/), pres)
    draw(plot)
    delete([/t2m_o,t2m_n,t2m_x,t_,wks,res,pres/])
  end do

end
