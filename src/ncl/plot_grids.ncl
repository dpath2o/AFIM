begin

  F_G_ERA5        = "/g/data/jk72/da1339/grids/ERA5_grid.nc"
  F_G_ERA5_SCRIP  = "/g/data/jk72/da1339/grids/ERA5_grid_rectilinear_SCRIP.nc"
  F_G_CICE5       = "/g/data/jk72/da1339/grids/CICE5_0p25_grid_with_t_deg.nc"
  F_G_CICE5_SCRIP = "/g/data/jk72/da1339/grids/CICE5_grid_curvilinear_SCRIP.nc"
  G               = addfile(F_G_ERA5,"r")
  lat             = G->lat
  lon             = G->lon
  ;; F_t2m   = addfile("/g/data/rt52/era5/single-levels/reanalysis/2t/2005/2t_era5_oper_sfc_20050101-20050131.nc","r")
  ;; t2m     = F_t2m->t2m(0,:,:)
  ;; t2m    := t2m - 273.15
  ;; t2m@units = "degC"
  ;; lat     = F_t2m->latitude
  ;; lon     = F_t2m->longitude
  ;; t2m!0   = "lat"
  ;; t2m!1   = "lon"
  ;; t2m&lat = lat
  ;; t2m&lon = lon
  ;; printVarSummary(t2m)

  ;-- initialise figure
  a4_height                 = 29.7
  a4_width                  = 23.0
  cm_per_inch               = 2.54
  wks                       = gsn_open_wks("png","ERA5_grid_rectilinear")
  res                       = True
  res@gsnMaximize           = True
  res@gsnPolar              = "NH"
  res@gsnDraw               = False   ; Don't draw plot
  res@gsnFrame              = False   ; Don't advance frame.
  res@mpDataBaseVersion     = "HighRes"
  res@mpOutlineBoundarySets = "Geophysical"   ; "AllBoundaries"
  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = 80
  res@tiMainString          = "ERA5 rectilinear grid"
  plot                      = gsn_csm_map_polar(wks,res)

;----------------------------------------------------------------------
; Draw lines at the lon/lon coordinate array values.
;----------------------------------------------------------------------
  pres                   = True
  pres@gsnCoordsAsLines  = True
  pres@gsnCoordsLat      = lat(::2,::2)
  pres@gsnCoordsLon      = lon(::2,::2)
  pres@gsLineThicknessF  = 1
  gsn_coordinates(wks,plot,lat,pres)
  delete(pres@gsnCoordsAsLines)

;----------------------------------------------------------------------
; Draw the lat/lon grid with markers, using black for grid locations
; where the data is not missing, and red for grid locations
; where the data is missing.
;---------------------------------------------------------------------
  pres@gsnCoordsNonMissingColor = "black"
  pres@gsnCoordsMissingColor    = "red"
  pres@gsMarkerSizeF            = 1     ; Default is a little large
  gsn_coordinates(wks,plot,lat,pres)
  delete([/wks,res,pres/])

end
