begin
  a       = addfile("/g/data/jk72/da1339/grids/CICE5_0p25_grid_in_degrees.nc","r")
  t       = a->tarea
  u       = a->uarea
  ulat    = a->ulat
  ulon    = a->ulon
  tlat    = a->tlat
  tlon    = a->tlon
  t@lat2d = tlat
  t@lon2d = tlon
  u@lat2d = ulat
  u@lon2d = ulon
  u_color = "NavyBlue"
  t_color = "orange"
  ; GRAPHICS
  wks                              = gsn_open_wks("png","CICE5_grid_northern_boundary")
  res                              = True
  res@gsnMaximize                  = True
  res@gsnDraw                      = False
  res@gsnFrame                     = False
  res@pmTickMarkDisplayMode        = "Always"   ; nicer map tickmarks
  res@mpDataBaseVersion            = "highres"
  res@mpMinLatF                    = 85
  res@mpMaxLatF                    = 90
  res@mpMinLonF                    = 175
  res@mpMaxLonF                    = 180
  res@tiMainString                 = "ACCESS-OM2 CICE5 Grid"
  res@gsnLeftString                = "U grid: " + u_color
  res@gsnRightString               = "T grid: " + t_color
  res@gsnStringFontHeightF         = 0.012
  res@gsnLeftStringOrthogonalPosF  = -0.01
  res@gsnLeftStringFontColor       = u_color
  res@gsnRightStringFontColor      = t_color
  res@gsnRightStringOrthogonalPosF = -0.01
  ; plot but don't draw
  plot                  = gsn_csm_map(wks,res)
  pres                  = True
  pres@gsnCoordsAttach  = True
  pres@gsLineThicknessF = 5.0
  pres@gsnCoordsAsLines = True 
  pres@gsLineColor      = t_color
  gsn_coordinates(wks,plot,t,pres)
  pres@gsLineColor = u_color
  gsn_coordinates(wks,plot,u,pres)
  draw(plot)
  frame(wks)
end