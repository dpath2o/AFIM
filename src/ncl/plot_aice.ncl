;*************************************************
; NCL Graphics: ice_4.ncl
;
; Concepts illustrated:
;    - Plotting ice data
;    - Filling gap between poles in CICE T-fold Tripole grid correctly
;     
; See http://oceans11.lanl.gov/trac/CICE; in particular the
; discussion of the T-fold Tripole GRID in the 
; CICE Documentation and User's Guide (PDF)
;
; In simplistic terms, the T-fold grid has an apparent gap that 
; manifests as a line between the two northern poles. The U-fold grid
; does not have the gap but since this data is aligned to the T-fold
; grid, it is displaced slightly from its correct position if drawn using the 
; U-fold grid. The solution is to add the top row of the U-fold grid
; to the T-grid prior to drawing.
; 
; Thanks to Petteri Uotila of CSIRO Marine & Atmospheric Research
; in Aspendale, Victoria, Australia for providing this solution.
;
;*************************************************
;
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

function fillArray(xi,xigap)
local dimxi, xo
begin
    dimxi = dimsizes(xi)
    ; add an extra row 
    xo = new((/dimxi(0)+1,dimxi(1)/),typeof(xi),getVarFillValue(xi))
    ; copy the original data
    xo(:dimxi(0)-1,:) = xi
    ; this line fills the extra cells with data in xigap
    xo(dimxi(0),:) = xigap(dimxi(0)-1,:)
    return(xo)
end

function fillGap(xi,ulat,ulon)
local xo, ulato, ulono
begin
    ; extend xi and fill with values next to extensions
    xo = fillArray(xi,xi)
    ; because 2d coordinates have been copied from xi and
    ; are of different dimsize than xo we need to delete them.
    ; Alternatively we could xo = (/' xi '/) but other 
    ; it is useful to copy all other attributes 
    delete(xo@lat2d)
    delete(xo@lon2d)
    ; locations in extensions are based on ulat and ulon
    ulato = fillArray(xi@lat2d,ulat)
    ulono = fillArray(xi@lon2d,ulon)
    xo@lat2d = ulato
    xo@lon2d = ulono
    return(xo)
end

function getContourLevels(varname,hemisphere)
begin
    if (varname .eq. "aice") then
       cntrs=(/5,10,15,20,30,40,50,60,70,80,85,90,95,99/) ; aice
    end if
    return(cntrs)
end

function setPolarPlotResources(varname,hemisphere)
    local NH_Pcntrs, SH_Pcntrs
begin
    ; add here a function that returns contours for different variables
    cntrs=getContourLevels(varname,hemisphere) ; aice
    ; resources settings good for sea-ice polar plots
    res = True
;
;   if panelling uncomment the next 2 lines and comment the line following
;   
;    res@gsnDraw  = False
;    res@gsnFrame = False
    res@gsnMaximize = True
;
;
    res@gsnPolar = hemisphere
    res@cnFillOn = True
    res@mpFillOn = False
    res@mpPerimDrawOrder = "PostDraw"
    res@gsnAddCyclic        = True
    res@cnLevelSelectionMode = "ExplicitLevels"     ; set manual contour levels
    res@cnLinesOn = False
    if( hemisphere .eq. "NH") then
        res@mpMinLatF = 50
    else
        res@mpMaxLatF = -50
    end if
    res@cnLevels       = cntrs
    res@cnMinLevelValF = cntrs(0)
    return(res) 
end

function maskIce(xi,hemisphere)
begin
    tlat2d = xi@lat2d
    size = dimsizes(tlat2d)
    if (size(0) .ne. dimsizes(xi(:,0))) then
      lat2d = tlat2d(:size(0)-2,:)
    else
      lat2d = tlat2d
    end if
    if( hemisphere .eq. "NH") then
        xi = mask(xi, lat2d.lt.50,False)
        xi = mask(xi, xi.lt.0.05,False)
    else
        xi = mask(xi, lat2d.gt.-50,False)
        xi = mask(xi, xi.lt.0.05,False)
    end if
    return(xi)
end

; main routine starts
begin
    fin = "/scratch/jk72/da1339/cice-dirs/runs/afim_era5/history/iceh.2005-01-31.nc"
    hemisphere = "SH"
    varname = "aice"
    fi = addfile(fin,"r")
    tlon2d = fi->TLON
    tlat2d = fi->TLAT
    ulon2d = fi->ULON
    ulat2d = fi->ULAT
    ; variable in t grid
    tice = 100*fi->aice(0,:,:)
    tice@lon2d = tlon2d
    tice@lat2d = tlat2d
    ; same variable but incorrectly in ugrid
    uice = tice
    uice@lon2d = ulon2d
    uice@lat2d = ulat2d
    ; extended tgrid by using ugrid in the gap
    eice = fillGap(tice,ulat2d,ulon2d)
    ulat = eice@lat2d
    ulon = eice@lon2d
    ; mask small values of ice concentration
    tice = maskIce(tice,hemisphere)
    uice = maskIce(uice,hemisphere)
    eice = maskIce(eice,hemisphere)
    printVarSummary(ulat)
    printVarSummary(ulon)
    printVarSummary(eice)

    ; colourmap
    wks = gsn_open_wks("png" ,"/g/data/jk72/da1339/GRAPHICAL/CICE/afim_era5/reG_patch/aice_2005-01-31")         ; send graphics to PNG file
;    cmap = "precip2_15lev"
;    gsn_define_colormap(wks,cmap)
    cmap = read_colormap_file("precip2_15lev")

    res = setPolarPlotResources(varname,hemisphere)
    plot = new(3,graphic)
    res@cnFillPalette = cmap(1:,:)
    res@cnFillMode = "RasterFill"
    ;res@cnRasterSmoothingOn = True          ; for smooth plot
    res@mpDataBaseVersion = "HighRes"
    res@tiMainString = "t-grid"
    ;
    ; uncomment the following lines to see the lat/lon and data rendered
    ; in index coordinate space (you can ignore warning messages)
    ;
    latplot = gsn_contour(wks,eice@lat2d,res)
    lonplot = gsn_contour(wks,eice@lon2d,res)
    dplot = gsn_contour(wks,eice,res)

    plot(0) = gsn_csm_contour_map_polar(wks,tice,res)
    res@tiMainString = "u-grid"
    plot(1) = gsn_csm_contour_map_polar(wks,uice,res)
    res@tiMainString = "t-grid extended with u-grid"
    plot(2) = gsn_csm_contour_map_polar(wks,eice,res)
;
;   if panelling uncomment the following lines 
;
    resP = True
    resP@gsnPaperOrientation = "portrait"
    resP@gsnMaximize      = True             ; maximize plot area
    gsn_panel(wks,plot,(/3,1/),resP)
;
end