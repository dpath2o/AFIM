begin
    F_t2m        = "/g/data/rt52/era5/single-levels/reanalysis/2t/2005/2t_era5_oper_sfc_20050101-20050131.nc"
    F_t2m_reG    = "/scratch/jk72/da1339/afim_input/ERA5/0p25/airtmp_01.nc"
    src_file     = addfile(F_t2m, "r")
    src_file_reG = addfile(F_t2m_reG, "r")
    t2m          = src_file->t2m(:, :, :)  ; Access all time steps
    time_units   = src_file->time@units

; Loop over each time step
    do time_step = 0, dimsizes(t2m(:, 0, 0)) - 1
        t2m_time     = t2m(time_step, :, :)   ; Extract data for a specific time step
        t2m_reG_time = src_file_reG->airtmp(time_step, :, :)  ; Corresponding regridded data
        t_           = src_file_reG->time(time_step)
        t_@units     = time_units
        time_date    = cd_calendar(t_, 0)
        wks          = gsn_open_wks("png", "plot_"+time_date(0,0)+"_"+time_date(0,1)+"_"+time_date(0,2)+"_"+time_date(0,3)+"00")  ; Output file name will be automatically generated
        plot         = new(4,graphic)

;--- Plot 1: Original T2M North Pole Stereographic ---
        res          = True
        res@gsnPolar = "NH"
        res@mpOutlineOn = True
        res@mpGridAndLimbOn = False
        res@cnFillOn = True
        res@cnLinesOn = False
        res@tiMainString = "Original T2M - North Pole"
        plot(0) = gsn_csm_contour_map_polar(wks, t2m_time, res)
        delete([/res@tiMainString/])

;--- Plot 2: Regridded T2M North Pole Stereographic ---
        res@tiMainString = "Regridded T2M - North Pole"
        plot(1) = gsn_csm_contour_map_polar(wks, t2m_reG_time, res)
        delete([/res@tiMainString,res@gsnPolar/])

;--- Plot 3: Original T2M South Pole Stereographic ---
        res@gsnPolar = "SH"
        res@tiMainString = "Original T2M - South Pole"
        plot(2) = gsn_csm_contour_map_polar(wks, t2m_time, res)
        delete([/res@tiMainString/])

;--- Plot 4: Regridded T2M South Pole Stereographic ---
        res@tiMainString = "Regridded T2M - South Pole"
        plot(3) = gsn_csm_contour_map_polar(wks, t2m_reG_time, res)

; Panel the plots
        gsn_panel(wks,plot,(/2,2/), False)

; Close the workstation
        ;; frame(wks)
        delete([/wks,res/])
        
    end do

end
