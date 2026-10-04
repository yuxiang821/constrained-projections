load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
f1=addfile("/mnt/e/north-future/ec/input/kj/all/rrv/ssp585/txx-rrv-ssp585.nc","r")
t1=f1->RRV
f2=addfile("/mnt/e/north-future/ec/input/kj/all/rrv/ssp585/tnn-rrv-ssp585.nc","r")
t2=f2->RRV
f3=addfile("/mnt/e/north-future/ec/input/kj/all/rrv/ssp585/rx5day-rrv-ssp585.nc","r")
t3=f3->RRV

r_all=new((/3,31,360/),float)
r_all(0,:,:)=t1
r_all(1,:,:)=t2
r_all(2,:,:)=t3


printMinMax(t1,False) 
printMinMax(t2,False) 
printMinMax(t3,False) 

wks = gsn_open_wks("pdf","/mnt/e/north-future/new-p/original/fig3/rrv/rrv-ssp585")
plot = new(3,graphic)     
res                      = True         
res@gsnDraw              = False        ; don't draw yet
res@gsnFrame             = False        ; don't advance frame yet
res@gsnAddCyclic          = True
res@gsnPolar             = "NH"

res@mpAreaMaskingOn     = True  ; open mask area function
res@mpOceanFillColor    = "white"  ; fill ocean with "white"
res@mpLandFillColor     = "transparent"  ; fill land with "dafult"
res@mpFillDrawOrder     = "PostDraw"  
; end mask
res@cnLinesOn     = False
res@cnInfoLabelOn = False

; res@mpFillOn             = False        ; turn off map fill
res@mpMaxLatF            = 90
res@mpMinLatF            = 60
res@mpCenterLonF         = 0

res@cnFillOn             = True         ; turn on color fill
;res@cnFillPalette        = "Hotres"  
res@cnLinesOn            = False        ; True is default
res@cnLineLabelsOn       = False      ; no longer automatic

res@gsnLeftStringFontHeightF=0.02
res@gsnRightStringFontHeightF=0.02
res@lbLabelBarOn        = True     ; no longer automatic

title=(/"SSP5-8.5_TXx_RRV (mm)","SSP5-8.5_TNn_RRV (mm)","SSP5-8.5_RX5day_RRV (mm)"/)

res1=res
res1@gsnLeftString  = title(0)
res1@gsnRightString  = ""
; res1@gsnRightString  = "~S~o~N~C/10yr"
; res1@lbLabelStride=2
res1@cnLevelSelectionMode =  "ExplicitLevels"   
res1@cnLevels = (/6,12,18,24,30,36,42,48,54,60/)
res1@cnFillColors = (/129,138,145,161,177,193,209,225,241,246,255/) 
plot(0)=gsn_csm_contour_map_polar(wks,r_all(0,:,:),res1)

res2                 = res
res2@gsnLeftString  = title(1)
res2@gsnRightString  = ""
; res2@gsnRightString  = "~S~o~N~C/10yr"
; res2@lbLabelStride=2
res2@cnLevelSelectionMode =  "ExplicitLevels"   
res2@cnLevels = (/4,8,12,16,20,24,28,32,36,40/)
res2@cnFillColors = (/129,138,145,161,177,193,209,225,241,246,255/) 
plot(1)=gsn_csm_contour_map_polar(wks,r_all(1,:,:),res2)

res3                 = res
res3@gsnLeftString  = title(2)
res3@gsnRightString  = ""
; res3@gsnRightString  = "mm/10yr"
; res3@lbLabelStride=2
res3@cnFillPalette = "MPL_BrBG"
res3@cnLevelSelectionMode =  "ExplicitLevels"  
res3@cnLevels = (/6,12,18,24,30,36,42,48/)
res3@cnFillColors =(/66,72,80,88,96,104,112,120,124/)  
; res3@cnLevelSpacingF      = 0.1 
plot(2)=gsn_csm_contour_map_polar(wks,r_all(2,:,:),res3)


resP                     = True         ; modify the panel plot
resP@gsnMaximize         = False        ; large format
resP@gsnPanelLabelBar    = False       ; add common colorbar
resP@gsnPaperOrientation = "portrait"   ; force portrait
;  resP@pmLabelBarOrthogonalPosF = 0.33
;  resP@lbLabelFontHeightF  = 0.007 
; draw panel without white space
;  resP@gsnPanelRowSpec = True 
;  drawNDCGrid(wks)
gsn_panel(wks,plot,(/1,3/),resP)


end