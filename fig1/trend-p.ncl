load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
f1=addfile("/mnt/e/north-future/trend/ssp585/txx_ssp585.nc","r")
t1=f1->txxETCCDI(:,{60:90},:)
; printVarSummary(t1)
f2=addfile("/mnt/e/north-future/trend/ssp585/tnn_ssp585.nc","r")
t2=f2->tnnETCCDI(:,{60:90},:)
; printVarSummary(t2)
f3=addfile("/mnt/e/north-future/trend/ssp585/rx5day_ssp585.nc","r")
t3=f3->rx5dayETCCDI(:,{60:90},:)

t_all=new((/3,151,31,360/), float)
t_all(0,:,:,:)=t1
t_all(1,:,:,:)=t2
t_all(2,:,:,:)=t3

t_all_anoma=new((/3,86,31,360/), float)

do i=0,dimsizes(t_all(:,0,0,0))-1
    base_avg=dim_avg_n_Wrap(t_all(i,45:64,:,:),0)
    base_avg_c=conform_dims((/86,31,360/),base_avg,(/1,2/))
    copy_VarCoords(t_all(i,65:150,:,:),base_avg_c)
    t_all_anoma(i,:,:,:)=t_all(i,65:150,:,:)-base_avg_c
end do
copy_VarCoords(t_all(:,65:150,:,:),t_all_anoma)

a    = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
lsdata = a->LSMASK
lsm  = landsea_mask(lsdata,t1&lat,t1&lon)

year=ispan(2015,2100,1)

trend1  = regCoef_n(year,t_all_anoma(0,:,:,:),0,0)                        ;trends(lat,lon) 
copy_VarCoords(t_all_anoma(0,0,:,:), trend1) 
trend1  = trend1*10 
printVarSummary(trend1)

trend2  = regCoef_n(year,t_all_anoma(1,:,:,:),0,0)                        ;trends(lat,lon)
copy_VarCoords(t_all_anoma(1,0,:,:), trend2) 
trend2  = trend2*10 

trend3  = regCoef_n(year,t_all_anoma(2,:,:,:),0,0)                        ;trends(lat,lon)
copy_VarCoords(t_all_anoma(2,0,:,:), trend3) 
trend3  = trend3*10 

tval1 = onedtond(trend1@tval,dimsizes(trend1))
df1 = onedtond(trend1@nptxy, dimsizes(trend1)) - 2
; print(df)
b = tval1;betainc(x, a, b)的x，a，b的维度必须一致
b = 0.5
b!0 = "lat"
b!1 = "lon"
b&lat = trend1&lat
b&lon = trend1&lon

prob1 = betainc(df1/(df1+tval1^2),df1/2.0,b);
prob1!0 = "lat"
prob1!1 = "lon"
prob1&lat = trend1&lat
prob1&lon = trend1&lon
prob1 = mask(prob1,lsm.eq.0,False)
prob1 = prob1 + random_uniform(0., 0.0001, dimsizes(prob1))
printMinMax(prob1,False)

tval2 = onedtond(trend2@tval,dimsizes(trend2))
df2 = onedtond(trend2@nptxy, dimsizes(trend2)) - 2
prob2 = betainc(df2/(df2+tval2^2),df2/2.0,b);
prob2!0 = "lat"
prob2!1 = "lon"
prob2&lat = trend2&lat
prob2&lon = trend2&lon
prob2 = mask(prob2,lsm.eq.0,False)
prob2 = prob2 + random_uniform(0., 0.0001, dimsizes(prob2))
printMinMax(prob2,False)

tval3 = onedtond(trend3@tval,dimsizes(trend3))
df3 = onedtond(trend3@nptxy, dimsizes(trend3)) - 2
prob3 = betainc(df3/(df3+tval3^2),df3/2.0,b);
prob3!0 = "lat"
prob3!1 = "lon"
prob3&lat = trend3&lat
prob3&lon = trend3&lon
prob3 = mask(prob3,lsm.eq.0,False)
prob3 = prob3 + random_uniform(0., 0.0001, dimsizes(prob3))
printMinMax(prob3,False)

trenda=mask(trend3,lsm.eq.0,False)
printMinMax(trenda,False) 

wks = gsn_open_wks("eps","/mnt/e/north-future/new-p/original/fig1/trend/trend-all")
plot= new(3,graphic)                         ; create graphic array
plot_prob = new(3,graphic)

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
res@lbLabelBarOn        = True

sres                 = True                                                         ; set up a second resource list               
sres@gsnDraw         = False                                                          ;         do not draw the plot               
sres@gsnFrame        = False                                                          ; do not advance the frame         
sres@cnLineLabelsOn  = False                                                          ; do not use line labels                     
sres@cnLinesOn       = False                                                          ; do not draw contour lines                  
sres@cnInfoLabelOn   = False                                                                                                        
sres@cnFillOn        = True                                                          ; color fill                                   
sres@lbLabelBarOn= False                                                                                                            
sres@cnLevelSelectionMode = "ExplicitLevels"                                                                                             
sres@cnLevels    = (/0.05/)                                                                                                
sres@cnMonoFillPattern    = False            ; want multiple patterns                                                               
sres@cnFillPatterns       = (/17,-1/)     ; the patterns                                                                        
sres@cnMonoFillScale      = False            ; want different densities                                                            
sres@cnFillScales         = (/1,1/)    ; change densities                                                                        
sres@cnMonoFillColor      =True                                                                                                     
sres@cnFillDotSizeF       = 0.003  
sres@cnFillOpacityF = 0.6 

res1=res
res1@gsnLeftString  = "(d)Trend_SSP5-8.5_TXx"
res1@gsnRightString  = "~S~o~N~C/10yr"
; res1@lbLabelStride=2
res1@cnLevelSelectionMode =  "ExplicitLevels"   
res1@cnLevels = (/0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9/)
res1@cnFillColors = (/135,145,161,177,193,209,225,243,255/) 
plot(0)=gsn_csm_contour_map_polar(wks,trend1,res1)
plot_prob(0)=gsn_csm_contour(wks,prob1,sres)
overlay(plot(0),plot_prob(0))

res2                 = res
res2@gsnLeftString  = "(e)Trend_SSP5-8.5_TNn"
res2@gsnRightString  = "~S~o~N~C/10yr"
; res2@lbLabelStride=2
res2@cnLevelSelectionMode =  "ExplicitLevels"   
res2@cnLevels = (/0.6,0.9,1.2,1.5,1.8,2.1,2.4,2.7/)
res2@cnFillColors = (/138,145,156,161,177,193,209,225,241/)   
plot(1)=gsn_csm_contour_map_polar(wks,trend2,res2)
plot_prob(1)=gsn_csm_contour(wks,prob2,sres)
overlay(plot(1),plot_prob(1))


res3                 = res
res3@gsnLeftString  = "(f)Trend_SSP5-8.5_RX5day"
res3@gsnRightString  = "mm/10yr"
; res3@lbLabelStride=2
res3@cnFillPalette = "MPL_BrBG"
res3@cnLevelSelectionMode =  "ExplicitLevels"   
res3@cnLevels = (/0.8,1.2,1.6,2.0,2.4,2.8,3.2,3.6/)
res3@cnFillColors =(/66,72,80,88,96,104,112,120,127/)  
; res3@cnLevelSpacingF      = 0.1 
plot(2)=gsn_csm_contour_map_polar(wks,trend3,res3)
plot_prob(2)=gsn_csm_contour(wks,prob3,sres)
overlay(plot(2),plot_prob(2))

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