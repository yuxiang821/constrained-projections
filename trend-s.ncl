load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
f1=addfile("/mnt/e/north-future/trend/ssp126/txx_ssp126.nc","r")
t1=f1->txxETCCDI(:,{60:90},:)
; printVarSummary(t1)
f2=addfile("/mnt/e/north-future/trend/ssp245/txx_ssp245.nc","r")
t2=f2->txxETCCDI(:,{60:90},:)
; printVarSummary(t2)
f3=addfile("/mnt/e/north-future/trend/ssp370/txx_ssp370.nc","r")
t3=f3->txxETCCDI(:,{60:90},:)

f5=addfile("/mnt/e/north-future/trend/ssp126/tnn_ssp126.nc","r")
t5=f5->tnnETCCDI(:,{60:90},:)
; printVarSummary(t5)
f6=addfile("/mnt/e/north-future/trend/ssp245/tnn_ssp245.nc","r")
t6=f6->tnnETCCDI(:,{60:90},:)
f7=addfile("/mnt/e/north-future/trend/ssp370/tnn_ssp370.nc","r")
t7=f7->tnnETCCDI(:,{60:90},:)

f9=addfile("/mnt/e/north-future/trend/ssp126/rx5day_ssp126.nc","r")
t9=f9->rx5dayETCCDI(:,{60:90},:)
f10=addfile("/mnt/e/north-future/trend/ssp245/rx5day_ssp245.nc","r")
t10=f10->rx5dayETCCDI(:,{60:90},:)
f11=addfile("/mnt/e/north-future/trend/ssp370/rx5day_ssp370.nc","r")
t11=f11->rx5dayETCCDI(:,{60:90},:)

t_all=new((/9,151,31,360/), float)
t_all(0,:,:,:)=t1
t_all(1,:,:,:)=t2
t_all(2,:,:,:)=t3
t_all(3,:,:,:)=t5
t_all(4,:,:,:)=t6
t_all(5,:,:,:)=t7
t_all(6,:,:,:)=t9
t_all(7,:,:,:)=t10
t_all(8,:,:,:)=t11

printVarSummary(t_all(0,65:150,:,:))

t_all_anoma=new((/9,86,31,360/), float)

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
; printVarSummary(trend1)

trend2  = regCoef_n(year,t_all_anoma(1,:,:,:),0,0)                        ;trends(lat,lon)
copy_VarCoords(t_all_anoma(1,0,:,:), trend2) 
trend2  = trend2*10 

trend3  = regCoef_n(year,t_all_anoma(2,:,:,:),0,0)                        ;trends(lat,lon)
copy_VarCoords(t_all_anoma(2,0,:,:), trend3) 
trend3  = trend3*10 

trend4  = regCoef_n(year,t_all_anoma(3,:,:,:),0,0)                        ;trends(lat,lon)
copy_VarCoords(t_all_anoma(3,0,:,:), trend4) 
trend4  = trend4*10 

trend5  = regCoef_n(year,t_all_anoma(4,:,:,:),0,0)                        ;trends(lat,lon)
copy_VarCoords(t_all_anoma(4,0,:,:), trend5) 
trend5  = trend5*10 
; printVarSummary(trend5)

trend6  = regCoef_n(year,t_all_anoma(5,:,:,:),0,0)                        ;trends(lat,lon)
copy_VarCoords(t_all_anoma(5,0,:,:), trend6) 
trend6  = trend6*10 

trend7  = regCoef_n(year,t_all_anoma(6,:,:,:),0,0)                        ;trends(lat,lon)
copy_VarCoords(t_all_anoma(6,0,:,:), trend7) 
trend7  = trend7*10 

trend8  = regCoef_n(year,t_all_anoma(7,:,:,:),0,0)                        ;trends(lat,lon)
copy_VarCoords(t_all_anoma(7,0,:,:), trend8) 
trend8  = trend8*10 

trend9  = regCoef_n(year,t_all_anoma(8,:,:,:),0,0)                        ;trends(lat,lon)
copy_VarCoords(t_all_anoma(8,0,:,:), trend9) 
trend9 = trend9*10 


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
printMinMax(prob1,False)

tval2 = onedtond(trend2@tval,dimsizes(trend2))
df2 = onedtond(trend2@nptxy, dimsizes(trend2)) - 2
prob2 = betainc(df2/(df2+tval2^2),df2/2.0,b);
prob2!0 = "lat"
prob2!1 = "lon"
prob2&lat = trend2&lat
prob2&lon = trend2&lon
prob2 = mask(prob2,lsm.eq.0,False)
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

tval4 = onedtond(trend4@tval,dimsizes(trend4))
df4 = onedtond(trend4@nptxy, dimsizes(trend4)) - 2
prob4 = betainc(df4/(df4+tval4^2),df4/2.0,b);
prob4!0 = "lat"
prob4!1 = "lon"
prob4&lat = trend4&lat
prob4&lon = trend4&lon
prob4 = mask(prob4,lsm.eq.0,False)
prob4 = prob4 + random_uniform(0., 0.0001, dimsizes(prob4))
printMinMax(prob4,False)

tval5 = onedtond(trend5@tval,dimsizes(trend5))
df5 = onedtond(trend5@nptxy, dimsizes(trend5)) - 2
prob5 = betainc(df5/(df5+tval5^2),df5/2.0,b);
prob5!0 = "lat"
prob5!1 = "lon"
prob5&lat = trend5&lat
prob5&lon = trend5&lon
prob5 = mask(prob5,lsm.eq.0,False)
printMinMax(prob5,False)

tval6 = onedtond(trend6@tval,dimsizes(trend6))
df6 = onedtond(trend6@nptxy, dimsizes(trend6)) - 2
prob6 = betainc(df6/(df6+tval6^2),df6/2.0,b);
prob6!0 = "lat"
prob6!1 = "lon"
prob6&lat = trend6&lat
prob6&lon = trend6&lon
prob6 = mask(prob6,lsm.eq.0,False)
prob6 = prob6 + random_uniform(0., 0.0001, dimsizes(prob6))
printMinMax(prob6,False)

tval7 = onedtond(trend7@tval,dimsizes(trend7))
df7 = onedtond(trend7@nptxy, dimsizes(trend7)) - 2
prob7 = betainc(df7/(df7+tval7^2),df7/2.0,b);
prob7!0 = "lat"
prob7!1 = "lon"
prob7&lat = trend7&lat
prob7&lon = trend7&lon
prob7 = mask(prob7,lsm.eq.0,False)
prob7 = prob7 + random_uniform(0., 0.0001, dimsizes(prob7))
printMinMax(prob7,False)

tval8 = onedtond(trend8@tval,dimsizes(trend8))
df8 = onedtond(trend8@nptxy, dimsizes(trend8)) - 2
prob8 = betainc(df8/(df8+tval8^2),df8/2.0,b);
prob8!0 = "lat"
prob8!1 = "lon"
prob8&lat = trend8&lat
prob8&lon = trend8&lon
prob8 = mask(prob8,lsm.eq.0,False)
prob8 = prob8 + random_uniform(0., 0.0001, dimsizes(prob8))
printMinMax(prob8,False)

tval9 = onedtond(trend9@tval,dimsizes(trend9))
df9 = onedtond(trend9@nptxy, dimsizes(trend9)) - 2
prob9 = betainc(df9/(df9+tval9^2),df9/2.0,b);
prob9!0 = "lat"
prob9!1 = "lon"
prob9&lat = trend9&lat
prob9&lon = trend9&lon
prob9 = mask(prob9,lsm.eq.0,False)
printMinMax(prob9,False)


; printMinMax(trend1,False) 
; printMinMax(trend2,False) 
; printMinMax(trend3,False) 
; printMinMax(trend4,False) 
; printMinMax(trend5,False) 
; printMinMax(trend6,False) 
; printMinMax(trend7,False) 
; printMinMax(trend8,False) 
; printMinMax(trend9,False) 

wks = gsn_open_wks("pdf","/mnt/e/north-future/new-p/original/supple/fig1/trend/trend-s")
plotA = new(3,graphic)                         ; create graphic array
plotB = new(3,graphic)
plotC = new(3,graphic)
plot_prob = new(9,graphic)

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
res@lbLabelBarOn        = False

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
res1@gsnLeftString  = "(a)SSP1-2.6_TXx"
res1@gsnRightString  = "~S~o~N~C/10yr"
; res1@lbLabelStride=2
res1@cnLevelSelectionMode =  "ExplicitLevels"   
res1@cnLevels = (/0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8/)
res1@cnFillColors = (/129,135,145,161,177,193,209,225,243/) 
plotA(0)=gsn_csm_contour_map_polar(wks,trend1,res1)
plot_prob(0)=gsn_csm_contour(wks,prob1,sres)
overlay(plotA(0),plot_prob(0))

res2                 = res
res2@gsnLeftString  = "(b)SSP2-4.5_TXx"
res2@gsnRightString  = "~S~o~N~C/10yr"
; res2@lbLabelStride=2
res2@cnLevelSelectionMode =  "ExplicitLevels"   
res2@cnLevels = (/0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8/)
res2@cnFillColors = (/129,135,145,161,177,193,209,225,243/) 
plotA(1)=gsn_csm_contour_map_polar(wks,trend2,res2)
plot_prob(1)=gsn_csm_contour(wks,prob2,sres)
overlay(plotA(1),plot_prob(1))


res3                 = res
res3@gsnLeftString  = "(c)SSP3-7.0_TXx"
res3@gsnRightString  = "~S~o~N~C/10yr"
; res3@lbLabelStride=2
res3@cnLevelSelectionMode =  "ExplicitLevels"   
res3@cnLevels = (/0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8/)
res3@cnFillColors = (/129,135,145,161,177,193,209,225,243/)  
; res3@cnLevelSpacingF      = 0.1 
plotA(2)=gsn_csm_contour_map_polar(wks,trend3,res3)
plot_prob(2)=gsn_csm_contour(wks,prob3,sres)
overlay(plotA(2),plot_prob(2))


res4                 = res
res4@gsnLeftString  = "(d)SSP1-2.6_TNn"
res4@gsnRightString  = "~S~o~N~C/10yr"
; res4@lbLabelStride=2
res4@cnLevelSelectionMode =  "ExplicitLevels"   
res4@cnLevels = (/0.3,0.6,0.9,1.2,1.5,1.8,2.1,2.4,2.7/)
res4@cnFillColors = (/129,138,145,156,161,177,193,209,225,241/)  
; res4@cnLevelSpacingF      = 0.1  
plotB(0)=gsn_csm_contour_map_polar(wks,trend4,res4) 
plot_prob(3)=gsn_csm_contour(wks,prob4,sres)
overlay(plotB(0),plot_prob(3))

res5                 = res
res5@gsnLeftString  = "(e)SSP2-4.5_TNn"
res5@gsnRightString  = "~S~o~N~C/10yr"
; res5@lbLabelStride=2
; res5@cnFillPalette = "MPL_YlOrRd"
res5@cnLevelSelectionMode =  "ExplicitLevels"   
res5@cnLevels = (/0.3,0.6,0.9,1.2,1.5,1.8,2.1,2.4,2.7/)
res5@cnFillColors = (/129,138,145,156,161,177,193,209,225,241/)  
plotB(1)=gsn_csm_contour_map_polar(wks,trend5,res5)
plot_prob(4)=gsn_csm_contour(wks,prob5,sres)
overlay(plotB(1),plot_prob(4))

res6=res
res6@gsnLeftString  = "(f)SSP3-7.0_TNn"
res6@gsnRightString  = "~S~o~N~C/10yr"
; res6@lbLabelStride=2
; res6@cnFillPalette = "MPL_YlOrRd"
res6@cnLevelSelectionMode =  "ExplicitLevels"   
res6@cnLevels = (/0.3,0.6,0.9,1.2,1.5,1.8,2.1,2.4,2.7/)
res6@cnFillColors = (/129,138,145,156,161,177,193,209,225,241/)  
plotB(2)=gsn_csm_contour_map_polar(wks,trend6,res6)
plot_prob(5)=gsn_csm_contour(wks,prob6,sres)
overlay(plotB(2),plot_prob(5))

res7=res
res7@gsnLeftString  = "(g)SSP1-2.6_RX5day"
res7@gsnRightString  = "mm/10yr"
; res1@lbLabelStride=2
res7@cnFillPalette = "MPL_BrBG"
res7@cnLevelSelectionMode =  "ExplicitLevels"   
res7@cnLevels = (/0.4,0.8,1.2,1.6,2.0,2.4,2.8,3.2/)
res7@cnFillColors =(/66,72,80,88,96,104,112,120,127/) 
plotC(0)=gsn_csm_contour_map_polar(wks,trend7,res7)
plot_prob(6)=gsn_csm_contour(wks,prob7,sres)
overlay(plotC(0),plot_prob(6))

res8                 = res
res8@gsnLeftString  = "(h)SSP2-4.5_RX5day"
res8@gsnRightString  = "mm/10yr"
res8@cnFillPalette = "MPL_BrBG"
res8@cnLevelSelectionMode =  "ExplicitLevels"   
res8@cnLevels = (/0.4,0.8,1.2,1.6,2.0,2.4,2.8,3.2/)
res8@cnFillColors =(/66,72,80,88,96,104,112,120,127/) 
plotC(1)=gsn_csm_contour_map_polar(wks,trend8,res8)
plot_prob(7)=gsn_csm_contour(wks,prob8,sres)
overlay(plotC(1),plot_prob(7))

res9                 = res
res9@gsnLeftString  = "(i)SSP3-7.0_RX5day"
res9@gsnRightString  = "mm/10yr"
; res3@lbLabelStride=2
res9@cnFillPalette = "MPL_BrBG"
res9@cnLevelSelectionMode =  "ExplicitLevels"   
res9@cnLevels = (/0.4,0.8,1.2,1.6,2.0,2.4,2.8,3.2/)
res9@cnFillColors =(/66,72,80,88,96,104,112,120,127/) 
plotC(2)=gsn_csm_contour_map_polar(wks,trend9,res9)
plot_prob(8)=gsn_csm_contour(wks,prob9,sres)
overlay(plotC(2),plot_prob(8))

pres1                  = True                ; mods desired
pres1@gsnFrame         = False               ; save panel until both ready
; pres1@gsnPanelRight    = 0.5                 ; draw from left edge to center
pres1@gsnPanelLabelBar = True                ; commom label bar
; pres1@lbLabelStride    = 2                   ; every other label
pres1@gsnPanelTop        = 0.9
pres1@gsnPanelBottom     =0.65
pres1@lbLabelFontHeightF   = 0.01


pres2                  = True                ; mods desired
pres2@gsnFrame         = False               ; save panel until both ready
; pres2@gsnPanelLeft     = 0.5		      ; draw from center to right edge 
pres2@gsnPanelLabelBar = True		      ; commom label bar
; pres2@lbLabelStride    = 4		      ; every 4th label
pres2@gsnPanelTop        = 0.65
pres2@gsnPanelBottom     =0.4
pres2@lbLabelFontHeightF   = 0.01

pres3                  = True                ; mods desired
pres3@gsnFrame         = False               ; save panel until both ready
; pres3@gsnPanelLeft     = 0.5		      ; draw from center to right edge 
pres3@gsnPanelLabelBar = True		      ; commom label bar
; pres3@lbLabelStride    = 4	
pres3@gsnPanelTop        = 0.4
pres3@gsnPanelBottom     =0.15
pres3@lbLabelFontHeightF   = 0.01

gsn_panel(wks,plotA,(/1,3/),pres1)          ; create first panel
gsn_panel(wks,plotB,(/1,3/),pres2)          ; create second panel
gsn_panel(wks,plotC,(/1,3/),pres3)
frame(wks)                        

end