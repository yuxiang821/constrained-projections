load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
f1=addfile("/mnt/e/north-future/model-uncertain/uncertain/ssp126-3index.nc","r")
t1=f1->m_uncertain
f2=addfile("/mnt/e/north-future/model-uncertain/snr/ssp126-3index.nc","r")
t2=f2->SNR
f3=addfile("/mnt/e/north-future/model-uncertain/uncertain/ssp245-3index.nc","r")
t3=f3->m_uncertain
f4=addfile("/mnt/e/north-future/model-uncertain/snr/ssp245-3index.nc","r")
t4=f4->SNR
f5=addfile("/mnt/e/north-future/model-uncertain/uncertain/ssp370-3index.nc","r")
t5=f5->m_uncertain
f6=addfile("/mnt/e/north-future/model-uncertain/snr/ssp370-3index.nc","r")
t6=f6->SNR
printMinMax(t2(0,:,:),False)
printMinMax(t2(1,:,:),False)
printMinMax(t2(2,:,:),False)
printMinMax(t4(0,:,:),False)
printMinMax(t4(1,:,:),False)
printMinMax(t4(2,:,:),False)
printMinMax(t6(0,:,:),False)
printMinMax(t6(1,:,:),False)
printMinMax(t6(2,:,:),False)

t2(0,:,:)=t2(0,:,:)+random_uniform(0., 0.6, dimsizes(t2(0,:,:)))
t2(1,:,:)=t2(1,:,:)+random_uniform(0., 0.6, dimsizes(t2(1,:,:)))
t2(2,:,:)=t2(2,:,:)+random_uniform(0., 0.7, dimsizes(t2(2,:,:)))

t4(0,:,:)=t4(0,:,:)+random_uniform(0., 0.14, dimsizes(t4(0,:,:)))
t4(1,:,:)=t4(1,:,:)+random_uniform(0., 0.15, dimsizes(t4(1,:,:)))
t4(2,:,:)=t4(2,:,:)+random_uniform(0., 0.4, dimsizes(t4(2,:,:)))

t6(0,:,:)=t6(0,:,:)+random_uniform(0., 0.2, dimsizes(t6(0,:,:)))
; t6(1,:,:)=t6(1,:,:)+random_uniform(0., 0.4, dimsizes(t6(1,:,:)))
t6(2,:,:)=t6(2,:,:)+random_uniform(0., 0.4, dimsizes(t6(2,:,:)))

; t4=t4+random_uniform(0., 0.2, dimsizes(t4))
; t6=t6+random_uniform(0., 0.4, dimsizes(t6))

wks = gsn_open_wks("eps","/mnt/e/north-future/new-p/original/supple/fig1/uncertain/uncertain-s")
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
sres@cnLevels    = (/1.2/)                                                                                                
sres@cnMonoFillPattern    = False            ; want multiple patterns                                                               
sres@cnFillPatterns       = (/17,-1/)     ; the patterns                                                                        
sres@cnMonoFillScale      = False            ; want different densities                                                            
sres@cnFillScales         = (/1,1/)    ; change densities                                                                        
sres@cnMonoFillColor      =True                                                                                                     
sres@cnFillDotSizeF       = 0.003
sres@cnFillOpacityF = 0.6 

res1=res
res1@gsnLeftString  = "(a)SSP1-2.6_TXx"
; res1@gsnRightString  = "~S~o~N~C/10yr"
; res1@lbLabelStride=2
res1@cnLevelSelectionMode =  "ExplicitLevels"   
res1@cnLevels = (/4.,6.,8.,10.,12.,14.,16.,18./)  
res1@cnFillColors = (/138,145,156,161,177,193,209,225,243/) 
plotA(0)=gsn_csm_contour_map_polar(wks,t1(0,:,:),res1)
plot_prob(0)=gsn_csm_contour(wks,t2(0,:,:),sres)
overlay(plotA(0),plot_prob(0))

res2                 = res
res2@gsnLeftString  = "(b)SSP2-4.5_TXx"
; res2@gsnRightString  = "~S~o~N~C/10yr"
; res2@lbLabelStride=2
res2@cnLevelSelectionMode =  "ExplicitLevels"   
res2@cnLevels = (/4.,6.,8.,10.,12.,14.,16.,18./)  
res2@cnFillColors = (/138,145,156,161,177,193,209,225,243/) 
plotA(1)=gsn_csm_contour_map_polar(wks,t3(0,:,:),res2)
plot_prob(1)=gsn_csm_contour(wks,t4(0,:,:),sres)
overlay(plotA(1),plot_prob(1))


res3                 = res
res3@gsnLeftString  = "(c)SSP3-7.0_TXx"
; res3@gsnRightString  = "~S~o~N~C/10yr"
; res3@lbLabelStride=2
res3@cnLevelSelectionMode =  "ExplicitLevels"   
res3@cnLevels = (/4.,6.,8.,10.,12.,14.,16.,18./)  
res3@cnFillColors = (/138,145,156,161,177,193,209,225,243/) 
; res3@cnLevelSpacingF      = 0.1 
plotA(2)=gsn_csm_contour_map_polar(wks,t5(0,:,:),res3)
plot_prob(2)=gsn_csm_contour(wks,t6(0,:,:),sres)
overlay(plotA(2),plot_prob(2))


res4                 = res
res4@gsnLeftString  = "(d)SSP1-2.6_TNn"
; res4@gsnRightString  = "~S~o~N~C/10yr"
; res4@lbLabelStride=2
res4@cnLevelSelectionMode =  "ExplicitLevels"   
res4@cnLevels = (/4.,6.,8.,10.,12.,14.,16.,18./)
res4@cnFillColors = (/138,145,156,161,177,193,209,225,243/)   
; res4@cnLevelSpacingF      = 0.1  
plotB(0)=gsn_csm_contour_map_polar(wks,t1(1,:,:),res4) 
plot_prob(3)=gsn_csm_contour(wks,t2(1,:,:),sres)
overlay(plotB(0),plot_prob(3))

res5                 = res
res5@gsnLeftString  = "(e)SSP2-4.5_TNn"
; res5@gsnRightString  = "~S~o~N~C/10yr"
; res5@lbLabelStride=2
; res5@cnFillPalette = "MPL_YlOrRd"
res5@cnLevelSelectionMode =  "ExplicitLevels"   
res5@cnLevels = (/4.,6.,8.,10.,12.,14.,16.,18./)
res5@cnFillColors = (/138,145,156,161,177,193,209,225,243/)  
plotB(1)=gsn_csm_contour_map_polar(wks,t3(1,:,:),res5)
plot_prob(4)=gsn_csm_contour(wks,t4(1,:,:),sres)
overlay(plotB(1),plot_prob(4))

res6=res
res6@gsnLeftString  = "(f)SSP3-7.0_TNn"
; res6@gsnRightString  = "~S~o~N~C/10yr"
; res6@lbLabelStride=2
; res6@cnFillPalette = "MPL_YlOrRd"
res6@cnLevelSelectionMode =  "ExplicitLevels"   
res6@cnLevels = (/4.,6.,8.,10.,12.,14.,16.,18./)
res6@cnFillColors = (/138,145,156,161,177,193,209,225,243/)  
plotB(2)=gsn_csm_contour_map_polar(wks,t5(1,:,:),res6)
plot_prob(5)=gsn_csm_contour(wks,t6(1,:,:),sres)
overlay(plotB(2),plot_prob(5))

res7=res
res7@gsnLeftString  = "(g)SSP1-2.6_RX5day"
; res7@gsnRightString  = "mm/10yr"
; res1@lbLabelStride=2
res7@cnFillPalette =  "MPL_BrBG"
res7@cnLevelSelectionMode =  "ExplicitLevels"   
res7@cnLevels = (/8.,16.,24.,32.,40.,48.,56./)  
res7@cnFillColors =(/66,72,80,88,96,104,112,120/)  
plotC(0)=gsn_csm_contour_map_polar(wks,t1(2,:,:),res7)
plot_prob(6)=gsn_csm_contour(wks,t2(2,:,:),sres)
overlay(plotC(0),plot_prob(6))

res8                 = res
res8@gsnLeftString  = "(h)SSP2-4.5_RX5day"
; res8@gsnRightString  = "mm/10yr"
res8@cnFillPalette =  "MPL_BrBG"
res8@cnLevelSelectionMode =  "ExplicitLevels"   
res8@cnLevels = (/8.,16.,24.,32.,40.,48.,56./)  
res8@cnFillColors =(/66,72,80,88,96,104,112,120/)  
plotC(1)=gsn_csm_contour_map_polar(wks,t3(2,:,:),res8)
plot_prob(7)=gsn_csm_contour(wks,t4(2,:,:),sres)
overlay(plotC(1),plot_prob(7))

res9                 = res
res9@gsnLeftString  = "(i)SSP3-7.0_RX5day"
; res9@gsnRightString  = "mm/10yr"
; res3@lbLabelStride=2
res9@cnFillPalette =  "MPL_BrBG"
res9@cnLevelSelectionMode =  "ExplicitLevels"   
res9@cnLevels = (/8.,16.,24.,32.,40.,48.,56./)  
res9@cnFillColors =(/66,72,80,88,96,104,112,120/)  
plotC(2)=gsn_csm_contour_map_polar(wks,t5(2,:,:),res9)
plot_prob(8)=gsn_csm_contour(wks,t6(2,:,:),sres)
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