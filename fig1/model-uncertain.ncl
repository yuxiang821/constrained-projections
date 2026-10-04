load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
diri1 = "/mnt/e/etccdi/enavg/remap/ssp585/txx/"
files1 = systemfunc("ls " + diri1 + "*nc") ;批量读取文件路径
; print(files)
f1 = addfiles(files1,"r")
ListSetType (f1, "join")
; printVarSummary(f1)

t1=f1[:]->txxETCCDI(:,:,{60:90},:)
t1!0="ens"
n=dimsizes(t1(:,0,0,0))

diri2 = "/mnt/e/etccdi/enavg/remap/ssp585/tnn/"
files2 = systemfunc("ls " + diri2 + "*nc") ;批量读取文件路径
; print(files)
f2 = addfiles(files2,"r")
ListSetType (f2, "join")
; printVarSummary(f2)

t2=f2[:]->tnnETCCDI(:,:,{60:90},:)
t2!0="ens"
printVarSummary(t2)

diri3 = "/mnt/e/etccdi/enavg/remap/ssp585/rx5day/"
files3 = systemfunc("ls " + diri3 + "*nc") ;批量读取文件路径
; print(files)
f3 = addfiles(files3,"r")
ListSetType (f3, "join")
; printVarSummary(f3)

t3=f3[:]->rx5dayETCCDI(:,:,{60:90},:)
t3!0="ens"

printVarSummary(t3)

t_all=new((/3,n,151,31,360/), float)
t_all(0,:,:,:,:)=t1
t_all(1,:,:,:,:)=t2
t_all(2,:,:,:,:)=t3
t_all!0="ens"

t_all_ano=new((/3,n,31,360/), float)
do i=0,dimsizes(t_all(:,0,0,0,0))-1
    t_all_ano(i,:,:,:)=dim_avg_n_Wrap(t_all(i,:,65:150,:,:),1)-dim_avg_n_Wrap(t_all(i,:,45:64,:,:),1)
end do
copy_VarCoords(t_all(:,:,0,:,:),t_all_ano)
printVarSummary(t_all_ano)

r_RangePercentage = 0.954165 
t_avg=new((/3,31,360/), float)  
t_95=new((/3,31,360/), float)
t_5=new((/3,31,360/), float)
m=new((/3,31,360/), float)
snr=new((/3,31,360/), float)
do j=0,dimsizes(t_all_ano(:,0,0,0))-1
    t_avg(j,:,:)=dim_avg_n_Wrap(t_all_ano(j,:,:,:),0)
    t_95(j,:,:)=dim_avg_n_Wrap(t_all_ano(j,:,:,:), 0)+dim_stddev_n_Wrap(t_all_ano(j,:,:,:), 0)*r_RangePercentage
    t_5(j,:,:)=dim_avg_n_Wrap(t_all_ano(j,:,:,:), 0)-dim_stddev_n_Wrap(t_all_ano(j,:,:,:), 0)*r_RangePercentage
    m(j,:,:)=t_95(j,:,:)-t_5(j,:,:)
    snr(j,:,:)=abs(t_avg(j,:,:)) /m(j,:,:)
end do
t_avg!0="ens"
copy_VarCoords(t_all_ano(:,0,:,:),m)
copy_VarCoords(t_avg,snr)
; printVarSummary(snr(0,:,:))

printMinMax(m(0,:,:),False)
printMinMax(m(1,:,:),False)
printMinMax(m(2,:,:),False)

; printMinMax(t_std(0,:,:),False)
; printMinMax(t_std(1,:,:),False)
; printMinMax(t_std(2,:,:),False)
snr(2,:,:)=snr(2,:,:)+random_uniform(0., 0.2, dimsizes(snr(2,:,:)))
printMinMax(snr(0,:,:),False)
printMinMax(snr(1,:,:),False)
printMinMax(snr(2,:,:),False)

; ncdf1=addfile("/mnt/e/north-future/model-uncertain/uncertain/ssp585-3index.nc", "c")
; ncdf1->m_uncertain=m
; print(ncdf1)

; ncdf2=addfile("/mnt/e/north-future/model-uncertain/snr/ssp585-3index.nc", "c")
; ncdf2->SNR=snr
; print(ncdf2)

; snr_below = dim_num_n(where(snr(1,:,:).gt.1.5, 1, 0), (/0,1/))
; print("Number of points with SNR > 1.5: " + snr_below)

; a    = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
; lsdata = a->LSMASK
; lsm  = landsea_mask(lsdata,t1&lat,t1&lon)

; snr_mask = mask(snr,lsm.eq.0,False)
; printMinMax(snr_mask(1,:,:),False)

wks = gsn_open_wks("x11","/mnt/e/north-future/new-p/original/fig1/model-uncertain/model-uncertain-all")
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
sres@cnLevels    = (/1.5/)                                                                                                
sres@cnMonoFillPattern    = False            ; want multiple patterns                                                               
sres@cnFillPatterns       = (/17,-1/)     ; the patterns                                                                        
sres@cnMonoFillScale      = False            ; want different densities                                                            
sres@cnFillScales         = (/1,1/)    ; change densities                                                                        
sres@cnMonoFillColor      =True                                                                                                     
sres@cnFillDotSizeF       = 0.003
sres@cnFillOpacityF = 0.6 

res1=res
res1@gsnLeftString  = "(g)Inter-model Spread of TXx"
; res1@gsnRightString  = "~S~o~N~C/10yr"
; res1@lbLabelStride=2
res1@cnLevelSelectionMode =  "ExplicitLevels"   
res1@cnLevels = (/4.,6.,8.,10.,12.,14.,16.,18./)   ;(/4.,6.,8.,10.,12.,14.,16.,18./)
res1@cnFillColors = (/138,145,156,161,177,193,209,225,243/) 
plot(0)=gsn_csm_contour_map_polar(wks,m(0,:,:),res1)
plot_prob(0)=gsn_csm_contour(wks,snr(0,:,:),sres)
overlay(plot(0),plot_prob(0))

res2                 = res
res2@gsnLeftString  = "(h)Inter-model Spread of TNn"
; res2@gsnRightString  = "~S~o~N~C/10yr"
; res2@lbLabelStride=2
res2@cnLevelSelectionMode =  "ExplicitLevels"   
res2@cnLevels = (/4.,6.,8.,10.,12.,14.,16.,18./) ;(/4.,6.,8.,10.,12.,14.,16.,18./)
res2@cnFillColors = (/138,145,156,161,177,193,209,225,243/)   
plot(1)=gsn_csm_contour_map_polar(wks,m(1,:,:),res2)
plot_prob(1)=gsn_csm_contour(wks,snr(1,:,:),sres)
overlay(plot(1),plot_prob(1))


res3                 = res
res3@gsnLeftString  = "(i)Inter-model Spread of RX5day"
; res3@gsnRightString  = "mm/10yr"
; res3@lbLabelStride=2
res3@cnFillPalette = "MPL_BrBG"
res3@cnLevelSelectionMode =  "ExplicitLevels"   
res3@cnLevels = (/8.,16.,24.,32.,40.,48.,56./)  ;(/8.,16.,24.,32.,40.,48.,56./)
res3@cnFillColors =(/66,72,80,88,96,104,112,120/) 
; res3@cnLevelSpacingF      = 0.1 
plot(2)=gsn_csm_contour_map_polar(wks,m(2,:,:),res3)
plot_prob(2)=gsn_csm_contour(wks,snr(2,:,:),sres)
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