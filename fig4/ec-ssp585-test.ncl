load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin

;=============================================================
; 0. 读取 Raw 数据  [member, lat, lon]
;=============================================================
f1 = addfile("/mnt/e/north-future/ec/input/kj-test/raw-relative1995-2014/ssp370/txx-raw-ssp370-all.nc","r")
t1 = f1->txxETCCDI
f2 = addfile("/mnt/e/north-future/ec/input/kj-test/raw-relative1995-2014/ssp370/tnn-raw-ssp370-all.nc","r")
t2 = f2->tnnETCCDI
f3 = addfile("/mnt/e/north-future/ec/input/kj-test/raw-relative1995-2014/ssp370/rx5day-raw-ssp370-all.nc","r")
t3 = f3->rx5dayETCCDI

printVarSummary(t1)

;=============================================================
; 1. 读取 Constrained 均值和标准差  [lat, lon]
;=============================================================
f4 = addfile("/mnt/e/north-future/ec/input/kj-test/c-relative1995-2014/ssp370/txx-c-ssp370.nc","r")
t4     = f4->txxETCCDI              ; EC_mean
t4_std = f4->std                    ; EC_std

f5 = addfile("/mnt/e/north-future/ec/input/kj-test/c-relative1995-2014/ssp370/tnn-c-ssp370.nc","r")
t5     = f5->tnnETCCDI
t5_std = f5->std

f6 = addfile("/mnt/e/north-future/ec/input/kj-test/c-relative1995-2014/ssp370/rx5day-c-ssp370.nc","r")
t6     = f6->rx5dayETCCDI
t6_std = f6->std

printVarSummary(t4)
printVarSummary(t4_std)

f_rrv1 = addfile("/mnt/e/north-future/ec/input/kj/all/rrv/ssp370/txx-rrv-ssp370.nc","r")
r1 = f_rrv1->RRV
f_rrv2 = addfile("/mnt/e/north-future/ec/input/kj/all/rrv/ssp370/tnn-rrv-ssp370.nc","r")
r2 = f_rrv2->RRV
f_rrv3 = addfile("/mnt/e/north-future/ec/input/kj/all/rrv/ssp370/rx5day-rrv-ssp370.nc","r")
r3 = f_rrv3->RRV

r=new((/3,31,360/), float)
r(0,:,:)=r1
r(1,:,:)=r2
r(2,:,:)=r3

a    = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
lsdata = a->LSMASK
lsm  = landsea_mask(lsdata,r1&lat,r1&lon)

t_std_mask=new((/3,31,360/), float)
do i=0,2
   t_std_mask(i,:,:) = mask(r(i,:,:),lsm.eq.0,False)
   printMinMax(t_std_mask(i,:,:),False)
end do

rad  = 4.0 * atan(1.0) / 180.0
wlat = cos(r1&lat * rad)

t_reg_mean = wgt_areaave_Wrap(t_std_mask, wlat, 1.0, 0)

print(t_reg_mean)

;=============================================================
; 2. 计算多模式平均（Unconstrained mean）
;=============================================================
t1_avg = dim_avg_n_Wrap(t1, 0)
t2_avg = dim_avg_n_Wrap(t2, 0)
t3_avg = dim_avg_n_Wrap(t3, 0)

;=============================================================
; 3. 计算差值场（d-f 面板）
;=============================================================
dif1_avg = t4 - t1_avg
dif2_avg = t5 - t2_avg
dif3_avg = t6 - t3_avg
copy_VarCoords(t4, dif1_avg)
copy_VarCoords(t5, dif2_avg)
copy_VarCoords(t6, dif3_avg)

printMinMax(t4, False)
printMinMax(t5, False)
printMinMax(t6, False)
printMinMax(dif1_avg, False)
printMinMax(dif2_avg, False)
printMinMax(dif3_avg, False)


;-------------------------------------------------------------
; 4b. (d-f) Difference 面板：单样本 t 检验
;-------------------------------------------------------------
dims  = dimsizes(t1)
n_ens = dims(0)
df    = tofloat(n_ens - 1)
n_flt = tofloat(n_ens)
print("Number of raw ensemble members = " + n_ens)

; --- TXx ---
d1_std  = dim_stddev_n(t1, 0)
d1_se   = d1_std / sqrt(n_flt)
d1_tval = where(d1_se .gt. 0, dif1_avg / d1_se, 0.)
prob_dif1 = dif1_avg
prob_dif1 = (/betainc(df/(df + d1_tval^2), \
              conform(d1_tval, 0.5*df, -1), \
              conform(d1_tval, 0.5,    -1) )/)

; --- TNn ---
d2_std  = dim_stddev_n(t2, 0)
d2_se   = d2_std / sqrt(n_flt)
d2_tval = where(d2_se .gt. 0, dif2_avg / d2_se, 0.)
prob_dif2 = dif2_avg
prob_dif2 = (/betainc(df/(df + d2_tval^2), \
              conform(d2_tval, 0.5*df, -1), \
              conform(d2_tval, 0.5,    -1) )/)

; --- Rx5day ---
d3_std  = dim_stddev_n(t3, 0)
d3_se   = d3_std / sqrt(n_flt)
d3_tval = where(d3_se .gt. 0, dif3_avg / d3_se, 0.)
prob_dif3 = dif3_avg
prob_dif3 = (/betainc(df/(df + d3_tval^2), \
              conform(d3_tval, 0.5*df, -1), \
              conform(d3_tval, 0.5,    -1) )/)


print("===== Difference panels (d-f): % grid points p<0.05 =====")
print("TXx:    " + sprintf("%5.1f", 100.0 * num(ndtooned(prob_dif1) .lt. 0.05) / \
      num(.not.ismissing(ndtooned(prob_dif1)))) + "%")
print("TNn:    " + sprintf("%5.1f", 100.0 * num(ndtooned(prob_dif2) .lt. 0.05) / \
      num(.not.ismissing(ndtooned(prob_dif2)))) + "%")
print("Rx5day: " + sprintf("%5.1f", 100.0 * num(ndtooned(prob_dif3) .lt. 0.05) / \
      num(.not.ismissing(ndtooned(prob_dif3)))) + "%")

;=============================================================
; 5. 画图
;=============================================================
wks = gsn_open_wks("x11","/mnt/e/north-future/new-p/original/fig3/all/ec-ssp370")
plot     = new(9, graphic)
plot_sig = new(3, graphic)

;--- 公共资源 ---
res                      = True         
res@gsnDraw              = False
res@gsnFrame             = False
res@gsnAddCyclic         = True
res@gsnPolar             = "NH"

res@mpAreaMaskingOn      = True
res@mpOceanFillColor     = "white"
res@mpLandFillColor      = "transparent"
res@mpFillDrawOrder      = "PostDraw"

res@cnLinesOn            = False
res@cnInfoLabelOn        = False
res@mpMaxLatF            = 90
res@mpMinLatF            = 60
res@mpCenterLonF         = 0
res@cnFillOn             = True
res@cnLineLabelsOn       = False

res@lbLabelFontHeightF   = 0.02
res@gsnLeftStringFontHeightF  = 0.02

;--- ★ 显著性打点资源 ★ ---
res_stip                 = True                                                         ; set up a second resource list               
res_stip@gsnDraw         = False                                                          ;         do not draw the plot               
res_stip@gsnFrame        = False                                                          ; do not advance the frame         
res_stip@cnLineLabelsOn  = False                                                          ; do not use line labels                     
res_stip@cnLinesOn       = False                                                          ; do not draw contour lines                  
res_stip@cnInfoLabelOn   = False                                                                                                        
res_stip@cnFillOn        = True                                                          ; color fill                                   
res_stip@lbLabelBarOn= False                                                                                                            
res_stip@cnLevelSelectionMode = "ExplicitLevels"                                                                                             
res_stip@cnLevels    = (/0.05/)                                                                                                
res_stip@cnMonoFillPattern    = False            ; want multiple patterns                                                               
res_stip@cnFillPatterns       = (/17,-1/)     ; the patterns                                                                        
res_stip@cnMonoFillScale      = False            ; want different densities                                                            
res_stip@cnFillScales         = (/1,1/)    ; change densities                                                                        
res_stip@cnMonoFillColor      =True                                                                                                     
res_stip@cnFillDotSizeF       = 0.003  
res_stip@cnFillOpacityF = 0.6 

;=============================================================
; (a) Constrained TXx
;=============================================================
res1 = res
res1@gsnLeftString  = "Constrained-TXx_SSP3-7.0 (~S~o~N~C)"
res1@gsnRightString = " "
res1@cnLevelSelectionMode = "ExplicitLevels"   
res1@cnLevels     = (/0.5,1.,1.5,2.0,2.5,3.0,3.5/)
res1@cnFillColors = (/129,145,161,177,193,209,225,241/) 
plot(0) = gsn_csm_contour_map_polar(wks, t4, res1)


;=============================================================
; (b) Constrained TNn
;=============================================================
res4 = res
res4@gsnLeftString  = "Constrained-TNn_SSP3-7.0 (~S~o~N~C)"
res4@gsnRightString = " "
res4@cnLevelSelectionMode = "ExplicitLevels"   
res4@cnLevels     = (/1,2,3,4,5,6,7,8,9/)
res4@cnFillColors = (/129,136,145,156,161,177,193,209,225,241/) 
plot(1) = gsn_csm_contour_map_polar(wks, t5, res4)


;=============================================================
; (c) Constrained Rx5day
;=============================================================
res7 = res
res7@gsnLeftString  = "Constrained-RX5day_SSP3-7.0 (mm)"
res7@gsnRightString = " "
res7@cnFillPalette = "MPL_BrBG"
res7@cnLevelSelectionMode = "ExplicitLevels"   
res7@cnLevels     = (/-6,-4,-2,0,2,4,6,8,10/)
res7@cnFillColors = (/18,26,34,42,72,80,88,96,104,112/) 
plot(2) = gsn_csm_contour_map_polar(wks, t6, res7)


;=============================================================
; (d) Difference TXx
;=============================================================
res2 = res
res2@gsnLeftString  = "Constrained-Unconstrained-TXx_SSP3-7.0 (~S~o~N~C)"
res2@gsnRightString = ""
res2@cnFillPalette = "cmp_b2r"
res2@cnLevelSelectionMode = "ExplicitLevels"   
res2@cnLevels     = (/-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0/)
res2@cnFillColors = (/2,6,10,12,14,16,20,34/) 
plot(3) = gsn_csm_contour_map_polar(wks, dif1_avg, res2)

plot_sig(0) = gsn_csm_contour(wks, prob_dif1, res_stip)
overlay(plot(3), plot_sig(0))

;=============================================================
; (e) Difference TNn
;=============================================================
res5 = res
res5@gsnLeftString  = "Constrained-Unconstrained-TNn_SSP3-7.0 (~S~o~N~C)"
res5@gsnRightString = ""
res5@cnFillPalette = "cmp_b2r"
res5@cnLevelSelectionMode = "ExplicitLevels"   
res5@cnLevels     = (/-2.4,-2.0,-1.6,-1.2,-0.8,-0.4,0/)
res5@cnFillColors = (/2,6,10,12,14,16,20,34/) 
plot(4) = gsn_csm_contour_map_polar(wks, dif2_avg, res5)

plot_sig(1) = gsn_csm_contour(wks, prob_dif2, res_stip)
overlay(plot(4), plot_sig(1))

;=============================================================
; (f) Difference Rx5day
;=============================================================
res8 = res
res8@gsnLeftString  = "Constrained-Unconstrained-RX5day_SSP3-7.0 (mm)"
res8@gsnRightString = ""
res8@cnFillPalette = "MPL_BrBG"
res8@cnLevelSelectionMode = "ExplicitLevels"   
res8@cnLevels     = (/-8,-7,-6,-5,-4,-3,-2,-1,0/)
res8@cnFillColors = (/2,8,14,20,26,32,38,44,50,80/) 
plot(5) = gsn_csm_contour_map_polar(wks, dif3_avg, res8)

plot_sig(2) = gsn_csm_contour(wks, prob_dif3, res_stip)
overlay(plot(5), plot_sig(2))

res_g = res
res_g@gsnLeftString        = "RRV-TXx_SSP3-7.0 (%)"
res_g@gsnRightString       = ""
res_g@cnLevelSelectionMode = "ExplicitLevels"
res_g@cnLevels             = (/6,12,18,24,30,36,42,48,54,60/)
res_g@cnFillColors         = (/129,138,145,161,177,193,209,225,241,246,255/)
plot(6) = gsn_csm_contour_map_polar(wks, r1, res_g)

res_h = res
res_h@gsnLeftString        = "RRV-TNn_SSP3-7.0 (%)"
res_h@gsnRightString       = ""
res_h@cnLevelSelectionMode = "ExplicitLevels"
res_h@cnLevels             = (/4,8,12,16,20,24,28,32,36,40/)
res_h@cnFillColors         = (/129,138,145,161,177,193,209,225,241,246,255/)
plot(7) = gsn_csm_contour_map_polar(wks, r2, res_h)

res_i = res
res_i@gsnLeftString        = "RRV-RX5day_SSP3-7.0 (%)"
res_i@gsnRightString       = ""
res_i@cnFillPalette        = "MPL_BrBG"
res_i@cnLevelSelectionMode = "ExplicitLevels"
res_i@cnLevels             = (/6,12,18,24,30,36,42,48/)
res_i@cnFillColors         = (/66,72,80,88,96,104,112,120,124/)
plot(8) = gsn_csm_contour_map_polar(wks, r3, res_i)

;=============================================================
; 6. Panel 输出
;=============================================================
resP                     = True
resP@gsnMaximize         = True
resP@gsnPanelLabelBar    = False
resP@gsnPaperOrientation = "portrait"
resP@pmLabelBarOrthogonalPosF = 0.33
gsn_panel(wks, plot, (/3,3/), resP)

end
