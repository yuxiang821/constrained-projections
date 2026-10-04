load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
tcvf=1.64485
tcvf66=0.954165
; print(tcvf)
; print(tcvf66)

f1 = systemfunc("ls /mnt/e/etccdi/enavg/remap_mask_time/ssp585/txx/*.nc")
; print(f1)
f2 = systemfunc("ls /mnt/e/etccdi/enavg/remap_mask_time/ssp585/tnn/*.nc")
f3 = systemfunc("ls /mnt/e/etccdi/enavg/remap_mask_time/ssp585/rx5day/*.nc")

fg1 = addfile("/mnt/e/north-future/ec/input/tas-new/tas-gmst/tas-gmst-ssp585.nc","r")
r_trend=fg1->trend
r_trend=r_trend*10.0
; printVarSummary(r_trend)
r_trend_avg=dim_avg_n(r_trend,0)
; print(r_trend_avg)
r_trend_std=stddev(r_trend)

r_trend_all=new((/21/),float)
r_trend_all(0)=r_trend_avg
r_trend_all(1:20)=r_trend
; print(r_trend_all)
fg_std=addfile("/mnt/e/north-future/ec/input/iv/iv-all.nc", "r")
internal=fg_std->std
internal=internal*10.0
internal_avg=avg(internal)
print(internal_avg)

e1=addfile("/mnt/d/cmip6-outdate/north-future/GMST/all/ERA5.nc","r")
te1=e1->t2m(252:659,:,:)
te1=te1-273.15
te1!0 = "time"
te1!1 = "lat"
te1!2 = "lon"

e2=addfile("/mnt/d/cmip6-outdate/north-future/GMST/all/NCEP.nc","r")
te2=e2->air(396:803,:,:)
te2=te2-273.15
te2!0 = "time"
te2!1 = "lat"
te2!2 = "lon"
; printVarSummary(te2)

e3=addfile("/mnt/d/cmip6-outdate/north-future/GMST/all/JRA55.nc","r")
te3=e3->GSMT(252:659,:,:)
te3=te3-273.15
te3!0 = "time"
te3!1 = "lat"
te3!2 = "lon"
; printVarSummary(te3)

;
te=new((/3, 408, 180,360/), float)
te(0,:,:,:)=te1
te(1,:,:,:)=te2
te(2,:,:,:)=te3
; te(3,:,:,:)=te4
te!0 = "ens"
; printVarSummary(te)

rad=4.*atan(1.)/180.
clat= te&lat           
weight=cos(rad*clat)

te_obs=new((/3, 1/), float)
do k = 0, dimsizes(te(:,0,0,0))-1
    t_mean=wgt_areaave(te(k,:,:,:), weight, 1.0, 1)
    t_mean!0="time"
    tmp_annual_mean=month_to_annual(t_mean, 1)
    te_trend=trend_manken(tmp_annual_mean, False, 0)
    te_obs(k,:)=te_trend(1)*10.0
end do
; print(te_obs)

e5=addfile("/mnt/d/cmip6-outdate/north-future/GMST/HadCRUT-1850-2021.nc","r")
te5=e5->GMST(131:164)
te5_trend=trend_manken(te5, False, 0)
te5_trend1=te5_trend(1)*10.0*1.074
; print(te5_trend1)
e6=addfile("/mnt/d/cmip6-outdate/north-future/GMST/GISTEMP-1950-2021.nc","r")
te6=e6->GMST(31:64)
te6_trend=trend_manken(te6, False, 0)
te6_trend1=(te6_trend(1)+0.014/10.0)*10.0
; print(te6_trend1)

te_obs_all=new((/5, 1/), float)
te_obs_all(0:2,:)=te_obs
te_obs_all(3,:)=te5_trend1
te_obs_all(4,:)=te6_trend1

obs_avg=dim_avg_n(te_obs_all, 0)
obs_avg!0="ens"
print(obs_avg)

; ncdf1=addfile("/mnt/e/north-future/ec/pre/kj/obs/obs_avg.nc","c")
; ncdf1->trend=obs_avg
; print(ncdf1)

obs_std=stddev(te_obs_all)
print(obs_std)

obs_1o=obs_avg-internal_avg
obs_2o=obs_avg+internal_avg
obs_3o=obs_avg-(internal_avg*2)
obs_4o=obs_avg+(internal_avg*2)

n=dimsizes(f1)
r_project_hot=new((/3,n/),float)

do i = 0, n-1
    fb1=addfile(f1(i), "r")
    t21=fb1->txxETCCDI(65:150,0,0)
    t21_base=fb1->txxETCCDI(45:64,0,0)
    t21=t21-dim_avg_n_Wrap(t21_base,0)
    r_project_hot(0,i)=avg(t21)

    fb2=addfile(f2(i), "r")
    t22=fb2->tnnETCCDI(65:150,0,0)
    t22_base=fb2->tnnETCCDI(45:64,0,0)
    t22=t22-dim_avg_n_Wrap(t22_base,0)
    r_project_hot(1,i)=avg(t22)

    fb3=addfile(f3(i), "r")
    t23=fb3->rx5dayETCCDI(65:150,0,0)
    t23_base=fb3->rx5dayETCCDI(45:64,0,0)
    t23=t23-dim_avg_n_Wrap(t23_base,0)
    r_project_hot(2,i)=avg(t23)
end do

; printMinMax(r_project_txx(3,:), False)
; printMinMax(r_project_tnn(3,:), False)
; printMinMax(r_project_rx5day(3,:), False)

r_project_hot_avg=dim_avg_n(r_project_hot,1)
r_project_hot_std=dim_stddev_n(r_project_hot, 1)

raw_all=new((/3,5/),float)
do i=0,2
    raw_all(i,0)=r_project_hot_avg(i)+r_project_hot_std(i)*(-tcvf)
    raw_all(i,1)=r_project_hot_avg(i)+r_project_hot_std(i)*(-tcvf66)
    raw_all(i,2)=r_project_hot_avg(i)
    raw_all(i,3)=r_project_hot_avg(i)+r_project_hot_std(i)*(tcvf66)
    raw_all(i,4)=r_project_hot_avg(i)+r_project_hot_std(i)*(tcvf)
end do
print(raw_all)

; ncdf2=addfile("/mnt/e/north-future/ec/new-output/box/raw/ssp585/hot-ssp585.nc","c")
; ncdf2->EC=raw_all
; print(ncdf2)

r_txx_all=new((/n+1/),float)
r_txx_all(0)=r_project_hot_avg(0)
r_txx_all(1:n)=r_project_hot(0,:)
printVarSummary(r_txx_all)
r_tnn_all=new((/n+1/),float)
r_tnn_all(0)=r_project_hot_avg(1)
r_tnn_all(1:n)=r_project_hot(1,:)

r_rx5day_all=new((/n+1/),float)
r_rx5day_all(0)=r_project_hot_avg(2)
r_rx5day_all(1:n)=r_project_hot(2,:)

mask_all=new((/3,n+1/),float)
mask_all(0,:)=r_txx_all
mask_all(1,:)=r_tnn_all
mask_all(2,:)=r_rx5day_all


coef1_hot=escorc(r_trend, r_project_hot(0,:))
coef2_hot=escorc(r_trend, r_project_hot(1,:))
coef3_hot=escorc(r_trend, r_project_hot(2,:))

coef_all=new((/3,1/),float)
coef_all(0,:)=coef1_hot
coef_all(1,:)=coef2_hot
coef_all(2,:)=coef3_hot

; 相关系数的显著性检验
prob1=rtest(coef1_hot, dimsizes(r_trend), 0)
prob2=rtest(coef2_hot, dimsizes(r_trend), 0)
prob3=rtest(coef3_hot, dimsizes(r_trend), 0)

prob=new((/3,1/),float)
prob(0,:)=prob1
prob(1,:)=prob2
prob(2,:)=prob3
print(prob)

; 原始回归线
yline_original=new((/3,51/),float)
; yline_95=new((/3,n+1/),float)
; yline_5=new((/3,n+1/),float)
do i=0,2
    re1_hot=regline_stats(r_trend_all, mask_all(i,:))
    xline=ispan(0,50,1)/10.0
    r_b=re1_hot@b 
    yline_original(i,:)=r_b(0)+re1_hot*xline   
end do

; re1_hot=regline_stats(r_trend_all, mask_all(i,:))
; y5=re1_hot@YPI025 
; print(y5)
; y95=re1_hot@YPI975
; q=dim_pqsort(r_trend_all, 2)
; yline_5=y5(q)
; yline_95=y95(q)

EC_mean=new((/3/),float)
EC_std=new((/3/),float)
k_EC=new((/3/),float)
std_EC=new((/3/),float)
EC_all=new((/3,5/),float)
EC=new((/3/),float)
RRV=new((/3/),float)

do i=0,2
    EC_mean(i)=r_project_hot_avg(i)+coef_all(i,:)*r_project_hot_std(i)*r_trend_std*(obs_avg-r_trend_avg)/(r_trend_std^2+internal_avg^2)
    EC_std(i)=r_project_hot_std(i)*sqrt(1-coef_all(i,:)^2/(1+internal_avg^2/r_trend_std^2))
    std_EC(i)=sqrt(1-coef_all(i,:)^2/(1+internal_avg^2/r_trend_std^2))
    k_EC(i)=coef_all(i,:)*r_project_hot_std(i)*r_trend_std/(r_trend_std^2+internal_avg^2)
    EC_all(i,0)=EC_mean(i)+EC_std(i)*(-tcvf)
    EC_all(i,1)=EC_mean(i)+EC_std(i)*(-tcvf66)
    EC_all(i,2)=EC_mean(i)
    EC_all(i,3)=EC_mean(i)+EC_std(i)*(tcvf66)
    EC_all(i,4)=EC_mean(i)+EC_std(i)*(tcvf)
    RRV(i)=(1-EC_std(i)^2/r_project_hot_std(i)^2)*100
end do
; line 
yline_EC=new((/3,51/),float)
do i=0,2
    b_EC=r_project_hot_avg(i)-k_EC(i)*r_trend_avg 
    yline_EC(i,:)=b_EC+k_EC(i)*xline
end do
EC_all!0="index"
EC_all!1="value"
print(EC_all)
; 单个模式
trend_dif=obs_avg-r_trend
trend_dif!0="ens"
k_EC!0="index"
printVarSummary(k_EC)

; ncdf1=addfile("/mnt/e/north-future/ec/pre/kj/trend-dif-ssp585.nc","c")
; ncdf1->trend_dif=trend_dif
; print(ncdf1)

; ncdf3=addfile("/mnt/e/north-future/ec/new-output/box/c/ssp585/hot-ssp585.nc","c")
; ncdf3->EC=EC_all
; print(ncdf3)


wks = gsn_open_wks("eps","/mnt/e/north-future/new-p/original/fig2/hec/hec-all")
plot = new(3,graphic)  
plot_l=new(3,graphic)
plot_lo=new(3,graphic)
plot_r1=new(3,graphic)
plot_r2=new(3,graphic)
l_plot=new(3,graphic)
c_plot=new(3,graphic)
dum1=new(3,graphic)
dum2=new(3,graphic)
lo_plot=new(3,graphic)
co_plot=new(3,graphic)


res                    = True              ; plot mods desired
res@gsnDraw            = False             ; don't draw yet
res@gsnFrame           = False  

res@vpHeightF 	 = 0.4            ; change aspect ratio of plot
res@vpWidthF 	   = 0.4

res@tmXBLabelFontHeightF    = 0.015     ;set the size of x-Axis words
res@tmYLLabelFontHeightF    = 0.015    ;set the size of y-Axis words

res@tiXAxisString	 = "Trend of global warming (~S~o~N~C pre decade)" 
; res@tiYAxisString =""
res@trXMinF	  = 0.
; res@trYMinF	  = 0

res1=res
res1@xyMarkLineMode     = "Markers"
res1@xyMarkerSizeF      = 0.015
res1@xyMarkerThicknessF = 2.
res1@xyMonoMarkerColor  = False          ; Allow different colors for markers
res1@xyMarkerColors     = (/"black","skyblue","skyblue","skyblue","skyblue",\
                            "darkorange","darkorange","darkorange","darkorange",\     ;"darkorange"
                            "lightpink","lightpink","lightpink","lightpink",\
                            "orchid","orchid","orchid","orchid",\
                            "green","green","green","green"/)  
res1@xyMarkers          = (/16,12, 7, 4, 3,\
                            12,7,4, 3,\  ; 7,
                            12, 7, 4, 3,\
                            12,7,4,3, \
                            12,7,4,3/)  

res2=res
res2@xyMarkLineMode   = "Lines" 
res2@xyDashPattern      = 0
res2@xyLineThicknessF  = 3       
res2@xyLineColor        =  "red"

res3=res
res3@xyMarkLineMode   = "Lines" 
res3@xyDashPattern      = 0
res3@xyLineThicknessF  = 3       
res3@xyLineColor        =  "black"

line_res = True
line_res@gsLineColor = "black"                       ; 设置竖线颜色
line_res@gsLineThicknessF = 2.0                    ; 设置竖线粗细
line_res@gsLineDashPattern = 0 
; line_res@PostDraw            = True   

; yx=(/2.1,4.5,9.0,18,30,70,50,100,270/)
; yn=(/0.3,1.0,2.0,0,6,20,0,10,60/)

; title=(/"(a)SSP1-2.6_TXx","(b)SSP1-2.6_tnn","(c)SSP1-2.6_rx5day"/)
; title=(/"(d)SSP2-4.5_TXx","(e)SSP2-4.5_tnn","(f)SSP2-4.5_rx5day"/)
; title=(/"(g)SSP3-7.0_TXx","(h)SSP3-7.0_tnn","(i)SSP3-7.0_rx5day"/)
title=(/"(a)SSP5-8.5_TXx","(b)SSP5-8.5_TNn","(c)SSP5-8.5_RX5day"/)   
v=(/"~F8~D~F~TXx(~S~o~N~C)","~F8~D~F~TNn(~S~o~N~C)","~F8~D~F~RX5day(mm)"/) 
p=(/" (p<0.01)"," (p<0.01)"," (p<0.01)"/) 

do j = 0, dimsizes(plot) - 1
    res1@gsnLeftString  = title(j)
    res1@gsnRightString  = "r=" +sprintf("%4.2f", coef_all(j,:))+ p(j)+"~C~" +"RRV=" + sprintf("%4.0f", RRV(j)) + "%"
    res1@tiYAxisString =v(j)
    ; res1@trYMinF	  = yn(j)
    ; res1@trYMaxF	  = yx(j)
    plot(j) = gsn_csm_xy (wks,transpose((/r_trend_all,r_trend_all/)), \
                         transpose((/mask_all(j,:),mask_all(j,:)/)),res1)
    plot_l(j)=gsn_csm_xy (wks,xline,yline_EC(j,:),res2) 
    plot_lo(j)=gsn_csm_xy (wks,xline,yline_original(j,:),res3) 
    ; plot_r1(j)=gsn_csm_xy (wks,r_trend_all(:),yline_5(j,:),res3) 
    ; plot_r2(j)=gsn_csm_xy (wks,r_trend_all(:),yline_95(j,:),res3) 
    overlay(plot(j),plot_l(j))
    overlay(plot(j),plot_lo(j))
    ; overlay(plot(j),plot_r1(j))
    ; overlay(plot(j),plot_r2(j))
end do

ymin = new(3,float)
ymax = new(3,float)
do l=0,2
    getvalues plot(l)
        "trYMinF": ymin(l)
        "trYMaxF": ymax(l)
    end getvalues
end do

xl = (/(/obs_avg,obs_avg/)/)
; print(xl)
yl = (/ (/ymin(0),EC_all(0,2)/), \
        (/ymin(1),EC_all(1,2)/), \
        (/ymin(2),EC_all(2,2)/)/)

xbars = (/ (/obs_1o, obs_2o, obs_2o, obs_1o, obs_1o/), \
            (/obs_3o, obs_4o, obs_4o, obs_3o, obs_3o/)/)

ybars = (/ (/ymin(0),ymin(0),ymax(0),ymax(0),ymin(0)/), \
            (/ymin(1),ymin(1),ymax(1),ymax(1),ymin(1)/), \
            (/ymin(2),ymin(2),ymax(2),ymax(2),ymin(2)/)/)

xc = (/ 0,obs_avg/)

yc = (/ (/EC_all(0,2),EC_all(0,2)/), \
        (/EC_all(1,2),EC_all(1,2)/), \
        (/EC_all(2,2),EC_all(2,2)/)/)

xbars1 = (/0.,0.5,0.5,0.,0./)

xo_h=(/(/r_trend_avg,r_trend_avg/)/)
yo_h=(/ (/ymin(0),r_project_hot_avg(0)/), \
        (/ymin(1),r_project_hot_avg(1)/), \
        (/ymin(2),r_project_hot_avg(2)/)/)

xo_l=(/(/0,r_trend_avg/)/)
yo_l=(/ (/r_project_hot_avg(0),r_project_hot_avg(0)/), \
        (/r_project_hot_avg(1),r_project_hot_avg(1)/), \
        (/r_project_hot_avg(2),r_project_hot_avg(2)/)/)


ybars_95 = (/ (/EC_all(0,0),EC_all(0,0),EC_all(0,4),EC_all(0,4),EC_all(0,0)/), \
            (/EC_all(1,0),EC_all(1,0),EC_all(1,4),EC_all(1,4),EC_all(1,0)/), \
            (/EC_all(2,0),EC_all(2,0),EC_all(2,4),EC_all(2,4),EC_all(2,0)/)/)

ybars_87 = (/ (/EC_all(0,1),EC_all(0,1),EC_all(0,3),EC_all(0,3),EC_all(0,1)/), \
            (/EC_all(1,1),EC_all(1,1),EC_all(1,3),EC_all(1,3),EC_all(1,1)/), \
            (/EC_all(2,1),EC_all(2,1),EC_all(2,3),EC_all(2,3),EC_all(2,1)/)/)

lres=True
lres@gsLineDashPattern=11
lres@gsLineColor="red"
lres@gsLineThicknessF=2
lres@PostDraw=True
; lres@gsnDrawOrder = "PreDraw" 

lres1=True
lres1@gsLineDashPattern=11
lres1@gsLineColor="black"
lres1@gsLineThicknessF=2
lres1@PostDraw=True


gres             = True
gres@gsFillColor = "gray80"
; gres@gsFillOpacityF = 0.4 
gres@tfPolyDrawOrder="Predraw"

gres1             = True
gres1@gsFillColor = "gray90"
; gres1@gsFillOpacityF = 0.4 
; gres1@PostDraw=True
gres1@tfPolyDrawOrder = "PreDraw" 

; gres2             = True
; gres2@gsFillColor = "lightpink"
; ; gres@gsFillOpacityF = 0.8 
; gres2@tfPolyDrawOrder="Predraw"

; gres3            = True
; gres3@gsFillColor = "lightpink"
; gres3@gsFillOpacityF = 0.4 
; ; gres3@PostDraw=True
; gres3@tfPolyDrawOrder = "PreDraw" 


; do i=0,1
;    do j=0,11
;       dum((i*12)+j) = gsn_add_polygon(wks,plot(j),xbars(i,:),ybars(j,:),gres)
;    end do
; end do


do m=0,2
   dum1(m) = gsn_add_polygon(wks,plot(m),xbars(1,:),ybars(m,:),gres1)
   dum2(m) = gsn_add_polygon(wks,plot(m),xbars(0,:),ybars(m,:),gres)
   l_plot(m)=gsn_add_polyline(wks, plot(m), xl(:), yl(m,:), lres)
   c_plot(m)=gsn_add_polyline(wks, plot(m), xc, yc(m,:), lres)
   lo_plot(m)=gsn_add_polyline(wks, plot(m), xo_h(:), yo_h(m,:), lres1)
   co_plot(m)=gsn_add_polyline(wks, plot(m), xo_l, yo_l(m,:), lres1)
end do

mkres                =  True
txres               =  True
txres@txFontHeightF =  0.008         ;-- default size is HUGE!
txres@txFontColor   = "black"
txres@txJust        = "CenterLeft"

label=(/"CMIP6 MME","ACCESS-CM2","ACCESS-ESM1-5","CMCC-ESM2","CanESM5","EC-Earth3-CC","EC-Earth3-Veg-LR",\
        "EC-Earth3-Veg","EC-Earth3","FGOALS-g3","INM-CM4-8","INM-CM5-0",\
        "IPSL-CM6A-LR","MIROC6","MPI-ESM1-2-HR","MPI-ESM1-2-LR","MRI-ESM2-0",\
        "NESM3","NorESM2-LM","NorESM2-MM","TaiESM1"/)

color=(/"black","skyblue","skyblue","skyblue","skyblue",\
        "darkorange","darkorange","darkorange","darkorange",\
        "lightpink","lightpink","lightpink","lightpink",\
        "orchid","orchid","orchid","orchid",\
        "green","green","green","green"/)  

markers1=(/16,12, 7, 4, 3,12, 7, 4, 3,12, 7, 4, 3,12,7,4,3, 12,7,4,3/)

; 计算每行每列的位置
cols_per_row = 7
nrows = ceil(tolong((n+2)/tofloat(cols_per_row))) ; 计算需要的行数

; 定义起始位置和间距
start_x = 0.06      ; 横向起始位置
start_y = 0.305      ; 纵向起始位置
x_space = 0.13     ; 横向间距
y_space = 0.018     ; 纵向间距

do i = 0, 20
  row = toint(i / cols_per_row)           ; 当前行
  col = mod(i, cols_per_row)              ; 当前列
  
  xpos = start_x + (col * x_space)        ; 计算x位置
  ypos = start_y - (row * y_space)        ; 计算y位置
  
  mkres@gsMarkerThicknessF = 2.
  mkres@gsMarkerSizeF      = 0.008
  mkres@gsMarkerIndex      = markers1(i)
  mkres@gsMarkerColor      = color(i)
  
  gsn_polymarker_ndc(wks, xpos, ypos, mkres)
  gsn_text_ndc(wks, label(i), xpos+0.02, ypos, txres)
end do


 resP                     = True         ; modify the panel plot
 resP@gsnMaximize         = True       ; large format
 resP@gsnPanelLabelBar    = False       ; add common colorbar
 resP@gsnPaperOrientation = "portrait"   ; force portrait

;  resP@pmLabelBarOrthogonalPosF = 0.33
;  resP@lbLabelFontHeightF  = 0.007 
; draw panel without white space
;  resP@gsnPanelRowSpec = True 
;  drawNDCGrid(wks)
 gsn_panel(wks,plot,(/1,3/),resP)

end