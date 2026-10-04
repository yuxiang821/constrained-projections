load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
f1=addfile("/mnt/e/north-future/time-all/ssp370/txx-1960-2100.nc","r")
t1=f1->txxETCCDI(:,:,0,0)
f2=addfile("/mnt/e/north-future/time-all/ssp370/tnn-1960-2100.nc","r")
t2=f2->tnnETCCDI(:,:,0,0)
f3=addfile("/mnt/e/north-future/time-all/ssp370/rx5day-1960-2100.nc","r")
t3=f3->rx5dayETCCDI(:,:,0,0)
printVarSummary(t3)

f4=addfile("/mnt/e/north-future/ec/input/time/c-relative1995-2014/ssp370/txx-c-ssp370.nc","r")
t4=f4->txxETCCDI
f5=addfile("/mnt/e/north-future/ec/input/time/c-relative1995-2014/ssp370/tnn-c-ssp370.nc","r")
t5=f5->tnnETCCDI
f6=addfile("/mnt/e/north-future/ec/input/time/c-relative1995-2014/ssp370/rx5day-c-ssp370.nc","r")
t6=f6->rx5dayETCCDI

t4_lower =f4->txxETCCDI_lower
t4_upper=f4->txxETCCDI_upper

t5_lower =f5->tnnETCCDI_lower
t5_upper=f5->tnnETCCDI_upper

t6_lower =f6->rx5dayETCCDI_lower
t6_upper=f6->rx5dayETCCDI_upper

size=dimsizes(t1(:,0))
t_all=new((/3,size,151/),float)
t_all(0,:,:)=t1
t_all(1,:,:)=t2
t_all(2,:,:)=t3
t_all!0="index"

c_all=new((/3,86/),float)
c_all(0,:)=t4
c_all(1,:)=t5
c_all(2,:)=t6
c_all!0="index"

c_all_b=new((/3,2,86/),float)
c_all_b(0,0,:)=t4_lower
c_all_b(0,1,:)=t4_upper
c_all_b(1,0,:)=t5_lower
c_all_b(1,1,:)=t5_upper
c_all_b(2,0,:)=t6_lower
c_all_b(2,1,:)=t6_upper

t_all_avg=new((/3,151/),float)
t_h=new((/3,55/),float)
t_all_h=new((/3,19,55/),float)
t_all_f=new((/3,19,86/),float)
t_f=new((/3,86/),float)

do i=0,2
   t_all_avg(i,:)=dim_avg_n_Wrap(t_all(i,:,:),0)
   t_h(i,:)=t_all_avg(i,10:64)-dim_avg_n_Wrap(t_all_avg(i,45:64),0)
   t_f(i,:)=t_all_avg(i,65:150)-dim_avg_n_Wrap(t_all_avg(i,45:64),0)
   do j=0,18
     t_all_h(i,j,:)=t_all(i,j,10:64)-dim_avg_n_Wrap(t_all(i,j,45:64),0)
     t_all_f(i,j,:)=t_all(i,j,65:150)-dim_avg_n_Wrap(t_all(i,j,45:64),0)
   end do 
end do
copy_VarCoords(t_all(:,0,10:64),t_h)
copy_VarCoords(t_all(:,0,65:150),t_f)
copy_VarCoords(t_all(:,:,45:64),t_all_h)
copy_VarCoords(t_all(:,:,65:150),t_all_f)
printVarSummary(t_h)

r_RangePercentage= 1.645   
t_max_min_h=new((/3,2,55/),float)
do k=0,2
   t_max_min_h(k,0,:)=dim_avg_n_Wrap(t_all_h(k,:,:), 0)+dim_stddev_n_Wrap(t_all_h(k,:,:), 0)*r_RangePercentage
   t_max_min_h(k,1,:)=dim_avg_n_Wrap(t_all_h(k,:,:), 0)-dim_stddev_n_Wrap(t_all_h(k,:,:), 0)*r_RangePercentage
end do

t_max_min_f=new((/3,2,86/),float)
do k=0,2
   t_max_min_f(k,0,:)=dim_avg_n_Wrap(t_all_f(k,:,:), 0)+dim_stddev_n_Wrap(t_all_f(k,:,:), 0)*r_RangePercentage
   t_max_min_f(k,1,:)=dim_avg_n_Wrap(t_all_f(k,:,:), 0)-dim_stddev_n_Wrap(t_all_f(k,:,:), 0)*r_RangePercentage
end do

tc1=dim_avg_n_Wrap(t_f(2,:),0)
tc2=dim_avg_n_Wrap(t_max_min_f(2,1,:),0)
tc3=dim_avg_n_Wrap(t_max_min_f(2,0,:),0)


year1=ispan(1960,2014,1)   ;66:85
year2=ispan(2015,2100,1)   ;66:85

y1=ispan(2021,2040,1)   ;66:85
y2=ispan(2041,2060,1)   ;66:85
y3=ispan(2081,2100,1)   ;66:85

ra1 = regline_stats(y1, t_f(0,66:85))                ; degC/year
ra1=ra1*10
print(ra1)

ra2 = regline_stats(y1, c_all(0,66:85))                ; degC/year
ra2=ra2*10
print(ra2)

ra3 = regline_stats(y1, t_f(1,66:85))                ; degC/year
ra3=ra3*10
print(ra3)

ra4 = regline_stats(y1, c_all(1,66:85))                ; degC/year
ra4=ra4*10
print(ra4)

ra5 = regline_stats(y1, t_f(2,66:85))                ; degC/year
ra5=ra5*10
print(ra5)

ra6 = regline_stats(y1, c_all(2,66:85))                ; degC/year
ra6=ra6*10
print(ra6)

df1=ra1@nptxy-2
a=ra1@tval
prob1=(1 -betainc(df1/(df1+ra1@tval(1)^2), df1/2.0, 0.5) ) 

df2=ra2@nptxy-2
prob2=(1 - betainc(df2/(df2+ra2@tval(1)^2), df2/2.0, 0.5) ) 

df3=ra3@nptxy-2
prob3=(1 - betainc(df3/(df3+ra3@tval(1)^2), df3/2.0, 0.5) ) 

df4=ra4@nptxy-2
prob4=(1 - betainc(df4/(df4+ra4@tval(1)^2), df4/2.0, 0.5) ) 

df5=ra5@nptxy-2
prob5=(1 - betainc(df5/(df5+ra5@tval(1)^2), df5/2.0, 0.5) ) 

df6=ra6@nptxy-2
prob6=(1 - betainc(df6/(df6+ra6@tval(1)^2), df6/2.0, 0.5) ) 

print(prob1)
print(prob2)
print(prob3)
print(prob4)
print(prob5)
print(prob6)

wks = gsn_open_wks("eps","/mnt/e/north-future/new-p/original/supple/fig1/time-range/ssp370")
plot = new(3,graphic)   
plot1 = new(3,graphic)   
plot2 = new(3,graphic)   
plot3 = new(3,graphic)   
plot4 = new(3,graphic) 
plot5 = new(3,graphic)     
poly = new(9,graphic)

res                    = True              ; plot mods desired
res@gsnDraw            = False             ; don't draw yet
res@gsnFrame           = False   
res@xyCurveDrawOrder = "PostDraw"

res@vpHeightF 	 = 0.3              ; change aspect ratio of plot
res@vpWidthF 	   = 0.5

res@tmXBLabelFontHeightF    = 0.015     ;set the size of x-Axis words
res@tmYLLabelFontHeightF    = 0.015    ;set the size of y-Axis words

; res@tiYAxisString	 = "Anomaly(~S~o~N~C)"       
res@tiXAxisString	 = "YEAR" 

res@trXMinF	  = 1960
res@trXMaxF	  = 2100
  ; res@tmYLMode        = "Manual"	
  ; res@tmYLTickSpacingF=0.5
res@tmXBLabelFontHeightF    = 0.02 
res@tmYLLabelFontHeightF    = 0.02 
res@tiYAxisFontThicknessF = 0.02
res@tiXAxisFontThicknessF = 0.02
res@gsnLeftStringFontHeightF=0.018

res@xyMonoLineColor    = False           
res@xyLineColor       = "grey40"
res@xyLineThicknesses	 = 3.      ; line thicknesses
res@xyDashPattern	 = 0.

res1=res
delete(res1@xyLineColor)
res1@gsnXYFillColors = "gray"
res1@xyLineColor     = -1  
res1@gsnXYFillOpacityF = 0.8

res2=res
delete(res2@xyLineColor)
res2@gsnXYFillColors = "papayawhip"
res2@xyLineColor     = -1  
res2@gsnXYFillOpacityF = 0.8

res3=res
delete(res3@xyLineColor)
res3@gsnXYFillColors = "lightblue"
res3@xyLineColor     = -1  
res3@gsnXYFillOpacityF = 0.6

res4=res
res4@xyLineColor       ="orangered"

res5=res
res5@xyLineColor       ="dodgerblue"

title1=(/"(g)Time Series of TXx under SSP3-7.0","(h)Time Series of TNn under SSP3-7.0","(i)Time Series of RX5day under SSP3-7.0"/)
title2=(/"Anomaly(~S~o~N~C)","Anomaly(~S~o~N~C)","Anomaly(mm)"/)
maxf=(/10,25,30/)
minf=(/-4,-10,-10/)

do l=0,dimsizes(plot)-1
   res@trYMaxF   =maxf(l)
   res@trYMinF   =minf(l)
   res@gsnLeftString  = title1(l)
   res@tiYAxisString  =title2(l)
   plot(l)=gsn_csm_xy (wks,year1,t_h(l,:),res)    
   plot1(l)=gsn_csm_xy(wks,year1,t_max_min_h(l,:,:),res1) 
   plot2(l)=gsn_csm_xy (wks,year2,t_max_min_f(l,:,:),res2)
   plot3(l)=gsn_csm_xy (wks,year2,c_all_b(l,:,:),res3)
   plot4(l)=gsn_csm_xy (wks,year2,t_f(l,:),res4)
   plot5(l)=gsn_csm_xy (wks,year2,c_all(l,:),res5)
   overlay(plot(l),plot1(l))
   overlay(plot(l),plot2(l))
   overlay(plot(l),plot3(l))
   overlay(plot(l),plot4(l))
   overlay(plot(l),plot5(l))
end do

xbars = (/ (/2021, 2040, 2040, 2021, 2021/), \
            (/2041, 2060, 2060, 2041, 2041/), \
            (/2081, 2100, 2100, 2081, 2081/)/)

ybars = (/ (/minf(0),minf(0),maxf(0),maxf(0),minf(0)/), \
            (/minf(1),minf(1),maxf(1),maxf(1),minf(1)/), \
            (/minf(2),minf(2),maxf(2),maxf(2),minf(2)/)/)

;  Create the bars and attach to appropriate plot
gres             = True
gres@gsFillColor = "lightgray"
gres@gsFillOpacityF = 0.8 

do i=0,2
   do j=0,2
      poly((i*3)+j) = gsn_add_polygon(wks,plot1(j),xbars(i,:),ybars(j,:),gres)
   end do
end do

plres1                        =  True
plres1@gsLineThicknessF       =  3.0        ;-- set line thickness
plres1@gsLineColor            = "grey40"       ;-- set line color
gsn_polyline_ndc(wks,(/0.055,0.075/),(/0.58,0.58/),plres1)

plres2                        =  True
plres2@gsLineThicknessF       =  3.0        ;-- set line thickness
plres2@gsLineColor            = "orangered"       ;-- set line color
gsn_polyline_ndc(wks,(/0.055,0.075/),(/0.565,0.565/),plres2)
  
plres3                        =  True
plres3@gsLineThicknessF       =  3.0        ;-- set line thickness
plres3@gsLineColor            = "dodgerblue"       ;-- set line color
gsn_polyline_ndc(wks,(/0.055,0.075/),(/0.55,0.55/),plres3)

txres               =  True
txres@txFontHeightF =  0.006        ;-- default size is HUGE!
txres@txFontColor   = "black"
gsn_text_ndc(wks,"Historical",0.095,0.58, txres)
gsn_text_ndc(wks,"SSP3-7.0",0.095,0.565, txres)
gsn_text_ndc(wks,"Constrained SSP3-7.0",0.119,0.55, txres)

; drawNDCGrid(wks)
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