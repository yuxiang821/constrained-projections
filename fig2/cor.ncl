
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

undef("add_filled_boxes")
procedure add_filled_boxes(wks, plot, data, levels, colors, prob, sig_level)
local dims, gnid, gnres, lnid, lnres, ii, jj, nrows, ncols, boxx, boxy, index, mkres, mkid
begin
    dims = dimsizes(data)
    nrows = dims(0)
    ncols = dims(1)
    gnid = new((/nrows, ncols/), "graphic")
    lnid = new((/nrows, ncols/), "graphic")  ; 用于存储轮廓线的图形对象
    mkid = new((/nrows, ncols/), "graphic")  ; 用于存储打点的图形对象
    
    ; 填充多边形的资源
    gnres = True
    
    ; 轮廓线的资源
    lnres = True
    lnres@gsLineColor = "black"     ; 边界线颜色
    lnres@gsLineThicknessF = 1.0    ; 边界线粗细
    
    ; 打点的资源
    mkres = True
    mkres@gsMarkerIndex = 3         ; 圆点
    mkres@gsMarkerSizeF = 0.015     ; 点的大小
    mkres@gsMarkerColor = "grey"    ; 点的颜色

    do jj = 0, nrows - 1
        do ii = 0, ncols - 1
            if (.not. ismissing(data(jj, ii))) then
                boxx = (/ii, ii, ii + 1, ii + 1, ii/)
                boxy = (/jj, jj + 1, jj + 1, jj, jj/)
                
                ; 添加填充色块
                index = closest_val(data(jj, ii), levels)
                gnres@gsFillColor = colors(index, :)
                gnid(jj, ii) = gsn_add_polygon(wks, plot, boxx, boxy, gnres)
                
                ; 添加格子边界线
                lnid(jj, ii) = gsn_add_polyline(wks, plot, boxx, boxy, lnres)

                ; 检查 p 值并添加打点
                if (prob(jj, ii) .lt. sig_level) then
                    mkx = ii + 0.5  ; 点的 x 坐标（单元中心）
                    mky = jj + 0.5  ; 点的 y 坐标（单元中心）
                    mkid(jj, ii) = gsn_add_polymarker(wks, plot, mkx, mky, mkres)
                end if
            end if
        end do
    end do
    
    ; 保存所有图形对象以避免内存回收
    plot@gnid = gnid    ; 填充区域
    plot@lnid = lnid    ; 边界线
    plot@mkid = mkid    ; 显著性打点
end


undef("add_labelbar")
procedure add_labelbar(wks,plot,colors,labels)
local nboxes, vph, vpw, nboxes, lbres, lbid, amres, annoid
begin
  nboxes = dimsizes(colors)

  getvalues plot                      ; Get plot size for use in
    "vpHeightF" : vph                 ; creating labelbar.
    "vpWidthF"  : vpw
  end getvalues

  lbres                    = True          ; labelbar only resources
  lbres@lbAutoManage       = False         ; Necessary to control sizes
  lbres@lbFillColors       = colors        ; labelbar colors
  lbres@lbMonoFillPattern  = True          ; Solid fill pattern
  lbres@lbLabelFontHeightF = 0.01          ; font height. Default is small
  lbres@lbLabelAlignment   = "BoxCenters"
  lbres@lbOrientation      = "Vertical"
  lbres@lbPerimOn          = False

  lbres@vpWidthF           = vpw * 0.2   ; Make labelbar slightly shorter than width of plot
  lbres@vpHeightF          = vph 

lbid = gsn_create_labelbar(wks,nboxes,labels,lbres)
amres                  = True
amres@amJust           = "BottomCenter"
amres@amParallelPosF   = 0.6     ; Centered about X axis
amres@amOrthogonalPosF = 0.5    ; Move labelbar down
plot@annoid = gsn_add_annotation(plot,lbid,amres)
end 

begin
f1 = systemfunc("ls /mnt/e/etccdi/enavg/remap_mask_time/ssp585/rx5day/*.nc")

n=dimsizes(f1)
r_project_hot=new((/n/),float)

do i = 0, n-1
    fb1=addfile(f1(i), "r")
    t21=fb1->rx5dayETCCDI(65:150,0,0)
    t21_base=fb1->rx5dayETCCDI(45:64,0,0)
    t21=t21-dim_avg_n_Wrap(t21_base,0)
    r_project_hot(i)=avg(t21)
end do

f2=addfile("/mnt/e/north-future/ec/input/tas-1850-2014/tas-ssp585.nc","r")
; time=f2->time
; time_str = cd_calendar(time, -3)
; print(time_str)
t=f2->tas
printVarSummary(t)

years=ispan(1850, 2014, 1)

start_years = ispan(1960, 2000, 10) 
end_years = ispan(1974, 2014, 10)     
; print(end_years)
n_start = dimsizes(start_years)
n_end = dimsizes(end_years)

corr_matrix = new((/n_end, n_start/), float)
corr_matrix@_FillValue=9.96921e+36
prob=new((/n_end, n_start/), float)
t_trend=new((/n, 2/), float)
do i = 0, n_start-1
    do j = 0, n_end-1
      if (start_years(i).lt.end_years(j)) then
        year_diff = end_years(j) - start_years(i)
        
        if (year_diff .le. 10) then
          corr_matrix(j,i) = corr_matrix@_FillValue
          prob(j,i) = prob@_FillValue
        else
          idx_start = ind(years.eq.start_years(i))
          idx_end = ind(years.eq.end_years(j))
          
          do m=0,n-1
            t_trend(m,:)=trend_manken(t(m,idx_start:idx_end),False,0)
          end do
          trend=t_trend(:,1)*10.0
          corr_matrix(j,i) = escorc(trend, r_project_hot)
          prob(j,i) =rtest(corr_matrix(j,i), dimsizes(trend), 0)
        end if
      else
        corr_matrix(j,i) = corr_matrix@_FillValue
        prob(j,i) = prob@_FillValue
      end if
    end do
end do
printVarSummary(corr_matrix)
print(corr_matrix(4,:))

c= corr_matrix(::-1,:)
printMinMax(c,False)
nrows   = dimsizes(c(:,0))
ncols   = dimsizes(c(0,:))

levels =(/-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7/)
rgba   = read_colormap_file("MPL_BrBG")
colors = (/rgba(40,:),rgba(56,:),rgba(68,:),rgba(72,:),rgba(82,:),rgba(92,:),rgba(100,:),rgba(108,:),rgba(116,:),rgba(124,:)/)

wks = gsn_open_wks("pdf","/mnt/e/north-future/new-p/original/fig2/cor/9")
res = True

res@gsnDraw     = False
res@gsnFrame    = False
res@gsnScale    = True
res@gsnMaximize = True

res@vpWidthF    = 0.4     ; Change aspect ratio such that
res@vpHeightF   = 0.4     ; boxes are square

res@trXMinF                 = 0 
res@trXMaxF                 = ncols
res@trYMinF                 = 0
res@trYMaxF                 = nrows

res@tmXBMajorOutwardLengthF = 0.
res@tmXBMajorLengthF        = 0.
res@tmYLMajorOutwardLengthF = 0.
res@tmYLMajorLengthF        = 0.

res@tmXTOn                  = False
res@tmYROn                  = False

res@tmXBMode                = "Explicit"
res@tmXBValues              = ispan(res@trXMinF,res@trXMaxF-1,1)
res@tmXBLabels              = (/"1960","1970","1980","1990","2000"/)
; res@tmXBLabelAngleF =45
res@tmYLMode               = "Explicit"
res@tmYLValues             = ispan(res@trYMinF,res@trYMaxF-1,1)
res@tmYLLabels             =  (/"1974","1984","1994","2004","2014"/)

res@tiMainString           = "Corr(trend of GMST,~F8~D~F~Rx5day_SSP3-7.0)"
res@tiMainFontHeightF      = 0.015

res@tmXBLabelFontHeightF   = 0.01
res@tmYLLabelFontHeightF   = 0.01
; res@gsnLeftString  = "(b)Corr(trend of GMST,~F8~D~F~rx5day_SSP5-8.5)"
res@tiXAxisString	 = "First year of segment" 
res@tiYAxisString  ="Last year of segment" 
res@tiXAxisFontHeightF=0.015
res@tiYAxisFontHeightF=0.015

plot = gsn_csm_blank_plot(wks,res)

sig_level = 0.05
add_filled_boxes(wks, plot, corr_matrix, levels, colors, prob, sig_level)
add_labelbar(wks, plot, colors, ""+levels)

draw(plot)
frame(wks)

end

