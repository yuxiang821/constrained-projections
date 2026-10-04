load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin
diri = "/mnt/e/north-future/uncertain/output/new/"
files = systemfunc("ls " + diri + "*nc")  ; 批量读取文件路径
print(files)

i_all= new((/9,151/), float)
m_all= new((/9,151/), float)
s_all= new((/9,151/), float)

do i = 0, 8
    files1 = files(i)
    f1 = addfile(files1, "r")
    i_all(i,:)= f1->I_variance
    m_all(i,:)= f1->M_variance
    s_all(i,:)= f1->S_variance
end do

i_all_s=new((/3,151/), float)
i_all_s(0,:)=i_all(7,:)
i_all_s(1,:)=i_all(5,:)
i_all_s(2,:)=i_all(2,:)

m_all_s=new((/3,151/), float)
m_all_s(0,:)=m_all(7,:)
m_all_s(1,:)=m_all(5,:)
m_all_s(2,:)=m_all(2,:)

s_all_s=new((/3,151/), float)
s_all_s(0,:)=s_all(7,:)
s_all_s(1,:)=s_all(5,:)
s_all_s(2,:)=s_all(2,:)


;2021-2040
i_e=dim_avg_n_Wrap(i_all_s(:,71:90),1)
m_e=dim_avg_n_Wrap(m_all_s(:,71:90),1)
s_e=dim_avg_n_Wrap(s_all_s(:,71:90),1)
print(s_e)
;2041-2060
i_m=dim_avg_n_Wrap(i_all_s(:,91:110),1)
m_m=dim_avg_n_Wrap(m_all_s(:,91:110),1)
s_m=dim_avg_n_Wrap(s_all_s(:,91:110),1)
print(s_m)
;2081-2100
i_l=dim_avg_n_Wrap(i_all_s(:,131:150),1)
m_l=dim_avg_n_Wrap(m_all_s(:,131:150),1)
s_l=dim_avg_n_Wrap(s_all_s(:,131:150),1)

input_all=new((/3,3,3/), float)

do i=0,2
    input_all(i,0,0)=i_e(i)*100
    input_all(i,0,1)=(i_e(i)+m_e(i))*100
    input_all(i,0,2)=(i_e(i)+m_e(i)+s_e(i))*100
    input_all(i,1,0)=i_m(i)*100
    input_all(i,1,1)=(i_m(i)+m_m(i))*100
    input_all(i,1,2)=(i_m(i)+m_m(i)+s_m(i))*100
    input_all(i,2,0)=i_l(i)*100
    input_all(i,2,1)=(i_l(i)+m_l(i))*100
    input_all(i,2,2)=(i_l(i)+m_l(i)+s_l(i))*100
end do

; print(input_all)
n=3
x= ispan(1,n,1)
titles = (/"(a)Uncertainty components of TXx", "(b)Uncertainty components of TNn", "(c)Uncertainty components of RX5day"/)

wks = gsn_open_wks("eps","/mnt/e/north-future/new-p/original/fig1/uncertain/uncertain-rx5day")
; plot = new(3,graphic)  
; bottom_plot= new(3,graphic)  
; mid_plot=new(3,graphic)
; top_plot=new(3,graphic)

res                    = True              ; plot mods desired
res@gsnDraw            = False             ; don't draw yet
res@gsnFrame           = False   

res@vpHeightF 	 = 0.4               ; change aspect ratio of plot
res@vpWidthF 	 = 0.6

res@trYMaxF   =100
res@trYMinF   =0
res@trXMinF   = 0
res@trXMaxF   = n+ 1

res@tmXBMode                 = "Explicit"
res@tmXBValues               = x
res@tmXBLabels               = (/"Near Term","Mid Term","Long Term"/)  ;,"SDII","RX5day","R95p
; res@gsnLeftString  = "(c)Long Term"
res@tiYAxisString	 = "Fraction of total variance(%)" 

res@tmXBLabelFontHeightF    = 0.012 
res@tmYLLabelFontHeightF    = 0.012 
res@tiYAxisFontThicknessF = 0.015
res@tiXAxisFontThicknessF = 0.015

res@gsnXYBarChart = True
res@gsnXYBarChartBarWidth = 0.4   

res1=res
res1@gsnXYBarChartColors = "orangered3"

res2=res
res2@gsnXYBarChartColors = "steelblue"

res3=res
res3@gsnXYBarChartColors = "green4"

; do i=0,2 
;     res@gsnLeftString  = titles(i)
;     plot(i)= gsn_csm_xy(wks,x,input_all(i,:,2),res)
;     bottom_plot(i)= gsn_csm_xy(wks,x,input_all(i,:,0),res1)
;     mid_plot(i)= gsn_csm_xy(wks,x,input_all(i,:,1),res2)
;     top_plot(i)= gsn_csm_xy(wks,x,input_all(i,:,2),res3)
;     overlay(plot(i),top_plot(i))    ; Long bars on outlined bars
;     overlay(plot(i),mid_plot(i))  
;     overlay(plot(i),bottom_plot(i))  
; end do 

res@gsnLeftString  = titles(2)
plot= gsn_csm_xy(wks,x,input_all(2,:,2),res)
bottom_plot= gsn_csm_xy(wks,x,input_all(2,:,0),res1)
mid_plot= gsn_csm_xy(wks,x,input_all(2,:,1),res2)
top_plot= gsn_csm_xy(wks,x,input_all(2,:,2),res3)
overlay(plot,top_plot)    ; Long bars on outlined bars
overlay(plot,mid_plot)  
overlay(plot,bottom_plot)

pgres                        =  True
pgres@gsFillColor            = "orangered3"      ;-- fill color
;pgres@gsFillOpacityF         =  0.3         ;-- set opacity of polygon
gsn_polygon_ndc(wks,(/0.2,0.24,0.24,0.2,0.2/),\
                    (/0.315,0.315,0.3,0.3,0.315/),\
                    pgres)

pgres1                        =  True
pgres1@gsFillColor            = "steelblue"      ;-- fill color
;pgres@gsFillOpacityF         =  0.3         ;-- set opacity of polygon
gsn_polygon_ndc(wks,(/0.4,0.44,0.44,0.4,0.4/),\
                    (/0.315,0.315,0.3,0.3,0.315/),\
                    pgres1)

pgres2                       =  True
pgres2@gsFillColor            = "green4"      ;-- fill color
;pgres@gsFillOpacityF         =  0.3         ;-- set opacity of polygon
gsn_polygon_ndc(wks,(/0.6,0.64,0.64,0.6,0.6/),\
                    (/0.315,0.315,0.3,0.3,0.315/),\
                    pgres2)

txres               =  True
txres@txFontHeightF =  0.01        ;-- default size is HUGE!
txres@txFontColor   = "black"
gsn_text_ndc(wks,"Internal Variability",0.3,0.305, txres)
gsn_text_ndc(wks,"Model Uncertainly",0.5,0.305, txres)
gsn_text_ndc(wks,"Scenario Uncertainly",0.705,0.305, txres)

; resP                     = True         ; modify the panel plot
; resP@gsnMaximize         = False        ; large format
; resP@gsnPanelLabelBar    = False       ; add common colorbar
; resP@gsnPaperOrientation = "portrait"   ; force portrait
; ;  resP@pmLabelBarOrthogonalPosF = 0.33
; ;  resP@lbLabelFontHeightF  = 0.007 
; ; draw panel without white space
; ;  resP@gsnPanelRowSpec = True 
; ;  drawNDCGrid(wks)
; gsn_panel(wks,plot,(/1,3/),resP)
draw(plot)
frame(wks)
end

