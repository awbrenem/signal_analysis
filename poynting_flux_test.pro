;Test poynting_flux.pro, which accepts an arbitrary direction for the "spinaxis" (less reliable)
;data, compared to rbsp_poynting_flux, which uses only XMGSE as that direction.

;The test will compare the following data:
;RBSPa at +20 mlat on 2014-08-27 near 07:42.
;Sees northward propagating chorus in B2 data


;NOTE: the results are identical and show northward propagating chorus.



timespan,'2014-08-27'
rbsp_load_efw_waveform,probe='a',datatype=['mscb2','eb2']


;Keep this > spinperiod for uvw-->mgse transformation
t0 = time_double('2014-08-27/07:00')
t1 = time_double('2014-08-27/08:00')


;reduce to times needed.
ebuvw = tsample('rbspa_efw_eb2',[t0,t1],times=tms)
mscbuvw = tsample('rbspa_efw_mscb2',[t0,t1],times=tms)
store_data,'rbspa_efw_eb2r',tms,ebuvw
store_data,'rbspa_efw_mscb2r',tms,mscbuvw



rbsp_uvw_to_mgse,'a','rbspa_efw_eb2r'
rbsp_uvw_to_mgse,'a','rbspa_efw_mscb2r'

tplot,['rbspa_efw_eb2r_mgse','rbspa_efw_mscb2r_mgse']


;-------------------------------------
;EMFISIS data

rbsp_load_emfisis,probe='a',coord='gse',cadence='1sec',level='l3'

;Some of the EMFISIS quicklook data extend beyond the day loaded.
;This messes things up later. Remove these data points now.
ttst = tnames('rbspa_emfisis_l3_1sec_gse_Mag',cnt)
if cnt eq 1 then time_clip,'rbspa_emfisis_l3_1sec_gse_Mag',t0,t1,replace=1,error=error
ttst = tnames('rbspa_emfisis_l3_1sec_gse_Mag',cnt)
if cnt eq 1 then time_clip,'rbspa_emfisis_l3_1sec_gse_Mag',t0,t1,replace=1,error=error


rbsp_load_spice_cdf_file,'a'


tinterpol_mxn,'rbspa_spinaxis_direction_gse','rbspa_emfisis_l3_1sec_gse_Mag'
get_data,'rbspa_spinaxis_direction_gse_interp',tt,wgse
rbsp_gse2mgse,'rbspa_emfisis_l3_1sec_gse_Mag',wgse,newname='rbspa_emfisis_l3_1sec_mgse_Mag'

;--------------------------------------


;For Pflux calculate further reduce data timerange to portion of single burst.
t0z = time_double('2014-08-27/07:42:02.600')
t1z = time_double('2014-08-27/07:42:06.300')


ebuvw = tsample('rbspa_efw_eb2r_mgse',[t0z,t1z],times=tms)
mscbuvw = tsample('rbspa_efw_mscb2r_mgse',[t0z,t1z],times=tms)
bo_reduced = tsample('rbspa_emfisis_l3_1sec_mgse_Mag',[t0z-20.,t1z+20.],times=tms2)

store_data,'rbspa_efw_eb2r_mgsez',tms,ebuvw
store_data,'rbspa_efw_mscb2r_mgsez',tms,mscbuvw
store_data,'bo_reduced_mgse',tms2,bo_reduced
tplot,['rbspa_efw_eb2r_mgsez','rbspa_efw_mscb2r_mgsez','bo_reduced_mgse']


Tshort = 0.001
Tlong = 0.01
;Tshort = 0.0001
;Tlong = 0.1


;poynting_flux,'rbspa_efw_mscb2r_mgsez','rbspa_efw_eb2r_mgsez',Tshort,Tlong,$
;  Bo='bo_reduced_mgse',baddirection=1

poynting_flux,'rbspa_efw_mscb2r_mgsez','rbspa_efw_eb2r_mgsez',Tshort,Tlong,$
  Bo='bo_reduced_mgse',baddirection=1



timespan,'2014-08-27/07:42:02',5,/sec

;Check to see if the filtered versions pick up the chorus
tplot,['rbspa_efw_eb2r_mgse','pflux_Ew_pfluxcoord']
tplot,['rbspa_efw_mscb2r_mgse','pflux_Bw_pfluxcoord']

fns = ['pflux_angle_pflux_Bo',$
 'pflux_angle_pflux_Bo_smoothed',$
 'pflux_angle_pflux_Bo_detrend',$
 'pflux_pfluxcoord',$
 'pflux_nospinaxis_pfluxcoord',$
 'pflux_Ew_pfluxcoord',$
 'pflux_Bw_pfluxcoord',$
 'pflux_inputcoord',$
 'pflux_nospinaxis_inputcoord',$
 'pflux_nospinaxis_inputcoord1',$
 'pflux_nospinaxis_pfluxcoord1']



for i=0,n_elements(fns)-1 do copy_data,fns[i],fns[i]+'_new'

split_vec,'pflux_pfluxcoord_new'
split_vec,'pflux_nospinaxis_pfluxcoord_new'

;plot component along Bo
tplot,['pflux_Ew_pfluxcoord_new','pflux_Bw_pfluxcoord_new','pflux_pfluxcoord_new_z','pflux_nospinaxis_pfluxcoord_new_z']




;plot component perp to Bo
tplot,['pflux_nospinaxis_pfluxcoord_new_y']


;compare full and nospinaxis perp components
tplot,['pflux_nospinaxis_pfluxcoord_new_y','pflux_pfluxcoord_new_y']



;--------------------------------------------------------
;Old method
rbsp_poynting_flux,'rbspa_efw_mscb2r_mgsez','rbspa_efw_eb2r_mgsez',Tshort,Tlong,$
  Bo='bo_reduced_mgse'

  ylim,'pflux_para',0,0
  ylim,'pflux_nospinaxis_para',0,0
  ylim,'pflux_nospinaxis_perp',0,0
  ylim,'pflux_perp?',0,0
  tplot,['pflux_Ew','pflux_Bw','pflux_para','pflux_nospinaxis_para']

  ;compare full and nospinaxis perp components
  tplot,['pflux_nospinaxis_perp','pflux_perp1']



;--------------------------------------------------------
  ;compare new and old

  tplot,['pflux_Ew_pfluxcoord_new','pflux_Ew']
  tplot,['pflux_Bw_pfluxcoord_new','pflux_Bw']
  tplot,['pflux_pfluxcoord_new_z','pflux_para']
  tplot,['pflux_nospinaxis_pfluxcoord_new_z','pflux_nospinaxis_para']
  tplot,['pflux_nospinaxis_pfluxcoord_new_y','pflux_nospinaxis_perp']
