;Divide spectral data into frequency bands. Return tplot variables of each band and properties of each band. 
;E.g. separate an RBSP spectrum into chorus upper and lower bands. 


;spec_tvar --> input spectrum tplot variable [ntimes, nvals, nfreqs]
;split_tvar --> input tplot variable with split lines  [ntimes, nlines]
;   e.g. to split a spectrum in two --> [ntimes, 1]. 
;        to split into 4 segments --> [ntimes, 3]  (e.g. 0.1*fce, 0.5*fce, 1*fce). 
;        Code will return nlines+1 spectra 
;        Split lines must be in ascending frequency order. 


pro spectrum_split_by_band, spec_tvar, split_tvar




;***************

timespan,'2013-01-01'
;cdf2tplot,'~/Downloads/rbsp-b_magnetometer_4sec-gse_emfisis-l3_20130101_v1.3.6.cdf'
;cdf2tplot,'~/Downloads/rbspb_efw-l2_spec_20130101_v02.cdf'


;rbsp_load_efw_spec_l2,probe='b'
;rbsp_load_emfisis,probe=probe,level='l3',coord='gsm'

get_data,'rbspb_efw_spec64_scmw',data=dd

fcevals = fltarr(n_elements(dd.x),2)
fcevals[*,0] = 100.
fcevals[*,1] = 1000.
store_data,'fces',dd.x,fcevals


spec_tvar = 'rbspb_efw_spec64_scmw'
split_tvar = 'fces'


;**************


get_data,spec_tvar,data=specv,dlim=dlim,lim=lim
get_data,split_tvar,data=splitv




;Make sure spec_tvar and split_tvar have the same cadence and number of data points
sz1 = n_elements(specv.x)
sz2 = n_elements(splitv.x)
if sz1 ne sz2 then tinterpol_mxn,split_tvar,spec_tvar,newname=spec_tvar+'_tmp' else copy_data,split_tvar,split_tvar+'_tmp'
get_data,split_tvar+'_tmp',data=splitv
splitv = splitv.y




;If split_tvar is of size [n] then recast it as [n,1]
nlines = size(splitv,/n_dimensions)
if nlines eq 1 then splitv = reform(splitv,n_elements(splitv),1)







for i=0,nlines do begin
  specv_tmp = specv.y
  specv_tmp[*] = !values.f_nan
  for qq=0,sz1-1 do begin
  
    if i eq 0. then goo = where((specv.v ge 0.) and (specv.v lt splitv[qq,i]))
    if i eq nlines then goo = where((specv.v ge splitv[qq,i-1]) and (specv.v lt max(specv.v)))
    if (i gt 0.) and (i lt nlines) then goo = where((specv.v ge splitv[qq,i-1]) and (specv.v lt splitv[qq,i]))

    if goo[0] ne -1 then specv_tmp[qq,goo] = 1    
  
  
  endfor

  store_data,'spec_'+strtrim(floor(i),2),specv.x,specv_tmp*specv.y,specv.v,dlim=dlim,lim=lim

endfor

rbsp_efw_init
options,['spec_0','spec_1'],'spec',1
ylim,[spec_tvar,'spec_0','spec_1'],1,7000,1
;zlim,'spec_?',0,1,0
;zlim,'spec_?',0,1,0
loadct,39
tplot,[spec_tvar,'spec_2','spec_1','spec_0']

stop

;
;;limit spectral data to chorus only (upper and lower band separately), and to
;;values outside of chorus range ("other")
;spectmpL = spectmp & spectmpU = spectmp & spectmpO = spectmp
;;Remove zero values that screw up average and median calculation
;goo = where(spectmp eq 0.)
;if goo[0] ne -1 then spectmpL[goo] = !values.f_nan
;if goo[0] ne -1 then spectmpU[goo] = !values.f_nan
;if goo[0] ne -1 then spectmpO[goo] = !values.f_nan
;
;if finite(fcetmp[0]) and is_struct(specv) then begin
;  for qq=0,n_elements(fcetmp)-1 do begin $
;    gooL = where((specv.v le fcetmp[qq]/10.) or (specv.v ge fcetmp[qq]/2.)) & $
;    gooU = where((specv.v lt fcetmp[qq]/2.)  or (specv.v gt fcetmp[qq])) & $
;    gooO = where((specv.v gt fcetmp[qq]/10.) or (specv.v lt 20.)) & $
;    if gooL[0] ne -1 then spectmpL[qq,gooL] = !values.f_nan & $
;    if gooU[0] ne -1 then spectmpU[qq,gooU] = !values.f_nan & $
;    if gooO[0] ne -1 then spectmpO[qq,gooO] = !values.f_nan
;endfor
;endif
;
;if finite(fcetmp[0]) then begin
;  store_data,'tmpL',tt,spectmpL,specv.v,dlim=dlim,lim=lim
;  store_data,'tmpU',tt,spectmpU,specv.v,dlim=dlim,lim=lim
;  store_data,'tmpO',tt,spectmpO,specv.v,dlim=dlim,lim=lim
;
;  ;vague idea of fill factor
;  totalchorusspecL_E = total(spectmpL,/nan)
;  totalchorusspecU_E = total(spectmpU,/nan)
;  totalnonchorusspec_E = total(spectmpO,/nan)
;
;  maxchorusspecL_E = max(spectmpL,/nan)
;  maxchorusspecU_E = max(spectmpU,/nan)
;  maxnonchorusspec_E = max(spectmpO,/nan)
;
;  avgchorusspecL_E = mean(spectmpL,/nan)
;  avgchorusspecU_E = mean(spectmpU,/nan)
;  avgnonchorusspec_E = mean(spectmpO,/nan)
;
;  medianchorusspecL_E = median(spectmpL)
;  medianchorusspecU_E = median(spectmpU)
;  mediannonchorusspec_E = median(spectmpO)
;
;  if totalchorusspecL_E eq 0. then totalchorusspecL_E = !values.f_nan
;  if totalchorusspecU_E eq 0. then totalchorusspecU_E = !values.f_nan
;  if totalnonchorusspec_E eq 0. then totalnonchorusspec_E = !values.f_nan
;  if maxchorusspecL_E eq 0. then maxchorusspecL_E = !values.f_nan
;  if maxchorusspecU_E eq 0. then maxchorusspecU_E = !values.f_nan
;  if maxnonchorusspec_E eq 0. then maxnonchorusspec_E = !values.f_nan
;  if avgchorusspecL_E eq 0. then avgchorusspecL_E = !values.f_nan
;  if avgchorusspecU_E eq 0. then avgchorusspecU_E = !values.f_nan
;  if avgnonchorusspec_E eq 0. then avgnonchorusspec_E = !values.f_nan
;  if medianchorusspecL_E eq 0. then medianchorusspecL_E = !values.f_nan
;  if medianchorusspecU_E eq 0. then medianchorusspecU_E = !values.f_nan
;  if mediannonchorusspec_E eq 0. then mediannonchorusspec_E = !values.f_nan
;
;
;  print,'Totals ',totalnonchorusspec_E,totalchorusspecL_E,totalchorusspecU_E
;  print,'Max ',maxnonchorusspec_E,maxchorusspecL_E,maxchorusspecU_E
;  print,'Avg ',avgnonchorusspec_E,avgchorusspecL_E,avgchorusspecU_E
;  print,'Median ',mediannonchorusspec_E,medianchorusspecL_E,medianchorusspecU_E
;
;endif else begin
;  totalchorusspecL_E = !values.f_nan
;  totalchorusspecU_E = !values.f_nan
;  totalnonchorusspec_E = !values.f_nan
;
;  maxchorusspecL_E = !values.f_nan
;  maxchorusspecU_E = !values.f_nan
;  maxnonchorusspec_E = !values.f_nan
;
;  avgchorusspecL_E = !values.f_nan
;  avgchorusspecU_E = !values.f_nan
;  avgnonchorusspec_E = !values.f_nan
;
;  medianchorusspecL_E = !values.f_nan
;  medianchorusspecU_E = !values.f_nan
;  mediannonchorusspec_E = !values.f_nan
;
;
;endelse
;



end
