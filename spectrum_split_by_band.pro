;Divide spectral data into frequency bands. Return tplot variables of each band and properties of each band. 
;E.g. separate an RBSP spectrum into chorus upper and lower bands. 





;spec_tvar --> input spectrum tplot variable [ntimes, nvals, nfreqs]
;split_tvar --> input tplot variable with split lines  [ntimes, nlines]
;   e.g. to split a spectrum in two --> [ntimes, 1]. 
;        to split into 4 segments --> [ntimes, 3]  (e.g. 0.1*fce, 0.5*fce, 1*fce). 
;        Code will return +1 spectra 
;        Split lines must be in ascending frequency order. 

;chnames --> keyword that can be set to help name the output tplot variables.
;Size should be nlines+2, where the extra elements are '0' and 'maxfreq' 
;   e.g. channelnames = ['0','0.1fce','0.5fce','fce','maxfreq']
;   output would be, for example, 'spec_0-0.1fce'
;
;   
;wv --> set to return wave values of final spectra. Size [nlines,4]. Includes:
;wave_vals[*,0] --> total spectral amp/power
;wave_vals[*,1] --> max value of spectral amp/power
;wave_vals[*,2] --> median value of spectral amp/power
;wave_vals[*,3] --> average value of spectral amp/power


;Example usage:
;spec_tvar = 'spectmp' --> input spectral data (x,y,v)
;channelnames = ['0','20','0.1fce','0.5fce','fce','7300Hz']
;split_tvar = 'fces'  --> tplot var containing all the fce lines (20, 0.1fce, 0.5fce, fce)
;spectrum_split_by_band,spec_tvar,split_tvar,chnames=channelnames,wv=wave_valsE




pro spectrum_split_by_band, spec_tvar, split_tvar, chnames=channelnames, wv=wave_vals, testing=testtmp


 
  get_data,spec_tvar,data=specv,dlim=dlim,lim=lim
  get_data,split_tvar,data=splitv
  
  
  
  ;-----------------------------------------------------------------------------------  
  ;Make sure spec_tvar and split_tvar have the same cadence and number of data points
  ;-----------------------------------------------------------------------------------

  sz1 = n_elements(specv.x)
  sz2 = n_elements(splitv.x)
  
  if sz1 ne sz2 then tinterpol_mxn,split_tvar,spec_tvar,newname=split_tvar+'_tmp' $
                else copy_data,split_tvar,split_tvar+'_tmp'
  
  get_data,split_tvar+'_tmp',data=splitv
  splitv = splitv.y
  
  

  
  
  ;-----------------------------------------------------------------------------------
  ;If split_tvar is of size [n] then recast it as [n,1]
  ;-----------------------------------------------------------------------------------

  tst = size(splitv,/n_dimensions)
  if tst eq 1 then splitv = reform(splitv,n_elements(splitv),1)
  
  nlines = size(splitv,/dimensions)
  nlines = nlines[1]
  
  
  wave_vals = fltarr(nlines+1,4)
  ;wave_vals[*,0] --> total spectral amp/power
  ;wave_vals[*,1] --> max value of spectral amp/power
  ;wave_vals[*,2] --> median value of spectral amp/power
  ;wave_vals[*,3] --> average value of spectral amp/power
  
  
  
  ;---------------------------------------------------------------------------------------------
  ;Main loop. For each combination of adjacent lines create the spectra that fall within them. 
  ;---------------------------------------------------------------------------------------------
  
  for i=0,nlines do begin
    specv_tmp = specv.y
    specv_tmp[*] = !values.f_nan
    
    
    for qq=0,sz1-1 do begin
  
      ;zero freq to first line freq  
      if i eq 0. then goo = where((specv.v ge 0.) and (specv.v lt splitv[qq,i]))
      ;last line freq to max spectral freq
      if i eq nlines then goo = where((specv.v ge splitv[qq,i-1]) and (specv.v lt max(specv.v)))
      ;freqs b/t adjacent lines
      if (i gt 0.) and (i lt nlines) then goo = where((specv.v ge splitv[qq,i-1]) and (specv.v lt splitv[qq,i]))
  
  
      if goo[0] ne -1 then specv_tmp[qq,goo] = 1    
    endfor
  
    specvresult = specv_tmp*specv.y
  
  
    ;get rid of zero spec values which screw up the below calculations.
    goo = where(specvresult eq 0.)
    if goo[0] ne -1 then specvresult[goo] = !values.f_nan
  
  
    wave_vals[i,0] = total(specvresult,/nan)
    wave_vals[i,1] = max(specvresult,/nan)
    wave_vals[i,2] = median(specvresult)
    wave_vals[i,3] = mean(specvresult,/nan)
   
  
    if not keyword_set(channelnames) then store_data,'spec_'+strtrim(floor(i),2),specv.x,specvresult,specv.v,dlim=dlim,lim=lim
    if keyword_set(channelnames) then store_data,'spec_'+channelnames[i]+'-'+channelnames[i+1],specv.x,specvresult,specv.v,dlim=dlim,lim=lim 
  
  endfor
  
  
  ;Set zero values from wave_vals array to NaNs.
  goo = where(wave_vals eq 0.)
  if goo[0] ne -1 then wave_vals[goo] = !values.f_nan
  
  
  
   
  
  ;----------------------------------------------------
  ;The following is for testing 
  ;----------------------------------------------------
  
  if keyword_set(testtmp) then begin 
    nametmp = strarr(n_elements(channelnames))
    for i=0,n_elements(channelnames)-2 do begin
      nametmp[i] = 'spec_'+channelnames[i]+'-'+channelnames[i+1]+'_comb'
      store_data,nametmp[i],data=['spec_'+channelnames[i]+'-'+channelnames[i+1],'fces']
    endfor
    
    rbsp_efw_init
    loadct,39
    tplot,nametmp
    stop
  endif
  
  ;----------------------------------------------------

end




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




