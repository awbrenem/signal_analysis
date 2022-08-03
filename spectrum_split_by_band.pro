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
;wave_vals[*,4] --> 0.25% quartile of spectral amp/power
;wave_vals[*,5] --> 0.50% quartile of spectral amp/power (same as median)
;wave_vals[*,6] --> 0.75% quartile of spectral amp/power


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
  
  
  wave_vals = fltarr(nlines+1,7)

  
  
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

    ;------------------------------------------------------------------------
    ;Quartiles
    ;------------------------------------------------------------------------

    ;Get rid of NaN values when calculating quartiles
    goo = where(finite(specvresult) eq 1)
    if goo[0] ne -1 then begin
      pttmp = cgPercentiles(specvresult[goo],Percentiles=[0.25,0.5,0.75])
      wave_vals[i,4] = pttmp[0]
      wave_vals[i,5] = pttmp[1]
      wave_vals[i,6] = pttmp[2]
    endif





   
  
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
      options,'fces','colors',250
      options,'fces','thick',2
    endfor
 
     
    rbsp_efw_init
    loadct,39
    ylim,nametmp,0.1,8000,1
    tplot,nametmp
    stop
  endif
  
  ;----------------------------------------------------

end


