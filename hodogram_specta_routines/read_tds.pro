;+
;*****************************************************************************************
;
;  FUNCTION :   read_tds
;  PURPOSE  :   Read multiple ascii files generated by fatds, sctds, swttds, mvtds...and perform
;				various wavelet calculations on the data
;				Returns a structure with TDS waveform data
;
;  CALLED BY:   
;               
;
;  CALLS:		wavelet.pro
;               
;
;  REQUIRES:    
;               
;
;  INPUT:
;               
;
;  EXAMPLES:    x = read_tds()
;               
;
;  KEYWORDS:    PATH : sets the root directory to search for ascii files. 
;					  Defaults to path='~/Desktop/swidl-0.1.3/tds_data/'
;				FILENAME : array of filenames to read (don't include path)
;				WAVELET_VALS : set this keyword to have function return the following values as a function of
;					time for each channel: 
;						bandwidth (Hz) 
;						max power (mV/m)^2 
;						freq of max power (Hz) 
;						deltaf/f		
;				SLIDESPEC_VALS : same as wavelet_vals but using a sliding spectrogram
;				ARRAY : pass array of data into program. If set, read_tds doesn't read in any files, it only
;					performs wavelet calculations on the inputted array. Setting the wavelet_vals keyword
;					is redundant in this case. Array has format [x,y]. x is the waveform and y is the 
;					number of different waveforms. 
;				DURATION : duration of waveform sample in (sec). Only needed if array keyword is set.
;
;   CHANGED:  1)  Added more accurate unix times [05/14/2010   v1.0.0]
;		2) (OBSOLETE) Added check to make sure all files used are of length 16384, raise warning if not [Paradise - 07/01/2011   v1.0.1] 
; 		3) Fixed error where when building date_time, msecs were added, but the seconds string never increased. [Paradise - 07/05/2011   v1.0.2]
;		4) Program can now handle arrays of different sizes. The size of each structure element now corresponds to the size of the largest
;			array [Breneman - 07/23/2011]
;
;   NOTES:      Bandwidth is calculated with the e-folding time. i.e. when the max power falls by 1/e^2
;               
;
;   CREATED:  11/19/2009
;   CREATED BY:  Aaron W. Breneman
;    LAST MODIFIED:  07/05/2011   v1.0.2
;    MODIFIED BY: Adiv Y. Paradise
;
;*****************************************************************************************
;-
;**************************************************************************

function read_tds,filename=filenameT,path=path,wavelet_vals=wavelet_vals,slidespec_vals=slidespec_vals,array=array,duration=duration


CATCH, Error_status    
IF Error_status NE 0 THEN BEGIN  
   print,'Error occurred, returning'
   CATCH, /CANCEL  
   return,-1
ENDIF  


if ~keyword_set(path) then path='~/Desktop/swidl-0.1.3/tds_data/'

if not keyword_set(array) then begin
if not keyword_set(filenameT) then filename = dialog_pickfile(path=path,filter='*.txt',/multiple_files) $
else filename = path + filenameT

nfiles = n_elements(filename)

;-------------------------------------------------------------------
;FIND HOW MANY LINES ARE IN THE LARGEST FILE. USED TO SET ARRAY SIZE
;-------------------------------------------------------------------

sz = 0.
junk = ''

good = intarr(n_elements(filename))
for i=0,n_elements(good)-1 do good[i] = query_ascii(filename[i])

tst = where(good eq 1)
if tst[0] eq -1 then begin
	print,'NONE OF THE REQUESTED FILES FOUND'
	return,-1
end


for i=0,n_elements(filename)-1 do begin
	tmp = query_ascii(filename[i],info)  ;DON'T DELETE

	if good[i] eq 1 then begin

		openr,lun,filename[i],/get_lun
		readf,lun,junk

		if (strmid(junk,0,22) eq 'STEREO/SWAVES TDS DATA') then offset = 23.	


		close,lun
		free_lun,lun
		if sz lt (info.lines - offset) then sz = info.lines - offset

	endif

endfor

;-------------------------------------------------------
;GRAB THE DATE AND INITIAL TIME OF CAPTURE FROM FILENAME
;-------------------------------------------------------

mths = ['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
date_time = strarr(nfiles)

for i=0,nfiles-1 do begin
	tmp = strpos(filename[i],'STEREO',/reverse_search)
	t2 = strmid(filename[i],tmp+15,20)  ;ex 08Apr2007_210450.617

	t3 = where(mths eq strlowcase(strmid(t2,2,3)))

	t3 = t3[0]
	t3 = t3+1

	if t3 lt 10. then t3 = '0' + strtrim(t3,2) else t3 = strtrim(t3,2)
	date_time[i] = strmid(t2,5,4) + '-' + t3 + '-' + strmid(t2,0,2) + '/' + strmid(t2,10,2) + ':' + strmid(t2,12,2) + ':' + strmid(t2,14,6)
endfor




;--------------------------------------------------------------------------------------------
;check to see if data are antenna coord. If so, then I have to include the pseudo dipole data
;--------------------------------------------------------------------------------------------

swtds = strmatch(filename[0],'*ANT*')

time = fltarr(sz)
Ex_raw = fltarr(sz)
Ey_raw = fltarr(sz)
Ez_raw = fltarr(sz)
Ex_cor = fltarr(sz)
Ey_cor = fltarr(sz)
Ez_cor = fltarr(sz)
if swtds eq 1. then pseudo_raw = fltarr(sz)
;----
endif else begin ;not keyword_set(array)
	tst = size(array,/dimensions)
	swtds = 0.
	sz = tst[0]
	wavelet_vals = 1.  ;set the keyword
	if n_elements(tst) eq 2 then nfiles = tst[1] else nfiles = 1.
endelse
;----

;------------------------
;DEFINE STRUCTURE
;------------------------

;----
if not keyword_set(wavelet_vals) or not keyword_set(slidespec_vals) then begin
	big_struct = {time:fltarr(sz),$
				date_time:strarr(sz),$
				Ex_raw:fltarr(sz),$
				Ey_raw:fltarr(sz),$
				Ez_raw:fltarr(sz),$
				Ex_cor:fltarr(sz),$
				Ey_cor:fltarr(sz),$
				Ez_cor:fltarr(sz),$
				filename:'' $
	}
	if swtds eq 1. then str_element,big_struct,'pseudo_raw',fltarr(sz),/ADD_REPLACE
endif
;----

if keyword_set(wavelet_vals) and not keyword_set(array) then begin

;dt_dummy = 0.1
;wavex = WAVELET(ex_cor,dt_dummy,PERIOD=period,COI=coi,/PAD,param=6)
wavex = wv_cwt(ex_cor,'Morlet',6,scale=period,/double,dscale=0.1)


nfreqs = n_elements(period)

	big_struct = {time:fltarr(sz),$
				date_time:strarr(sz),$
				Ex_raw:fltarr(sz),$
				Ey_raw:fltarr(sz),$
				Ez_raw:fltarr(sz),$
				Ex_cor:fltarr(sz),$
				Ey_cor:fltarr(sz),$
				Ez_cor:fltarr(sz),$
				Ex_max_power:fltarr(sz),$
				Ey_max_power:fltarr(sz),$
				Ez_max_power:fltarr(sz),$
				Ex_bandwidth:fltarr(sz),$
				Ey_bandwidth:fltarr(sz),$
				Ez_bandwidth:fltarr(sz),$
				Ex_timewidth:fltarr(nfreqs),$
				Ey_timewidth:fltarr(nfreqs),$
				Ez_timewidth:fltarr(nfreqs),$
				Ex_freq_of_max_power:fltarr(sz),$
				Ey_freq_of_max_power:fltarr(sz),$
				Ez_freq_of_max_power:fltarr(sz),$
				Ex_deltaf_to_f:fltarr(sz),$
				Ey_deltaf_to_f:fltarr(sz),$	
				Ez_deltaf_to_f:fltarr(sz),$				
				filename:'' $
}
	if swtds eq 1. then begin
		str_element,big_struct,'pseudo_raw',fltarr(sz),/ADD_REPLACE
		str_element,big_struct,'pseudo_max_power',fltarr(sz),/ADD_REPLACE
		str_element,big_struct,'pseudo_bandwidth',fltarr(sz),/ADD_REPLACE
		str_element,big_struct,'pseudo_freq_of_max_power',fltarr(sz),/ADD_REPLACE
		str_element,big_struct,'pseudo_deltaf_to_f',fltarr(sz),/ADD_REPLACE
	endif
endif
;----
if keyword_set(slidespec_vals) and not keyword_set(array) then begin

	step = 0.1   ;0.1-->good freq res, 0.001-->good time res
	step_overlap = 0.8 ;--> decreasing this increases the timestep
	foo = slide_spec(time,ex_cor,step,step_overlap,/db,iwindow=1,time_index=time_index,freq_bins=freqs) 
	nfreqs = n_elements(freqs)
	ntimes = n_elements(time_index)

	big_struct = {time:fltarr(sz),$
				date_time:strarr(sz),$
				Ex_raw:fltarr(sz),$
				Ey_raw:fltarr(sz),$
				Ez_raw:fltarr(sz),$
				Ex_cor:fltarr(sz),$
				Ey_cor:fltarr(sz),$
				Ez_cor:fltarr(sz),$
				Ex_max_power:fltarr(ntimes),$
				Ey_max_power:fltarr(ntimes),$
				Ez_max_power:fltarr(ntimes),$
				Ex_bandwidth:fltarr(ntimes),$
				Ey_bandwidth:fltarr(ntimes),$
				Ez_bandwidth:fltarr(ntimes),$
				Ex_timewidth:fltarr(nfreqs),$
				Ey_timewidth:fltarr(nfreqs),$
				Ez_timewidth:fltarr(nfreqs),$
				Ex_freq_of_max_power:fltarr(ntimes),$
				Ey_freq_of_max_power:fltarr(ntimes),$
				Ez_freq_of_max_power:fltarr(ntimes),$
				Ex_deltaf_to_f:fltarr(ntimes),$
				Ey_deltaf_to_f:fltarr(ntimes),$	
				Ez_deltaf_to_f:fltarr(ntimes),$				
				filename:'' $
}
	if swtds eq 1. then begin
		str_element,big_struct,'pseudo_raw',fltarr(sz),/ADD_REPLACE
		str_element,big_struct,'pseudo_max_power',fltarr(sz),/ADD_REPLACE
		str_element,big_struct,'pseudo_bandwidth',fltarr(sz),/ADD_REPLACE
		str_element,big_struct,'pseudo_freq_of_max_power',fltarr(sz),/ADD_REPLACE
		str_element,big_struct,'pseudo_deltaf_to_f',fltarr(sz),/ADD_REPLACE
	endif
endif
;----

if keyword_set(array) then begin
	big_struct = {timeseries:fltarr(sz),$
				max_power:fltarr(sz),$
				bandwidth:fltarr(sz),$
				freq_of_max_power:fltarr(sz),$
				deltaf_to_f:fltarr(sz) $
	}
endif
;----

big_struct=replicate(big_struct,nfiles)



;------------------------------------------------------------------------------------------
;read data for each file
;------------------------------------------------------------------------------------------



for q=0,nfiles-1 do begin

		if not keyword_set(array) then begin
	
			if good[q] eq 1 then begin
					openr,lun,filename[q],/get_lun;,error=err
					print,'READ_TDS.PRO : FILE ' + filename[q] + ' FOUND'
				
					junk=''
					for i=0,offset-1 do begin
						readf,lun,junk
					endfor
		
					str=''
					i=0
					while not eof(lun) do begin
						readf,lun,str
						time[i] = double(strmid(str,0,17))
						Ex_raw[i] = float(strmid(str,20,17))
						Ey_raw[i] = float(strmid(str,40,17))
						Ez_raw[i] = float(strmid(str,60,17))
						Ex_cor[i] = float(strmid(str,80,17))
						Ey_cor[i] = float(strmid(str,100,17))
						Ez_cor[i] = float(strmid(str,120,17))
						if swtds eq 1. then pseudo_raw[i] = float(strmid(str,160,17))
						i=i+1
					endwhile
				;endelse
		
			endif else print,'READ_TDS.PRO : FILE ' + filename[q] + ' NOT FOUND'
	
		endif else timeseries = array

	
	
	;--------------------------------------------------------------------------------------------
	;CHECK TO MAKE SURE CURRENTLY LOADED IN FILE HAS SAME SIZE AS "SZ". IF NOT, THEN ARTIFICIALLY
	;ENLARGE IT. I HAVE TO MAKE SURE THAT EVERY ARRAY I TACK ONTO THE STRUCTURE IS THE SAME SIZE
	;--------------------------------------------------------------------------------------------


	szdiff = sz - n_elements(time)
	if szdiff ne 0 then begin
		arr2 = replicate(!values.f_nan,szdiff)
		
		time = [time,arr2]
		Ex_raw = [Ex_raw,arr2]
		Ey_raw = [Ey_raw,arr2]
		Ez_raw = [Ez_raw,arr2]
		Ex_cor = [Ex_cor,arr2]
		Ey_cor = [Ey_cor,arr2]
		Ez_cor = [Ez_cor,arr2]
		if swtds eq 1. then pseudo_raw = [pseudo_raw,arr2]
	
	endif

	
	;for q=0,nfiles-1 do begin
	;
	;if not keyword_set(array) then begin
	;
	;	if good[q] eq 1 then begin
	;		tmp2 = query_ascii(filename[q],info2)
	;		if info2.lines-23 ne 16384 then begin
	;			print,'------------------------------------------'
	;			print,'FILE ' + filename[q] + ' INCORRECT LENGTH'
	;			print,'FILE ' + filename[q] + ' NOT INCLUDED'
	;			print,'------------------------------------------'
	;		endif else begin
	;			openr,lun,filename[q],/get_lun;,error=err
	;	;		if (err ne 0) then begin
	;	;			print,'--------------------------------'
	;	;			print,'FILE ' + filename[q] + ' NOT FOUND'
	;	;			print,'--------------------------------'
	;	;			return,-1
	;	;		endif else print,'FILE ' + filename[q] + ' FOUND'
	;			print,'READ_TDS.PRO : FILE ' + filename[q] + ' FOUND'
	;		
	;			junk=''
	;			for i=0,offset-1 do begin
	;				readf,lun,junk
	;			endfor
	;
	;			str=''
	;			i=0
	;			while not eof(lun) do begin
	;				readf,lun,str
	;				time[i] = float(strmid(str,0,17))
	;				Ex_raw[i] = float(strmid(str,20,17))
	;				Ey_raw[i] = float(strmid(str,40,17))
	;				Ez_raw[i] = float(strmid(str,60,17))
	;				Ex_cor[i] = float(strmid(str,80,17))
	;				Ey_cor[i] = float(strmid(str,100,17))
	;				Ez_cor[i] = float(strmid(str,120,17))
	;				if swtds eq 1. then pseudo_raw[i] = float(strmid(str,160,17))
	;				i=i+1
	;			endwhile
	;		endelse
	;
	;	endif else print,'READ_TDS.PRO : FILE ' + filename[q] + ' NOT FOUND'
	;
	;endif else timeseries = array
	
	
	;-------------------------
	;PLACE DATA INTO STRUCTURE
	;-------------------------
	
	if not keyword_set(array) then begin
	
	big_struct[q].time = time
	big_struct[q].ex_raw = ex_raw
	big_struct[q].ey_raw = ey_raw
	big_struct[q].ez_raw = ez_raw
	big_struct[q].ex_cor = ex_cor
	big_struct[q].ey_cor = ey_cor
	big_struct[q].ez_cor = ez_cor
	if keyword_set(filename) then big_struct[q].filename = filename[q] else big_struct[q].filename = ''
	if swtds eq 1. then big_struct[q].pseudo_raw = pseudo_raw
	
	endif else big_struct[q].timeseries = timeseries
	
	
	;----------------------------------------------
	;CALCULATE VARIOUS PROPERTIES FROM SLIDING SPEC
	;----------------------------------------------
	
	if keyword_set(slidespec_vals) then begin
	e = 2.71828
	
	step = 0.1   ;0.1-->good freq res, 0.001-->good time res
	step_overlap = 0.8 ;--> decreasing this increases the timestep
	powerx = slide_spec(time,ex_cor,step,step_overlap,iwindow=1,time_index=time_index,freq_bins=freqs) 
	powery = slide_spec(time,ey_cor,step,step_overlap,iwindow=1)
	powerz = slide_spec(time,ez_cor,step,step_overlap,iwindow=1)
	
	if swtds eq 1. then powerp = slide_spec(time,pseudo_raw,step,step_overlap,iwindow=1)
	
	freqs = 1000.*freqs
	
	tms = (max(time,/nan)-min(time,/nan))*indgen(n_elements(time_index))/(n_elements(time_index)-1.)
	
	for i=0,n_elements(time_index)-1 do begin
		big_struct[q].ex_max_power[i] = max(powerx[i,*],subscrx,/nan)
		big_struct[q].ey_max_power[i] = max(powery[i,*],subscry,/nan)
		big_struct[q].ez_max_power[i] = max(powerz[i,*],subscrz,/nan)
		if swtds eq 1. then	big_struct[q].pseudo_max_power[i] = max(powerp[i,*],subscrp,/nan)
		
		big_struct[q].ex_freq_of_max_power[i] = freqs[subscrx]
		big_struct[q].ey_freq_of_max_power[i] = freqs[subscry]
		big_struct[q].ez_freq_of_max_power[i] = freqs[subscrz]
		if swtds eq 1. then big_struct[q].pseudo_freq_of_max_power[i] = freqs[subscrp]/1000.
		
		minpow = max(powerx[i,*])/e^2
		tmp = where(powerx[i,*] ge minpow,cnt)
		if tmp[0] ne -1 then big_struct[q].ex_bandwidth[i] = abs(freqs[tmp[0]] - freqs[tmp[cnt-1]]) else big_struct[q].ex_bandwidth[i] = !values.f_nan
		if tmp[0] ne -1 then big_struct[q].ex_deltaf_to_f[i] = big_struct[q].ex_bandwidth[i]/big_struct[q].ex_freq_of_max_power[i]
	
		minpow = max(powery[i,*])/e^2
		tmp = where(powery[i,*] ge minpow,cnt)
		if tmp[0] ne -1 then big_struct[q].ey_bandwidth[i] = abs(freqs[tmp[0]] - freqs[tmp[cnt-1]]) else big_struct[q].ey_bandwidth[i] = !values.f_nan
		if tmp[0] ne -1 then big_struct[q].ey_deltaf_to_f[i] = big_struct[q].ey_bandwidth[i]/big_struct[q].ey_freq_of_max_power[i]
	
		minpow = max(powerz[i,*])/e^2
		tmp = where(powerz[i,*] ge minpow,cnt)
		if tmp[0] ne -1 then big_struct[q].ez_bandwidth[i] = abs(freqs[tmp[0]] - freqs[tmp[cnt-1]]) else big_struct[q].ez_bandwidth[i] = !values.f_nan
		if tmp[0] ne -1 then big_struct[q].ez_deltaf_to_f[i] = big_struct[q].ez_bandwidth[i]/big_struct[q].ez_freq_of_max_power[i]
	
		if swtds eq 1. then begin
			minpow = max(powerp[i,*])/e^2
			tmp = where(powerp[i,*] ge minpow,cnt)
			if tmp[0] ne -1 then big_struct[q].pseudo_bandwidth[i] = abs(freqs[tmp[0]] - freqs[tmp[cnt-1]]) else big_struct[q].pseudo_bandwidth[i] = !values.f_nan
			if tmp[0] ne -1 then big_struct[q].pseudo_deltaf_to_f[i] = big_struct[q].pseudo_bandwidth[i]/big_struct[q].pseudo_freq_of_max_power[i]
		endif
	
	endfor ;each time
	
	for i=0,n_elements(freqs)-1 do begin
	
		minpow = max(powerx[*,i])/e^2
		tmp = where(powerx[*,i] ge minpow,cnt)
		if tmp[0] ne -1 then big_struct[q].ex_timewidth[i] = abs(tms[tmp[0]] - tms[tmp[cnt-1]]) else big_struct[q].ex_timewidth[i] = !values.f_nan
	
		minpow = max(powery[*,i])/e^2
		tmp = where(powery[*,i] ge minpow,cnt)
		if tmp[0] ne -1 then big_struct[q].ey_timewidth[i] = abs(tms[tmp[0]] - tms[tmp[cnt-1]]) else big_struct[q].ey_timewidth[i] = !values.f_nan
	
		minpow = max(powerz[*,i])/e^2
		tmp = where(powerz[*,i] ge minpow,cnt)
		if tmp[0] ne -1 then big_struct[q].ez_timewidth[i] = abs(tms[tmp[0]] - tms[tmp[cnt-1]]) else big_struct[q].ez_timewidth[i] = !values.f_nan
	
	endfor
	
	endif
	
	
	;----------------------------------------------
	;CALCULATE VARIOUS PROPERTIES FROM WAVELET SPEC
	;----------------------------------------------
	
	if keyword_set(wavelet_vals) then begin
	e = 2.71828
	
	if keyword_set(time) then begin
		ntime = n_elements(time)
	;	dt = (max(time,/nan)-min(time,/nan))/ntime
	endif else begin
		ntime = n_elements(timeseries)
	;	if not keyword_set(duration) then begin
	;		duration = max(time,/nan) - min(time,/nan)
	;	endif
	;	dt = duration/ntime
	endelse
	
	if not keyword_set(array) then begin
	
	;wavex = WAVELET(ex_cor,dt,PERIOD=period,COI=coi,/PAD,SIGNIF=signifx,fft_theor=spec_slopex,param=6)
	;wavey = WAVELET(ey_cor,dt,PERIOD=period,COI=coi,/PAD,SIGNIF=signify,fft_theor=spec_slopey,param=6)
	;wavez = WAVELET(ez_cor,dt,PERIOD=period,COI=coi,/PAD,SIGNIF=signifz,fft_theor=spec_slopez,param=6)
	
	wavex = wv_cwt(ex_cor,'Morlet',6,scale=period,/double,dscale=0.1)
	wavey = wv_cwt(ey_cor,'Morlet',6,scale=period,/double,dscale=0.1)
	wavez = wv_cwt(ez_cor,'Morlet',6,scale=period,/double,dscale=0.1)
	
	;if swtds eq 1. then wavep = WAVELET(pseudo_raw,dt,PERIOD=period,COI=coi,/PAD,SIGNIF=signifp,fft_theor=spec_slopep,param=6)
	if swtds eq 1. then wavep = wv_cwt(pseudo_raw,'Morlet',6,scale=period,/double,dscale=0.1)
	
	
	
	
	nscale = N_ELEMENTS(period)
	freqs = 1/period
	powerx = abs(wavex)^2
	powery = abs(wavey)^2
	powerz = abs(wavez)^2
	if swtds eq 1. then powerp = abs(wavep)^2
	
	subscr=0.
	
	for i=0,n_elements(time)-1 do begin
		big_struct[q].ex_max_power[i] = max(powerx[i,*],subscrx,/nan)
		big_struct[q].ey_max_power[i] = max(powery[i,*],subscry,/nan)
		big_struct[q].ez_max_power[i] = max(powerz[i,*],subscrz,/nan)
		if swtds eq 1. then	big_struct[q].pseudo_max_power[i] = max(powerp[i,*],subscrp,/nan)
		
		big_struct[q].ex_freq_of_max_power[i] = freqs[subscrx]
		big_struct[q].ey_freq_of_max_power[i] = freqs[subscry]
		big_struct[q].ez_freq_of_max_power[i] = freqs[subscrz]
		if swtds eq 1. then big_struct[q].pseudo_freq_of_max_power[i] = freqs[subscrp]
		
		minpow = max(powerx[i,*])/e^2
		tmp = where(powerx[i,*] ge minpow,cnt)
		if tmp[0] ne -1 then big_struct[q].ex_bandwidth[i] = abs(freqs[tmp[0]] - freqs[tmp[cnt-1]]) else big_struct[q].ex_bandwidth[i] = -1.
		if tmp[0] ne -1 then big_struct[q].ex_deltaf_to_f[i] = big_struct[q].ex_bandwidth[i]/big_struct[q].ex_freq_of_max_power[i]
	
		minpow = max(powery[i,*])/e^2
		tmp = where(powery[i,*] ge minpow,cnt)
		if tmp[0] ne -1 then big_struct[q].ey_bandwidth[i] = abs(freqs[tmp[0]] - freqs[tmp[cnt-1]]) else big_struct[q].ey_bandwidth[i] = -1.
		if tmp[0] ne -1 then big_struct[q].ey_deltaf_to_f[i] = big_struct[q].ey_bandwidth[i]/big_struct[q].ey_freq_of_max_power[i]
	
		minpow = max(powerz[i,*])/e^2
		tmp = where(powerz[i,*] ge minpow,cnt)
		if tmp[0] ne -1 then big_struct[q].ez_bandwidth[i] = abs(freqs[tmp[0]] - freqs[tmp[cnt-1]]) else big_struct[q].ez_bandwidth[i] = -1.
		if tmp[0] ne -1 then big_struct[q].ez_deltaf_to_f[i] = big_struct[q].ez_bandwidth[i]/big_struct[q].ez_freq_of_max_power[i]
	
		if swtds eq 1. then begin
			minpow = max(powerp[i,*])/e^2
			tmp = where(powerp[i,*] ge minpow,cnt)
			if tmp[0] ne -1 then big_struct[q].pseudo_bandwidth[i] = abs(freqs[tmp[0]] - freqs[tmp[cnt-1]]) else big_struct[q].pseudo_bandwidth[i] = -1.
			if tmp[0] ne -1 then big_struct[q].pseudo_deltaf_to_f[i] = big_struct[q].pseudo_bandwidth[i]/big_struct[q].pseudo_freq_of_max_power[i]
		endif
	
	endfor ;each time
	
	for i=0,n_elements(freqs)-1 do begin
	
		minpow = max(powerx[*,i])/e^2
		tmp = where(powerx[*,i] ge minpow,cnt)
		if tmp[0] ne -1 then big_struct[q].ex_timewidth[i] = abs(time[tmp[0]] - time[tmp[cnt-1]]) else big_struct[q].ex_timewidth[i] = !values.f_nan
	
		minpow = max(powery[*,i])/e^2
		tmp = where(powery[*,i] ge minpow,cnt)
		if tmp[0] ne -1 then big_struct[q].ey_timewidth[i] = abs(time[tmp[0]] - time[tmp[cnt-1]]) else big_struct[q].ey_timewidth[i] = !values.f_nan
	
		minpow = max(powerz[*,i])/e^2
		tmp = where(powerz[*,i] ge minpow,cnt)
		if tmp[0] ne -1 then big_struct[q].ez_timewidth[i] = abs(time[tmp[0]] - time[tmp[cnt-1]]) else big_struct[q].ez_timewidth[i] = !values.f_nan
	
	endfor
	
	
	
	endif  ;not keyword_set(array)
	
	;----
	if keyword_set(array) then begin
	
	;	wave = WAVELET(timeseries,dt,PERIOD=period,COI=coi,/PAD,SIGNIF=signif,fft_theor=spec_slope,param=6)
		wave = wv_cwt(timeseries,'Morlet',6,scale=period,/double,dscale=0.1)
	
		nscale = N_ELEMENTS(period)
		freqs = 1/period
		power = abs(wave)^2
		subscr=0.
	
		for i=0,n_elements(timeseries)-1 do begin
			big_struct[q].max_power[i] = max(power[i,*],subscr,/nan)
			big_struct[q].freq_of_max_power[i] = freqs[subscr]
			minpow = max(power[i,*])/e^2
			tmp = where(power[i,*] ge minpow,cnt)
			if tmp[0] ne -1 then big_struct[q].bandwidth[i] = abs(freqs[tmp[0]] - freqs[tmp[cnt-1]]) else big_struct[q].bandwidth[i] = -1.
			if tmp[0] ne -1 then big_struct[q].deltaf_to_f[i] = big_struct[q].bandwidth[i]/big_struct[q].freq_of_max_power[i]
		endfor
	endif ;keyword_set(array)
	;----
	
	endif ;wavelet_vals ne 0.
	
	
	
	if not keyword_set(array) then begin
		close,lun
		free_lun,lun
	endif

endfor ;for each file


;----------------------
;append the unix times
;----------------------

if not keyword_set(array) then begin

	t0 = strmid(date_time,20,3)

	for i=0,nfiles-1 do begin
		msecs = 1000d*double(big_struct[i].time) + double(t0[i])
		foo = strtrim(msecs/1000.,2)
		foosecs=fix(strmid(foo,0,1))
		foomsecs = strmid(foo,2,100)
		secs=fix(strmid(date_time[i],17,2))
		big_struct[i].date_time = strmid(date_time[i],0,17) + strmid(string(secs+foosecs),6,2) + '.' + foomsecs
	endfor

endif

return,big_struct
end

