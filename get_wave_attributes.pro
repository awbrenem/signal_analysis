;+
;*****************************************************************************************
;
;  FUNCTION :   get_wave_attributes.pro
;  PURPOSE  :   Calculate various attibutes for the input waveform such as:
;					-freq of peak power
;					-bandwidth, df/f
;					-bandwidth standard deviation (only if wssn ne 'neither')
;					-etc...
;
;  CALLED BY:
;
;
;  CALLS:
;       wv_cwt.pro
;		slide_spec.pro
;
;  INPUT:
;		wave = [n] array of wave amplitude values
;		time = [n] array of the time values
;		step = the step size for the sliding spec or wavelet transform
;		step_overlap = the amount of overlap between steps for sliding spec transform
;		wssn -> 'wavelet', 'slidespec' or 'neither'
;				neither -> only return single waveform values. The FFT performed here included
;							the entire waveform to maximize freq resolution
;				wavelet -> in addition to single waveform values, use the wavelet transform
;							to return various values for each time bin
;				slidespec -> same as wavelet but uses slide_spec.pro
;		mVm -> set to return power in mV/m rather than dB.
;
;
;  EXAMPLES:    x = get_wave_attributes(time,wave,wssn='wavelet',step=0.05,step_overlap=.9)
;               x = get_wave_attributes(time,wave,step=0.075,step_overlap=.85)
;
;
;   CHANGED:
;	*Bandwidth is calculated as twice the distance from the location of max power to the
;		furthest point (farthest right in the waveform plot) in the waveform above the
;		cutoff power ('minpow').
;		This helps the program avoid calculating high bandwidths at the max power
;		due to some points above the power cutoff at very low frequencies.
;
;
;	NOTES: This program ignores freqs lt 10 Hz when figuring out the freq of peak power.
;			These aren't reliable in the STEREO burst captures.
;
;		Sometimes the bandwidth is zero b/c the power drops off to below the 1/e^2
;		value within the frequency resolution. In this case, I'll just define the bandwidth
;		as the frequency resolution.
;
;
;
;   CREATED:  08/02/2011
;   CREATED BY:  Sam Schreiner
;    LAST MODIFIED:
;    MODIFIED BY: Aaron W Breneman 09/08/2011
;				  Aaron W Breneman 12/02/2011 - fixed Bw = 0 issue (see notes)
;
;*****************************************************************************************
;-


;wa = get_wave_attributes2(time,wave,wssn='neither')


;--------------------------------------------------------------------------------------
;--------------------------------------------------------------------------------------

function BW_cntr_right,freqs,power,subscr,plotpow=plotpow,mVm=mVm
;Takes a plot of frequency vs. power (i.e. from an FFT) and calculates
;the bandwidth (BW) as double the distance from the frequency at
;maximum power to the highest frequency above the power cutoff.
;Note that I'm not counting any power at freqs lower than the peak freq
;because there is often a peak at low freqs that causes the bandwidth to be huge.

	e = 2.71828d0

	minpow = power[subscr]/e^2		;power is on a linear scale

	power_tmp = power[subscr:n_elements(power)-1]  ;only data at and to the right of the peak
	tmp = where(power_tmp ge minpow) + subscr
	if tmp[0] ne -1 then begin
		BW = (max(freqs[tmp]) - freqs[subscr])*2

		;Sometimes the bandwidth is zero b/c the power drops off to below the 1/e^2
		;value within the frequency resolution. In this case, I'll just define the bandwidth
		;as the frequency resolution.
		if BW eq 0 then BW = freqs[1] - freqs[0]
	endif else BW = !values.f_nan

	;option to plot the power vs. frequency and which data points are above the cutoff
	if keyword_Set(plotpow) then begin
		  plot,freqs,power-power[0],xrange=[0,5*freqs[subscr]]
		  oplot,freqs[tmp],power[tmp]-power[0],psym=5,color=250
		  wait,.5
	endif

	return,BW
end

;--------------------------------------------------------------------------------------
;--------------------------------------------------------------------------------------



function get_wave_attributes,time,wave,wssn=wssn,$
	  pets=step,step_overlap=step_overlap

if ~keyword_set(step) then step = 0.1
if ~keyword_set(step_overlap) then step_overlap = 0.8
if ~keyword_set(wssn) then wssn = 'slidespec'

e = 2.71828


notes = ["time_bins -> the time bins returned from slidespec or wavelet (secs)",$
		 "freq_bins -> the freq bins returned from slidespec or wavelet (Hz)",$
		 "max_power -> array of the maximum power for each time column",$
		 "bandwidth -> array of bandwidth (Hz) as a function of time",$
		 "timewidth -> array of the time of max power for each frequency row",$
		 "freq_of_max_power -> array of the frequency of the max power for each time",$
		 "deltaf_to_f -> the bandwidth divided by the frequency for each time"]

notes2 = ["max_power -> max power for the entire burst capture",$
		 "bandwidth -> bandwidth for the time of maximum power in burst capture",$
		 "freq_of_max_power -> the freq of the max power for the entire TDS capture",$
		 "deltaf_to_f -> the bandwidth divided by the frequency for the time of maximum power"]

if ~keyword_set(wssn) then wssn = 'slidespec'

sz = n_elements(wave)
dt = (max(time,/nan)-min(time,/nan))/n_elements(time)

;-------------------------------------------------------------------------------------------------
;CALCULATE QUANTITIES FROM SLIDING SPECTROGRAM
;-------------------------------------------------------------------------------------------------

if wssn eq 'slidespec' then begin

	;SET UP STRUCTURE
	power = slide_spec(time,wave,step,step_overlap,iwindow=1,time_index=time_index,freq_bins=freqs)

	freqs = 1000.*freqs
	nfreqs = n_elements(freqs)
	ntimes = n_elements(time_index)

	tms = (max(time,/nan)-min(time,/nan))*indgen(n_elements(time_index))/(n_elements(time_index)-1.)

	x = {time:fltarr(sz),$
		wave:fltarr(sz),$
		time_bins:fltarr(ntimes),$
		freq_bins:fltarr(nfreqs),$
		max_power:fltarr(ntimes),$
		bandwidth:fltarr(ntimes),$
		timewidth:fltarr(nfreqs),$
		freq_of_max_power:fltarr(ntimes),$
		deltaf_to_f:fltarr(ntimes),$
		notes:notes $
	}

	;CALCULATE STUFF
	x.time_bins = tms    ;s
	x.freq_bins = freqs  ;Hz

	for i=0,n_elements(time_index)-1 do begin
		x.max_power[i] = max(power[i,*],subscr,/nan)
		x.freq_of_max_power[i] = freqs[subscr]
		x.bandwidth[i] = BW_cntr_right(freqs,power[i,*],subscr);,/plotpow)
		x.deltaf_to_f[i] = x.bandwidth[i]/x.freq_of_max_power[i]
	endfor

	for i=0,n_elements(freqs)-1 do begin
		minpow = max(power[*,i])/e^2
		tmp = where(power[*,i] ge minpow,cnt)
		if tmp[0] ne -1 then x.timewidth[i] = abs(tms[tmp[0]] - tms[tmp[cnt-1]]) else x.timewidth[i] = !values.f_nan
	endfor

;	str_element,x,'tms',tms,/ADD_REPLACE	;special time elements corresponding to slide_spec quantities
	x.time = time
	x.wave = wave


endif



;-------------------------------------------------------------------------------------------------
;CALCULATE QUANTITIES FROM WAVELET SPECTROGRAM
;-------------------------------------------------------------------------------------------------

if wssn eq 'wavelet' then begin

	;SET UP STRUCTURE
 	waveWV = wv_cwt(wave,'Morlet',6,scale=period,/double,dscale=step,/pad)

	nfreqs = n_elements(period)
	period2 = period*dt
	freqs = 1/period2

	power = abs(waveWV)^2

	x = {		time:fltarr(sz),$
				wave:fltarr(sz),$
				time_bins:fltarr(n_elements(period2)),$
				freq_bins:fltarr(nfreqs),$
				max_power:fltarr(sz),$
				bandwidth:fltarr(sz),$
				timewidth:fltarr(nfreqs),$
				freq_of_max_power:fltarr(sz),$
				deltaf_to_f:fltarr(sz),$
				f_resolution:!values.f_nan,$
				notes:notes $
	}

	;CALCULATE VARIOUS QUANTITIES
	x.time_bins = period2
	x.freq_bins = reverse(freqs)
	x.f_resolution = freqs[1] - freqs[0]
	for i=0,n_elements(time)-1 do begin
		x.max_power[i] = max(power[i,*],subscr,/nan)
		x.freq_of_max_power[i] = freqs[subscr]
		x.bandwidth[i] = BW_cntr_right(freqs,power[i,*],subscr)
		x.deltaf_to_f[i] = x.bandwidth[i]/x.freq_of_max_power[i]
	endfor

	for i=0,n_elements(freqs)-1 do begin
		minpow = max(power[*,i])/e^2
		tmp = where(power[*,i] ge minpow,cnt)
		if tmp[0] ne -1 then x.timewidth[i] = abs(time[tmp[0]] - time[tmp[cnt-1]]) else x.timewidth[i] = !values.f_nan
	endfor

	x.time = time
	x.wave = wave


endif



;-----------------------------------------------------------------------
;STDEV OF BANDWIDTH (MIDDLE 3/5 OF SIGNAL) - HELPS TO SEPARATE WHISTLERS FROM LANGMUIR WAVES
;THE REALLY NICE NARROWBAND WHISTLERS HAVE A STDEV OF ~< 20
;Note that I don't want to use entire signal b/c sometimes a nice sine wave trails off into
;noise at the end. This takes away the favoritism for long-duration signals
;-----------------------------------------------------------------------


if wssn ne 'neither' then begin
	fraction = n_elements(time_index)/8.
	tmpp = moment(x.bandwidth[floor(fraction):floor(n_elements(x.bandwidth)-fraction)])

;	print,tms[fraction]
;	print,tms[n_elements(tms)-fraction]
;	plot,tms,x.bandwidth

	bandwidth_stdev = sqrt(tmpp[1]) ;* 1000.


endif else bandwidth_stdev = !values.f_nan


;------------------------------------------------------------------------------------
;FIND THE VARIOUS PROPERTIES OF THE WAVE OVERALL (FFT perfomed on the entire waveform
;------------------------------------------------------------------------------------

y = {time:!values.f_nan,$
	amp:!values.f_nan,$
	max_power:!values.f_nan,$
	min_power:!values.f_nan,$
	max_powerdB:!values.f_nan,$
	min_powerdB:!values.f_nan,$
	bandwidth:!values.f_nan,$
	bandwidth_stdev:bandwidth_stdev,$
	freq_of_max_power:!values.f_nan,$
	deltaf_to_f:!values.f_nan,$
	f_resolution:!values.f_nan,$
	notes2:notes2 $
}


fft_freqs = lindgen(sz/2+1)/(sz*dt)  ;Donnelly and Rust vol1 eqn 1

;various forms of power - probably doesn't matter which one is used
pow = 2*abs(fft(wave,/double))^2
;pow = pow[0:sz/2-1]
pow = pow[0:sz/2]
powdb = 10*alog(pow)
if keyword_set(mVm) then pow = sqrt(pow/2.)  ;mV/m (not sure if this is correct)


;make sure that we don't take into account frequencies less than 10 Hz
;when determining the frequency of max power. These freqs can be artificially introduced
;by the gain correction.

goo = where(fft_freqs le 10)
if goo[0] ne -1 then elem = goo[n_elements(goo)-1] else elem = 0
y.max_power = max(pow[elem:n_elements(pow)-1],subscr,/nan)
subscr = subscr + elem


y.min_power = y.max_power/e^2       ;make sure this is the linear version
y.freq_of_max_power = fft_freqs[subscr]
y.bandwidth = BW_cntr_right(fft_freqs,pow,subscr)
y.deltaf_to_f = y.bandwidth/y.freq_of_max_power
y.f_resolution = fft_freqs[1] - fft_freqs[0]


tmpp = where(fft_freqs ge (y.freq_of_max_power + y.bandwidth/2.))
if tmpp[0] eq n_elements(powdB) then tmpp[0] = tmpp[0] - 1

y.max_powerdB = powdB[subscr]
y.min_powerdB = powdB[tmpp[0]]       ;make sure this is the linear version


;########################################################################

y.amp = max(wave,ss)
y.time = time[ss]



if wssn ne 'neither' then begin
	window,0,retain=1,xsize=550,ysize=710
	!p.multi=[0,0,5,0,0]
	!p.charsize=1.8

	plot,x.time,x.wave,ytitle='Ew'
	plot,x.time_bins,x.max_power,ytitle='max power'
	plot,x.time_bins,x.bandwidth,ytitle='BW',yrange=[0,10000]
	plot,x.time_bins,x.freq_of_max_power,ytitle='F at pmax'
	plot,x.time_bins,x.deltaf_to_f,ytitle='df/f',yrange=[0,1]
endif



if wssn eq 'neither' then array = y else array = {x:x,y:y}

;!p.multi = [0,0,3]
;plot,fft_freqs,powdB,xrange=[y.freq_of_max_power/2.,y.freq_of_max_power*2]
;plot,time,wave,xrange=[0,max(tms)],xstyle=1
;plot,tms[fraction:n_elements(tms)-fraction],x.bandwidth[fraction:n_elements(tms)-fraction],xrange=[0,max(tms)],xstyle=1


;stop

return,array
end
