;+
;FUNCTION:   slide_spec.pro
;  Calculates a sliding power spectrum for time series data 
;
;							
;  version 0.1, 20080827
;
;ARGUMENTS:
;		TIME		-> time array (sec)
;		DATA		-> time series data array
;		STEP		-> length of each slice, in percent of total length
;		STEP_OVERLAP	-> slice overlap, given in percent of slice length
;
;KEYWORDS:
;       TIME_INDEX=time_index	<- integer array containing the index for each slice time.
;								   i.e., time[time_index] is the array of x-axis times for the spec
;		FREQ_BINS=freq_bins		<- array of frequency bins (kHz)
;		/DB						<> if set, calculates power in dB: power_temp=10*alog10(2*padded_step_length*(abs(fft(temp_time_series)))^2)
;								if not set then: power_temp=2*padded_step_length*(abs(fft(temp_time_series)))^2
;		mVm_arr 			    <- identical to fft_array but in mV/m. Summing this array
;									over all the freq bins for a single time will give you roughly the waveform
;									amplitude at this time
;		iwindow					-> window the data (1:Hamming, 2:Hanning, 3:Gaussian)
;		zero_pad				-> set to have program zero pad the data. Can run MUCH faster if you set this. 
;
;RETURNS:
;		FFT_ARRAY	<- 2D sliding FFT array (# time steps by # freq_bins)
;
;CALLING SEQUENCE:
;       IDL> 
;
;NOTES:
;       (none)
;-
;CREATED BY:    Kris Kersten, August 2008
;
;HISTORY:
;   08/30/2008  - KK - created, AWB added mVm_arr option
;   01/08/2011  - AWB - fixed the end_pad_index variable...at least I think.
;						Changed all the "lindgen" to "dindgen". This was occasionally an issue
;						when plotting with dB
;						
;
;INCLUDED MODULES:
;       slide_spec
;
;LIBS USED:
;       (none)
;
;DEPENDENCIES
;       
;-


function slide_spec,time,time_series,step,step_overlap,time_index=time_index,freq_bins=freq_bins,db=db,iwindow=iwindow,zero_pad=zero_pad,mVm_arr=mVm_arr

	npoints=0L		; number of points in the time series
	step_points=0L	; number of points in each "slice"
	step_overlap_points=0L
	step_start=0L	; starting index of sliding "slice"
	step_end=0L		; ending index of sliding "slice"
	step_count=0L	; counter for stepping through the time series

	; figure out number of points per step
	if step gt 1. then step=1.
	if step_overlap gt 1. then step_overlap=1.
	if step le 0. then step=.25
	npoints=size(time,/n_elements)
	if npoints ne size(time_series,/n_elements) then message,'size mismatch in time series.'
	step_points=long(npoints*step)
	step_overlap_points=long(step_points*step_overlap)
	; make sure we have a working overlap (i.e., in the range 0,step_points-1)
	if step_overlap_points lt 0 then step_overlap_points=0
	if step_overlap_points eq step_points then step_overlap_points=step_points-1
	
	; figure out time step (assuming input times in seconds)
	total_time=time[npoints-1]-time[0]
	time_step=total_time/npoints
	step_length=step_points*time_step
	
	step_start=0
	step_end=step_points-1

	; zero pad each slice for the FFT if desired
	padded_step_points=0L
	if ~keyword_set(zero_pad) then begin
		padded_step_points=step_points
		padded_step_length=step_length	
	endif else begin
		while padded_step_points le step_points do begin
			padded_step_points=2L^zero_pad
			zero_pad+=1
		endwhile
		padded_step_length=padded_step_points*time_step
		front_zero_pad=make_array((padded_step_points-step_points)/2,value=0.,/float)
		back_zero_pad=make_array(padded_step_points-step_points-(padded_step_points-step_points)/2,value=0.,/float) ; this looks silly, but is necessary when the total number of pad points is odd
	endelse

	; create the y-axis freq bins [kHz] (dependent on step, again assuming times in seconds)
	freq_bins=(dindgen(padded_step_points/2-1)+1)*npoints/(padded_step_points*total_time*1000.)

	; set up windowing, if desired
	if ~keyword_set(iwindow) then iwindow=0
	case iwindow of
		0: window_array=make_array(step_points,value=1.0,/float)
		1: window_array=(25./46.)-(21./46.)*cos(2.*!pi*dindgen(step_points)/step_points) ;Hammimg
		2: window_array=0.5-0.5*cos(2.*!pi*dindgen(step_points)/step_points) ;Hanning
		3: window_array=exp(-0.5*(3.*(dindgen(step_points)-step_points/2)/(step_points/2.))^2) ;Gaussian
		else:
	endcase
	sumsq=total(window_array^2)

	rms=sqrt(sumsq/step_points)
	window_array=window_array/rms
	
	; calculate the sliding fft
	fft_array=fltarr((npoints-step_overlap_points)/(step_points-step_overlap_points),padded_step_points/2-1)

	while step_end le npoints-1 do begin
		if keyword_set(zero_pad) then temp_time_series=[front_zero_pad,window_array*time_series[step_start:step_end],back_zero_pad] $
			else temp_time_series=window_array*time_series[step_start:step_end]
	    if keyword_set(db) then $
			power_temp=10*alog10(2*padded_step_length*(abs(fft(temp_time_series)))^2) $
		else $
			power_temp=2*padded_step_length*(abs(fft(temp_time_series)))^2
			
		fft_array[step_count,*]=power_temp[1:padded_step_points/2-1]
		step_count+=1
		step_start+=(step_points-step_overlap_points)
		step_end+=(step_points-step_overlap_points)
	endwhile

	; set up the time index with padding
	time_index=dindgen(step_count)*(step_points-step_overlap_points)+step_points/2
	start_pad_index=(dindgen(round(float(time_index[0])/(step_points-step_overlap_points))))*(step_points-step_overlap_points)


;###########################################
;This is a change made by AWB on 01/08/2011
;I haven't examined why prog would fail here on occasion, but the fix seems to work.
;	end_pad_index=(dindgen((npoints-time_index[step_count-1])/(step_points-step_overlap_points)))*(step_points-step_overlap_points)+time_index[step_count-1]


if (npoints-time_index[step_count-1])/(step_points-step_overlap_points) ge 1 then $
	end_pad_index=(dindgen((npoints-time_index[step_count-1])/(step_points-step_overlap_points)))*(step_points-step_overlap_points)+time_index[step_count-1] $
else end_pad_index = time_index[step_count-1]
;##########################################




	time_index=[start_pad_index,time_index[0],time_index,end_pad_index,npoints-1]

	; pad the fft array
	start_pad_data=make_array(size(start_pad_index,/n_elements)+1,padded_step_points/2-1,value=min(fft_array))
	end_pad_data=make_array(size(end_pad_index,/n_elements)+1,padded_step_points/2-1,value=min(fft_array))
	fft_array=[start_pad_data,fft_array,end_pad_data]


	if keyword_set(db) then mVm_arr = sqrt(10^(fft_array/10.)/(2.*padded_step_length))  $
	else mVm_arr = sqrt(fft_array/2./padded_step_length)


	return,fft_array 
end


	
	
