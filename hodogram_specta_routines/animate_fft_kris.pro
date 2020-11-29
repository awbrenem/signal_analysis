;+
;PROGRAM: animate_fft_kris.pro
;
;PURPOSE: Outputs a MPEG movie (.mov), a sliding FFT of timeseries data.
;
;ARGUMENTS:
;       TIMESERIES  ->  input timeseries data
;       TIMESTEP    ->  sample period (time per point) of the input timeseries
;       FFTLENGTH   ->  number of points of timeseries data to be used in each
;                       FFT
;       FFTSTEP     ->  number of points to shift for each consecutive FFT
;						(values less than 0.02*FFTLENGTH give smooth animation)
;
;KEYWORDS:
;       MOVIENAME   ->  moviename='my_movie'  Defines the output filename for
;                                             the animated FFT movie (input
;                                             basename only, extension is
;                                             appended automatically).
;       SCREEN      /   /screen               By default, the movie is generated
;                                             from PNG captures of Z buffer
;                                             plots.  The /screen keyword forces
;                                             plots to be drawn to the X window
;                                             display (slightly slower than Z
;                                             buffered plots).
;       FULLSERIES  /   /fullseries
;       XSIZE       ->  xsize=some_xsize
;       YSIZE       ->  ysize=some_ysize
;
;RETURNS:     N/A
;
;CALLING SEQUENCE:
;
;       IDL> 
;
;NOTES:
;
;
;-
;CREATED BY:    Kris Kersten, Nov 2009
;
;MODIFICATION HISTORY:
;       11/20/2009  created, KK
;
;INCLUDED MODULES:
;       animate_fft
;
;LIBS USED:
;
;
;DEPENDENCIES:
;
;-
pro animate_fft_kris,timeseries,timestep,fftlength,fftstep,moviename=moviename, $
				screen=screen,fullseries=fullseries,xsize=xsize,ysize=ysize, $
				fft_array=fft_array,no_interp=no_interp


	start_time=systime(1)
	
	; set up a randomly named image directory
	imgdir='shot.'+strcompress(string(randomu(systime(1),/long)),/remove_all)
	spawn,'mkdir '+imgdir
	
	; figure out what to call the movie
	if ~keyword_set(moviename) then begin
		moviename=''
		print,'Movie filename not set on input (set with moviename="example" keyword.)'
		print,'Please enter the base filename.' 
		read,moviename,prompt='(.mov extension will be added): '
	endif
	moviename=moviename+'.mov'
	
	; set up the window or Z buffer
	if keyword_set(fullseries) then begin
		!p.multi=[0,1,2]
		if ~keyword_set(xsize) then xsize=800
		if ~keyword_set(ysize) then ysize=600
	endif else begin
		!p.multi=[0,2,1]
		if ~keyword_set(xsize) then xsize=900
		if ~keyword_set(ysize) then ysize=450
	endelse
	
	; plot to screen or z buffer (default)
	if keyword_set(screen) then begin
		set_plot,'X'
		!p.charsize=0.8
		window,xsize=xsize,ysize=ysize
	endif else begin
		set_plot,'Z'
		device,set_pixel_depth=24,set_resolution=[xsize,ysize]
		!p.charsize=0.65
	endelse

	device,decomposed=0
	loadct,0,/silent
	

	npoints=long(size(timeseries,/n_elements))
	; nsteps is the number of individual FFTs
	nsteps=(npoints-fftlength)/fftstep

	; initialize the sub series vars
	subseries=fltarr(fftlength)
	subtimes=fltarr(fftlength)
	subpower=fltarr(nsteps,fftlength/2)
	
    freq_bins=lindgen(fftlength/2)/(timestep*fftlength)
	; convert to MHz
	freq_bins=freq_bins/1.e6
	fmin_index=max(where(freq_bins lt 10))
	fmin=freq_bins[fmin_index]

	times=dindgen(npoints)*timestep
	; convert to microseconds
	times=times*1.e6

	; set up the timeseries plot range
	ymax=max(abs(timeseries),/nan)
	plotyrange=[-1.1*ymax,1.1*ymax]

	; set up the start and end points for the sub timeseries
	startpoint=long(0)
	endpoint=long(startpoint+fftlength-1)

	
	length=npoints*timestep
	print,'length (s)',length
	sublength=fftlength*timestep
	print,'sublength (s)',sublength

	subseries=timeseries[startpoint:endpoint]
	subtimes=times[startpoint:endpoint]

	subpower[0,*]=((abs(fft(timeseries[startpoint:endpoint])))^2)[0:fftlength/2-1]

	startpoint+=fftstep
	endpoint+=fftstep

	for stepcount=1,nsteps-1 do begin

		subpower[stepcount,*]=((abs(fft(timeseries[startpoint:endpoint])))^2)[0:fftlength/2-1]

		; now move to the next interval		
		startpoint+=fftstep
		endpoint+=fftstep
		; spit out a little status update so we know how fast we're moving through the ffts/
		if stepcount mod 100 eq 0 then print,'stepcount:',stepcount
	endfor

	fft_array=subpower


	readmov='y' ;force movie
	while readmov ne 'y' and readmov ne 'n' do read,readmov,prompt="Make the movie? (y/n): "
	
	if readmov eq 'y' then begin

		powermax=max(subpower,/nan)
		powermin=min(subpower,/nan)

		startpoint=long(0)
		endpoint=long(startpoint+fftlength-1)
	
		; for the sake of speedy movie output, we're going to plot an interpolated
		; array when we're working with the full series
	
		if keyword_set(no_interp) then begin
			scale=1
			small_times=times
			small_timeseries=timeseries
		endif else begin
			scale=10
			small_times=interpol(times,npoints/scale)
			small_timeseries=interpol(timeseries,npoints/scale)
		endelse
		
		for stepcount=1,nsteps-1 do begin
			; get the current starting time in the overall timeseries for plot label
			stringtime=string(startpoint*timestep*1.e6,format='(f6.2)')
			stringcount=string(stepcount,format='(i4.4)')
			if keyword_set(fullseries) then begin
				plot,times,timeseries,yrange=plotyrange,/xstyle,/ystyle, $
					ytitle='channel voltage [V]', xtitle='time [microsec]', $
						title=stringtime+' microsec',/nodata

				; interpolated plot of the full timeseries, sub series highlighted (fast)
				oplot,small_times[0:startpoint/scale],small_timeseries[0:startpoint/scale],color=195
				oplot,small_times[startpoint/scale:endpoint/scale],small_timeseries[startpoint/scale:endpoint/scale],color=255
				oplot,small_times[endpoint/scale+1:npoints/scale-1],small_timeseries[endpoint/scale+1:npoints/scale-1],color=195

			endif else begin
				; single plot of the sub series
				plot,times[startpoint:endpoint],timeseries[startpoint:endpoint],$
					yrange=plotyrange,/ystyle,/xstyle,ytitle='channel voltage [V]',$
					xtitle='time [microsec]',title=stringtime+' microsec'
			endelse
		
			plot,freq_bins,subpower[stepcount,*],yrange=[powermin,powermax],$
				/xstyle,/ystyle,/ylog,/xlog,xrange=[fmin,freq_bins[fftlength/2-1]],$
				xtitle='f [MHz]',ytitle='power [arb]'
			write_png,imgdir+'/'+stringcount+'.png',tvrd(true=1)

			if stepcount eq 1 then begin
				plot,times,timeseries,yrange=plotyrange,/xstyle,/ystyle, $
					ytitle='channel voltage [V]', xtitle='time [microsec]', $
						title=stringtime+' microsec',/nodata
				plot,freq_bins,subpower[stepcount,*],yrange=[powermin,powermax],$
					/xstyle,/ystyle,/ylog,/xlog,xrange=[fmin,freq_bins[fftlength/2-1]],$
					xtitle='f [MHz]',ytitle='power [arb]',/nodata	
				write_png,'blank_plot.png',tvrd(true=1)
			endif

			startpoint+=fftstep
			endpoint+=fftstep
			if stepcount mod 100 eq 0 then print,'stepcount:',stepcount
		
		endfor

		print,'Encoding movie: '+moviename
		encode_command='ffmpeg -r 60 -i '+imgdir+'/%04d.png -qscale 2 '+moviename
		
		print,encode_command
		stop
		spawn,encode_command
		
		print,''
		print,'Movie complete: '+moviename
	endif
	
	stop
	
	
	print,'Cleaning up frame images...'
	rmcommand='rm -rf '+imgdir
	spawn,rmcommand
	print,'Done!'

	end_time=systime(1)
	
	execution_time=end_time-start_time
	print,''
	print,'EXECUTION TIME: '+strcompress(string(execution_time),/remove_all)+' sec'
	if execution_time gt 60. then print,'                      ( '+strcompress(string(long(execution_time)/60),/remove_all)+' min '+strcompress(string(long(execution_time) mod 60),/remove_all)+' sec )'
	print,''


end
