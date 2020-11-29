; simple routine for calculating a power spectrum using ctime

pro get_fft, tvar=tvar, power=power, windowed=windowed, win=win, exact=exact,$
	flog=flog


	ctime,time,npoints=2,vname=vname,exact=exact
	
	; make sure we have a usable time range and phase variable

	if time[0] eq time[1] then begin

		message,'Invalid time range.  Returning...'
		phcurve=-1
		return
		
	endif

	if ~keyword_set(tvar) then tvar=vname[0]
	
	get_data,tvar,data=d

	t=d.x
	
	if time[0] gt time[1] then time=reverse(time)
	
	tstart=time[0]
	tend=time[1]
	evlength=tend-tstart
		
	istart=min(where(t gt tstart))
	iend=max(where(t lt tend))
	npoints=iend-istart+1
	
	times=d.x[istart:iend]
	timeseries=d.y[istart:iend]
	dts=times[1:npoints-1]-times[0:npoints-2]
	deltat=median(dts)
	
	print,''
	print,'median dt(s):', deltat
	print,'average dt(s):', mean(dts)
	print,'sample rate (kHz):',1./(1000.*deltat),1./(1000.*mean(dts))
	
	freq_bins=findgen(npoints/2)/(deltat*1000.*npoints)
	
	if keyword_set(windowed) then $
		fft_win=(25./46.)-(21./46.)*cos(2.*!pi*lindgen(npoints)/npoints) $ ;Hammimg
	else fft_win=make_array(npoints,value=1.0,/float)
	
	power=2*evlength*(ABS(FFT(timeseries*fft_win)))^2
	power=10*alog10(power[0:npoints/2-1])

	if ~keyword_set(win) then win=2
	window,win
	
	if keyword_set(flog) then frange=[.01,freq_bins[npoints/2-1]] $
		else frange=[freq_bins[0],freq_bins[npoints/2-1]]
	plot,freq_bins,power,xtitle='f [kHz]',ytitle='Power [dB]',$
		/xstyle,/ystyle,xrange=frange, xlog=flog
	
end