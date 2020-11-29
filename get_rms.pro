; simple routine for calculating RMS of a timeseries

pro get_rms, tvar=tvar, rms=rms, exact=exact


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
		
	istart=min(where(t gt tstart))
	iend=max(where(t lt tend))
	npoints=iend-istart+1
	
	v=double(d.y[istart:iend])
	v2=v*v
	rms=sqrt(total(v2)/npoints)
	
	print,''
	print,'RMS:', rms
	print,''
		
end