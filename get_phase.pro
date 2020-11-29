; simple routine for calculating phase curve using ctime

pro get_phase, tvar=tvar, phcurve=phcurve, deg=deg, width=width, win=win


	ctime,time,npoints=2,vname=vname
	
	; make sure we have a usable time range and phase variable

	if time[0] eq time[1] then begin

		message,'Invalid time range.  Returning...'
		phcurve=-1
		return
		
	endif

	if ~keyword_set(tvar) then begin
		if (vname[0] ne vname[1]) $
			or (strpos(vname[0],'pha') eq -1) $
			or (strpos(vname[1],'pha') eq -1) then begin
			
			message,'No valid phase variable found for selected times (time points must be selected from within the phase panel).  Returning...'
			phcurve=-1
			return
	
		endif else tvar=vname[0]

	endif
	
	get_data,tvar,data=d

	t=d.x
	
	if time[0] gt time[1] then time=reverse(time)
	
	tstart=time[0]
	tend=time[1]
		
	istart=min(where(t gt tstart))
	iend=max(where(t lt tend))
	
	ph=d.y
	phy=d.v
	physize=size(phy,/n_elements)
	
	phcurve=fltarr(physize)
	
	for i=0,physize-1 do phcurve[i]=mean(ph[istart:iend,i],/nan)
	
	if ~keyword_set(width) then width=20
	sphcurve=smooth(phcurve,width,/nan)
	
	if ~keyword_set(deg) then deg=3
	pfit=poly_fit(phy,phcurve,deg,chisq=chisq,covar=covar,measure_errors=errors,yfit=yfit)
	
	if ~keyword_set(win) then win=3
	window,win
	plot,phy,phcurve,psym=-4,/xstyle,yrange=[-180,180],/ystyle, $
		xtitle='f [Hz]', ytitle='phase [deg]'
	oplot,phy,sphcurve
	oplot,phy,yfit
	
end