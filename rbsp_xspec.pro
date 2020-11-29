;+
; PROCEDURE rbsp_xspec
;
; PURPOSE: Calculates cross spectrum coherence and phase for RBSP waveform data
;
; CALLING SEQUENCE: rbsp_xspec,tplot_var1, tplot_var2
;
; INPUTS:
;		tplot_var1 - TPLOT variable name or number
;		tplot_var2 - TPLOT variable name or number
;
; OUTPUTS:
;
; KEYWORDS:
;		tplot_var_xspec = string name for new TPLOT XSPEC variable
;		npts = number of points per FFT
;		n_ave = number of specs to average for each spectrum
;		tspec = return variable for spectrum times
;		xspec_coh = return variable for coherence spectrum
;		xspec_pha = return variable for phase spectrum
;		freq = return variable for frequency bins
;		df = return variable for frequency steps
;		/nan_fill_gaps : fill all gaps in spectrum with NaNs
;		/median_subract : subtract median from each signal
;		median_width = width of neighborhood used for median filtering
;		median_val1 = manually specified median value for var 1
;		median_val2 = manually specified median value for var 2
;		mingap = minimum gap in timeseries used to determine streaks in data
;		cthresh = coherence threshold.  anything with coherence below cthresh
;				is set to NaN in both the coherence and phase spectrums
;		/verbose
;
;
;	HISTORY (fields_xspec)
;	v. 1.0, John Bonnell, UCBSSL, 22 Feb 2011.
;	-- initial version, based on FIELDS_SPEC, v. 2.1.
;
;
;	HISTORY (rbsp_xspec)
;   v. 1.0, Kris Kersten, UMN, June 2012
;	-- branched from JB's fields_xspec, modified for RBSP burst captures
;	-- added usage info
;	-- now depends on match.pro from NASA's IDLastro library at:
;		http://idlastro.gsfc.nasa.gov/
;		(used to find timeseries intersections)
;
;-
;
; DEPENDENCIES:
;	cross_spec.pro
;	defined.pro (dependency of cross_spec.pro)
;	match.pro from NASA's IDLastro library at: http://idlastro.gsfc.nasa.gov/
;
;-

pro rbsp_xspec, tplot_var1, tplot_var2, $
	tplot_var_xspec=tplot_var_xspec, $
	npts=npts, n_ave=n_ave, $
	tspec=tspec, xspec_coh=xspec_coh, xspec_pha=xspec_pha, freq=freq, df=df, $
	nan_fill_gaps=nan_fill_gaps, $
	median_subtract=median_subtract, median_width=median_width, $
	median_val1=median_val1, median_val2=median_val2, $
	mingap=mingap, cthresh=cthresh, $
	verbose=verbose

	; check for valid input.
	if not keyword_set( tplot_var1) or not keyword_set( tplot_var2) then begin
		if keyword_set( verbose) then $
			message, /info, 'TPLOT_VAR1 and TPLOT_VAR2 must be set; exiting.'
		return
	endif

	; setup default keyword values.
	if not keyword_set( npts) then $
		npts = 512L

	if not keyword_set( n_ave) then $
		n_ave = 2L

	; determine number of points consumed per spectrum calculated.
	pts_per_spec = npts*n_ave

	; determine number of frequency bins.
	nf = npts/2L + 1L

	; retrieve the time series from their respective TPLOT variables.
	; assumes that TPLOT_VAR1 and _VAR2 have the same sample rates and times.
	get_data, tplot_var1, data=d1, lim=lim1, dlim=dlim1
	get_data, tplot_var2, data=d2, lim=lim2, dlim=dlim2

	t1temp = d1.x
	t2temp = d2.x
	
	x1temp = d1.y
	x2temp = d2.y

	; We need to find intervals where the timetags in t1 and t2 line up.
	; It seems like NASA idlastro match.pro routine should work for this.
	; i1 and i2 return the indices of matching elements in t1temp and t2temp
	match,t1temp,t2temp,i1,i2

	x1=x1temp[i1]
	x2=x2temp[i2]
	t1=t1temp[i1]
	t2=t2temp[i2] ; not necessary.  after matching, t2=t1

	; subtract median from x1, x2?
	if keyword_set( median_subtract) then begin
		if keyword_set( median_width) then begin
			median_val1 = median( x1, median_width)
			median_val2 = median( x2, median_width)
		endif else begin
			median_val1 = median( x1)
			median_val2 = median( x2)
		endelse
		x1 = x1 - median_val1
		x2 = x2 - median_val2
	endif


	; find streaks in the data
	nt=size(t1,/n_elements)
	dts=t1[1L:nt-1L]-t1[0L:nt-2L]
	mediandt=median(dts)

	; set a sensible minimum gap size for finding breaks in the timeseries
	if not(keyword_set(mingap)) then mingap=100.*mediandt
	gaps=where(dts gt mingap, ngaps)

	; set up streak start/stop indices
	nstreaks=ngaps+1
	streak_start = [ 0L, gaps+1]
	streak_stop = [ gaps, nt-1L ]
	streak_len = [ streak_stop - streak_start + 1L ]

	; count the number of specs required
	nspec = 0L
	for i=0L,nstreaks-1L do begin
		nspec = nspec + streak_len[ i]/pts_per_spec
	endfor

	if nspec eq 0L then begin
		message, /info, 'No streak was long enough to accomodate a spectrum; exiting.'
		spec = !values.f_nan
		tspec = !values.d_nan
		freq = !values.f_nan
		df = !values.f_nan

		return
	endif

	; allocate space for XSPEC coherence and phase arrays, along with XSPEC time tag array.
	xspec_coh = fltarr( nspec, nf) + !values.f_nan
	xspec_pha = fltarr( nspec, nf) + !values.f_nan
	tspec = dblarr( nspec) + !values.d_nan


	; now calculate cross spec for each streak
	ispec = 0L
	for i=0L,nstreaks-1L do begin

		ns = streak_len[ i]/pts_per_spec

		for j=0L,ns-1L do begin

			j1 = streak_start[ i] + j*pts_per_spec
			j2 = streak_start[ i] + (j+1L)*pts_per_spec - 1L
			dt = t1[ j1 + 1L] - t1[ j1]

			power_spec, x1[ j1:j2], $
				analytic_sig=analytic_sig, $
	    		n_ave=n_ave, npts=npts, $
	    		sample = dt, ff, ss, /over

            cross_spec, x1[ j1:j2], x2[ j1:j2], $
              n_ave=n_ave, $
              npts=npts, sample = dt, coh, phase, freq, $
              /overlap, $
              verbose=verbose

	    	xspec_coh[ ispec, *] = coh[ *]
	    	xspec_pha[ ispec, *] = phase[ *]
	    	tspec[ ispec] = 0.5*( t1[ j1] + t1[ j2])

			ispec = ispec + 1L

		endfor

	endfor

	freq = ff
	df = freq[ 1L] - freq[ 0L]


	; kill phase info for bins where coherenece is less than some threshold
	if not(keyword_set(cthresh)) then begin
		cthresh=0.8
	endif else begin
		if cthresh ge 1. then begin
			message,'Coherence threshold must be less than 1.  Ignoring threshold.'
			cthresh=0.
		endif
	endelse
	
	ilow=where(xspec_coh lt cthresh, nlow)
	if nlow ne 0 then begin
		xspec_pha[ilow]=!values.f_nan
		xspec_coh[ilow]=!values.f_nan
	endif


	; fill gaps in the spectrum with NaNs?
	if keyword_set(nan_fill_gaps) then begin

		; find streaks in the spec so we can sandwich them in NaNs
		nt=size(tspec,/n_elements)
		dts=tspec[1L:nt-1L] - tspec[0L:nt-2L]
		mediandt=median(dts)
		
		; find gaps, ignoring any gaps less than 100x the normal spec cadence
		gaps=where(dts gt 100.*mediandt, ngaps)
	
		NaNs=fltarr(1,nf)
		NaNs[0:*]=!values.f_nan
		
		; sandwich arrays in NaNs
		xspec_coh=[NaNs, xspec_coh, NaNs]
		xspec_pha=[NaNs, xspec_pha, NaNs]
		tspec=[tspec[0]-0.1*mediandt, tspec, tspec[nt-1]+0.1*mediandt]
		
		; shift the gaps index by 1 for the prefixed NaN
		gaps=gaps+1L
		
	
		if ngaps ne 0 then begin
		
			new_coh=[xspec_coh[0L:gaps[0L],*], NaNs]
			new_pha=[xspec_pha[0L:gaps[0L],*], NaNs]
			newt=[tspec[0L:gaps[0L]], tspec[gaps[0L]]+0.1*mediandt]
	
			for gapcount=1L,ngaps-1L do begin
				new_coh=[new_coh, NaNs, xspec_coh[gaps[gapcount-1L]+1L:gaps[gapcount],*], NaNs]
				new_pha=[new_pha, NaNs, xspec_pha[gaps[gapcount-1L]+1L:gaps[gapcount],*], NaNs]
				newt=[newt, tspec[gaps[gapcount-1L]+1L]-0.1*mediandt, $
						tspec[gaps[gapcount-1L]+1L:gaps[gapcount]], $
						tspec[gaps[gapcount]]+0.1*mediandt]
			endfor
			new_coh=[new_coh, NaNs, xspec_coh[gaps[ngaps-1L]+1L:nt,*]]
			new_pha=[new_pha, NaNs, xspec_pha[gaps[ngaps-1L]+1L:nt,*]]
			newt=[newt, tspec[gaps[ngaps-1L]+1L]-0.1*mediandt, tspec[gaps[ngaps-1L]+1L:nt]]
	
		endif
	
		xspec_coh=new_coh
		xspec_pha=new_pha
		tspec=newt

	endif

	; Let's do this in degrees
	xspec_pha=xspec_pha/!dtor

	; Now pack spec, tspec, and freq into a tplot variable.
	str_element, lim, 'spec', 1, /add
	str_element, lim, 'zlog', 0, /add
	str_element, lim, 'ylog', 0, /add
	str_element, lim, 'yrange', [ 0., max( freq)], /add
	str_element, lim, 'ystyle', 1, /add
	str_element, lim, 'x_no_interp', 1, /add
	str_element, lim, 'y_no_interp', 1, /add


	if not keyword_set( tplot_var_xspec) then $
		tplot_var_xspec = string( tplot_var1, tplot_var2, format='(A,"_X_",A)')

	; options for Coherence:
	str_element, lim, 'ztitle', $
		string( 'COH', tplot_var1, tplot_var2, format='(A,"!C!C[",A,"x",A,"]")'), /add

	tplot_var_xspec_coh = tplot_var_xspec + '_coh'

	store_data, tplot_var_xspec_coh, $
		data={ x:tspec, y:xspec_coh, v:freq, df:df, npts:npts, n_ave:n_ave }, $
		lim=lim, dlim=dlim
	options, tplot_var_xspec_coh, 'zrange', [ cthresh, 1.]

	; options for Phase:
	str_element, lim, 'ztitle', $
		string( 'PHA', tplot_var1, tplot_var2, 'deg', format='(A,"!C!C[",A,"x",A,X,A,"]")'), /add

	tplot_var_xspec_pha = tplot_var_xspec + '_pha'

	store_data, tplot_var_xspec_pha, $
		data={ x:tspec, y:xspec_pha, v:freq, df:df, npts:npts, n_ave:n_ave }, $
		lim=lim, dlim=dlim
	options, tplot_var_xspec_pha, 'zrange', [ -180., 180.]

return
end
