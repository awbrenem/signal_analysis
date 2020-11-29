;+
; PROCEDURE rbsp_spec
;
; PURPOSE: Calculates power spectrum for RBSP waveform data
;
; CALLING SEQUENCE: rbsp_spec,tplot_var
;
; INPUTS:
;		tplot_var - TPLOT variable name or number
;
; OUTPUTS:
;
; KEYWORDS:
;		tplot_var_spec = string name for new TPLOT XSPEC variable
;		npts = number of points per FFT
;		n_ave = number of specs to average for each spectrum
;		tspec = return variable for spectrum times
;		spec = return variable for power spectrum
;		freq = return variable for frequency bins
;		df = return variable for frequency steps
;		/nan_fill_gaps : fill all gaps in spectrum with NaNs
;		/median_subract : subtract median from each signal
;		median_width = width of neighborhood used for median filtering
;		median_val = manually specified median value tplot_var
;		mingap = minimum gap in timeseries used to determine streaks in data
;		/verbose
;
;	HISTORY (fields_spec)
;	v. 1.0, John Bonnell, UCBSSL, 2 Sept 2005.
;	-- initial version.
;	v. 1.1, John Bonnell, UCBSSL, 2 Sept 2005.
;	-- added streak handling code.
;	v. 2.0, JWB, UCBSSL, 10 July 2008.
;	-- added analytical signal handling code.
;	v. 2.1, JWB, UCB SSL, 16 Oct 2009.
;	-- added MEDIAN_SUBTRACT option.
;
;
;	HISTORY (rbsp_spec)
;   v. 1.0, Kris Kersten, UMN, June 2012
;	-- branched from JB's fields_spec, modified for RBSP burst captures
;	-- removed incomplete analytic_sig code
;	-- now depends on match.pro from NASA's IDLastro library at:
;		http://idlastro.gsfc.nasa.gov/
;		(used to find timeseries intersections)
;
;	TO DO:
;	-- UNITS!!!!!!!!!!!!!!!
;
;-
;
; DEPENDENCIES:
;	power_spec.pro
;	defined.pro (dependency of power_spec.pro)
;	match.pro from NASA's IDLastro library at: http://idlastro.gsfc.nasa.gov/
;
;
;
	;AARON'S NOTE *******************************************************
	;SEE NOTE BELOW
	;THERE'S A PROBLEM WITH ispec. ispec CAN BECOME LARGER THAN THE tspec and spec ARRAYS
	;BECAUSE OF THE i=0,nstreaks-1 LOOP
;
	;**********************************************************************
;
;
;-

pro rbsp_spec, tplot_var, $
	tplot_var_spec=tplot_var_spec, $
	npts=npts, n_ave=n_ave, $
	tspec=tspec, spec=spec, freq=freq, df=df, $
	nan_fill_gaps=nan_fill_gaps, $
	median_subtract=median_subtract, median_width=median_width, median_val=median_val, $
	mingap=mingap, $
	verbose=verbose

	; check for valid input.
	if not keyword_set( tplot_var) then begin
		if keyword_set( verbose) then $
			message, /info, 'no TPLOT_VAR set; exiting.'
		return
	endif

	; setup default keyword values.
	if not keyword_set(npts) then $
		npts = 512L

	if not keyword_set(n_ave) then $
		n_ave = 2L

	; determine number of points consumed per spectrum calculated.
	pts_per_spec = double(npts)*double(n_ave)

	; determine number of frequency bins.
	nf = npts/2L + 1L

	; retrieve time series data from TPLOT var
	get_data, tplot_var, data=d, lim=lim, dlim=dlim
	t = d.x
	x = d.y

	if keyword_set(median_subtract) then begin
		if keyword_set(median_width) then $
			median_val = median(x, median_width) else $
			median_val = median(x)
		x = x - median_val
	endif


;******This doesn't seem to work if there are no streaks********
;	; find streaks in the data
	nt=size(t,/n_elements)


	;; dts=t[1L:nt-1L]-t[0L:nt-2L]
	;; mediandt=median(dts)

	;; ; set a sensible minimum gap size for finding breaks in the timeseries
	;; if not(keyword_set(mingap)) then mingap=100.*mediandt


;****Aaron's testing*****
;dts[5:8] = mingap + 0.1
;*************************

	;; gaps=where(dts gt mingap, ngaps)

	;; ; set up streak start/stop indices (times b/t gaps in data..i.e. good data)
	;; nstreaks=ngaps+1
	;; streak_start = [0L,gaps+1]
	;; streak_stop = [gaps,nt-1L]


;****AARON CHANGED THESE
	streak_start = 0
	streak_stop = nt-1L
;***************************

	streak_len =  streak_stop - streak_start + 1L




	; count the number of specs required
	nspec = 0L
;	for i=0L,nstreaks-1L do nspec = nspec + streak_len[i,*]/pts_per_spec
;	for i=0L,nstreaks-1L do print,streak_len[i,*]/long(pts_per_spec)
;	for i=0L,nstreaks-1L do nspec = nspec + streak_len[i]/pts_per_spec
        nspec = streak_len/pts_per_spec
        if nspec eq 0 then nspec = 1

;	if nspec eq 0L then begin
;		message, /info, 'No streak was long enough to accomodate a spectrum; exiting.'
;		spec = !values.f_nan
;		tspec = !values.d_nan
;		freq = !values.f_nan
;		df = !values.f_nan
;
;		return
;	endif

	; allocate space for XSPEC coherence and phase arrays, along with XSPEC time tag array.
	spec = fltarr( nspec, nf) + !values.f_nan
	tspec = dblarr( nspec) + !values.d_nan


;	; now calculate spectrum for each streak
	ispec = 0L
;	for i=0L,nstreaks-1L do begin

;		ns = streak_len[i,*]/pts_per_spec
		ns = streak_len/pts_per_spec


;.compile /Users/aaronbreneman/Desktop/code/Aaron/rbsp/cribs/TPLOT_RBSP_XSPEC_SPEC


;*****Aaron changed this
if ns eq 0 then ns = 1

		for j=0L,ns-1L do begin
			j1 = streak_start + j*pts_per_spec
			j2 = streak_start + (j+1L)*pts_per_spec - 1L
			dt = t[ j1 + 1L] - t[ j1]

			power_spec, x[ j1:j2], $
				analytic_sig=analytic_sig, $
				n_ave=n_ave, npts=npts, $
				sample = dt, ff, ss, /over
			spec[ ispec, *] = ss[ *]
			tspec[ ispec] = 0.5*( t[ j1] + t[ j2])

			ispec = ispec + 1L

		endfor


;	STOP
	;AARON'S NOTE *******************************************************

	;THERE'S A PROBLEM HERE. ispec CAN BECOME LARGER THAN THE tspec and spec ARRAYS
	;BECAUSE OF THE i=0,nstreaks-1 LOOP

	;**********************************************************************

	;endfor

	freq = ff
	df = freq[ 1L] - freq[ 0L]

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
		spec=[NaNs, spec, NaNs]
		tspec=[tspec[0]-0.1*mediandt, tspec, tspec[nt-1]+0.1*mediandt]

		; shift the gaps index by 1 for the prefixed NaN
		gaps=gaps+1L


		if ngaps ne 0 then begin

			new_spec=[spec[0L:gaps[0L],*], NaNs]
			newt=[tspec[0L:gaps[0L]], tspec[gaps[0L]]+0.1*mediandt]

			for gapcount=1L,ngaps-1L do begin
				new_spec=[new_spec, NaNs, spec[gaps[gapcount-1L]+1L:gaps[gapcount],*], NaNs]
				newt=[newt, tspec[gaps[gapcount-1L]+1L]-0.1*mediandt, $
						tspec[gaps[gapcount-1L]+1L:gaps[gapcount]], $
						tspec[gaps[gapcount]]+0.1*mediandt]
			endfor
			new_spec=[new_spec, NaNs, spec[gaps[ngaps-1L]+1L:nt,*]]
			newt=[newt, tspec[gaps[ngaps-1L]+1L]-0.1*mediandt, tspec[gaps[ngaps-1L]+1L:nt]]

			spec=new_spec
			tspec=newt

		endif


	endif


	; pack spec, tspec, and freq into a tplot variable.
	str_element, lim, 'spec', 1, /add
	str_element, lim, 'zlog', 1, /add
	str_element, lim, 'ylog', 0, /add
	str_element, lim, 'yrange', [ 0., max( freq)], /add
	str_element, lim, 'ystyle', 1, /add
	str_element, lim, 'x_no_interp', 1, /add
	str_element, lim, 'y_no_interp', 1, /add


; STILL NEED TO WORK OUT UNITS FOR RBSP
;	units = thm_get_units( tplot_var, default_val=default_val)
;	spec_units = string( units, format='("(",A,")!U2!N/Hz")')
;	str_element, lim, 'ztitle', $
;		string( tplot_var, spec_units, format='(A,"!C!C[",A,"]")'), /add
;	str_element, dlim.data_att.units, spec_units, /add

	str_element, lim, 'ytitle', string( tplot_var, format='(A,"!C!C[Hz]")'), /add

;	store_data, tplot_var+'_SPEC', $
;		data={ x:tspec, y:spec, v:freq, df:df, npts:npts, n_ave:n_ave, units:spec_units }, $
;		lim=lim, dlim=dlim

	if not(keyword_set(tplot_var_spec)) then tplot_var_spec=tplot_var+'_SPEC'
	store_data, tplot_var_spec, $
		data={ x:tspec, y:spec, v:freq, df:df, npts:npts, n_ave:n_ave }, $
		lim=lim, dlim=dlim

	return

end
