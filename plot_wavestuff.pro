;+
;*****************************************************************************************
;
;  PROCEDURE :  plot_wavestuff
;  PURPOSE  :   Creates tplot variables of and plot really nice waveforms,
;				spectrograms, wavelets, hodograms, power spectra
;				Times in sec and freqs in Hz
;
;
;  REQUIRES:
;
;
;  INPUT:		wave -> tplot variable with waveform (ex. 'wave')
;
;
;  KEYWORDS:
;
;				PLOT TOGGLE OPTIONS
;					hod --> set to plot hodograms
;					spec --> set to plot freq/time spectrograms (freq in Hz)
;					wavelet --> set to plot freq/time wavelets
;					psd --> set to plot power spectral density
;					wf --> set to plot 3 components of waveform. If none of the above plot options
;						are set then this defaults to 'y'
;
;				GENERAL PLOT OPTIONS FOR TIMESERIES TYPE PLOTS
;					freqs --> frequencies to overlay onto plots (Hz)
;					freq_labels --> the names of the "freqs"  (NOT YET IMPLEMENTED)
;					vline -> vertical lines to overlay on plots (unix time) [Probs should set vvarname too]
;					hline -> NOT WORKING YET...horizontal lines to overlay on plots [Must also set hvarname]
;					hvarname -> name of the tplot variable in which to plot the hline
;					vvarname -> ''      ''			''			''			''  vline
;
;				SLIDING SPEC OPTIONS
;					npts -> number of points per FFT (defaults to 512)
;					n_ave -> number of specs to average for each spectrum
;
;					**************
;					OBSOLETE....slide_spec.pro options
;					pets -> stepsize for slide_spec.pro
;					step_overlap -> step_overlap from slide_spec.pro
;					mincomp, maxcomp -> Use these to "squeeze" out more color in sliding spectrograms.
;						These are used as: 	csx=bytscl(sldx,min=min(sldx*mincomp,/nan),max=max(sldx/maxcomp,/nan))
;						Setting these can be useful if you want to see more color contrast in certain
;						frequency bins that may not have as much power as others.
;						Maybe start by trying mincomp=10, maxcomp=10 ??
;					**************
;
;				PSD OPTIONS
;					nowindow -> Don't use windowing in PSD
;					Epsd -> _extra structure used for the call to fft_power_calc
;					vals_PSD -> array of returned values for the PSD
;
;				WAVELET OPTIONS
;					samescale -> 'x','y', or 'z' Use the same power scaling on all three
;							 wavelets (defined by the input component)
;					fac -> scale factor for the wavelet power (minval/fac,maxval*fac)
;
;				HODOGRAM OPTIONS
;					extra1vec -> an [3] vector to overplot on hodograms (BLUE)
;					extra2vec -> same as vec (GREEN)
;					extra2vec -> same as vec (ORANGE)
;					extra3vec -> same as vec (RED)
;
;
;				FILE OUTPUT OPTIONS
;					postscript --> set to only plot the .ps version
;					root --> root directory to plot in. Defaults to current directory
;					filename --> name of file to be saved (don't include path)
;
;
;				WINDOW DELETE OPTIONS
;					You can call plot_wavestuff with the intention of only deleting a particular
;					window.
;					del_wf -> set to delete waveform plot
;					del_hod
;					del_spec
;					del_wavelet
;					del_psd
;
;				_EXTRA STRUCTURES FOR PLOT OPTIONS
;					extra_wf
;					extra_hod
;					extra_spec
;					extra_wavelet
;					extra_psd
;
;				OTHER KEYWORDS
;					dscale -> change the dscale of the wavelet. Defaults to 1 for speed.
;							  0.1 gives much nicer values.
;					plotcg -> plot using Dfanning's cgwindow function. (WORK IN PROGRESS)
;
;					nodelete -> don't delete the tplot variables at end.
;					force_N -> pass via an extra structure to zero pad
;
;  EXAMPLES:    Since all these print out as tplot variables, use options to change
;				plot parameters. This program prints out a list of tplot variables that
;				it has created so you know what variable to change.
;					Ex. options,'slide_specx','yrange',[0,1]
;
;
;   NOTES:		ADD IN RBSP_SPEC AS AN OPTION
;
;
;
;    CREATED:  02/02/2011
;    CREATED BY:  Aaron W. Breneman
;    LAST MODIFIED:     v1.0.0
;    MODIFIED BY: Aaron W. Breneman -08/15/2011 -> all input times in sec and freqs in kHz
;									-12/29/2011 -> added nowindow option
;									-01/24/2011 -> Added "vec" keyword
;									------------------------------------------------------
;									-06/26/2012 -> Major update: quantities are now turned
;												   into tplot variables and plotted using
;												   the tplot routines to take advantage of
;												   the nice time formatting
;									-04/09/2013 -> Can now overplot vlines on the tplotxy PSD plots
;									-05/06/2013 -> Fixed ps plotting. Other small changes
;
;*****************************************************************************************
;-
;*****************************************************************************************


;pro psd_plots,sz=sz

;	tplotxy,'psd_x',multi='3 1',/noisotropic
;	if sz gt 1 then tplotxy,'psd_y',/addpanel,ymargin=[0.2,0.2],/noisotropic;,xlog=xlog,ylog=ylog,xrange=xr,yrange=yr
;	if sz gt 2 then tplotxy,'psd_z',/addpanel,ymargin=[0.2,0.2],/noisotropic;,xlog=xlog,ylog=ylog,xrange=xr,yrange=yr

;	tplotxy,'psd_x',multi='3 1',ymargin=[0.2,0.2],xsize=1100,ysize=300,/noisotropic,$
;		xstyle=1
;	if sz gt 1 then tplotxy,'psd_y',/addpanel,ymargin=[0.2,0.2],/noisotropic
;	if sz gt 2 then tplotxy,'psd_z',/addpanel,ymargin=[0.2,0.2],/noisotropic
;end


pro plot_wavestuff,wave,$
	hod=hod,$
	spec=spec,$
	wavelet=wavelet,$
	psd=psd,$
	wf=wf,$
	freqs=freqs,$
	freq_labels=freq_labels,$
	vline=vline,$
	hline=hline,$
	hvarname=hvarname,$
	vvarname=vvarname,$
	extra1vec=vec1,$    ;Blue
	extra2vec=vec2,$    ;Green
	extra3vec=vec3,$    ;Orange
	extra4vec=vec4,$    ;Red
	root=root,$
	filename=filename,$
	npts=npts,$
	n_ave=n_ave,$
	other=other,$
	var_other=var_other,$
	dscale=ds,$
	nocb=nocb,$
	postscript=postscript,$
	samescale=samescale,$
	fac=fac,$
	nowindow=nowindow,$
	vals_PSD=vals_PSD,$
	del_wf=del_wf,del_spec=del_spec,del_wavelet=del_wavelet,del_hod=del_hod,del_psd=del_psd,del_other=del_other,$
	extra_wf=extra_wf,extra_hod=extra_hod,extra_spec=extra_spec,extra_psd=extra_psd,extra_wavelet=extra_wavelet,extra_other=extra_other,$
	plotcg=cg,$
	nodelete=nodelete,$
	Epsd=Epsd,$
	waveform_ytitles=wfyt,$
	waveform_title=wft,$
	title_hodogram=hodt,$
	zlim_spec=zlim_spec
	;,Ewf=Ewf,Ehod=Ehod,EspecE=Espec,Ewavelet=Ewavelet


;Set so that I can restore defaults upon returning
pinit = !p


get_data,wave,data=wave2


;if ~keyword_set(filename) then filename = ''
if ~keyword_set(root) then root = '~/Desktop/'
if ~keyword_set(nocb) then nocb=0 else nocb=1
if ~keyword_set(ds) then ds = 1.

if ~keyword_set(postscript) then ps = 0 else ps = 1
if ~keyword_set(samescale) then samescale = 'no'
if ~keyword_set(fac) then fac = 1

;just plot the waveform if no plots are explicitly requested
if ~keyword_set(wf) $
and ~keyword_set(spec) $
and ~keyword_set(wavelet) $
and ~keyword_set(other) $
and ~keyword_set(hod) $
and ~keyword_set(psd) then plot_wf = 1

if ~keyword_set(plot_wf) and ~keyword_set(wf) then plot_wf = 0
if keyword_set(wf) then plot_wf = 1

if ~keyword_set(hod) then plot_hod = 0 else plot_hod = 1
if ~keyword_set(spec) then plot_spec = 0 else plot_spec = 1
if ~keyword_set(wavelet) then plot_wavelet = 0 else plot_wavelet = 1
if ~keyword_set(other) then plot_other = 0 else plot_other = 1
if ~keyword_set(psd) then plot_psd = 0 else plot_psd = 1

if ~keyword_set(del_wf) then del_wf = 0
if ~keyword_set(del_hod) then del_hod = 0
if ~keyword_set(del_spec) then del_spec = 0
if ~keyword_set(del_wavelet) then del_wavelet = 0
if ~keyword_set(del_psd) then del_psd = 0

if keyword_set(cg) then cgdelete,/all


;make sure "root" has a trailing slash
goo = strpos(root,'/',/reverse_search)
len = strlen(root)
if goo+1 ne len then root = root + '/'


if keyword_set(ps) then !p.charsize = 0.8 else !p.charsize = 1.2
device,decomposed=0
loadct,39


tmp = size(wave2.y)

if tmp[0] eq 1 then sz = 1
if tmp[0] eq 2 then sz = tmp[2]


tplot_options,'xmargin',[20.,15.]
tplot_options,'ymargin',[3,6]
tplot_options,'xticklen',0.08
tplot_options,'yticklen',0.02
tplot_options,'xthick',1
tplot_options,'ythick',1
tplot_options,'labflag',-1
tplot_options,'title',''

;Dummy plot to define !tplotxy variable
;tplot,wave


span = wave2.x[n_elements(wave2.x)-1] - wave2.x[0]
timespan,wave2.x[0],span,/seconds




;------------------------------------
;SET UP WINDOWS IF NOT PLOTTING TO PS
;------------------------------------

;Determine the size and location of the windows based on screen size
ss = get_screen_size()

	win_wf = 1
	win_psd = 2
	win_spec = 3
	win_wavelet = 4
	win_hod = 5
	win_other = 6



if ~keyword_set(cg) then begin



	;WAVEFORM PLOTS
	if plot_wf then begin
		xs=floor(ss[0]/2.5)
		ys=floor(ss[1]/2.3)
		xp = ss[0]/5.
		yp=floor(ss[1])/25
		if ~ps then window,win_wf,xsize=xs,ysize=ys,xpos=xp,ypos=yp,retain=2,title='waveform'

	endif

	;PSD PLOTS
	if plot_psd then begin
		xs=ss[0]/2.
		ys=ss[1]/4.
;		xp = 250.
		if ~ps then window,win_psd,xsize=xs,ysize=ys,retain=2,title='PSD'
	endif

	;SPEC PLOTS
	if plot_spec then begin
		xs=floor(ss[0]/2.5)
		ys=floor(ss[1]/2.3)
		xp=floor(ss[0]-xs)
		yp=floor(ss[1])/25
		if ~ps then window,win_spec,xsize=xs,ysize=ys,xpos=xp,ypos=yp,retain=2,title='spec'
	endif

	;WAVELET PLOTS
	if plot_wavelet then begin
		xs=floor(ss[0]/2.5)
		ys=floor(ss[1]/2.3)
		xp=floor(ss[0]-xs)
		yp = ss[1]/1.7
		if ~ps then window,win_wavelet,xsize=xs,ysize=ys,xpos=xp,ypos=yp,retain=2,title='wavelet'
	endif

	;HOD PLOTS
	if plot_hod then begin
		xs=ss[0]/7.5
		ys=ss[1]/1.7
		xp = ss[0]/40.
		yp = ss[1]/15
;		xp = 0  ;ss[0]-xsize
;		yp = ss[1]
		if ~ps then window,win_hod,xsize=xs,ysize=ys,xpos=xp,ypos=yp,retain=2,title='hodograms'
	endif

	;OTHER PLOTS
	if plot_other then begin
		xs=floor(ss[0]/2.5)
		ys=floor(ss[1]/2.3)
		xp=floor(ss[0]-xs)
		yp=floor(ss[1])/1.7
		if ~ps then window,win_other,xsize=xs,ysize=ys,xpos=xp,ypos=yp,retain=2,title='other'
	endif


endif




;******should this be wave2.x or wave2.y??????????*********
nelem = n_elements(wave2.x)


windowarr=lonarr(n_elements(wave2.y[*,0]))
my_windowf,nelem-1,2,windowarr
;1 = Hamming
;2 = Hanning



;---------------
;WAVEFORM PLOTS
;---------------

if plot_wf then begin

	;put into tplot variables
	store_data,'fieldx',data={x:wave2.x,y:wave2.y[*,0]}
	if sz gt 1 then store_data,'fieldy',data={x:wave2.x,y:wave2.y[*,1]}
	if sz gt 2 then store_data,'fieldz',data={x:wave2.x,y:wave2.y[*,2]}



	;Set up default options
	tplot_options,'title','Waveform(s)'

	if keyword_set(wfyt) then begin
		options,'fieldx','ytitle',wfyt[0]
		options,'fieldy','ytitle',wfyt[1]
		options,'fieldz','ytitle',wfyt[2]
	endif else begin
		options,'fieldx','ytitle','waveform!Cx'
		options,'fieldy','ytitle','waveform!Cy'
		options,'fieldz','ytitle','waveform!Cz'
	endelse

	if KEYWORD_SET(wft) then tplot_options,'title',wft
	;Overwrite options with _extra structure
;	if is_struct(extra_wf) then options,['fieldx','fieldy','fieldz'],_extra=extra_wf
;	if is_struct(extra_wf) then tplot_options,['fieldx','fieldy','fieldz'],_extra=extra_wf



	if keyword_set(cg) then begin
		xs=ss[0]/2.5
		ys=ss[1]/2.
		xp=0
		yp=ss[1]-1.7*ss[1]/3.

		cgwindow,wxsize=xs,wysize=ys,wxpos=xp,wypos=yp


		if sz eq 1 then cgwindow,'tplot','fieldx',window=-1,/addcmd
		if sz eq 2 then cgwindow,'tplot',['fieldx','fieldy'],window=-1,/addcmd
		if sz eq 3 then cgwindow,'tplot',['fieldx','fieldy','fieldz'],window=-1,/addcmd
	endif else begin


		;figure out .ps stuff
		if ~ps then wset,win_wf
		if ps then !p.font = 0 else !p.font = 1
		if ~ps then device,set_font = '7X13'  ;built in Unix font

		if ps and keyword_set(filename) then popen,root + 'wf_' + filename
		if ps and ~keyword_set(filename) then popen,root + 'wf'

		if sz eq 1 then tplot,'fieldx',window=-1
		if sz eq 2 then tplot,['fieldx','fieldy'],window=-1
		if sz eq 3 then tplot,['fieldx','fieldy','fieldz'],window=-1

		if ps then pclose


	endelse

		tplot_options,'title',''

		print,'WAVEFORM TPLOT VARIABLES PLOTTED:'
		print,'...fieldx'
		if sz ge 2 then print,'...fieldy'
		if sz eq 3 then print,'...fieldz'


endif



;--------------------
;PSD PLOTS
;--------------------


if plot_psd then begin


	if keyword_set(nowindow) then begin
		power_x = fft_power_calc(wave2.x,wave2.y[*,0],_extra=Epsd,FORCE_N=force_n)
		if sz gt 1 then power_y = fft_power_calc(wave2.x,wave2.y[*,1],_extra=Epsd,FORCE_N=force_n)
		if sz gt 2 then power_z = fft_power_calc(wave2.x,wave2.y[*,2],_extra=Epsd,FORCE_N=force_n)
	endif else begin
			power_x = fft_power_calc(wave2.x,windowarr*wave2.y[*,0],_extra=Epsd,FORCE_N=force_n)
			if sz gt 1 then power_y = fft_power_calc(wave2.x,windowarr*wave2.y[*,1],_extra=Epsd,FORCE_N=force_n)
			if sz gt 2 then power_z = fft_power_calc(wave2.x,windowarr*wave2.y[*,2],_extra=Epsd,FORCE_N=force_n)
	endelse

	if keyword_set(nowindow) then begin
		power_x = fft_power_calc(wave2.x,wave2.y[*,0],FORCE_N=force_n)
		if sz gt 1 then power_y = fft_power_calc(wave2.x,wave2.y[*,1],FORCE_N=force_n)
		if sz gt 2 then power_z = fft_power_calc(wave2.x,wave2.y[*,2],FORCE_N=force_n)
	endif else begin
		power_x = fft_power_calc(wave2.x,windowarr*wave2.y[*,0],FORCE_N=force_n)
		if sz gt 1 then power_y = fft_power_calc(wave2.x,windowarr*wave2.y[*,1],FORCE_N=force_n)
		if sz gt 2 then power_z = fft_power_calc(wave2.x,windowarr*wave2.y[*,2],FORCE_N=force_n)
	endelse




	;Remove vlines outside of current plot range
	if keyword_set(vline) then begin
		xr = [0,0]
		if is_struct(extra_psd) then if keyword_set(extra_psd.xrange) then xr = extra_psd.xrange $
								else xr = [min(power_x.freq),max(power_x.freq)]

		goo = where(vline le xr[0])
		if goo[0] ne -1 then vline[goo] = !values.f_nan
		goo = where(vline ge xr[1])
		if goo[0] ne -1 then vline[goo] = !values.f_nan


	endif


;	if keyword_set(cg) then begin
;
;		!p.charsize = 1.6
;		xs=floor(ss[0]/2.5)
;		ys=floor(ss[1]/1.5)
;		xp=floor(ss[0]-xs)
;		yp=floor(ss[1])/20
;
;		cgwindow,wxsize=xs,wxpos=xp,wypos=yp,wmulti=[0,3,0],waspect=1/2.8
;
;		cgwindow,'plot',power_x.freq,power_x.power_x,/addcmd,title='x'
;		if sz ge 2 then cgwindow,'plot',power_y.freq,power_y.power_x,/addcmd
;		if sz eq 3 then cgwindow,'plot',power_z.freq,power_z.power_x,/addcmd
;	endif else begin



		;Create a fake array of times so I can use the tplotxy routine
		times2 = interpol(wave2.x,n_elements(power_x.freq))


		;powvals = {x:times2,y:[[power_x.freq],[power_x.power_x]]}
		store_data,'psd_x',data={x:times2,y:[[power_x.freq],[power_x.power_x]]}
		if sz gt 1 then store_data,'psd_y',data={x:times2,y:[[power_y.freq],[power_y.power_x]]}
		if sz gt 2 then store_data,'psd_z',data={x:times2,y:[[power_z.freq],[power_z.power_x]]}

		fmin = 0.1
		fmax = max(power_x.freq)



		secs = times2[n_elements(times2)-1] - times2[0]
		timespan,times2[0],secs,/seconds

		;Set default options
		tplot_options,'title','Frequency power spectrum'
		options,'psd_x','ytitle','PSD x!CdB vs freq (Hz)'
		options,'psd_y','ytitle','PSD y!CdB vs freq (Hz)'
		options,'psd_z','ytitle','PSD z!CdB vs freq (Hz)'
		xlog = 0
		ylog = 0
;		xr = [0,max(power_x.freq)]
		xr = [0,0.004]
		yr = [-20,60]

		;Set options with _extra structure
		if is_struct(extra_psd) then options,['psd_x','psd_y','psd_z'],_extra=extra_psd
		if is_struct(extra_psd) and tag_exist(extra_psd,'ylog') then ylog = extra_psd.ylog
		if is_struct(extra_psd) and tag_exist(extra_psd,'xlog') then xlog = extra_psd.xlog
		if is_struct(extra_psd) and tag_exist(extra_psd,'yrange') then yr = extra_psd.yrange
		if is_struct(extra_psd) and tag_exist(extra_psd,'xrange') then xr = extra_psd.xrange


		if ~ps and ~keyword_set(cg) then wset,win_psd
		if ps then !p.font = 0
		if ps and keyword_set(filename) then popen,root + 'psd_' + filename,/landscape
		if ps and ~keyword_set(filename) then popen,root + 'psd',/landscape
		if ps then multiv = '3 2' else multiv = '3 1'


		if keyword_set(cg) then begin
		    cgwindow,'tplotxy','psd_x',wmulti=[0,3,1,0,0],wxsize=1100,wysize=300,_extra={multi:'3 1',ymargin:[0.2,0.2],$
												  noisotropic:1,xstyle:1,xlog:xlog,$
												  ylog:ylog,xrange:xr,yrange:yr}
		    if sz gt 1 then cgwindow,'tplotxy','psd_y',_extra={addpanel:1,ymargin:[0.2,0.2],noisotropic:1,xlog:xlog,ylog:ylog,$
									xrange:xr,yrange:yr},/addcmd
		    if sz gt 2 then cgwindow,'tplotxy','psd_z',_extra={addpanel:1,ymargin:[0.2,0.2],noisotropic:1,xlog:xlog,ylog:ylog,$
									xrange:xr,yrange:yr},/addcmd
		endif else begin
		    tplotxy,'psd_x',multi=multiv,ymargin=[0.2,0.2],xsize=1100,ysize=300,/noisotropic,$
			    xstyle=1,xlog=xlog,ylog=ylog,xrange=xr,yrange=yr
				if keyword_set(vline) then for qq=0,n_elements(vline)-1 do plots,[[vline[qq],yr[0]],[vline[qq],yr[1]]]
		    if sz gt 1 then begin
		    	tplotxy,'psd_y',/addpanel,ymargin=[0.2,0.2],/noisotropic,xlog=xlog,ylog=ylog,xrange=xr,yrange=yr
				if keyword_set(vline) then for qq=0,n_elements(vline)-1 do plots,[[vline[qq],yr[0]],[vline[qq],yr[1]]]
			endif
		    if sz gt 2 then begin
		    	tplotxy,'psd_z',/addpanel,ymargin=[0.2,0.2],/noisotropic,xlog=xlog,ylog=ylog,xrange=xr,yrange=yr
				if keyword_set(vline) then for qq=0,n_elements(vline)-1 do plots,[[vline[qq],yr[0]],[vline[qq],yr[1]]]
			endif

		endelse

	;	tplotxy,'psd_x',multi='3 1',ymargin=[0.2,0.2],xsize=1100,ysize=300,/noisotropic,$
	;		yrange=[-80,10],xtitle='freq (Hz)',ytitle='X dB',$
	;		title='PSD for ' + filename,window=2,/xlog,xrange=[fmin,fmax],xstyle=1
	;	if sz gt 1 then tplotxy,'psd_y',/addpanel,ymargin=[0.2,0.2],/noisotropic,$
	;		yrange=[-80,10],xtitle='freq (Hz)',ytitle='Y dB',$
	;		title='PSD for ' + filename,/xlog;,xrange=[fmin,fmax],xstyle=1
	;	if sz gt 2 then tplotxy,'psd_z',/addpanel,ymargin=[0.2,0.2],/noisotropic,$
	;		yrange=[-80,10],xtitle='freq (Hz)',ytitle='Z dB',$
	;		title='PSD for ' + filename,/xlog;,xrange=[fmin,fmax],xstyle=1


	;**************
	;Vline and Hline don't work at the moment b/c I'm using tplotxy
	;	if keyword_set(vline) then timebar,vline,varname='psd_x'
	;	if keyword_set(hline) then timebar,hline,/databar,varname=hvarname


		if ps then pclose

		if sz eq 1 then vals_PSD = power_x
		if sz eq 2 then vals_PSD = [power_x,power_y]
		if sz eq 3 then vals_PSD = [power_x,power_y,power_z]


;	endelse
;

	tplot_options,'title',''

	print,'PSD TPLOT VARIABLES PLOTTED:'
	print,'...psd_x'
	if sz eq 2 then print,'...psd_y'
	if sz eq 3 then print,'...psd_z'

endif

;--------------------
;SPECTROGRAM PLOTS
;--------------------

if plot_spec then begin



	get_data,wave,data=dat

	store_data,['x','y','z'],/delete

	if sz ge 1 then store_data,'x',data={x:dat.x,y:reform(dat.y[*,0])}
	if sz ge 2 then store_data,'y',data={x:dat.x,y:reform(dat.y[*,1])}
	if sz ge 3 then store_data,'z',data={x:dat.x,y:reform(dat.y[*,2])}


	if sz ge 1 then rbsp_spec,'x',npts=npts,n_ave=n_ave,/nan_fill_gaps;,mingap=7d-5;,/median_subtract
	if sz ge 2 then rbsp_spec,'y',npts=npts,n_ave=n_ave,/nan_fill_gaps
	if sz ge 3 then rbsp_spec,'z',npts=npts,n_ave=n_ave,/nan_fill_gaps



	tplot_options,'title','Sliding spec(s)'

	;set default options
	tplot_options,'title','Sliding spectrum'
	options,'x_SPEC','ytitle','spec x!C(Hz)'
	options,'y_SPEC','ytitle','spec y!C(Hz)'
	options,'z_SPEC','ytitle','spec z!C(Hz)'


	;Set options with _extra structure
	if is_struct(extra_spec) then options,['x_SPEC','y_SPEC','z_SPEC'],_extra=extra_spec



	if keyword_set(cg) then begin
		xs=floor(ss[0]/2.5)
		ys=floor(ss[1]/2.)
		xp=floor(ss[0]-xs)
		yp=floor(ss[1])/20

		cgwindow,wxsize=xs,wysize=ys,wxpos=xp,wypos=yp

		if sz eq 1 then cgwindow,'tplot','x_SPEC',window=-1,/addcmd
		if sz eq 2 then cgwindow,'tplot',['x_SPEC','y_SPEC'],window=-1,/addcmd
		if sz eq 3 then cgwindow,'tplot',['x_SPEC','y_SPEC','z_SPEC'],window=-1,/addcmd
	endif else begin


		if ~ps then wset,win_spec
		if ps then !p.font = 0
;		if ps then popen,root + 'spec_' + filename + '.eps'
		if ps and keyword_set(filename) then popen,root + 'spec_' + filename
		if ps and ~keyword_set(filename) then popen,root + 'spec'

		if KEYWORD_SET(zlim_spec) then zlim,['x_','y_','z_']+'SPEC',zlim_spec[0],zlim_spec[1],1

		if sz eq 1 then tplot,'x_SPEC',window=-1
		if sz eq 2 then tplot,['x_SPEC','y_SPEC'],window=-1
		if sz eq 3 then tplot,['x_SPEC','y_SPEC','z_SPEC'],window=-1


		if ps then pclose

	endelse

	tplot_options,'title',''

	print,'SPECTROGRAM TPLOT VARIABLES PLOTTED:'
	print,'...' + wave + '_SPEC'
	if sz eq 2 then print,'...slide_specy'
	if sz eq 3 then print,'...slide_specz'


endif


;--------------------
;WAVELET PLOTS
;--------------------


if plot_wavelet then begin




	wavelet_to_tplot,wave2.x,wave2.y[*,0],new_name='waveletx',dscale=ds
	if sz ge 2 then wavelet_to_tplot,wave2.x,wave2.y[*,1],new_name='wavelety',dscale=ds
	if sz eq 3 then wavelet_to_tplot,wave2.x,wave2.y[*,2],new_name='waveletz',dscale=ds

	options,['waveletx','wavelety','waveletz'],'ytitle','Freq!C[Hz]'




;	if keyword_set(cg) then begin
;		xs=floor(ss[0]/2.5)
;		ys=floor(ss[1]/2.)
;		xp=floor(ss[0]-xs)
;		yp=floor(ss[1])/20
;
;		cgwindow,wxsize=xs,wysize=ys,wxpos=xp,wypos=yp
;
;		if sz eq 1 then cgwindow,'tplot','waveletx',window=-1,/addcmd
;		if sz eq 2 then cgwindow,'tplot',['waveletx','wavelety'],window=-1,/addcmd
;		if sz eq 3 then cgwindow,'tplot',['waveletx','wavelety','waveletz'],window=-1,/addcmd
;
;	endif else begin

		if ps then !p.font = 0
		if ps and keyword_set(filename) then popen,root + 'wavelet_' + filename
		if ps and ~keyword_set(filename) then popen,root + 'wavelet'

		tplot_options,'title','Wavelets(s)'

		if keyword_set(cg) then begin
		  if sz eq 1 then cgwindow,'tplot','waveletx',wxsize=floor(ss[0]/2.5),wysize=floor(ss[1]/2.3),wxpos=floor(ss[0]-floor(ss[0]/2.5)),wypos=floor(ss[1])/25,wtitle='spec'
		  if sz eq 2 then cgwindow,'tplot',['waveletx','wavelety'],wxsize=floor(ss[0]/2.5),wysize=floor(ss[1]/2.3),wxpos=floor(ss[0]-floor(ss[0]/2.5)),wypos=floor(ss[1])/25,wtitle='spec'
		  if sz eq 3 then cgwindow,'tplot',['waveletx','wavelety','waveletz'],wxsize=floor(ss[0]/2.5),wysize=floor(ss[1]/2.3),wxpos=floor(ss[0]-floor(ss[0]/2.5)),wypos=floor(ss[1])/25,wtitle='spec'
		endif else begin
		  if ~ps then wset,win_wavelet
		  if sz eq 1 then tplot,'waveletx'
		  if sz eq 2 then tplot,['waveletx','wavelety']
		  if sz eq 3 then tplot,['waveletx','wavelety','waveletz']
		endelse
		if keyword_set(hline) then timebar,hline,/databar,varname=hvarname

		if ps then pclose

;	endelse

	tplot_options,'title',''

	print,'WAVELET TPLOT VARIABLES PLOTTED:'
	print,'...waveletx'
	if sz ge 2 then print,'...wavelety'
	if sz eq 3 then print,'...waveletz'


endif


;-----------------------------------------------------
;HODOGRAMS
;-----------------------------------------------------

if plot_hod then begin


	get_data,wave,data=goob,dlimits=dgoob,limits=lgoob
	store_data,'wave_tmp',data=goob,dlimits=dgoob,limits=lgoob

	maxv = max(wave2.y)



	;make sure the vector to be overplotted (vec) is scaled to the max value of the hodogram
	;so that it plots nicely
	if keyword_set(vec1) then begin
		vec1tmp = double(vec1)
		vec1tmp = maxv*vec1tmp/sqrt(vec1tmp[0]^2 + vec1tmp[1]^2 + vec1tmp[2]^2)
	endif
	if keyword_set(vec2) then begin
		vec2tmp = double(vec2)
		vec2tmp = maxv*vec2tmp/sqrt(vec2tmp[0]^2 + vec2tmp[1]^2 + vec2tmp[2]^2)
	endif
	if keyword_set(vec3) then begin
		vec3tmp = double(vec3)
		vec3tmp = maxv*vec3tmp/sqrt(vec3tmp[0]^2 + vec3tmp[1]^2 + vec3tmp[2]^2)
	endif
	if keyword_set(vec4) then begin
		vec4tmp = double(vec4)
		vec4tmp = maxv*vec4tmp/sqrt(vec4tmp[0]^2 + vec4tmp[1]^2 + vec4tmp[2]^2)
	endif



	;set default options
	if keyword_set(hodt) then mtitle=hodt else mtitle='Hodograms'
	xr = [-1*maxv,maxv]
	yr = xr


	;Set options with _extra structure
	if is_struct(extra_hod) then options,'wave_tmp',_extra=extra_hod
	if is_struct(extra_hod) and tag_exist(extra_hod,'yrange') then yr = extra_hod.yrange
	if is_struct(extra_hod) and tag_exist(extra_hod,'xrange') then xr = extra_hod.xrange






	if sz lt 3 then print,'NEED TO ENTER A [N,3] WAVE QUANTITY TO PLOT HODOGRAMS' else begin


;		if keyword_set(cg) then begin
;
;
;			!p.charsize = 1.6
;			xs=floor(ss[0]/4.5)
;			ys=floor(ss[1]/1.5)
;			xp=floor(ss[0]-xs)
;			yp=floor(ss[1])/20
;
;			cgwindow,wxsize=xs,wysize=ys,wxpos=xp,wypos=yp,wmulti=[0,0,3],waspect=2.8
;
;			cgwindow,'plot',wave2.y[*,0],wave2.y[*,1],/addcmd
;			if sz ge 2 then cgwindow,'plot',wave2.y[*,0],wave2.y[*,2],/addcmd
;			if sz eq 3 then cgwindow,'plot',wave2.y[*,1],wave2.y[*,2],/addcmd
;
;		endif else begin

			tplot_options,'title','Hodograms'


; 			if ~ps then wset,win_hod
			if ps then !p.font = 0
			if ps and keyword_set(filename) then popen,root + 'hod_' + filename
			if ps and ~keyword_set(filename) then popen,root + 'hod'

			if keyword_set(cg) then begin
			    cgwindow,'tplotxy','wave_tmp',wmulti=[0,1,3,0,0],$
				    _extra={multi:'1,3',versus:'xy',xrange:xr,yrange:xr},$
				    wxsize=ss[0]/7.5,wysize=ss[1]/1.7,wxpos=ss[0]/40.,wypos=ss[1]/15
			    cgwindow,'tplotxy','wave_tmp',_extra={versus:'xz',addpanel:1,xrange:xr,yrange:xr},/addcmd
			    cgwindow,'tplotxy','wave_tmp',_extra={versus:'yz',addpanel:1,xrange:xr,yrange:xr},/addcmd
			endif else begin
			    if ~ps then wset,win_hod
					if ps then multfac = 5. else multfac = 2.
			    tplotxy,'wave_tmp',multi='1 3',versus='xy',xrange=xr,yrange=yr,mtitle=mtitle
					if keyword_set(vec1tmp) then oplot,[0,vec1tmp[0]],[0,vec1tmp[1]],color=50,thick=multfac*2.5
					if keyword_set(vec2tmp) then oplot,[0,vec2tmp[0]],[0,vec2tmp[1]],color=120,thick=multfac*2
					if keyword_set(vec3tmp) then oplot,[0,vec3tmp[0]],[0,vec3tmp[1]],color=200,thick=multfac*1.5
					if keyword_set(vec4tmp) then oplot,[0,vec4tmp[0]],[0,vec4tmp[1]],color=250,thick=multfac*1.
			    tplotxy,'wave_tmp',/addpanel,versus='xz',xrange=xr,yrange=xr
					if keyword_set(vec1tmp) then oplot,[0,vec1tmp[0]],[0,vec1tmp[2]],color=50,thick=multfac*2.5
					if keyword_set(vec2tmp) then oplot,[0,vec2tmp[0]],[0,vec2tmp[2]],color=120,thick=multfac*2
					if keyword_set(vec3tmp) then oplot,[0,vec3tmp[0]],[0,vec3tmp[2]],color=200,thick=multfac*1.5
					if keyword_set(vec4tmp) then oplot,[0,vec4tmp[0]],[0,vec4tmp[2]],color=250,thick=multfac*1.
			    tplotxy,'wave_tmp',/addpanel,versus='yz',xrange=xr,yrange=xr
					if keyword_set(vec1tmp) then oplot,[0,vec1tmp[1]],[0,vec1tmp[2]],color=50,thick=multfac*2.5
					if keyword_set(vec2tmp) then oplot,[0,vec2tmp[1]],[0,vec2tmp[2]],color=120,thick=multfac*2
					if keyword_set(vec3tmp) then oplot,[0,vec3tmp[1]],[0,vec3tmp[2]],color=200,thick=multfac*1.5
					if keyword_set(vec4tmp) then oplot,[0,vec4tmp[1]],[0,vec4tmp[2]],color=250,thick=multfac*1.
			endelse
; 			tplotxy,'wave_tmp',versus='xz',/addpanel,xrange=xr,yrange=xr
; 			tplotxy,'wave_tmp',versus='yz',/addpanel,xrange=xr,yrange=xr



			if ps then pclose

;		endelse

	endelse


	tplot_options,'title',''

	print,'HODOGRAM TPLOT VARIABLES PLOTTED:'
	if sz ge 3 then begin
		print,'...fieldx'
		print,'...fieldy'
		print,'...fieldz'
	endif

	if ~keyword_set(cg) then store_data,'wave_tmp',/delete

endif


;------------------------------------------
;PLOT OTHER QUANITIY
;------------------------------------------





if plot_other then begin


	;Set up default options
	tplot_options,'title','Additional quantity'


	;Overwrite options with _extra structure
	if is_struct(extra_other) then options,var_other,_extra=extra_other


	if keyword_set(cg) then begin
		xs=ss[0]/2.5
		ys=ss[1]/2.
		xp=0
		yp=ss[1]-1.7*ss[1]/3.

		cgwindow,wxsize=xs,wysize=ys,wxpos=xp,wypos=yp
		cgwindow,'tplot',var_other,window=-1,/addcmd
; 		stop
	endif else begin

		;figure out .ps stuff
		if ~ps and ~cg then wset,win_other
		if ps then !p.font = 0
		if ps and keyword_set(filename) then popen,root + 'other_' + filename
		if ps and ~keyword_set(filename) then popen,root + 'other'

		tplot,var_other,window=-1



		if ps then pclose


	endelse

		tplot_options,'title',''
		print,'OTHER TPLOT VARIABLES PLOTTED:'

; stop
endif


;------------------------------------------------------------------------------------









; stop

if del_wf and keyword_set(win_wf) then wdelete,win_wf
if del_hod and keyword_set(win_hod) then wdelete,win_hod
if del_spec and keyword_set(win_spec) then wdelete,win_spec
if del_wavelet and keyword_set(win_wavelet) then wdelete,win_wavelet
if del_psd and keyword_set(win_psd) then wdelete,win_psd

; stop

;Delete all the variables temporarily created unless we are using Dfanning's
;cgwindow method. Keeping the variables around then allows you to resize them.

if ~keyword_set(nodelete) then begin
	if ~keyword_set(cg) then begin
; 		stop
		store_data,['fieldx','fieldy','fieldz','psd_x','psd_y','psd_z','x_SPEC','y_SPEC',$
			'z_SPEC','waveletx','wavelety','waveletz','wave_tmp'],/delete
	endif ;else store_data,['waveletx','wavelety','waveletz'],/delete
endif




tplot_options,'title',''

;restore defaults
!p = pinit


end







;	dt = (max(wave2.x)-min(wave2.x))/n_elements(wave2.x)   ;time b/t samples
;
;	wavex = WV_CWT(wave2.y[*,0],'Morlet',6,scale=period,/double,dscale=ds,/pad)
;	if sz gt 1 then wavey = WV_CWT(wave2.y[*,1],'Morlet',6,scale=period,/double,dscale=ds,/pad)
;	if sz gt 2 then wavez = WV_CWT(wave2.y[*,2],'Morlet',6,scale=period,/double,dscale=ds,/pad)
;
;	period2 = period*dt
;	freqs2 = 1/period2
;	powerx = abs(wavex^2)
;	if sz gt 1 then powery = abs(wavey^2)
;	if sz gt 2 then powerz = abs(wavez^2)
;
;	ntime = n_elements(wave2.x)
;	nscale = n_elements(period)
;
;	if samescale eq 'no' then begin
;		csx = bytscl(alog10(powerx),min=min(alog10(powerx)/fac,/nan),max=max(alog10(powerx)*fac,/nan))
;		if sz gt 1 then csy = bytscl(alog10(powery),min=min(alog10(powery)/fac,/nan),max=max(alog10(powery)*fac,/nan))
;		if sz gt 2 then csz = bytscl(alog10(powerz),min=min(alog10(powerz)/fac,/nan),max=max(alog10(powerz)*fac,/nan))
;	endif else begin
;		if samescale eq 'x' then whichpow = powerx
;		if samescale eq 'y' then whichpow = powery
;		if samescale eq 'z' then whichpow = powerz
;		csx = bytscl(alog10(powerx),min=min(alog10(whichpow)/fac,/nan),max=max(alog10(whichpow)*fac,/nan))
;		if sz gt 1 then csy = bytscl(alog10(powery),min=min(alog10(whichpow)/fac,/nan),max=max(alog10(whichpow)*fac,/nan))
;		if sz gt 2 then csz = bytscl(alog10(powerz),min=min(alog10(whichpow)/fac,/nan),max=max(alog10(whichpow)*fac,/nan))
;	endelse
;
;
;	lim1   = {YLOG:1,ZLOG:0,YSTYLE:1,PANEL_SIZE:2,XMINOR:5,XTICKLEN:0.04,$
;			  ZSTYLE:1,spec:1}
;	store_data,'waveletx',data={X:wave2.x,Y:csx,V:freqs2/1000.},LIMIT=lim1
;	if sz gt 1 then store_data,'wavelety',data={X:wave2.x,Y:csy,V:freqs2/1000.},LIMIT=lim1
;	if sz gt 2 then store_data,'waveletz',data={X:wave2.x,Y:csz,V:freqs2/1000.},LIMIT=lim1
;














;OBSOLETE - OLD slide_spec.pro CODE

;	if ~keyword_set(step) then step = 0.1
;	if ~keyword_set(step_overlap) then step_overlap = 0.8

;	sldx = slide_spec(wave2.x,wave2.y[*,0],step,step_overlap,/db,iwindow=1,time_index=time_index,freq_bins=freq_binsx,/zero_pad)
;	if sz gt 1 then sldy = slide_spec(wave2.x,wave2.y[*,1],step,step_overlap,/db,iwindow=1,freq_bins=freq_binsy,/zero_pad)
;	if sz gt 2 then sldz = slide_spec(wave2.x,wave2.y[*,2],step,step_overlap,/db,iwindow=1,freq_bins=freq_binsz,/zero_pad)


;	if ~keyword_set(mincomp) then mincomp = 1.
;	if ~keyword_set(maxcomp) then maxcomp = 1.

;	csx=bytscl(sldx,min=min(sldx*mincomp,/nan),max=max(sldx/maxcomp,/nan))
;	if sz gt 1 then csy=bytscl(sldy,min=min(sldy*mincomp,/nan),max=max(sldy/maxcomp,/nan))
;	if sz gt 2 then csz=bytscl(sldz,min=min(sldz*mincomp,/nan),max=max(sldz/maxcomp,/nan))


;	;PUT INTO TPLOT VARIABLES
;	;lim1   = {YLOG:0,ZLOG:0,YSTYLE:1,PANEL_SIZE:2,XMINOR:5,XTICKLEN:0.04,$
;	;		  zstyle:1,spec:1,zrange:[0,255]}
;	store_data,'slide_specx',data={x:wave2.x[time_index],y:csx,v:freq_binsx};,limit=lim1
;	if sz gt 1 then store_data,'slide_specy',data={x:wave2.x[time_index],y:csy,v:freq_binsx};,limit=lim1
;	if sz gt 2 then store_data,'slide_specz',data={x:wave2.x[time_index],y:csz,v:freq_binsx};,limit=lim1
;
;	wset,3


;	if sz eq 1 then tplot,'slide_specx',window=-1
;	if sz eq 2 then tplot,['slide_specx','slide_specy'],window=-1
;	if sz eq 3 then tplot,['slide_specx','slide_specy','slide_specz'],window=-1


;	if keyword_set(vline) then time_bar,vline,varname=vvarname
;	if keyword_set(hline) then timebar,hline,/databar,varname=hvarname

;	print,'SPECTROGRAM TPLOT VARIABLES PLOTTED:'
;	print,'...slide_specx'
;	if sz eq 2 then print,'...slide_specy'
;	if sz eq 3 then print,'...slide_specz'
