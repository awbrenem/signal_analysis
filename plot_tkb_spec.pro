;+
;*****************************************************************************************
;
;  FUNCTION  : plot_tkb_spec.pro  
;  PURPOSE : Plots the results from get_tkb_spec.pro including theta_kb and polarization
;
;  CALLED BY:  
;               
;
;  CALLS:
;               
;
;  REQUIRES:    
;               
;
;  INPUT:		tkb -> array returned from get_tkb_spec.pro
;				wf -> FA waveform to be plotted
;				times -> times for the FA waveform (make sure first time is zero)
;				title - to be appended to plots
;				minpow -> minimum power to plot (dB)
;				maxepol -> max value for polarization ratio plot.
;				fn -> title (along with directory structure) for filename (include .ps)
;				ps_only - don't plot to window
;				foplot - freqs to overplot in kHz
;				vline  - overplot vertical line
;				_wfE=wfE,_modE=modE,_tkbE=tkbE,_polE=polE,_specE=specE
;
; 
;  EXAMPLES:    
;               
;
;  KEYWORDS:    
;               
;
;   CHANGED:  1)  NA [MM/DD/YYYY   v1.0.0]
;
;   NOTES:    
;
;			  This program also outputs a .sav file for each TDS captures that contains
;			  the theta_kb and amplitude spec info.
;               
;			  See also theta_kb_plots.pro for plotting other quantities such as the res
;			  energies and theta_kb_hists.pro for histograms. Both use the .sav files 
;			  created from this program.
;
;   CREATED:  07/08/2010
;   CREATED BY:  Aaron W. Breneman
;    LAST MODIFIED:  09/29/2011   v1.0.0
;			  Added minpow, vline keywords
;    MODIFIED BY: Aaron W. Breneman
;
;*****************************************************************************************
;-

;Possible methods to test goodness of epol
;1. consistency of zero-crossings (like my wave classification method)
;2. Angle of semimajor axis to x-axis and semiminor axis to y-axis. 
;		-are the angles small
;		-are the two angles at 90 degrees? 


pro plot_tkb_spec,tkb,wf,times,title=title,fn=fn,$
	ps_only=ps_only,minpow=minpow,maxepol=maxepol,foplot=foplot,vline=vline,$
	_wfE=wfE,_modE=modE,_tkbE=tkbE,_polE=polE,_specE=specE



tkb2 = tkb

if not keyword_set(fn) then fn = '~/Desktop/code/Aaron/tkb_spec.ps'
if not keyword_set(maxepol) then maxepol=10

window_ps = {xsize:7.5,ysize:9.5,xoffset:0.425,yoffset:1.68,inches:1,portrait:1,landscape:0}


efieldx = reform(wf[*,0])
efieldy = reform(wf[*,1])
efieldz = reform(wf[*,2])

epol = fix(tkb.epol)
polarization = fix(tkb2.polarization)
theta_kb = fix(tkb2.theta_kb)
time_bins = tkb2.time_bins
freq_bins= tkb2.freq_bins
ex_spec = tkb2.ex_spec
time_bins = tkb2.time_bins
freq_bins = tkb2.freq_bins

device,decomposed=0
!p.charsize=2.
xmargin_left=6 
xmargin_right=11


maxrange = max(efieldx) > max(efieldy) > max(efieldz)

if ~keyword_set(title) then title2 = '' else title2 = title


;--------------------------------------
;Adjust things to pretty up the plots
;--------------------------------------

;set the color values for polarization for color table 39
;90 = black for theta_kb
;255 = white for polarization


goo = where(epol gt maxepol)
if goo[0] ne -1 then epol[goo] = maxepol-1


goo = where(tkb2.polarization eq 2)
if goo[0] ne -1 then polarization[goo] = 2

if keyword_set(minpow) then begin
	goo = where(ex_spec le minpow)
	if goo[0] ne -1 then begin
		epol[goo] = 90
		theta_kb[goo] = 90
		polarization[goo] = 255
	endif
endif


;Check for NaN values. 
goo = where((finite(theta_kb) eq 0) or (finite(polarization) eq 0))
if goo[0] ne -1 then theta_kb[goo] = 90
if goo[0] ne -1 then epol[goo] = 90
if goo[0] ne -1 then polarization[goo] = 255

goo1 = where(polarization eq 0)
if goo1[0] ne -1 then polarization[goo1] = 255 ;white
goo2 = where(polarization eq -1)
if goo2[0] ne -1 then polarization[goo2] = 75   ;LH = blue
goo3 = where(polarization eq 1)
if goo3[0] ne -1 then polarization[goo3] = 250  ;RH = red

goo = where(theta_kb eq 0.)
if goo[0] ne -1 then theta_kb[goo] = 90 
if goo[0] ne -1 then epol[goo] = 90 
if goo[0] ne -1 then polarization[goo] = 255



loadct,39   ;last element is white


;set the range from 0 to the max value of theta_kb to stretch to 0-255
color_shades = bytscl(theta_kb,min=0,max=90,/nan)
color_shades2 = bytscl(polarization,min=0,max=255,/nan)
color_shades3 = bytscl(epol,min=0,max=maxepol,/nan)



if ~keyword_set(ps_only) then begin
	!p.multi = [0,0,6,0,0]
	plot,times,efieldx,xrange=[0,max(times)],xstyle=1,xmargin=[xmargin_left,xmargin_right],$
		yrange=[-maxrange,maxrange],title='Ex_stix (blue), Ey (green), Ez (orange) (mV/m)',$
		xtitle='time (msec)',/nodata,_extra=wfE
	oplot,times,efieldx,color=50
	oplot,times,efieldy,color=130
	oplot,times,efieldz,color=220
	
	
	amp = sqrt(efieldx^2 + efieldy^2 + efieldz[*,0]^2)
	plot,times,amp,xrange=[0,max(times)],xstyle=1,xmargin=[xmargin_left,xmargin_right-1],yrange=[0,$
		max(amp,/nan)],title='Total amplitude mod (all 3 components) (mV/m)',xtitle='time (msec)',_extra=ampE
	plot_spec,theta_kb,time_bins,freq_bins,/no_interp,ytitle='freq (kHz)',xtitle='time (msec)',$
		title='theta_kb from cold plasma dispersion',zrange=[0,90],ztitle='theta_kb',xstyle=1,$
		shades=color_shades,_extra=tkbE,hline=foplot,vline=vline
	plot_spec,epol,time_bins,freq_bins,/no_interp,ytitle='freq (kHz)',xtitle='time (msec)',$
		title='epol (Ex/Ey Stix)',zrange=[0,maxepol],ztitle='epol',xstyle=1,$
		shades=color_shades3,_extra=tkbE,hline=foplot,vline=vline

	plot_spec,polarization,time_bins,freq_bins,/no_interp,ytitle='freq (kHz)',xtitle='time (msec)',$
		title='Handedness (Red=right handed)',ztitle='Polarization',xstyle=1,$
		_extra=polE,hline=foplot,vline=vline,shades=color_shades2
	plot_spec,ex_spec,time_bins,freq_bins,/no_interp,ytitle='freq (kHz)',xtitle='time (msec)',/db,$
		xstyle=1,title=title2,_extra=specE,hline=foplot,vline=vline


endif

device,decomposed=0
loadct,39

;plot to .ps  (Note that I'm not plotting the color bar here b/c this doesn't allow me to
;load the file into Adobe Illustrator.
set_plot,'ps'
!p.font = 0
device,filename=fn,decomposed=0,encapsulated=1,font_size=8,/helvetica,/color

!p.multi = [0,0,6,0,0]

plot,times,efieldx,xrange=[0,max(times)],xstyle=1,xmargin=[xmargin_left,xmargin_right],$
	yrange=[-maxrange,maxrange],title='Ex_stix (blue), Ey (green), Ez (orange) (mV/m)',$
	xtitle='time (msec)',/nodata,_extra=wfE
oplot,times,efieldx,color=50
oplot,times,efieldy,color=130
oplot,times,efieldz,color=220

amp = sqrt(efieldx^2 + efieldy^2 + efieldz[*,0]^2)
plot,times,amp,xrange=[0,max(times)],xstyle=1,xmargin=[xmargin_left,xmargin_right-1],yrange=[0,$
	max(amp,/nan)],title='Total amplitude mod (all 3 components) (mV/m)',xtitle='time (msec)',_extra=empE
plot_spec,theta_kb,time_bins,freq_bins,/no_interp,ytitle='freq (kHz)',xtitle='time (msec)',$
	title='theta_kb from cold plasma dispersion',zrange=[0,90],ztitle='theta_kb',xstyle=1,$
	shades=color_shades,_extra=tkbE,hline=foplot,vline=vline,/nocb
plot_spec,epol,time_bins,freq_bins,/no_interp,ytitle='freq (kHz)',xtitle='time (msec)',$
	title='epol (Ex/Ey Stix)',zrange=[0,maxepol],ztitle='epol',xstyle=1,$
	shades=color_shades3,_extra=tkbE,hline=foplot,vline=vline,/nocb
plot_spec,polarization,time_bins,freq_bins,/no_interp,ytitle='freq (kHz)',xtitle='time (msec)',$
	title='Polarization (Red=right handed)',ztitle='Polarization',xstyle=1,_extra=polE,$
	hline=foplot,vline=vline,shades=color_shades2,/nocb
plot_spec,ex_spec,time_bins,freq_bins,/no_interp,ytitle='freq (kHz)',xtitle='time (msec)',/db,$
	xstyle=1,title=title2,_extra=specE,hline=foplot,vline=vline,/nocb


device,/close
set_plot,'x'

end









