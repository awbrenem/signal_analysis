;+
;*****************************************************************************************
;
;  FUNCTION :   get_tkb_spec
;  PURPOSE  :   Returns a structure with the following spectrograms from Ew component:
;					power spec in the same components as inputted
;					wavenormal spec
;					kmag spec
;					polarization spec
;
;
;  CALLED BY:
;
;
;  CALLS:		cold_dispersion_function
;				slide_spec
;				power_of_2
;				fft_bandpass
;				rotate_field_2_vec
;				fft_power_calc
;
;
;  REQUIRES:
;
;
;  INPUT:		time -> [n] array of time values (sec)
;				efield -> [n,3] array of electric field values in same coord as Bo
;				Bo -> [3] element vector of magnetic field to align waves to (nT)
;				dens -> plasma density (cm^-3)
;               pets -> variable for slide_spec.pro (I'm not calling this "step" b/c of the
;						stupid IDL rule where it doesn't look at the entire keyword name. Thus
;						it confuses "step" with "step_overlap" and thinks they are both called
;						"step" and returns a redundant keyword error.
;				step_overlap -> variable for slide_spec.pro
;				minamp -> wave amplitude in each f/t box must be at least this value
;						for the wave normal angle to be calculated. Defaults to 0 in whatever
;						units are input (e.g. mV/m). Note that it is probably best to leave
;						this at zero because you can always adjust the minimum amplitude plotted
;						in plot_tkb_spec.pro
;				suppressplot -> Don't plot hodogram for each f/t box
;				identifier -> name of the TDS capture.
;				loadsv -> Filename of the .sav version if it exists. Need the "identifier" to properly save
;				epolmax -> max value of polarization ratio (Ex/Ey) that is acceptable (defaults to 4)
;
;  EXAMPLES:    x = read_tds(filenames='test.txt')
;				Bo = [1,2,3]
;				dens = 5000.
;				r = get_tkb_spec(x.time,[[x.ex_cor],[x.ey_cor],[x.ez_cor]],Bo,dens)
;
;
;  KEYWORDS:
;
;
;   CHANGED:  1)  NA [MM/DD/YYYY   v1.0.0]
;
;   NOTES:
;
;
;   CREATED:  June 2010
;   CREATED BY:  Aaron W. Breneman
;    LAST MODIFIED:  MM/DD/YYYY   v1.0.0
;    MODIFIED BY: Aaron W. Breneman
;
;*****************************************************************************************
;-

function get_tkb_spec,time,efield,Bo,dens,pets=step,step_overlap=step_overlap,$
	minamp=minamp,suppressplot=suppressplot,H_plus=H_plus,He_plus=He_plus,O_plus=O_plus,$
	identifier=identifier,loadsv=loadsv,epolmax=epolmax

init_path
directory = !data.general + 'idlsave/theta_kb/'

if ~keyword_set(epolmax) then epolmax=4

if keyword_set(loadsv) then begin
;see if the .sav file has already been created. If so, then just load it and skip all the computation
filenm = loadsv
file_exists = file_info(filenm)
if file_exists.exists eq 1 then fileyes = 1. else fileyes = 0.
endif else fileyes = 0.

if fileyes eq 0. then begin

	if ~keyword_set(minamp) then minamp = 0.0
	if ~keyword_set(step) then step=0.1
	if ~keyword_set(step_overlap) then step_overlap=0.0

	efieldx = reform(efield[*,0])
	efieldy = reform(efield[*,1])
	efieldz = reform(efield[*,2])
	bmag = sqrt(Bo[0]^2 + Bo[1]^2 + Bo[2]^2)

	;create arrays
	ex_spec = slide_spec(time,efieldx,step,step_overlap,iwindow=1,time_index=time_bins,freq_bins=freq_bins,/db,/zero_pad)
	ey_spec = slide_spec(time,efieldy,step,step_overlap,iwindow=1,time_index=time_bins,freq_bins=freq_bins,/db,/zero_pad)
	ez_spec = slide_spec(time,efieldz,step,step_overlap,iwindow=1,time_index=time_bins,freq_bins=freq_bins,/db,/zero_pad)

	freq_bins = abs(freq_bins)
	epol = fltarr(n_elements(ex_spec[*,0]),n_elements(ex_spec[0,*]))
	epol_slice = fltarr(n_elements(time_bins))  ;polarization ratio for each time slice (all freqs)
	fmax = epol
	theta_kb = epol
	theta_res = epol
	kmag = epol
	polarization = epol  ;takes on values of -1,2,1 for LH, linear, RH polarization

	fstep = (max(freq_bins)-min(freq_bins))/(n_elements(freq_bins)-1)
	tstep = (max(time_bins)-min(time_bins))/(n_elements(time_bins)-1)

	;get window
	nelem = n_elements(efieldx)
	windowarr = replicate(0,nelem)
	my_windowf,nelem-1,2,windowarr

	srt = n_elements(time)/(max(time,/nan)-min(time,/nan))	;sample rate
	;srt = 4096/(max(time,/nan)-min(time,/nan))	;sample rate

;	PRINT,'************'
;	print,'HEY DIPSHIT, DONT FORGET TO ADJUST SRT'
;	PRINT,'************'

	window,1,xsize=600,ysize=600

	for i=0.,n_elements(time_bins)-2 do begin

		print,'i = ' + strtrim(i,2) + '/' + strtrim(n_elements(time_bins)-2,2)
	;	tc = time[i*tstep:(i+1)*(tstep-1)]
		tc = time[i*tstep:(i+1)*(tstep)]


		;pad the array to the correct size
	;	if nelem le 4096. then force = 2L^12 else force = 2L^14

		force = 2L^4
		kk = 5
		while force lt nelem do begin
			force = 2L^kk
			kk = kk+1
		endwhile

		ecx = power_of_2(efieldx[i*tstep:(i+1)*(tstep)],FORCE_N=force)
		ecy = power_of_2(efieldy[i*tstep:(i+1)*(tstep)],FORCE_N=force)
		ecz = power_of_2(efieldz[i*tstep:(i+1)*(tstep)],FORCE_N=force)

		epol_slice[i] = max(ecx[tc],/nan)/max(ecy[tc],/nan)


		;now bandpass each time chunk into the separate freq bins
		!p.multi = [0,4,5,0,0]

		for j=0.,n_elements(freq_bins)-2 do begin

			vbp = vector_bandpass([[ecx],[ecy],[ecz]],srt,1000*j*fstep,1000*(j+1)*fstep)

			wf = [[vbp[*,0]],[vbp[*,1]],[vbp[*,2]]]


			;align to field
			wf_FA = rotate_field_2_vec(Bo,wf,times=time)

			;Find the polarization of perp component w/r to magnetic field
			Etots = sqrt(wf_FA.vecFA[*,0]^2 + wf_FA.vecFA[*,1]^2)
			phi = acos(wf_FA.vecFA[*,0]/Etots)/!dtor
			goobar = where(wf_FA.vecFA[*,1] lt 0.)
			if goobar[0] ne -1 then phi[goobar] = 360.-phi[goobar]    ;dot product only goes from 0-180. Extend to full range

			dphi = deriv(phi)

			;Median filter to remove spike values
			filterv = 10
			dphi2 = median(dphi,filterv)
			;Remove fluctuating end elements
			dphi2[0:filterv-1] = 0.
			dphi2[n_elements(dphi2)-filterv:n_elements(dphi2)-1] = 0.

			if total(dphi2,/nan) ge 0. then polarization[i,j] = 1. else polarization[i,j] = -1.

			range = max(wf_FA.vecFA,/nan)

			;find epol and peak freq (use value from x_stix-direction)

			epol[i,j] = max(wf_FA.vecFA[*,0])/max(wf_FA.vecFA[*,1])

	;		if ((range lt minamp) or (epol[i,j] gt epolmax)) then epol[i,j] = 1
			if range lt minamp then epol[i,j] = 1

			;IDENTIFY THE LINEARLY POLARIZED WAVES
			if ((epol[i,j] gt epolmax) and (range ge minamp)) then polarization[i,j] = 2

			power_x = fft_power_calc(wf_FA.times,windowarr*wf_FA.vecfa[*,0])
			tmp = max(power_x.power_x,loc)
			fmax[i,j] = power_x.freq[loc]

			if ~keyword_set(suppressplot) then begin
				if ((range ge minamp) and (epol[i,j] le epolmax)) then plot,wf_FA.vecFA[*,0],wf_FA.vecFA[*,1],xrange=[-1*range,range],yrange=[-1*range,range],title='freq=' + strtrim(freq_bins[j],2)+'kHz'
	;	print,epol[i,j]
	;	print,range
	;	print,'______'
			endif

		endfor ;j
;		!p.multi = [0,0,1]
;		window,2,xsize=300,ysize=300
;		plot,freq_bins,ex_spec[i,*],xrange=[0,360],xtitle='freq kHz',ytitle='ex_spec (dB)'

	endfor ;i


 	for i=0,n_elements(time_bins)-1 do begin
		;st = cold_dispersion_function(epol=reform(epol[i,*]),freq=reform(fmax[i,*]),dens=dens,bo=bmag,H_plus=H_plus,He_plus=He_plus,O_plus=O_plus)
		st = cold_plasma_dispersion(epol=reform(epol[i,*]),freq=reform(fmax[i,*]),dens=dens,bo=bmag,H_plus=H_plus,He_plus=He_plus,O_plus=O_plus)
		theta_kb[i,*] = st.theta_kb
		theta_res[i,*] = st.resangle
		kmag[i,*] = st.kmag

		goob = where(reform(epol[i,*]) eq 1)
		if goob[0] ne -1 then begin
			theta_kb[i,goob] = !values.f_nan
			theta_res[i,goob] = !values.f_nan
			kmag[i,goob] = !values.f_nan
			polarization[i,goob] = !values.f_nan
		endif
	endfor

	struct = {ex_spec:ex_spec,ey_spec:ey_spec,ez_spec:ez_spec,epol_slice:epol_slice,theta_kb:theta_kb,epol:epol,$
		theta_res:theta_res,polarization:polarization,kmag:kmag,time_bins:time_bins,$
		freq_bins:freq_bins}
	if keyword_set(identifier) then save,filename=directory + identifier + '_theta_kb.sav'

endif else begin
	print,'.sav file found for ' + filenm
	restore,filenm
	struct = {ex_spec:ex_spec,ey_spec:ey_spec,ez_spec:ez_spec,epol_slice:epol_slice,theta_kb:theta_kb,epol:epol,$
	theta_res:theta_res,polarization:polarization,kmag:kmag,time_bins:time_bins,$
	freq_bins:freq_bins,epolmax:epolmax}
endelse



return,struct
end
