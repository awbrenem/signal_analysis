;+
;*****************************************************************************************
;
;  PROCEDURE :   poynting_flux_vs_freq.pro
;  PURPOSE  :   Calculates the freq space Poynting flux from electric and
;               magnetic field inputs (in whatever coord). Returns the Poynting
;               flux spectra and the total output power (summed over all freqs and
;               multiplied by bandwidth) as function of time
;               Output flux is in µW m^(-2).
;
;  CALLED BY:
;               NA
;
;  CALLS:
;               tplot_struct_format_test.pro
;               sample_rate.pro
;
;  INPUT:
;               Ewvar           :  IDL TPLOT structure containing the E-Field data with
;                                  the format:
;                                  X  :  [N]-Element array of Unix times
;                                  Y  :  [N,3]-Element array of E-field vectors (mV/m)
;               Bwvar           :  IDL TPLOT structure containing the B-Field data with
;                                  the format:
;                                  X  :  [N]-Element array of Unix times
;                                  Y  :  [N,3]-Element array of B-field vectors (nT)
;                          *******************************************
;                          ** timestamps for the Bw must match EW **
;                          *******************************************
;
;  EXAMPLES:
;
;
;  KEYWORDS:
;               FLOW          :  Scalar [float/double] defining the low frequency [Hz]
;                                  cutoff to use for the Poynting flux analysis
;                                  [Default = 0.0]
;               NFFT          :  Scalar [long] defining the # of frequency bins in each
;                                  FFT
;                                  [Default = 128]
;               NSHFT         :  Scalar [long] defining the # of points to shift between
;                                  each FFT
;                                  [Default = 32]
;               nophase       ;  Set if you don't trust the phase calibration b/t E and B.
;                                This will return the maximum possible Poynting flux by taking
;                                abs(Ew) and abs(Bw) (= sqrt(Re^2 + Im^2)) and then adding
;                                the Poynting flux components instead of subtracting them.
;
;
;   CHANGED:  1)  NA [MM/DD/YYYY   v1.0.0]
;
;   NOTES:
;
;
;   CREATED:  2017-05-23
;   CREATED BY:  Aaron W Breneman (modified from Lynn Wilson's code calc_poynting_flux_freq.pro)
;*****************************************************************************************
;-

pro poynting_flux_vs_freq,Ewvar,Bwvar,$
  FLOW=flow,NFFT=nfft,NSHFT=nshft,erg=erg,nophase=nophase



  ;;----------------------------------------------------------------------------------------
  ;; => Define some parameters
  ;;----------------------------------------------------------------------------------------

  muo         = (4d0*!DPI)*1d-7            ;; Permeability of free space (N/A^2 or H/m)
  sfac        = 1d-3*1d-9*1d6/muo[0]       ;; mV->V, nT->T, W->µW, divide by µ_o


  IF (N_ELEMENTS(flow) NE 1) THEN  flow  = 0d0  ELSE flow  = DOUBLE(flow[0])
  IF (N_ELEMENTS(nfft) NE 1) THEN  nfft  = 128L ELSE nfft  = LONG(nfft)
  IF (N_ELEMENTS(nshft) NE 1) THEN nshft = 32L  ELSE nshft = LONG(nshft[0])


  get_data,Bwvar,data=Bw
  get_data,Ewvar,data=Ew

  times = Bw.X
  evecs = Ew.Y
  bvecs = Bw.Y
  ntimes = N_ELEMENTS(times)


  srt = DOUBLE(ROUND(sample_rate(times,GAP_THRESH=2d0,/AVE)))  ;sample rate


  ;;----------------------------------------------------------------------------------------
  ;; => Define some dummy parameters
  ;;----------------------------------------------------------------------------------------
  ntimes_a      = ntimes - 1L                    ;; # of vectors in each interval
  n_dft       = ntimes_a - nfft             ;; # of timestamps per interval
  ;;  Create a Hanning window
  win         = HANNING(nfft,/DOUBLE)*2d0
  winfull     = HANNING(ntimes,/DOUBLE)*2d0
  ;;  Define the bandwidth [Hz] of each FFT frequency bin
  bandw       = srt/nfft
  freq        = DINDGEN(nfft/2d0)*bandw   ;; FFT frequency bin values [Hz]
  ind_sf      = WHERE(freq GE flow,gdsf)




  ;;----------------------------------------------------------------------------------------
  ;; => Convert to Frequency Domain
  ;;----------------------------------------------------------------------------------------


  ;; Define # of timestamps for FFTs
  n_ts = ntimes/nshft + 1L
  tt = DBLARR(n_ts)  ;dummy timestamps

  ;; Define dummy arrays for the complex FFTs of the two fields
  Ew_fft = DCOMPLEXARR(n_ts,nfft,3L)  ;;  [K,W,3]-Element array
  Bw_fft = DCOMPLEXARR(n_ts,nfft,3L)


  i = 0L
  for j=0L, n_dft-1L, nshft do begin
    upj = j + nfft - 1L
    for k=0L, 2L do Ew_fft[i,*,k] = FFT(evecs[j[0]:upj[0],k]*win)
    for k=0L, 2L do Bw_fft[i,*,k] = FFT(bvecs[j[0]:upj[0],k]*win)
    i++
  endfor



  ;; Keep only relevant values
  ind = (i - 1L)

  ;; Remove extra spectra on end
  Ew_fft = Ew_fft[0L:(ind - 1L),*,*]
  Bw_fft = Bw_fft[0L:(ind - 1L),*,*]
  ;; Define timestamps
  tt = times[0] + (DINDGEN(ind)*nshft + nfft/2d0)/srt



  ;;----------------------------------------------------------------------------------------
  ;; => Calculate Poynting Flux [Frequency Domain]
  ;;----------------------------------------------------------------------------------------



  ;Poynting flux of spectral FFT ;[µW m^(-2)/Hz]

  if ~KEYWORD_SET(nophase) then begin
    sx_f = DOUBLE(Ew_fft[*,*,1]*CONJ(Bw_fft[*,*,2]) - Ew_fft[*,*,2]*CONJ(Bw_fft[*,*,1]))*sfac[0]/2d0
    sy_f = DOUBLE(Ew_fft[*,*,2]*CONJ(Bw_fft[*,*,0]) - Ew_fft[*,*,0]*CONJ(Bw_fft[*,*,2]))*sfac[0]/2d0
    sz_f = DOUBLE(Ew_fft[*,*,0]*CONJ(Bw_fft[*,*,1]) - Ew_fft[*,*,1]*CONJ(Bw_fft[*,*,0]))*sfac[0]/2d0
  endif

  if KEYWORD_SET(nophase) then begin
    Ew_fftR = abs(Ew_fft)
    Bw_fftR = abs(Bw_fft)

    sx_f = DOUBLE(Ew_fftR[*,*,1]*Bw_fftR[*,*,2] + Ew_fftR[*,*,2]*Bw_fftR[*,*,1])*sfac[0]/2d0
    sy_f = DOUBLE(Ew_fftR[*,*,2]*Bw_fftR[*,*,0] + Ew_fftR[*,*,0]*Bw_fftR[*,*,2])*sfac[0]/2d0
    sz_f = DOUBLE(Ew_fftR[*,*,0]*Bw_fftR[*,*,1] + Ew_fftR[*,*,1]*Bw_fftR[*,*,0])*sfac[0]/2d0
  endif

  s_f  = [[[sx_f]],[[sy_f]],[[sz_f]]]
  power_sf = SQRT(TOTAL(TOTAL(s_f[*,ind_sf,*]^2,3L,/NAN),2L))*srt/nfft  ;[µW m^(-2)]




  if KEYWORD_SET(erg) then begin

    ;10d6 erg/s = 1 W
    ;(uW/m2)/1d6 = W/m2
    ;(W/m2)/100^2 = W/cm2
    ;(W/cm2)*10d6 = erg/s/cm2

    sx_f *= 0.001
    sy_f *= 0.001
    sz_f *= 0.001
    power_sf *= 0.001

  endif


  ;;  remove negative or zero values from POWER_SF
  bad = WHERE(power_sf LE 0,bd)
  IF (bd GT 0) THEN power_sf[bad] = !VALUES.D_NAN


  if ~KEYWORD_SET(erg) then store_data,'power',tt,power_sf,dlim={yrange:[0,30],ytitle:'Power(integrated)!CuW/m2'}
  if KEYWORD_SET(erg)  then store_data,'power',tt,power_sf,dlim={yrange:[0,30],ytitle:'Power(integrated)!Cerg/s/cm2'}
  options,'power','psym',-4


  ;;----------------------------------------------------------------------------------------
  ;; => Return to user
  ;;----------------------------------------------------------------------------------------

  if ~KEYWORD_SET(erg) then begin
    store_data,'Sx',tt,sx_f[*,ind_sf],freq,$
    dlim={ylog:1,yrange:[1d-5,500],spec:1,zlog:1,zrange:[1d-5,1],ytitle:'PFlux-Xhat!CuW/m2/Hz'}
    store_data,'Sy',tt,sy_f[*,ind_sf],freq,$
    dlim={ylog:1,yrange:[1d-5,500],spec:1,zlog:1,zrange:[1d-5,1],ytitle:'PFlux-Yhat!CuW/m2/Hz'}
    store_data,'Sz',tt,sz_f[*,ind_sf],freq,$
    dlim={ylog:1,yrange:[1d-5,500],spec:1,zlog:1,zrange:[1d-5,1],ytitle:'PFlux-Zhat!CuW/m2/Hz'}
  endif else begin
    store_data,'Sx',tt,sx_f[*,ind_sf],freq,$
    dlim={ylog:1,yrange:[1d-5,500],spec:1,zlog:1,zrange:[1d-5,1],ytitle:'PFlux-Xhat!Cerg/s/cm2/Hz'}
    store_data,'Sy',tt,sy_f[*,ind_sf],freq,$
    dlim={ylog:1,yrange:[1d-5,500],spec:1,zlog:1,zrange:[1d-5,1],ytitle:'PFlux-Yhat!Cerg/s/cm2/Hz'}
    store_data,'Sz',tt,sz_f[*,ind_sf],freq,$
    dlim={ylog:1,yrange:[1d-5,500],spec:1,zlog:1,zrange:[1d-5,1],ytitle:'PFlux-Zhat!Cerg/s/cm2/Hz'}
  endelse

  store_data,'Ex',tt,Ew_fft[*,ind_sf,0],freq,$
  dlim={ylog:1,yrange:[1d0,1d3],zrange:[1d0,1d3],spec:1,zlog:1,ytitle:'Ew spec-Xhat!CmV/m/sqrt(Hz)'}
  store_data,'Ey',tt,Ew_fft[*,ind_sf,1],freq,$
  dlim={ylog:1,yrange:[1d0,1d3],zrange:[1d0,1d3],spec:1,zlog:1,ytitle:'Ew spec-Yhat!CmV/m/sqrt(Hz)'}
  store_data,'Ez',tt,Ew_fft[*,ind_sf,2],freq,$
  dlim={ylog:1,yrange:[1d0,1d3],zrange:[1d0,1d3],spec:1,zlog:1,ytitle:'Ew spec-Zhat!CmV/m/sqrt(Hz)'}


  store_data,'Bx',tt,Bw_fft[*,ind_sf,0],freq,$
  dlim={ylog:1,yrange:[1d0,1d3],zrange:[1d-5,1d3],spec:1,zlog:1,ytitle:'Bw spec-Xhat!CnT/sqrt(Hz)'}
  store_data,'By',tt,Bw_fft[*,ind_sf,1],freq,$
  dlim={ylog:1,yrange:[1d0,1d3],zrange:[1d-5,1d3],spec:1,zlog:1,ytitle:'Bw spec-Yhat!CnT/sqrt(Hz)'}
  store_data,'Bz',tt,Bw_fft[*,ind_sf,2],freq,$
  dlim={ylog:1,yrange:[1d0,1d3],zrange:[1d-5,1d3],spec:1,zlog:1,ytitle:'Bw spec-Zhat!CnT/sqrt(Hz)'}


END
