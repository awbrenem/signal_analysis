;+
; NAME: POWER_SPEC
;
; PURPOSE: Transform time series data into averaged frequency power
;          spectra using FFT's. Result is the two-sided Power Spectral
;          Density (PSD) given in (input units)^2/Hz.
;
; CALLING SEQUENCE: power_spec, data, n_ave=n_ave, npnts=npts,
;                   sample=sample, freq, ave_spec, overlap=overlap
;
; INPUTS:
;         DATA -     A one-dimensional array containing the values of
;                    the time series data for which an average power
;                    spectra is desired.
;
; OUTPUTS:
;         FREQ     - An array returned by the routine containing the
;                    frequency axis of the power spectra.
;
;         AVE_SPEC -An array returned by the routine containing the
;                   amplitudes of the power spectra. The units
;                   returned are (input units)^2/Hz. Thus EACH
;                   AMPLITUDE in the FFT is divided by the bandwidth
;                   in Hz of each frequency interval.
;
; KEYWORDS:
;         OVERLAP - Set this keyword to slide the FFT interval by
;                   one-half of an interval instead of averaging
;                   together each separate interval. For most data
;                   types this will yield a higher number of averages
;                   and less error per data point than straight
;                   sequential averaging (see below!). Note that N_AVE
;                   specifies how many sequential averages are taken
;                   without overlap. Thus N_AVE=4 with /OVERLAP yields
;                   a total of 7 averages, etc.
;
;         N_AVE -   The number of FFT's to average together within
;                   the data series provided. This specifies how many
;                   sequential segments in which to divide the
;                   data. Setting the /OVERLAP keyword increases the
;                   number of segments by 2*N_AVE-1.
;
;         NPNTS -   The number of time series points that each
;                   interval will be for the discrete FFT's.
;
;         SAMPLE -   The sample time (in secs) of each time series
;                    point. If unspecified, an arbitrary sample time
;                    of 1 is assigned to the data.
;
; EXAMPLE: You have 8192 continuous time series points. You want to
;          average together 8 FFT's to obtain a single average spectra
;          for this data. Set n_ave = 8, npnts=1024. Result is
;          a power_spectra obtained by taking a single FFT of each
;          1024 pt sequential interval in the data and averaging them
;          together. For the same parameters, if you set /overlap,
;          then the interval slides 1/2 and then averages, i.e. an FFT
;          is taken of the first 1024 points, then moves forward 512
;          points, another FFT, then slide again by 512 points. Thus
;          for num_ave = 8 and /overlap, you average together 8+7 = 15
;          FFT's for the interval of 8192 points.
;
;    A note about errors:
;          For n sequential (non-sliding) averages FFT's have a
;          statistical error of 1/n^0.5. For the same time interval,
;          if a sliding average is used of 1/2 a window overlap, n
;          goes up by a factor of 2 but the errors are no longer
;          completely independent -- so the new error per point is
;          1/(0.81*N)^0.5 -- but N=2n so there is less error per point
;          than in the non-sliding average case. See any of the
;          numerical recipes books by Press et al for a complete
;          description of this and other aspects of computational FFT's.
;
;    A note about the data:
;          Results may not be valid if the number of points per FFT
;          (npts) is not 2^j where j = any integer.
;          This routine assumes that the input data is continuous and
;          that there are no data gaps or bad time points.
;
; REVISION HISTORY
;  ver 1.0 11/15/96 G.T.Delory
;  ver 1.1 03/26/97 G.T.Delory
;  ver 1.2 04/10/97 G.T.Delory
;  ver 1.21 John Bonnell, UCBSSL, April 18, 2000
;    moved window calculation out of loop over data segments.
;-

pro power_spec, data, $
  N_AVE = n_ave, $
  NPTS = npts, $
  SAMPLE = sample, $
  freq, $
  ave_spec, $
  OVERLAP = overlap, $
  verbose=verbose, $
  analytic_sig=analytic_sig

  if not keyword_set(data) then begin
    print,'No data provided -- check calling sequence and try again'
    data = replicate(!values.f_nan,256)
  endif

  if not keyword_set(n_ave) then begin
    n_ave=4
    npts = n_elements(data)/n_ave
    if keyword_set(verbose) then $
    message, /info, $
    string( n_ave, format='("using default value for N_AVE (",I,").")')
  endif

  fft_out = fltarr(npts)
  ave_spec = fltarr(npts/2+1)

  ww = hanning( npts)
  W = npts*total( ww*ww)

  if not keyword_set(sample) then begin
    sample = 1.0
    if keyword_set(verbose) then $
    message, /info, $
    string( sample, format='("using default sampling interval for SAMPLE (",E10.3,").")')
  endif
  


  ; If overlap keyword is set, slide half of a segment. If not, slide a
  ; complete segment.
  if keyword_set(overlap) then begin
    n_loop=2*n_ave-1
    seg_step=npts/2
  endif else begin
    seg_step=npts
    n_loop=n_ave
  endelse

  ; Initialize array indexing
  istart = 0L
  istop = npts-1L

  ; Begin loop for averaging segments through the time series.
  for i =0L,n_loop-1L do begin
    ; window the data using a Hanning envelope,
    ; and then take the FFT of the windowed segment.
    fft_out = fft( data[ istart:istop]*ww, 1)
    ave_spec[0L] = ave_spec[0L] + (1./W)*fft_out[0L]*conj(fft_out[0L])
    ave_spec[1:(npts/2-1)] = ave_spec[1:(npts/2-1)] + $
    (1./W)*(fft_out[1:(npts/2-1)]*conj(fft_out[1:(npts/2-1)]) + $
    reverse(fft_out[npts/2+1:npts-1]) * $
    conj(reverse(fft_out[npts/2+1:npts-1])))
    ave_spec[npts/2] = ave_spec[npts/2] + $
    (1./W)*(fft_out[npts/2]*conj(fft_out[npts/2]))
    istop = istop + seg_step
    istart = istart + seg_step
  endfor

  ; Define frequency array and delta frequency between each point.
  freq = findgen(npts/2+1)/(npts*sample)
  delta_freq = freq[1]-freq[0]

  ; Normalize the average power spectra.
  ave_spec = ave_spec/n_loop

  ; Convert spectral power (amplitude^2) to (amplitude^2/Hz)
  ave_spec[1:npts/2-1] = ave_spec[1:npts/2-1]/delta_freq
  ; The zeroeth and NPTS/2th PSD bins are only half as wide as the rest
  ; of the frequency bins -- add factor of 2 for (amplitude^2/Hz)
  ave_spec[0] = 2*ave_spec[0]/delta_freq
  ave_spec[npts/2] = 2*ave_spec[npts/2]/delta_freq

end
