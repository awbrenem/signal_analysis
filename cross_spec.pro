;+    
; NAME: CROSS_SPEC
;
; PURPOSE: Routine to calculate the cross-spectra of two input time
;          series sig1 and sig2. Program returns coherence and phase
;          delay vs. frequency for the two input signals.
;
; CALLING SEQUENCE: cross_spec, sig1, sig2, n_ave, npts, sample, coh,
;                   frequency, overlap=overlap
;
; INPUTS:
;         sig1  - A one dimensional array containing time series
;                 waveform amplitudes.
;
;         sig2  - A second time series signal similar to sig1.
;
;  npts
;    long, number of points in each FFT.
;
;         N_AVE - The number of segments to divide up the time series
;                 data in sig1 and sig2 for averaging.
;         
;         SAMPLE -The sample time for each point in the time
;                       series data in sig1 and sig2.
;
;  do_xcorr
;    if set, the cross-correlation in time of the two signals will
;    also be computed using the averaged cross- and auto-spectral
;    estimates.
;
; OUTPUTS -
;
;         COH   - The resulting coherence data for sig1 and sig2. This
;                 represents gamma^2 in the popular literature.
;
;         PHASE - The resulting phase delay between the two input
;                 signals sig1 and sig2.
;
;         FREQ  - The frequency axis for the coherence and phase data.
;
;  xcorr
;    float[ M], cross-correlation in time between sig1 and sig2.
;
;  lag
;    float[ M], lag times corresponding to the elements of xcorr.
;    Note that lag is not in monotonically increasing order.  It is
;    instead in the typical, negative lags following positive lags.
;
; KEYWORDS:
;
;         OVERLAP - Set this keyword to slide the cross-spectral
;                   interval by one-half of an interval instead of
;                   averaging together each separate interval. For
;                   most data types this will yield a higher number of
;                   averages and less error per data point than straight
;                   sequential averaging. Note that N_AVE
;                   specifies how many sequential averages are taken
;                   without overlap. Thus N_AVE=4 with /OVERLAP yields
;                   a total of 7 averages, etc. 
;
;
; EXAMPLE:   You have two time series signals sig1 and sig2 that you
;            want to derive a cross-spectra from. Each of these time
;            series data must be of the same length and should be a
;            multiple of 2^n, say 8192 points. To get decent
;            statistics at least 4-8 sequential or sliding averages
;            should be used -- so pick n_ave = 8, npts=1024, and
;            sample = the sample time (in seconds) for the time series
;            data. Alternatively, to get a higher frequency resolution
;            set n_ave=4, npts=2048, and set /overlap. This yields 4+3
;            = 7 sliding averages for the data.
;
;    A note about errors:
;            For n sequential (non-sliding) averages FFT's have a
;            statistical error of 1/n^0.5. For the same time interval,
;            if a sliding average is used of 1/2 a window overlap, n
;            goes up by a factor of 2 but the errors are no longer
;            completely independent -- so the new error per point is
;            1/(0.81*N)^0.5 -- but N=2n so there is less error per point
;            than in the non-sliding average case. For more
;            information about errors and correlation analyses in
;            general, consult Bendat & Piersol, "Engineering
;            Applications of Correlation and Spectral Analysis", 1993.
;
;    A note about the data:
;            Results may not be valid if the number of points per FFT
;            (npts) is not 2^j where j = any integer.
;            This routine assumes that the input data is continuous and
;            that there are no data gaps or bad time points.
;
;
; REVISION HISTORY
;  ver 1.0 04/03/97 G.T.Delory
;  ver 1.1 04/11/97 G.T.Delory
;  ver 1.2 04/24/00 John Bonnell UCBSSL
;    moved HANNING call out of loop over data segments.
;    changed all array indexing from parentheses to square brackets.
;    simplified phase quadrant code (single call to atan( y, x)).
;  ver 2.0 John Bonnell, UCBSSL, September 25, 2000
;    removed n_ave normalizing factor from corr and norm
;      (not needed for cross-spectrum
;       and wrong for overlap-averaged case).
;    added computation of cross-correlation in time
;      and lag time array upon request.
;    removed explicit allocation and normalization 
;      of coh and phase arrays.
;-

pro cross_spec, $
         sig1, $
         sig2, $
         N_AVE=n_ave, $
         NPTS=npts, $
         SAMPLE=sample, $
         coh, $
         phase, $
         freq, $
         OVERLAP=overlap, $
         do_xcorr=do_xcorr, xcorr=xcorr, lag=lag, $
         verbose=verbose
         

; Allocate and Initialize accumulation arrays.
corr = complexarr(1,npts)
norm = fltarr(2,npts)

; Basic error checking and setting defaults for unspecified keywords.

if n_elements(sig1) NE n_elements(sig2) then begin
    print,'CROSS_SPEC: Input signals must have the same number of ' + $
      'elements'
    sig1 = replicate(!values.f_nan,256)
    sig2 = replicate(!values.f_nan,256)
endif

if not (defined(sig1) or defined(sig2)) then begin
    print,'No data provided for one or both of the input channels -- ' + $
      'check calling sequence and try again'
    sig1 = replicate(!values.f_nan,256)
    sig2 = replicate(!values.f_nan,256)
endif

if not (defined(npts) and defined(n_ave)) then begin
    print,'CROSS_SPEC: Number of points per average interval ' + $
      'unspecified.'
    npts = n_elements(sig1)/n_ave
    print,'Assuming '+strcompress(string(npts), /remove_all) + ' ' + $
      'points per interval.'
endif

if not defined(n_ave) then begin
    n_ave=4
    npts = n_elements(sig1)/n_ave
    print,'CROSS_SPEC: Number of averages undefined; assuming n_ave = 4, '
    print,'There are '+strcompress(string(npts), /remove_all)+' ' + $
      'points per averaging interval' 
endif

if keyword_set( verbose) then $
  message, /info, $
  string( npts, n_ave, format='("npts:",X,I6,",",X,"n_ave:",X,I3)')

; If overlap keyword is set, slide half of a segment. If not, slide a
; complete segment.
if keyword_set(overlap) then begin
    n_loop=2*n_ave-1
    seg_step=npts/2
endif else begin
    seg_step=npts
    n_loop = n_ave
endelse

; Initialize indexes for averaging windows.
start=long(0)
stop=long(npts-1)

; compute window coefficients.
ww = hanning( npts)

; Averaging loop - slide window through the data.
for i = long(0),long(n_loop-1) do begin
  sig1_window = sig1[start:stop]*ww
  sig2_window = sig2[start:stop]*ww
  fftsig1 = fft(sig1_window,+1)
  fftsig2 = fft(sig2_window,+1)
  corr = corr+(fftsig1*conj(fftsig2))
  norm[0,*] = norm[0,*]+(fftsig1*conj(fftsig1))
  norm[1,*] = norm[1,*]+(fftsig2*conj(fftsig2))
  stop = stop+seg_step
  start = start+seg_step 
endfor

; calculate cross-correlation in time, if desired.
; note that the cross-correlations come out in the usual,
; wrapped-around, negative-lags following positive lags.
if keyword_set( do_xcorr) then begin
    xcorr = float( fft( corr, -1))
    var = total( norm, 2)/(npts + 0.)
    
    xcorr = xcorr/(sqrt( var[ 0]*var[1]))
    
    n2 = npts/2
    lag = findgen( npts)
    lag[ n2:npts-1L] = lag[ n2:npts-1L] - npts
    lag = sample*lag
endif

; Form cross-spectral function.
cross_spec = corr/((norm[0,*]*norm[1,*])^0.5)

; Derive coherence quantity.
coh = float(cross_spec*CONJ(cross_spec))

; Take first half of coherence function and ignore mirrored part.
coh = coh[0:npts/2]

y = imaginary(cross_spec)       ;Extract imaginary part of coherence
x = float(cross_spec)     ;Extract real part of coherence

; Determine phase angle.
phase = atan( y, x)

;Ignore mirrored part of phase angle data.
phase = phase[0:npts/2]

;Construct frequency axis.
freq = findgen(npts/2+1)/(npts*sample)

end
