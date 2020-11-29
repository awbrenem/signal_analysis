;+
;*****************************************************************************************
;
;  FUNCTION :   fft_power_calc.pro
;  PURPOSE  :    Returns the FFT power spectrum from an input array of vectors.  
;                  The power spectrum can be calculated using windowing by setting
;                  the keyword READ_WIN.  The frequency bins (Hz) associated with
;                  the power spectrum are returned too, along with a Y-scale 
;                  estimate for plotting routines.
;
;  CALLED BY:   
;               NA
;
;  CALLS:
;               my_windowf.pro
;               power_of_2.pro
;
;  REQUIRES:    
;               1)  UMN Modified Wind/3DP IDL Libraries
;
;  INPUT:
;               TT        :  times (s) associated with vector data [DBLARR]
;               DAT       :  [n] array of vectors
;
;  EXAMPLES:    
;               NA
;
;  KEYWORDS:    
;               READ_WIN  :   If set, program uses windowing for FFT calculation
;               RANGE     :   2-element array defining the start and end point elements
;                               to use for plotting the data
;               FORCE_N   :  Set to a scalar (best if power of 2) to force the program
;                               my_power_of_2.pro return an array with this desired
;                               number of elements [e.g.  FORCE_N = 2L^12]
;               SAMP_RA   :  Set to a scalar defining the sample rate of data if odd time
;                               series is sent in (make sure units are consistent)
;
;   CHANGED:  1)  Forced 1D data entry                             [09/25/2008   v1.2.0]
;             2)  Changed logic statement                          [09/29/2008   v1.2.1]
;             3)  Fixed prime # issue near 2^15                    [11/04/2008   v1.2.2]
;             4)  Fixed indexing issue in 3)                       [11/05/2008   v1.2.3]
;             5)  Changed output to return ONLY arrays that have powers of 2 for 
;                  # of elements                                   [11/20/2008   v3.0.0]
;             6)  Fixed scaling issue in 5)                        [11/21/2008   v3.0.1]
;             7)  Fixed windowing issue in 6)                      [11/21/2008   v3.0.2]
;             8)  Fixed scaling issue in 7)                        [11/22/2008   v3.0.3]
;             9)  Changed output structure                         [11/22/2008   v3.1.0]
;            10)  Added KEYWORD:   FORCE_N                         [11/22/2008   v3.1.1]
;            11)  Changed some syntax                              [12/03/2008   v3.1.2]
;            12)  Added KEYWORD:   SAMP_RA                         [04/06/2009   v3.1.3]
;            13)  Changed program my_power_of_2.pro to power_of_2.pro
;                                                                  [08/10/2009   v3.2.0]
;            14)  Fixed units of FFT power spectrum and updated man page
;                                                                  [06/15/2011   v3.3.0]
;            15)  Fixed typo:  had not previously included windowed response in FFT
;                   and corrected issue when zero-padded result had N -> too large
;                                                                  [09/26/2011   v3.4.0]
;            16)  Fixed order of applying window and zero-padding  [04/12/2012   v3.5.0]
;            17)  Fixed indexing typo                              [05/22/2012   v3.5.1]
;
;   NOTES:      
;               1)  The power spectrum units are now correctly in [units^2/Hz]
;
;  REFERENCES:  
;               1)  Paschmann, G. and P.W. Daly (1998), "Analysis Methods for Multi-
;                      Spacecraft Data," ISSI Scientific Report, Noordwijk, 
;                      The Netherlands., Int. Space Sci. Inst.
;
;   CREATED:  08/26/2008
;   CREATED BY:  Lynn B. Wilson III
;    LAST MODIFIED:  05/22/2012   v3.5.1
;    MODIFIED BY: Lynn B. Wilson III
;
;*****************************************************************************************
;-

FUNCTION fft_power_calc,tt,dat,READ_WIN=read_win,RANGE=range,FORCE_N=force_n,SAMP_RA=sra

;-----------------------------------------------------------------------------------------
; => Define dummy variables
;-----------------------------------------------------------------------------------------
nt  = N_ELEMENTS(tt)
d   = REFORM(dat)               ; -Prevent any [1,n] or [n,1] array from going on
px  = DBLARR(nt)                ; -dummy power associated w/ X-comp. if dat is a vector
py  = DBLARR(nt)                ; -dummy power associated w/ Y-comp. if dat is a vector
pz  = DBLARR(nt)                ; -dummy power associated w/ Z-comp. if dat is a vector
pr  = [0d0,0d0]                 ; -dummy power range
fr  = LINDGEN(nt)*1d0           ; -dummy freq. bin array
dum = CREATE_STRUCT('POWER_X',px,'POWER_Y',py,'POWER_Z',pz,'POWER_RA',pr,'FREQ',fr)
;-----------------------------------------------------------------------------------------
; => Define data range
;-----------------------------------------------------------------------------------------
IF KEYWORD_SET(range) THEN BEGIN
  myn = range[1] - range[0]  ; -number of elements used for min. var. calc.
  IF (myn LE nt AND range[0] GE 0 AND range[1] LE nt) THEN BEGIN
    s = range[0]
    e = range[1]
  ENDIF ELSE BEGIN
    print,'Too many elements demanded in keyword: RANGE'
    s   = 0
    e   = nt - 1L
    myn = nt
  ENDELSE
ENDIF ELSE BEGIN
  s   = 0
  e   = nt - 1L
  myn = nt
ENDELSE
;-----------------------------------------------------------------------------------------
; => Set up vars for power spectrum
;-----------------------------------------------------------------------------------------
evlength = tt[e] - tt[s]                        ; => length of data (s)
nfbins   = 1L + e - s
IF KEYWORD_SET(sra) THEN BEGIN
  nsps = (DOUBLE(sra))[0]                       ; => Avg. sample rate (# points/time)
ENDIF ELSE BEGIN
  nsps = (nfbins - 1L)/evlength[0]              ; => Avg. sample rate (# points/time)
ENDELSE
;-----------------------------------------------------------------------------------------
; => Pad array with zeros to force number of elements to a power of 2
;  [Note:  This also increases the frequency resolution of your power spectrum]
;-----------------------------------------------------------------------------------------
; => LBW III  [04/12/2012   v3.5.0]
;      the zero-padding should actually be done AFTER the window is applied [see ref.]
fft_win  = FLTARR(nfbins)                       ; => var. used for FFT windowing routine
IF KEYWORD_SET(read_win) THEN BEGIN
  my_windowf,nfbins - 1L,2,fft_win              ; => Uses a Hanning Window first
ENDIF ELSE BEGIN
  fft_win = REPLICATE(1.,nfbins)
ENDELSE
wsign     = d[s:e]*fft_win

IF KEYWORD_SET(force_n) THEN BEGIN
  focn = force_n
ENDIF ELSE BEGIN
  focn = nt
ENDELSE
signal   = power_of_2(wsign,FORCE_N=focn)   ; => now pad with zeroes
nfbins2  = N_ELEMENTS(signal)                    ; => Redefine # of frequency bins, but NOT sample rate
;-----------------------------------------------------------------------------------------
; => Calculate Power Spectral Density (PSD)
;-----------------------------------------------------------------------------------------
;wsign     = signal*fft_win
freqbins  = DINDGEN(nfbins2/2L)*(nsps/(nfbins2 - 1L))         ; => use Hz scale on FFT power spec.
nevlength = evlength + (nfbins2*1d0 - nfbins*1d0)/nsps        ; => New effective event length

; => LBW III  [04/12/2012   v3.5.0]
t_fftx    = ABS(FFT(signal))^2
;t_fftx    = ABS(FFT(wsign))^2     ;  LBW III 09/26/2011   v3.4.0
t_powx    = (2d0*(1d0*nfbins2)^2)/(1d0*nfbins*nsps[0])*t_fftx ; => [units^2/Hz]
tpx       = t_powx[0L:(nfbins2/2L - 1L)]                      ; => PSD
;-----------------------------------------------------------------------------------------
; => Convert power to dB above background
;-----------------------------------------------------------------------------------------
px = 1d1*ALOG10(tpx)
rx = 1.1*MAX(ABS(px),/NAN)
rn = MIN(ABS(px),/NAN)/1.1
ry = rx + 1d-1*ABS(rx)       ; -MAX power = max power + 10% of max power
rz = rn - 1d-1*ABS(rn)       ; -MIN power = min power - 10% of min power
pr = [rz,ry]                 ; -Power range [MIN,MAX]

p_str = CREATE_STRUCT('POWER_X',px,'POWER_RA',pr,'FREQ',freqbins,'POWER_A',tpx)
RETURN,p_str
END