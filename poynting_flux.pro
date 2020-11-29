;+
;  PROCEDURE :  poynting_flux    (generalized and simplified from rbsp_poynting_flux)
;
;  PURPOSE  : calculates Poynting flux (ergs/s/cm^2). Separates into field-aligned and perp components
;			  Returns results as tplot variables.
;       ------------------------------------------------------------------------
;       NOTE: Includes values with and w/o spinaxis component. This can be generalized
;       to say that if there is one component that is less reliable, then have it be
;       the X component (on RBSP EFW this would correspond to X-MGSE spinaxis direction)
;       ------------------------------------------------------------------------
;
;
;  REQUIRES:  tplot library
;
;  KEYWORDS:
;       Bw -> tplot name of the [n,3] magnetic field
;             waveform (nT).
;			  Ew -> tplot name of the [n,3] electric field
;			        waveform (mV/m)
;			  Tshort, Tlong -> short and long period of waveform to use.
;
;       Bo -> (optional keyword) array of DC magnetic field directions.
;					    Use, for example, if Bw is from AC-coupled data.
;             If not included then Bw is downsampled and used as background field
;
;       method1 -> (default) uses smoothing instead of a bandpass_filter which
;                  results in a sharper freq rolloff than smoothing method.
;                  Both methods give very similar results for test chorus waves
;       method2 -> The waveform data are first downsampled to 1/Tshort to avoid the wasted
;			             computation of having to run program with data at unnecessarily high
;			             sample rate. The waveform is then run through
;			             bandpass_filter to the frequency range flow=1/Tlong...fhigh=1/Tshort
;
;       mapiono -> Set to additionally map variables to the ionosphere
;
;       baddirection -> set to 1,2, or 3. Used to specify which of the input
;               components is less reliable (e.g. spinaxis direction). Also used
;               to define the Pflux coordinate system (P1,P2,P3). Defaults to
;               x-hat. If there is no bad direction then just ignore the output
;               variables with the "nospinaxis" tag.
;
;
;   NOTES:     DO NOT INPUT DATA WITH SIGNIFICANT GAPS
;
;********************************************************************
;           TESTING
;
;  1) Rotation testing: rotating the Evec = [[Ep1],[Ep2],[Ep3]]
;     back to the inputcoord gives the exact same result as Ew input.
;  2) Tested variations of "baddirection" keyword.
;       a) The output values in "pfluxcoord" are identical, as they should be
;          since this direction is used to define the Pflux coordinate system.
;          This is true both of the full 3D version and the "nospinaxis" versions.
;       b) The output values of full 3D pflux transformed back to inputcoord are also identical, as they
;          should be b/c all 3 components are used.
;       c) The output values of pflux_nospinaxis_inputcoord should be different.
;          b/c we're selectively ignoring a different component depending on
;          what baddirection is set to. NOTE: NOT CURRENTLY WORKING. THE RESULTS
;           ARE ALWAYS THE SAME!!!!!!!!!!!!
;  3) Unit vectors are orthogonal
;  4) Tested on Van Allen Probes chorus from EFW's B1
;             data on 2014-08-27 at ~07:42 on probe A. A is at +20
;             mlat and the Pflux indicates propagation away from eq
;             with magnitude values consistent with those in Li et
;             al., 2013 and Santolik et al., 2010
;  5) "nospinaxis" results nearly identical if all 3 Ew components are input vs
;       only inputting the two good components.
;
;
;********************************************************************
;
;			 Poynting flux coord system
;   		 	P1 = B x unitvec  (unitvec defaults to input xhat. Can also set it to
;                            whatever direction contains less reliable data, like
;                            the spin axis component)
;				  P2 = B x P1
;  		   	P3 = B
;
;
;			 The output tplot variables are:
;
;			 	These three output variables contain a mix of spin axis and spin plane components:
;		 		  pflux_perp1  -> Poynting flux in perp1 (to Bo) direction
;			 		pflux_perp2  -> Poynting flux in perp2 (to Bo) direction
; 			 	pflux_para   -> Poynting flux along Bo
;
;			 	These partial Poynting flux calculations contain only values in plane
;       perpendicular to "unitvec":
;			 		pflux_nospinaxis_perp
;			 		pflux_nospinaxis_para
;
;
;   CREATED:  11/28/2012
;   CREATED BY:  Aaron W. Breneman
;    LAST MODIFIED:  MM/DD/YYYY   v1.0.0
;    MODIFIED BY:
;
;-


pro poynting_flux,Bw,Ew,Tshort,Tlong,$
  Bo=Bo,$
  method2=method2,$
  mapiono=mapiono,$
  baddirection=badir


  ;Need to define unit vector used for rotation to Pflux coord.
  ;By default we'll use input x-direction. However, if one of the three input
  ;coord is less trustworthy (e.g. short spin axis direction) then we'll use
  ;that direction
  if ~KEYWORD_SET(badir) then $
    unitvec = [1.,0.,0.] $
  else begin
    if badir eq 1 then unitvec = [1.,0.,0.]
    if badir eq 2 then unitvec = [0.,1.,0.]
    if badir eq 3 then unitvec = [0.,0.,1.]
  endelse


  get_data,Bw,data=Bw_test
  get_data,Ew,data=Ew_test

  if ~is_struct(Bw_test) or ~is_struct(Ew_test) then begin
    print,'**************************************************'
    print,'ONE OR BOTH OF THE TPLOT VARIABLES CONTAIN NO DATA'
    print,'....RETURNING....'
    print,'**************************************************'
    return
  endif



  ;Get DC magnetic field and use to define P1,P2,P3 directions
  ;If DC Bo is not input separately, then downsample Bw to use as Bo
  ;Otherwise just use inputted Bo
  if ~keyword_set(Bo) then begin
    rbsp_downsample,Bw,suffix='_DC',1/40.
    Bdc = Bw + '_DC'
  endif else begin
    tinterpol_mxn,Bo,Bw,newname='Mag_DC',/spline
    Bdc = 'Mag_DC'
  endelse



  ;Interpolate to get MagDC and Ew data to be on the same times as the Bw data
  get_data,Bw,data=goo
  times = goo.x
  tinterpol_mxn,Ew,times,/spline
  Ew = Ew + '_interp'
  tinterpol_mxn,Bdc,times,/spline
  Bdc = Bdc + '_interp'






  ;Define new coordinate system
  nelem = n_elements(times)

  ;P1 unit vector
  p1 = double([[replicate(0,nelem)],[replicate(0,nelem)],[replicate(0,nelem)]])
  get_data,Bdc,data=Bo_dc
  for i=0L,nelem-1 do p1[i,*] = crossp(Bo_dc.y[i,*],unitvec)
  ;normalize p1
  p1mag = fltarr(nelem)
  for i=0L,nelem-1 do p1mag[i] = sqrt(p1[i,0]^2 + p1[i,1]^2 + p1[i,2]^2)
  for i=0L,nelem-1 do p1[i,*] = p1[i,*]/p1mag[i]


  ;P2 unit vector
  p2 = p1
  for i=0L,nelem-1 do p2[i,*] = crossp(Bo_dc.y[i,*],p1[i,*])
  ;normalize p2
  p2mag = fltarr(nelem)
  for i=0L,nelem-1 do p2mag[i] = sqrt(p2[i,0]^2 + p2[i,1]^2 + p2[i,2]^2)
  for i=0L,nelem-1 do p2[i,*] = p2[i,*]/p2mag[i]


  ;P3 unit vector: Background magnetic field direction
  Bmag_dc = sqrt(Bo_dc.y[*,0]^2 + Bo_dc.y[*,1]^2 + Bo_dc.y[*,2]^2)
  Bo_dc_uvec = Bo_dc.y
  Bo_dc_uvec[*,0] = Bo_dc.y[*,0]/Bmag_dc
  Bo_dc_uvec[*,1] = Bo_dc.y[*,1]/Bmag_dc
  Bo_dc_uvec[*,2] = Bo_dc.y[*,2]/Bmag_dc
  p3 = Bo_dc_uvec



  ;*********************************************
  ;Test to make sure unit vectors are orthogonal
  ;*********************************************

  ;for i=0,3000 do print,acos(total(p1[i,*]*p2[i,*])/(p1mag[i]*p2mag[i]))/!dtor   ;perp!
  ;for i=0,3000 do print,acos(total(p1[i,*]*Bmgse.y[i,*])/(p1mag[i]*Bmag[i]))/!dtor   ;perp!
  ;for i=0,3000 do print,acos(total(p2[i,*]*Bmgse.y[i,*])/(p2mag[i]*Bmag[i]))/!dtor   ;perp!




  ;----------------------------------------------------------------------------
  ;Now we've defined our Poynting flux unit vectors as p1,p2,
  ;P3=B_uvec. Project the Ew, Bw and Bo data into these three directions
  ;----------------------------------------------------------------------------




  ;Calculate the time-varying rotation matrix
  rotmat = replicate(0.,3,3,nelem)
  for i=0L,nelem-1 do rotmat[*,*,i] = rot_mat(reform(P3[i,*]),reform(P1[i,*]))




  get_data,Ew,data=Egoo
  get_data,Bw,data=Bgoo

  Ewave = Egoo.y
  Bwave = Bgoo.y

  Ep = rotmat ## Ewave
  Bp = rotmat ## Bwave



  ;Determine rates, etc.
  goo = rbsp_sample_rate(times,out_med_avg=medavg)
  rate = medavg[0]


;--------------------------------------------------
;Method 1 - use IDL's smooth function
;--------------------------------------------------

  if ~keyword_set(method2) then begin


    ;Find lower and upper periods for smoothing
    detren = floor(Tlong * rate)
    smoo = floor(Tshort * rate)

    dE=[[Ep[*,0]-smooth(Ep[*,0],detren,/nan)],[Ep[*,1]-smooth(Ep[*,1],detren,/nan)],[Ep[*,2]-smooth(Ep[*,2],detren,/nan)]]
    dB=[[Bp[*,0]-smooth(Bp[*,0],detren,/nan)],[Bp[*,1]-smooth(Bp[*,1],detren,/nan)],[Bp[*,2]-smooth(Bp[*,2],detren,/nan)]]

    fE1=smooth(dE[*,0],smoo,/nan)
    fE2=smooth(dE[*,1],smoo,/nan)
    fE3=smooth(dE[*,2],smoo,/nan)
    fB1=smooth(dB[*,0],smoo,/nan)
    fB2=smooth(dB[*,1],smoo,/nan)
    fB3=smooth(dB[*,2],smoo,/nan)

    B1bg=smooth((Bo_dc_uvec[*,0]),smoo,/nan)
    B2bg=smooth((Bo_dc_uvec[*,1]),smoo,/nan)
    B3bg=smooth((Bo_dc_uvec[*,2]),smoo,/nan)

    fE1 = fE1/1000. & fE2 = fE2/1000. & fE3 = fE3/1000.  ;V/m
    fB1 = fB1/1d9 & fB2 = fB2/1d9 & fB3 = fB3/1d9        ;Tesla


  endif else begin

  ;--------------------------------------------------
  ;Method 2 - use bandpass filter
  ;--------------------------------------------------

    ;Downsample both the Bw and Ew based on Tshort.
    ;No need to have the cadence at a higher rate than 1/Tshort
    sr = 2/Tshort
    nyquist = sr/2.

    rbsp_downsample,[Bw,Ew],suffix='_DS_tmp',sr

    Bw = Bw +  '_DS_tmp'
    Ew = Ew +  '_DS_tmp'



    ;--------------------------------------------------
    ;At this point Ep, Bp and Bdc have a sample rate of
    ;2/Tshort and are sampled at the same times.
    ;We want to bandpass so that the lowest possible
    ;frequency is 1/Tlong Samples/sec
    ;--------------------------------------------------


    ;Define frequencies as a fraction of Nyquist
    flow = (1/Tlong)/nyquist
    fhigh = (1/Tshort)/nyquist


    ;Zero-pad these arrays to speed up FFT
    fac = 1
    nelem = n_elements(Ep[*,0])
    while 2L^fac lt n_elements(Ep[*,0]) do fac++

    tst = 2L^fac - nelem
    if tst ne 0. then begin
      addarr = fltarr(2L^fac - nelem) ;array of zeros
      Ep2 = [Ep,[[addarr],[addarr],[addarr]]]
      Bp2 = [Bp,[[addarr],[addarr],[addarr]]]
    endif else begin
      Ep2 = Ep
      Bp2 = Bp
    endelse

    Epf = BANDPASS_FILTER(Ep2,flow,fhigh) ;,/gaussian)
    Bpf = BANDPASS_FILTER(Bp2,flow,fhigh) ;,/gaussian)


    smoo = floor(Tshort * nyquist)
    B1bg=smooth((Bo_dc_uvec[*,0]),smoo,/nan)
    B2bg=smooth((Bo_dc_uvec[*,1]),smoo,/nan)
    B3bg=smooth((Bo_dc_uvec[*,2]),smoo,/nan)


    ;Remove the padded zeros
    fE = Epf[0:nelem-1,*]/1000.   ;V/m
    fB = Bpf[0:nelem-1,*]/1d9     ;Tesla

    fE1 = fE[*,0] & fE2 = fE[*,1] & fE3 = fE[*,2]
    fB1 = fB[*,0] & fB2 = fB[*,1] & fB3 = fB[*,2]

  endelse







;--------------------------------------------------
;Calculate Poynting flux
;--------------------------------------------------

  muo = 4d0*!DPI*1d-7        ; -Permeability of free space (N/A^2)

  ;J/m^2/s
  S1=(fE2*fB3-fE3*fB2)/muo
  S2=(fE3*fB1-fE1*fB3)/muo
  S3=(fE1*fB2-fE2*fB1)/muo



  ;THE P1 DIRECTION IS THE ONLY COMPONENT GUARANTEED TO NOT HAVE ANY CONTAMINATED
  ;SPINAXIS DATA. SAFE COMPONENTS ARE FE1, AND FB1,FB2,FB3. THEREFORE THE
  ;PROJECTION OF THIS ONTO BO IS FE1*FB2 AND THE PERP PROJECTION IS FE1*FB3

  ;Sp below stands for spinplane. Note that Sperp_clean-hat is a different direction
  Sperp_clean = (fE1*fB3)/muo  ;(in -1*p2-hat direction (same as original unitvec direction))
  Spar_clean = (fE1*fB2)/muo  ;(in p3-hat direction)


  ;erg/s/cm2
  S1 = S1*(1d7/1d4) & S2 = S2*(1d7/1d4) & S3 = S3*(1d7/1d4)
  Sperp_clean = Sperp_clean*(1d7/1d4) & Spar_clean = Spar_clean*(1d7/1d4)

  Svec_pfluxcoord = [[S1],[S2],[S3]]
  S1tmp = S1
  S1tmp[*] = 0.

  Svec_nospinaxis_pfluxcoord = [[S1tmp],[-1*Sperp_clean],[Spar_clean]]

;--------------------------------------------------
;Find angle b/t Poynting flux and Bo
;--------------------------------------------------

  Bbkgnd = [[B1bg],[B2bg],[B3bg]]

  ;Bo_dc_uvec[*,0]
  S_mag = sqrt(S1^2 + S2^2 + S3^2)

  ;Bo defined to be along the third component
  angle_SB = acos((S3*Bbkgnd[*,2])/S_mag)/!dtor
  store_data,'pflux_angle_pflux_Bo',data={x:times,y:angle_SB}


  ;Smooth the wave normal angle calculation
  stime = 100.*(1/rate)
  if stime lt (times[n_elements(times)-1]-times[0]) then rbsp_detrend,'pflux_angle_pflux_Bo',stime




  ;Define the P unit vectors in terms of input coord system (xhat,yhat,zhat)
  Svec_inputcoord = Svec_pfluxcoord
  Svec_inputcoord[*] = 0.
  for i=0,nelem-1 do Svec_inputcoord[i,*] = reform(rotmat[*,*,i]) # reform(Svec_pfluxcoord[i,*])


  ;Same as above but with the "nospinaxis" version. From above, Sperp_clean is in the
  ;-1*p2-hat direction and Spar_clean is in the p3-hat direction
  Svec_nospinaxis_inputcoord = Svec_nospinaxis_pfluxcoord
  Svec_nospinaxis_inputcoord[*] = 0.
  for i=0,nelem-1 do Svec_nospinaxis_inputcoord[i,*] = reform(rotmat[*,*,i]) # reform(Svec_nospinaxis_pfluxcoord[i,*])



    ;  ;Test rotate Ew back to input coord (WORKS!!!!)
    ;  Evecx = fltarr(n_elements(S1)) & Evecy = Evecx & Evecz = Evecx
    ;  for i=0L,n_elements(S1)-1 do Evecx[i] = Ep1[i]*p1[i,0] + Ep2[i]*p2[i,0] + Ep3[i]*p3[i,0]
    ;  for i=0L,n_elements(S1)-1 do Evecy[i] = Ep1[i]*p1[i,1] + Ep2[i]*p2[i,1] + Ep3[i]*p3[i,1]
    ;  for i=0L,n_elements(S1)-1 do Evecz[i] = Ep1[i]*p1[i,2] + Ep2[i]*p2[i,2] + Ep3[i]*p3[i,2]
    ;  Evec_inputcoord = [[Evecx],[Evecy],[Evecz]]



;------------------------------------
;Estimate mapped Poynting flux
;From flux tube conservation B1*A1 = B2*A2  (A=cross sectional area of flux tube)
;B2/B1 = A1/A2
;P = EB/A ~ 1/A
;P2/P1 = A1/A2 = B2/B1
;Assume an ionospheric magnetic field of 45000 nT at 100km. This value shouldn't change too much
;and is good enough for a rough mapping estimate.
;------------------------------------

  if KEYWORD_SET(mapiono) then begin
    S1_ion = 45000d * S1/Bmag_dc
    S2_ion = 45000d * S2/Bmag_dc
    S3_ion = 45000d * S3/Bmag_dc
    Sperp_clean_ion = 45000d * Sperp_clean/Bmag_dc
    Spar_clean_ion = 45000d * Spar_clean/Bmag_dc
  endif



  ;Store all as tplot variables
  store_data,'pflux_pfluxcoord',data={x:times,y:Svec_pfluxcoord}
  store_data,'pflux_nospinaxis_pfluxcoord',data={x:times,y:Svec_nospinaxis_pfluxcoord}
  store_data,'pflux_Ew_pfluxcoord',data={x:times,y:1000.*[[fE1],[fE2],[fE3]]} ;change back to mV/m
  store_data,'pflux_Bw_pfluxcoord',data={x:times,y:1d9*[[fB1],[fB2],[fB3]]}   ;change back to nT
  store_data,'pflux_angle_pflux_Bo',data={x:times,y:angle_SB}
  store_data,'pflux_inputcoord',data={x:times,y:Svec_inputcoord}
  store_data,'pflux_nospinaxis_inputcoord',data={x:times,y:Svec_nospinaxis_inputcoord}



  if KEYWORD_SET(mapiono) then begin
    store_data,'pflux_para_iono_pfluxcoord',data={x:times,y:S3_ion}
    store_data,'pflux_perp1_iono_pfluxcoord',data={x:times,y:S1_ion}
    store_data,'pflux_perp2_iono_pfluxcoord',data={x:times,y:S2_ion}
    store_data,'pflux_nospinaxis_para_iono_pfluxcoord',data={x:times,y:Spar_clean_ion}
    store_data,'pflux_nospinaxis_perp_iono_pfluxcoord',data={x:times,y:Sperp_clean_ion}

    options,'pflux_nospinaxis_perp_iono_pfluxcoord','ytitle','pflux!Cmapped!Cto ionosphere!Cperp to Bo!C[erg/cm^2/s]'
    options,'pflux_nospinaxis_para_iono_pfluxcoord','ytitle','pflux!Cmapped!Cto ionosphere!Cpara to Bo!C[erg/cm^2/s]'
    options,'pflux_perp1_iono_pfluxcoord','ytitle','pflux!Cmapped!Cperp1 to Bo!C[erg/cm^2/s]'
    options,'pflux_perp2_iono_pfluxcoord','ytitle','pflux!Cmapped!Cperp2 to Bo!C[erg/cm^2/s]'
    options,'pflux_para_iono_pfluxcoord','ytitle','pflux!Cmapped!Cparallel to Bo!C[erg/cm^2/s]'

    options,'pflux_nospinaxis_perp_iono_pfluxcoord','labels','No spin axis!C comp'
    options,'pflux_nospinaxis_para_iono_pfluxcoord','labels','No spin axis!C comp!C+ along Bo'
    options,'pflux_nospinaxis_para_iono_pfluxcoord','colors',2
    options,'pflux_nospinaxis_perp_iono_pfluxcoord','colors',1
    ylim,['pflux_nospinaxis_perp_iono_pfluxcoord','pflux_nospinaxis_para_iono_pfluxcoord'],-10,10
    ylim,['pflux_perp1_iono_pfluxcoord','pflux_perp2_iono_pfluxcoord','pflux_para_iono_pfluxcoord'],0,0
  endif


  options,'pflux_pfluxcoord','ytitle','Poynting flux in P1,P2,P3 coord system!C[erg/cm^2/s]'
  options,'pflux_nospinaxis_pfluxcoord','ytitle','Poynting flux in P1,P2,P3 coord system!Ccalculated without spinaxis component!C[erg/cm^2/s]'
  options,'pflux_inputcoord','ytitle','Poynting flux in input coord system!C[erg/cm^2/s]'
  options,'pflux_nospinaxis_inputcoord','ytitle','Poynting flux in input coord system!Ccalculated without spinaxis component!C[erg/cm^2/s]'
  options,'pflux_Ew_pfluxcoord','ytitle','Ew!Cpflux coord!C[mV/m]'
  options,'pflux_Bw_pfluxcoord','ytitle','Bw!Cpflux coord!C[nT]'
  options,'pflux_angle_pflux_Bo','ytitle','Angle (deg) b/t!CBo and Pflux'
  options,'pflux_nospinaxis_pfluxcoord','labels','No spin axis!C comp'



end
