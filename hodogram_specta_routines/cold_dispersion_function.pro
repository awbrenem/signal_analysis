;+
;*****************************************************************************************
;
;  FUNCTION : cold_dispersion.pro  
;  PURPOSE  : Returns useful cold plasma parameters from full cold plasma dispersion relation.
;			  Requires input freq, Emax/Eint ratio, dens, Bo. 
;
;  CALLED BY: N/A 
;               
;  CALLS:	         
;
;  REQUIRES:    
;               
;
;  INPUT: 	epol -> Ex/Ey ratio (plane parallel to z-hat (Bo) in Stix coord. 
;			freq -> array of freqs (Hz) 
;			dens -> density (cm-3) 
;			Bo   -> background magnetic field (nT).
;			mr -> mass ratio of proton to electron. Default is 1836. Useful to change if
;				 you'd like to compare results to CMA diagrams with artificial mr. 
;			H_plus, He_plus, O_plus -> fractions of H+, He+ and O+. If not set 
;				program assumes H+ plasma. 
;
;  NOTES: -Valid for any cold plasma with no density or magnetic field gradients on the scale
;		  of wavelength.
;		  -Whistlers, ion cyclotron, X-mode, O-mode, hydromagnetic modes....
;		  -On the nightside the plasma should be completely
;			H+ by about 1400 km. 
;	 	  -Hot plasma effects not important for beta<<1 and for freqs above but not too close
;			to the ion cyclotron freq.
;
;
;
;		  All equation references are from Stix 1992.
;
;		  If we know a-priori the wave normal angle then we'd just use eqn 34 to calculate
;		  the index of refraction. This would be the way to proceed when you have Bw
;		  measurements and can determine theta_kb directly. However, with Ew measurements
;		  we have to find the index of refraction via eqn 42. Then we can figure out
;		  theta_kb, etc. Both methods are exact for a cold plasma. 
;
;  ASSUMPTIONS: Assumes plane waves. i.e. that each value of w has a single k and that
;			    there is no spatial dependence of density or magnetic field. 
;
;  OUTPUT:  Structure that contains arrays of |k| (1/km), theta_kb (deg), n, wavelength (km), 
;			resonance cone angle, phase velocity (km/s), cyclotron res energies (eV)
;
;   CHANGED:  1)  NA [MM/DD/YYYY   v1.0.0]
;
;   CREATED:  11/02/2009
;   CREATED BY:  Aaron W. Breneman
;    LAST MODIFIED:    v1.1.0
;    MODIFIED BY: Hacked from J. Wygant's cold plasma dispersion solver. 
;				07/14/2011 - AWB - the speed of light variable is now defined as "c_ms". This
;									was being overwritten by the cold plasma param "C" causing
;									the resonance energies to be incorrect. The cyclotron 
;									resonance energies correspond well to Fig2 in Abel and Thorne98
;
;
;*****************************************************************************************
;-
;******************************************************************************************

function cold_dispersion,epol=epol,freq=freq,dens=dens,Bo=Bo,H_plus=pH,He_plus=pHe,O_plus=pO


;-----------------------------------------------------------------------
;SETUP AND OTHER STUFF
;-----------------------------------------------------------------------


notes = ['freq in Hz','kmag in 1/km','theta_kb in deg','n=index of refraction','wavelength in km',$
		'resangle in deg','phasevel in km/sec','bpol is magnetic polarization ratio','energy in eV for nonrelativistic energies',$
		'resonance and cutoff values in Hz, calculated without ion contributions']

nfreqs = n_elements(freq)

c_ms      = 2.99792458d8      ; -Speed of light in vacuum (m/s)
me     = 9.1093897d-31     ; -Electron mass (kg)
mp     = 1.6726231d-27     ; -Proton mass (kg)
e1eV = 1.6e-19
erg2joule = 1e-7           ; -ergs to joules
e1eV = 1.6e-19              ;joules in 1 eV
BnT2BG = 1e-5               ;nT to Gauss (multiply by this number)

mr1 = 1836.  ;mass ratio of H+ (proton) to electron
mr2 = 1836.*4.    ;              He+ to electron     
mr3 = 1836.*16.    ;              O+  to electron

if ~keyword_set(pHe) and n_elements(pHe) eq 0. then pHe = 0.
if ~keyword_set(pO) and n_elements(pO) eq 0. then pO = 0.
if ~keyword_set(pH) and n_elements(pH) eq 0. then pH = 1.

;density ratios
nr1 = pH
nr2 = pHe
nr3 = pO

dtmp = ''
Botmp = ''
epoltmp = ''
freqtmp = ''

if not keyword_set(dens) then begin
	read,dtmp,prompt='Enter density in cm-3: '
	dens=dtmp
endif	
if not keyword_set(Bo) then begin
	read,Botmp,prompt='Enter Magnitude of B-field (nT)' 
	Bo=Botmp
endif	
if not keyword_set(epol) then begin
	read,epoltmp,prompt='Enter Ex/Ey ratio (Stix coord)' 
	epol=epoltmp
endif
if not keyword_set(freq) then begin
	read,freqtmp,prompt='Enter wave freq in plasma frame (Hz)'
	freq=freqtmp
endif

w = 2.*!pi*freq

;consider just H+ (which doesn't have a neutron -- deuterium does)
w_ce = -28.*Bo*2.*!pi
w_cH = -1*w_ce/mr1
w_cHe = -1*w_ce/mr2
w_cO = -1*w_ce/mr3

w_pe = 8980.*sqrt(dens)*2.*!pi
w_pH = w_pe/sqrt(mr1/nr1)
w_pHe = w_pe/sqrt(mr2/nr2)
w_pO = w_pe/sqrt(mr3/nr3)

fce = abs(w_ce/2./!pi)
fpe = w_pe/2./!pi
    
;------------------------------------------------------------------------------------------------
; We now have all the angular  frequencies in place
; (electron gyrofreq is negative)
; -----------------------------------------------------------------------------------------------
; The next step is to calculate R, L, and P as defined in Stix page 7 (1992 version) equations
; 20, 21, 22 . There is an ion and electron contribution for each. Remember that in 
; Stix the gyrofrequency has a negative sign for the electrons and a positive
; sign for the ions (eqn 12). By using this convention, the contributions to 
; R and L from the ions and electrons have the same algebraic form. We use
; this fact when we sum over the ion and electron contributions to R   
; the ion term is R_i_tmp and the electron term is R_e_tmp
; the total expression is R= 1-(R_i_tmp + R_e_tmp)
;------------------------------------------------------------------------------------------------

R_e_tmp= (w_pe^2)/(w*(w+w_ce))
R_H_tmp= (w_pH^2)/(w*(w+w_cH))
R_He_tmp= (w_pHe^2)/(w*(w+w_cHe))
R_O_tmp= (w_pO^2)/(w*(w+w_cO))
R = 1 - (R_e_tmp + R_H_tmp + R_He_tmp + R_O_tmp)

;we now do the same thing for L where L_e_tmp and L_i_tmp are the electron and
;ion contributions
L_e_tmp = (w_pe^2)/(w*(w-w_ce))
L_H_tmp = (w_pH^2)/(w*(w-w_cH))
L_He_tmp = (w_pHe^2)/(w*(w-w_cHe))
L_O_tmp = (w_pO^2)/(w*(w-w_cO))
L = 1-(L_e_tmp + L_H_tmp + L_He_tmp + L_O_tmp)


; Now for P which is the along B part
P_e_tmp = (w_pe^2)/(w^2)
P_H_tmp = (w_pH^2)/(w^2)
P_He_tmp = (w_pHe^2)/(w^2)
P_O_tmp = (w_pO^2)/(w^2)
P = 1 - (P_e_tmp + P_H_tmp + P_He_tmp + P_O_tmp)


; we know have R,L, and P
;Now following Stix we define S and D 
S=0.5*(R+L)
D=0.5*(R-L)


;------------------------------------------------------------------------------------------------
; Now in this first use of these relation we know the polarization of E. 
; For Stix, the Bfield is along z. k is in the xz plane, and there is a
; relation for the polarization i(Ex/Ey)=(n**2-S)/D (eqn 42 page 10).
; Now we know what (Ex/Ey)
; is and so we can figure out what n = kc/w is. Then we can go back and
; calculate other stuff like the k vector etc.
;------------------------------------------------------------------------------------------------

n2= epol*D+S  ;eqn 42 Stix
n=sqrt(n2)
kvect=n*w/(3.0d5)       ;1/km
wavelength=2.*!pi/kvect  ;km

; Calculate angle for resonance cone Stix page 12 eq 1-45
tanthetasq=-P/S
tantheta=sqrt(tanthetasq)
resangle= (360./6.2830)* atan(tantheta)

tmp = where(finite(resangle) eq 0.)
if tmp[0] ne -1 then resangle[tmp] = 90.

phasevel=3d5/n  ;km/sec


;------------------------------------------------------------------------------------------------
; Now that we know n2=n**2 (the index of refraction) we can use equation 36
; in Stix to get the tan(theta) were theta is the angle between k and B.
; Pretty thrilling.
;------------------------------------------------------------------------------------------------

tan2=-1*P*(n2-R)*(n2-L)/((S*n2-R*L)*(n2-P))
tanx=sqrt(tan2)
theta=atan(tanx)
thetadeg=(360./2./!pi)*theta
kvalue = 1000.*sqrt(4.*!pi^2*fpe^2*freq/(c_ms^2*(fce*cos(thetadeg*!dtor)-freq)))  ;alternative calculation
	;it gives the same value as kvect

A = S*sin(thetadeg*!dtor)^2 + P*cos(thetadeg*!dtor)^2
B = R*L*sin(thetadeg*!dtor)^2 + P*S*(1 + cos(thetadeg*!dtor)^2)
C = P*R*L
F2 = (R*L - P*S)^2 * sin(thetadeg*!dtor)^4 + 4*P^2*D^2*cos(thetadeg*!dtor)^2
F = sqrt(F2)


;------------------------------------------------------------------------------------------------
;Frequencies of cutoffs and resonances
;------------------------------------------------------------------------------------------------

;R=0 cutoff --> no ion term
acoeff = 1.
bcoeff = w_ce
ccoeff = -1*w_pe^2
w_ro_e = (-bcoeff + sqrt(bcoeff^2 - 4*acoeff*ccoeff))/(2*acoeff)
f_ro_e = w_ro_e/2./!pi
;L=0 cutoff --> no ion term, but otherwise no approximations
acoeff = 1.
bcoeff = -1*w_ce
ccoeff = -1*w_pe^2 
w_lo_e = (-bcoeff + sqrt(bcoeff^2 - 4*acoeff*ccoeff))/(2*acoeff)
f_lo_e = w_lo_e/2./!pi


;Gurnett's vesion....gives the same value as mine
;w_lo_e2 = -1*abs(w_ce)/2. +  sqrt((w_ce/2.)^2 + w_pe^2)
;f_lo_e2 = w_lo_e2/2./!pi





s_o = sqrt(w_ce^2 + w_pe^2)/2./!pi

cp_params = {R:R,L:L,D:D,S:S,P:P,A:A,B:B,C:C,F:F,R_cutoff:f_ro_e,L_cutoff:f_lo_e,P_cutoff:w_pe/2./!pi,S_resonance:s_o,$
	R_resonance:abs(w_ce)/2./!pi}


;--------------------------------------------------------------------------
;Find relativistic cyclotron resonance energies for a range of pitch angles
;--------------------------------------------------------------------------

;test pitch angles
pa = [0.,10.,20.,30.,40.,50.,60.,70.,80.,89.]

vz_cycl = fltarr(nfreqs,10) & vz_anom = fltarr(nfreqs,10) & vz_landau = fltarr(nfreqs,10)
vtots_cycl = fltarr(nfreqs,10) & vtots_anom = fltarr(nfreqs,10) & vtots_landau = fltarr(nfreqs,10)
vc_ratio_cycl = fltarr(nfreqs,10) & vc_ratio_anom = fltarr(nfreqs,10) & vc_ratio_landau = fltarr(nfreqs,10)
relativistic_cycl = strarr(nfreqs,10) & relativistic_anom = strarr(nfreqs,10) & relativistic_landau = strarr(nfreqs,10)


kz = abs((kvect/1000.)*cos(thetadeg*!dtor))

for qb = 0,n_elements(pa)-1 do begin

	;define test velocities for solving transcendental equation
	vz = indgen(10000)*c_ms/9999. * cos(pa[qb]*!dtor)

	
	;cyclotron resonance energy
	gama = 1/sqrt(1-(vz/c_ms/cos(pa[qb]*!dtor))^2)  ;relativistic gamma factor
	f1cycl = vz
	f2cycl = (2*!pi/kz)*(1*fce/gama - freq)
	diff = abs(f1cycl-f2cycl)
	tmp = min(diff,val)	
	vz_cycl[*,qb] = vz[val]   ;m/s
	;vpar2 = (2*!pi/kz)*[fce*sqrt(1-(vz[val]/c_ms/cos(pa[qb]*!dtor))^2) - f]  ;m/s
	vtots_cycl[*,qb] = vz_cycl[*,qb]/cos(pa[qb]*!dtor)  ;electron velocity
	vc_ratio_cycl[*,qb] = vtots_cycl[*,qb]/c_ms
	tpp = where(vc_ratio_cycl[*,qb] ge 0.1)
	;------------------------------------------------
	;anomalous resonance energy
	f1anom = -1*vz
	f2anom = (2*!pi/kz)*(-1*fce*sqrt(1-(vz/c_ms/cos(pa[qb]*!dtor))^2) - freq)
	diff = abs(f1anom-f2anom)
	tmp = min(diff,val)
	vz_anom[*,qb] = vz[val]   ;m/s
	;vpar2 = (2*!pi/kz)*[fce*sqrt(1-(-1*vz[val]/c_ms/cos(pa[*,qb]*!dtor))^2) - f]  ;m/s
	vtots_anom[*,qb] = vz_anom[*,qb]/cos(pa[qb]*!dtor)  ;ion velocity
	vc_ratio_anom[*,qb] = vtots_anom[*,qb]/c_ms
	;-----------------------------------------------
	;landau resonance energy
	vz_landau[*,qb] = 2*!pi*freq/kz
	vtots_landau[*,qb] = vz_landau[*,qb]/cos(pa[qb]*!dtor)
	vc_ratio_landau[*,qb] = vtots_landau[*,qb]/c_ms

endfor ;for each pitch angle

;Relativistic energy in eV (e.g. p37 in "Modern Physics, 2nd edition")
Etots_cycl = 0.511d6/sqrt(1-(vtots_cycl^2/c_ms^2)) - 0.511d6
Etots_anom = 0.511d6/sqrt(1-(vtots_anom^2/c_ms^2)) - 0.511d6
Etots_landau = 0.511d6/sqrt(1-(vtots_landau^2/c_ms^2)) - 0.511d6

Ez_cycl = 0.511d6/sqrt(1-(vz_cycl^2/c_ms^2)) - 0.511d6
Ez_anom = 0.511d6/sqrt(1-(vz_anom^2/c_ms^2)) - 0.511d6
Ez_landau = 0.511d6/sqrt(1-(vz_landau^2/c_ms^2)) - 0.511d6



notes2 = ['values are for each pitch angle','resonance energies in keV','velocities in km/s','vc_ratio_xxx is the ratio v/c']
cyclo_res = {vz:vz_cycl/1000.,vtots:vtots_cycl/1000.,vc_ratio:vc_ratio_cycl,Ez:Ez_cycl/1000.,Etots:Etots_cycl/1000.,pitch_angles:pa,notes:notes2}
anoma_res = {vz:vz_anom/1000.,vtots:vtots_anom/1000.,vc_ratio:vc_ratio_anom,Ez:Ez_anom/1000.,Etots:Etots_anom/1000.,pitch_angles:pa,notes:notes2}
landau_res = {vz:vz_landau/1000.,vtots:vtots_landau/1000.,vc_ratio:vc_ratio_landau,Ez:Ez_landau/1000.,Etots:Etots_landau/1000.,pitch_angles:pa,notes:notes2}

;------------------------------------------------------------------------------------------------
;VARIOUS PLASMA QUANTITIES
;------------------------------------------------------------------------------------------------


fpe = w_pe/2./!pi
fpH = w_pH/2./!pi
fpHe = w_pHe/2./!pi
fpO= w_pO/2./!pi

fce = abs(w_ce/2./!pi)
fcH = w_cH/2./!pi
fcHe = w_cHe/2./!pi
fcO = w_cO/2./!pi

meff = 1/(pH/1. + pHe/4. + pO/16.)    
flhr = sqrt(fpe^2*fce^2/(fpe^2 + fce^2)*(1/1836./meff))


e_inertial = 5.31e5/sqrt(dens)/100./1000.               ;e- inertial length (skin depth) (km)
ion_inertial = 2.28e7/sqrt(dens)/100./1000.             ;H+ inertial length (km)
VA = 2.18e11*Bo*BnT2BG/sqrt(dens)/100./1000.             ;H+ Alfven vel (km/s)
VA2 = fcH/fpH*3e5                                    ;Alternate H+ Alfven vel (km/s)
VAe = VA * sqrt(1836.)                               ;electron Alfven vel (km/sec)
mepp = (Bo*BnT2BG)^2/8./!pi/dens * erg2joule /e1ev       ;characteristic magnetic energy per particle (eV)


;------------------------------------------------------------------------------------------------
; now we calculate the magnetic field polarization in the plane perpendicular
; to the k vector. We use the fact that the plane of the electric 
; field variations is perpendicular to B and we rotate to a plane perpendicular
; to K (by theta) and use the phase velocity measurement in the relation between; E and B.
; This give a simple relation between electric and magnetic field polarization:
;------------------------------------------------------------------------------------------------

bpol=epol*cos(theta)


;Calculate Ex,Ey,Ez
;Since I only know epol=Ex/Ey I'll assume that Ey=1 and Ex=epol. Ez will then be related to these.
Ex2Ey = epol
Ez2Ex = (n^2*cos(theta)*sin(theta))/(n^2*sin(theta)^2-P)


;------------------------------------------------------------------------------------------------

struct =   {freq:freq, $
			kmag:kvect, $
			theta_kb:thetadeg, $
			n:n, $
			dens:dens,$
			Bo:Bo,$
			wavelength:wavelength, $
			resangle:resangle, $
			phasevel:phasevel, $
			bpol:bpol, $
			cyclo_res:cyclo_res,$
			anoma_res:anoma_res,$
			landau_res:landau_res,$
			energy_notes:notes2,$
			fpe:fpe,$
			fpH:fpH,$
			fpHe:fpHe,$			
			fpO:fpO,$
			fce:fce,$
			fcH:fcH,$
			fcHe:fcHe,$
			fcO:fcO,$
			flhr:flhr,$
			e_inertial:e_inertial,$
			ion_inertial:ion_inertial,$
			VA_ion:VA,$
			VA_electron:VAe,$
			mepp:mepp,$
			Ex2Ey:epol,$
			Ez2Ex:Ez2Ex,$
			cp_params:cp_params,$
			notes:notes}

return,struct
end


