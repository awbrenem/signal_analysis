;+
;FUNCTION: 	rotate_field_2_vec.pro (compare to Lynn's field_rot.pro)
;
;PURPOSE: 	Returns the tansformation matrix rotated to the input vector. Rotated input vector is 
;			defined to be along new z-hat axis.
;			The new x-hat and z-hat plane will contain the max variance E-field
;			Use the transformation matrix as: vecFA = transpose(m2##vec)
;
;ARGUMENTS:
;       waveform -> [n,3] input vector of waveform data (can get this from read_tds.pro)
;		vec   -> [3] element vector to represent the z-hat direction. Ex. Bo in the coord system of "waveform".
;					Note that the coords of "vec" and "waveform" must be the same!
;		vec2 -> [3] element vector of another vector that you'd like rotated to the coord system of vec. Must be in the same coord as "vec" and "waveform"
;		times --> the times associated w/ waveform data (not used - just returned as part of structure)
;
;KEYWORDS: N/A
;
;RETURNS:  waveform in FA coordinates
;
;CALLING SEQUENCE:
;       	IDL> rotate_FAcoord,vec,waveform
;
;NOTES: 	z_fac is along input vector
;			y_fac is given by vec cross x_max, where x_max is the maximum variance eigenvector
;			x_fac - max variance eigenvector always lies in x-z plane
;
;NOTES2:	I've tested this with the following sample input:
;				1. Bo = [1,0,0], waveform=[[stuff],[0],[0]]
;					Bo_rot = [0,0,1], waveform_rot = [[0],[0],[stuff]]  ;along new z-hat as it should be
;				2. Bo = [1,0,0], waveform=[[0],[stuff],[0]]
;					Bo_rot = [0,0,1], waveform_rot = [[stuff],[0],[0]] ;along Emax as it should be.
;				3. Comparison with fatds.pro (z component will be the same)
;-
;CREATED BY:    Aaron Breneman, 03/16/2010
;
;MODIFICATION HISTORY:
;
;INCLUDED MODULES:
;
;LIBS USED: (none)
;
;DEPENDENCIES: my_min_var_rot.pro
;
;-----------------------------------------------------------------------------------------------------------



function rotate_field_2_vec,vec,waveform,vec2=vec2,times=times

if not keyword_set(vec) then begin
	print,'NO B-FIELD DATA INPUTTED. CANNOT ROTATE'
	return,1
endif
if not keyword_set(waveform) then begin
	print,'NO WAVEFORM DATA INPUTTED. CANNOT ROTATE'
	return,1
endif else begin
	tmp = size(waveform)
	if tmp[2] ne 3 then begin
		print,'INCORRECT SIZE ON INPUT WAVEFORM VECTOR'
		return,1
	endif
endelse
if not keyword_set(times) then begin
	times = indgen(n_elements(waveform[*,0]))
endif

Vecmag = sqrt(vec[0]^2 + vec[1]^2 + vec[2]^2)

;-----------------------------------------------
;find minimum variance field components
;-----------------------------------------------

vals = my_min_var_rot(waveform,bkg_field=vec,/nomssg)
Emax = vals.eigenvectors[*,2]*vals.eigenvalues[0]
Eint = vals.eigenvectors[*,1]*vals.eigenvalues[1]
Emin = vals.eigenvectors[*,0]*vals.eigenvalues[2]

Emax_mag = sqrt(total(Emax^2))
Eint_mag = sqrt(total(Eint^2))
Emin_mag = sqrt(total(Emin^2))
Emax_hat = Emax/Emax_mag

;------------------------------
;PERFORM ROTATION TO STIX COORD
;------------------------------

;zs = [vec[0]/Vecmag,vec[1]/Vecmag,vec[2]/Vecmag]
zs = vec/Vecmag
ys = crossp(zs,Emax_hat)
ysmag = sqrt(total(ys[0]^2 + ys[1]^2 + ys[2]^2))
ys = ys/ysmag
xs = crossp(ys,zs)
xsmag = sqrt(total(xs[0]^2 + xs[1]^2 + xs[2]^2))
xs = xs/xsmag



;redefine vectors in Stix coord -----------
vec_stix = [0,0,Vecmag]
if keyword_set(vec2) then vec2_stix = [total(vec2*xs),total(vec2*ys),total(vec2*zs)] else vec2_stix = replicate(!values.f_nan,3)
Emax_stix = [total(Emax*xs),total(Emax*ys),total(Emax*zs)]
Eint_stix = [total(Eint*xs),total(Eint*ys),total(Eint*zs)]
Emin_stix = [total(Emin*xs),total(Emin*ys),total(Emin*zs)]


vecFAx = waveform[*,0]*xs[0] + waveform[*,1]*xs[1] + waveform[*,2]*xs[2]
vecFAy = waveform[*,0]*ys[0] + waveform[*,1]*ys[1] + waveform[*,2]*ys[2]
vecFAz = waveform[*,0]*zs[0] + waveform[*,1]*zs[1] + waveform[*,2]*zs[2]

struct = {vecFA:[[vecFAx],[vecFAy],[vecFAz]],$
		  times:times,$
		  Emax:Emax,Eint:Eint,Emin:Emin,$
		  Emax_stix:Emax_stix,Eint_stix:Eint_stix,Emin_stix:Emin_stix,$
		  Vec_input:vec,$
		  vec_stix:vec_stix,$
		  vec2_stix:vec2_stix,$
		  theta_kb_BW_ONLY:vals.theta_kb,$
		  dtheta_kb_BW_ONLY:vals.dtheta,$
		  eigenvalues:vals.eigenvalues,$
		  deig_vals:vals.deig_vals,$
		  eigenvectors:vals.eigenvectors,$
		  dmin_vec:vals.dmin_vec,$
		  notes:['all coord given are in terms of the input coord','dtheta_kb is the uncertainty in the estimate of theta_kb: note that this is only the wave normal angle if the input vec is the magnetic field']}


return,struct


end