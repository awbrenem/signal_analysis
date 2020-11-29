;+
;*****************************************************************************************
;
;  FUNCTION : helicity.pro  
;  PURPOSE  : Finds the helicity as a function of time for two waveform components (perp to Bo)
;			 
;
;  CALLED BY: N/A 
;               
;  CALLS:	         
;
;  REQUIRES:    
;               
;
;  INPUT: 	Ex -> Stix perp efield component (x-z plane contains k-vector)
;			Ey -> other perp direction
;			maxang -> the maximum angle two points can be separated by before it is ignored.
;					Defaults to 300 degrees		
;
;
;
;  NOTES: 
;
;   CREATED:  08/15/2011
;   CREATED BY:  Aaron W. Breneman
;    LAST MODIFIED:    v1.1.0
;    MODIFIED BY: 
;
;*****************************************************************************************
;-
;******************************************************************************************

pro helicity,Ex,Ey,time,maxang=maxang;,filterv=filterv

if ~keyword_set(maxang) then maxang=300.
;if ~keyword_set(filterv) then filterv = 3


Etots = sqrt(Ex^2 + Ey^2)
phi = acos(Ex/Etots)/!dtor
goobar = where(Ey lt 0.)
phi[goobar] = 360.-phi[goobar] ;dot product only goes from 0-180. Extend to full range


diff = fltarr(n_elements(time))
for i=0L,n_elements(time)-2 do diff[i] = phi[i+1]-phi[i]

goo = where(abs(diff) gt maxang)
if goo[0] ne -1 then diff[goo] = !values.f_nan


diffN = replicate(!values.f_nan,n_elements(diff))


goobp = where(diff gt 0)
goobn = where(diff lt 0)

if goobp[0] ne -1 then diffN[goobp] = 1
if goobn[0] ne -1 then diffN[goobn] = -1

!p.font=0
;set_plot,'ps'
;device,filename='~/Desktop/helicity.ps'

;window,1,xsize=800,ysize=500
!p.multi = [0,0,2]
!p.charsize = 1.2
plot,time,diff,title='Angle diff b/t successive elements (positive = RH)',xtitle='time',xstyle=1,xrange=[0,0.032];,yrange=[-1,1]
plot,time,diffN,title='Handedness (positive = RH)',xtitle='time',yrange=[-1.2,1.2],xstyle=1,xrange=[0,0.032]

;device,/close
;set_plot,'x'

end




 

;-----------------------------
;OLD METHOD USING DERIVATIVES
;-----------------------------

;dphi = deriv(phi)

;;Median filter to remove spike values
;dphi2 = median(dphi,filterv)

;;Remove fluctuating end elements
;dphi2[0:filterv-1] = 0.
;dphi2[n_elements(dphi2)-filterv:n_elements(dphi2)-1] = 0.


;dphi2 = dphi2/max(dphi2,/nan)
;
;tots = total(dphi2,/nan);/abs(total(dphi2,/nan))


;!p.charsize = 2.
;plot,time*1000.,dphi2,title='Normalized helicity (positive = RH) ; total = ' + strtrim(tots,2),xtitle='time',yrange=[-1,1]



;stop
