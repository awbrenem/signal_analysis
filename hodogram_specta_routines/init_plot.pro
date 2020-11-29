;+
;*****************************************************************************************
;
;  FUNCTION :   init_plot
;  PURPOSE  :   creates the !plot system variable, useful for setting up plots
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
;  INPUT:
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
;
;   CREATED:  01/21/2011
;   CREATED BY:  Aaron W. Breneman
;    LAST MODIFIED:  MM/DD/YYYY   v1.0.0
;    MODIFIED BY: Aaron W. Breneman
;
;*****************************************************************************************
;-

pro init_plot

;for plotting...
;postscript _extra structure for line plots
	ps1 = {linestyle:0,thick:5.8,xthick:5.8,ythick:5.8,zthick:5.8,charsize:1.5,charthick:5.8,xstyle:1,ystyle:1,font:0,font_size:12}
;ps _extra structure for contour plots
	psc1 = {thick:5.8,xthick:5.8,ythick:5.8,zthick:5.8,charsize:1.5,charthick:5.8,xstyle:1,ystyle:1,font:0,font_size:12}
;x-window _extra structure for line plots
	xwin1 = {linestyle:0,thick:1.8,xthick:1.8,ythick:1.8,zthick:1.8,charsize:1.7,charthick:1.2,xstyle:1,ystyle:1,font:1}
;x-win _extra for contour plots
	xwinc1 = {thick:1.8,xthick:1.8,ythick:1.8,zthick:1.8,charsize:1.7,charthick:1.2,xstyle:1,ystyle:1,font:1}


struct = {ps1:ps1,$
		  psc1:psc1,$
		  xwin1:xwin1,$
		  xwinc1:xwinc1,$
		  }

defsysv,'!plot',exists=exists
if exists ne 0 then print,'!plot already set!!' 


if exists eq 0 then begin
   defsysv,'!plot', struct
   print,'!plot system variable set'
endif




end