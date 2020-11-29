;+
;*****************************************************************************************
;
;  PROCEDURE : plot_spec
;  PURPOSE  :  plots nice looking spectrograms
;
;               
;
;  CALLS:
;               
;  CALLING SEQUENCE:
;       IDL> PLOT_SPEC, Z, X, Y
;
;
;
;  KEYWORDS:    
;				tz -> Z values
;				tx -> time bins for spec (get from slide_spec.pro)
;				ty -> freq bins for spec (get from slide_spec.pro)	
;
;				/nocb -> no color bar. Note that dfanning's color bar program can fail
;					if you have data that is spaced over too many decades. Probably a 
;					floating point accuracy issue. Also, note that the plots won't 
;					be openable in adobe illustrator unless this keyword is set. 
;				/no_interp -> don't interpolate the colors
;				shades -> the color shades to plot. Defaults to: 
;					if (zlog) then cs=bytscl(alog10(z),min=min(alog10(z),/nan),max=max(alog10(z),/nan)) else cs=bytscl(z)
;				julian -> plot julian dates
;				hline  -> oplot horizontal lines
;				vline  -> oplot vertical lines
;				_extra -> for plots except for colorbar
;				_extracb -> for colorbar
;
;
;  NOTES:
;       The auto-shading won't work if the input array tz has NaNs. Program searches for these
;		and changes them to zeros. 
;		The shade_surf function with the /t3d option leaves a z-axis relic that can't (as far
;		as we can tell) be removed. This is an error in IDL. 
;
;		Can get data from slide_spec.pro
;
;  CREATED BY:    Kris Kersten, Sept. 2008
;
;  MODIFICATION HISTORY:
;       09/28/2008  created, KK
;       01/28/2009  minor modifications by Aaron Breneman to be compatible with morlet specs 
;		04/07/2011  AWB - added nicer time labels via Julian dates
;		09/29/2011  AWB - added hline and vline keywords
;		02/12/2012  AWB - added _extracb for colorbar
;
;  INCLUDED MODULES:
;		plot_spec
;
;  LIBS USED:
;
;
;  DEPENDENCIES:
;	  From the Coyote Library (http://www.dfanning.com):
;		colorbar.pro
;		error_message.pro
;		fpufix.pro
;		logscl.pro
;		scale_vector.pro
;
;               
;
;  EXAMPLES:  
;			Example for data that spans many orders of magnitude:
;			In this case most of the power is in the 0-50 kHz range, but the total freq range extends up to
;			MHz. Therefore, to avoid having the higher freq stuff determine the plot range I'm manually setting
;			the keyword "shades" with a minimum cutoff.
;
;		plot_spec,data,times,freq_bins,/no_interp,shades=bytscl(alog10(sld1),min=-5.),/zlog
;  
;               
;
;               
;
;*****************************************************************************************
;-
;*****************************************************************************************



pro plot_spec,tz,tx,ty,xlog=xlog,ylog=ylog,zlog=zlog,xtitle=xtitle,$
			  ytitle=ytitle,ztitle=ztitle,title=title,zrange=zrange,julian=julian,$
			  yrange=yrange,no_interp=no_interp,shades=color_shades,_extra=extra,_extracb=extracb,$
			  nocb=nocb,hline=hline,vline=vline
			

	;test for Nan values
	tmp = where(finite(tz) eq 0.)
	if tmp[0] ne -1 then begin
		tz[tmp] = 0.
		print,'**********************************************************************************************'
		print,'FROM PLOT_SPEC.PRO : NAN VALUES HAVE BEEN SET TO ZERO B/C THE COLOR_SHADES WONT WORK OTHERWISE'
		print,'**********************************************************************************************'
	endif	
	
			
			
	; don't clobber the original data!
	z=tz
	x=tx
	y=ty
	if keyword_set(color_shades) then cs = color_shades	
	
	xlog=keyword_set(xlog)
    ylog=keyword_set(ylog)
    zlog=keyword_set(zlog)    
	
	
	
	if ~keyword_set(color_shades) then begin
		if (zlog) then cs=bytscl(alog10(z),min=min(alog10(z),/nan),max=max(alog10(z),/nan)) else cs=bytscl(z)
	endif 
	
	cs = reform(cs)

	if ~keyword_set(zrange) then begin
		zrange = [min(cs,/nan),max(cs,/nan)]
	endif

		
	; the zeroth bit of the binary representation of !d.flags indicates the the current
	; graphics device's scalable pixel capability, so test the base-10 representation
	; for evenness to determine support (borrowed from D. Fanning's colorbar.pro)
	IF (!D.Flags AND 1) NE 0 THEN scalablePixels = 1 ELSE scalablePixels = 0

    
    
	; see how many plots we're working with and size columns/rows accordingly (mainly for colorbar)
	if !p.multi[1] eq 0 then plot_column_width=1. else plot_column_width=1./!p.multi[1]
	if !p.multi[2] eq 0 then plot_row_width=1. else plot_row_width=1./!p.multi[2] ; not used atm
	
	; kludge the plot x-margins based on window size (default is [6,10] in units of character width)
	xmargin_left=6 
	xmargin_right=10
	if (~scalablePixels) then begin
		if !p.charsize eq 0 then !p.charsize=2
		if !d.x_size ne 0 then begin
			xmargin_right=long((!d.x_size^.51)/(2.*sqrt(!p.charsize))) ; yeah, this is somewhat ridiculous... but it works for "reasonable" plot window and character sizes
		endif
	endif

	; set up the y range	
	if keyword_set(yrange) then begin
	  good_range=where(y ge yrange[0] and y le yrange[1])
	  if size(good_range,/n_elements) ge 2 then begin
	    z=z[*,good_range]
	    y=y[good_range]
        cs=cs[*,good_range]
	  endif else print,"Ignoring spec YRANGE: No intersection with data!"
	endif

	; note: !x.window and !y.window contain the normalized coordinates of the last plot pane drawn.
	; these coordinates will be used after calling shade_surf to position the axes, image (for PostScript only), and colorbar

	
	; draw the plot
	if ~keyword_set(no_interp) then begin
		; make a "smooth" spectrogram

		if keyword_set(julian) then shade_surf,z,julian,y,/t3d,shades=cs,xstyle=1,ystyle=1,zstyle=1,xmargin=[xmargin_left,xmargin_right],ymargin=[3,2], $
		           xlog=xlog,ylog=ylog,image=image,/nodata,_extra=extra,xtickformat='LABEL_DATE',xtickunits=['Time']
		if ~keyword_set(julian) then shade_surf,z,x,y,/t3d,shades=cs,xstyle=1,ystyle=1,zstyle=1,xmargin=[xmargin_left,xmargin_right],ymargin=[3,2], $
		           xlog=xlog,ylog=ylog,image=image,yrange=yrange,zrange=zrange,_extra=extra
		if keyword_set(title) then xyouts,(!x.window[1]-!x.window[0])/2.+!x.window[0],!y.window[1]+.025*(!y.window[1]-!y.window[0]),title,/normal,alignment=0.5,charsize=!p.charsize/1.2
		; check scalablePixels for PS output
		if scalablePixels then tv,image,!x.window[0],!y.window[0],xsize=!x.window[1]-!x.window[0],ysize=!y.window[1]-!y.window[0],/normal,/t3d
		
		if keyword_set(hline) then for b=0,n_elements(hline)-1 do oplot,[x[0],x[n_elements(x)-1]],[hline[b],hline[b]],color=1
		if keyword_set(vline) then for b=0,n_elements(vline)-1 do vline,vline[i]


	endif else begin
		; make a block spectrogram using shade_surf to set up the axes
		

		;------------------------
		;SET UP NICER TIME LABELS
		;------------------------

		if keyword_set(julian) then begin			
			dummy = LABEL_DATE(DATE_FORMAT=['%H:%I:%S'])  	
			goo = time_string(x)
			julian = JULDAY(strmid(goo,5,2),strmid(goo,8,2),strmid(goo,0,4),strmid(goo,11,2),strmid(goo,14,2),strmid(goo,17,2))
		endif

		if keyword_set(julian) then shade_surf,z,julian,y,/t3d,shades=cs,xstyle=1,ystyle=1,zstyle=1,xmargin=[xmargin_left,xmargin_right],ymargin=[3,2], $
		           xlog=xlog,ylog=ylog,image=image,/nodata,_extra=extra,xtickformat='LABEL_DATE',xtickunits=['Time']
		if ~keyword_set(julian) then shade_surf,z,x,y,/t3d,shades=cs,xstyle=1,ystyle=1,zstyle=1,xmargin=[xmargin_left,xmargin_right],ymargin=[3,2], $
		           xlog=xlog,ylog=ylog,image=image,/nodata,_extra=extra;,xtickformat="(A1)"


		if keyword_set(title) then xyouts,(!x.window[1]-!x.window[0])/2.+!x.window[0],!y.window[1]+.025*(!y.window[1]-!y.window[0]),title,/normal,alignment=0.5,charsize=!p.charsize/1.2
		nx=size(x,/n_elements)
		ny=size(y,/n_elements)
		for xcount=0.,nx-2 do begin
			for ycount=0.,ny-2 do begin
				if keyword_set(julian) then polyfill,[julian[xcount],julian[xcount+1],julian[xcount+1],julian[xcount]],$
						 [y[ycount],y[ycount],y[ycount+1],y[ycount+1]],$
						 color=cs[xcount,ycount]			

				if ~keyword_set(julian) then polyfill,[x[xcount],x[xcount+1],x[xcount+1],x[xcount]],$
						 [y[ycount],y[ycount],y[ycount+1],y[ycount+1]],$
						 color=cs[xcount,ycount]		
				endfor
		endfor

		if keyword_set(hline) then for b=0,n_elements(hline)-1 do oplot,[julian[0],julian[n_elements(julian)-1]],[hline[b],hline[b]],color=1
		if keyword_set(vline) then for b=0,n_elements(vline)-1 do vline,vline[i]
		
		if keyword_set(hline) then for b=0,n_elements(hline)-1 do oplot,[x[0],x[n_elements(x)-1]],[hline[b],hline[b]],color=1
		if keyword_set(vline) then for b=0,n_elements(vline)-1 do vline,vline[i]
		
	endelse		
	
		
	; draw the axes
	axis,xstyle=1,ystyle=1,xtitle=xtitle,xtickformat='(A1)'
;	axis,xstyle=1,ystyle=1,xtitle=xtitle,xtickformat='LABEL_DATE',xtickunits=['Time']
	axis,yaxis=0,xstyle=1,ystyle=1,ytitle=ytitle
	axis,xaxis=1,xtickname=replicate(' ',30),xstyle=1,ystyle=1                                                 
	axis,yaxis=1,ytickname=replicate(' ',30),xstyle=1,ystyle=1


	; set up vertical color bar positions using !x.window and !y.window
	y0=!y.window[0]
	y1=!y.window[1]
	x0=!x.window[1]+0.01*plot_column_width
	x1=x0+0.025*plot_column_width
	
	
	color_bar_position=[x0,y0,x1,y1]
;	if keyword_set(zrange) then color_bar_range=zrange else color_bar_range=[min(z,/nan),max(z,/nan)]
	if keyword_set(zrange) then color_bar_range=zrange else color_bar_range=[min(cs),max(cs)]
	
	
	if keyword_set(zlog) then begin
		color_bar_range[0]=max([color_bar_range[0],1.e-3])
;		zlogticks=long(color_bar_range[1]-color_bar_range[0])/100.+1
		zlogticks=long(color_bar_range[1]-color_bar_range[0])/10.+1
		zlogticknames=0.
		for zlogtickcount=0,zlogticks-1 do begin
			zlogticknames=[zlogticknames,10^zlogtickcount]
		endfor
	endif
	; draw the colorbar - REQUIRES D. Fanning's colorbar.pro routine!
	; note: the spectrogram ztitle corresponds to the ytitle of the vertical color bar


if ~keyword_set(nocb) then begin
	if ~keyword_set(ztitle) then ztitle=''
	if (zlog) then colorbar,position=[x0,y0,x1,y1],range=color_bar_range,/vertical,/right,ylog=1,divisions=zlogticks,ticknames=zlogticknames,ncolors=255,ytitle=ztitle,minrange=0,maxrange=255.,_extra=extracb $   ;, $
			  else colorbar,position=[x0,y0,x1,y1],range=color_bar_range,/vertical,/right,ncolors=256,ytitle=ztitle,minrange=0,maxrange=255,_extra=extracb  
endif			  
			  
end

