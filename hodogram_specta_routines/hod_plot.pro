;+
;*****************************************************************************************
;
;  PROCEDURE : hod_plot.pro
;  PURPOSE  :  Plots hodograms and "slow" hodograms along w/ a movie of a waveform
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
;  INPUT:		wf -> [*,3] input waveform
;				Bo -> three element magnetic field vector to overplot
;				slow -> keyword to plot the slow hodograms
;				wait -> the wait time b/t subsequent plots (0.01 secs by default).
;						Set to zero for minimum wait time.
;				start -> The initial array element of waveform to plot
;				final -> final array element of waveform to plot
;               ps -> create ps file (saved as hod_plot.ps to desktop)
;				chunksz -> the number of points to keep on the screen while
;					plotting the slow hodograms. Defaults to 10
;				plotnth -> plot every nth point
;
;  EXAMPLES:
;
;
;
;   CHANGED:  1)  NA [MM/DD/YYYY   v1.0.0]
;
;   NOTES:
;
;
;   CREATED:  08/02/2010
;   CREATED BY:  Aaron W. Breneman
;    LAST MODIFIED:  MM/DD/YYYY   v1.0.0
;    MODIFIED BY: Aaron W. Breneman  07/27/2011 -> added chunksz functionality.
;										-> added a movie of the location of hodogram over waveform
;									 08/09/2011 -> added the plotnth keyword
;												   added * symbol for current hod point
;
;*****************************************************************************************
;-

pro hod_plot,wf,Bo=Bo,slow=slow,wait=wait,start=s,final=e,ps=ps,chunksz=chunksz,plotnth=plotnth

	;	window,1,xsize=250,ysize=1500
	;	window,2,xsize=1200,ysize=350
	;	window,3,xsize=1200,ysize=350

	!p.charsize=1.2

	xl = 0.25
	xr = 0.95
	y1 = 0.
	y2 = 0.33
	y3 = 0.66
	y4 = 1.
	offset = 0.05



	if ~keyword_set(chunksz) then chunksz=10.
	if ~keyword_set(s) then s=0.
	if ~keyword_set(e) then e = n_elements(wf[*,0])-1.
	if ~keyword_set(wait) then wait = 0.01
	if ~keyword_set(plotnth) then plotnth=1   ;plot every point
	if keyword_set(slow) then print,"****ENTER 's' TO STOP SLOW HODOGRAM****"


	if keyword_set(ps) then begin
		window_hod = {xsize:3.0,ysize:8.,xoffset:0.425,yoffset:1.68,inches:1,portrait:1,landscape:0}
		set_plot,'ps'
		!p.font=0
		device,filename='/hod_plot.ps',decomposed=0,encapsulated=2,font_size=12,/helvetica,/color,_extra=window_hod
	endif


	device,decomposed=0
	loadct,39

	nelem = e-s

	;-----------------------
	;x vs y
	;-----------------------
	hrange1 = max(wf[s:e,0],/nan) > max(wf[s:e,1],/nan)

	;plot hodogram skeleton
	plot1 = PLOT([0,0],xtitle='X',ytitle='Y',xstyle=1,ystyle=1,xrange=[-hrange1,hrange1],yrange=[-hrange1,hrange1],/NODATA,$
	position=[xl,y3+offset,xr,y4-offset])




	;plot waveform
	xvals = reform(indgen(n_elements(wf[*,0])))
	plot2 = PLOT(xvals,wf[*,0],color=[0,0,50],TITLE="Wave")
	plot2 = PLOT(xvals,wf[*,1], /overplot,color=[0,150,0])

;	;plot wave magnitude
;	plot3 = plot(xvals,wf[*,0]^2 + wf[*,1]^2 + wf[*,2]^2)

	if keyword_set(slow) then begin
		stop_slow_plot = ''
		if keyword_set(Bo) then begin
			;oplot,[0,Bo[0]],[0,Bo[1]],color=0,thick=2
			;oplot,[0],[0],color=0,thick=2,psym=7
		endif


		;		plots,wf[s,0],wf[s,1] ;,color=chcolors[0]

		scc=s+1
		while scc le (nelem+s)-1 and stop_slow_plot ne 's' do begin
			plot1.window.SetCurrent

			;			if scc le chunksz+1 then plots,wf[scc,0],wf[scc,1],/continue else begin
	;		if scc le chunksz+1 then print, 'test' else begin ;plot1 = plot([wf[scc,0],wf[scc,1]],/overplot)
;			endif else begin

				if scc mod plotnth eq 0 then begin
					;plot an entire hodogram block of size "chunksz"
					plot1.window.SetCurrent

					;erase previous chunk
					if scc gt 2*chunksz then plot1 = plot(wf[scc-2*chunksz:scc-chunksz,0],wf[scc-2*chunksz:scc-chunksz,1],/overplot,color=[255,255,255])
					if scc gt 2*plotnth then plot1 = plot([wf[scc-2*plotnth,0],0],[wf[scc-2*plotnth,1],0],/overplot,symbol="*",color=[255,255,255])
					if scc gt 2*plotnth then plot1 = plot([wf[0,0],0],[wf[0,1],0],/overplot,symbol="*",color=[255,255,255])
					;plot new chunk
					if scc gt chunksz then plot1 = plot(wf[scc-chunksz:scc,0],wf[scc-chunksz:scc,1],/overplot)
					if scc gt chunksz then plot1 = plot([wf[scc,0],0],[wf[scc,1],0],/overplot,symbol="*")

;stop
					plot2.window.SetCurrent

					;replot previous chunk to be original color
					if scc gt 2*chunksz then plot2 = PLOT(xvals[scc-2*chunksz:scc-chunksz],wf[scc-2*chunksz:scc-chunksz,0],/overplot,color=[0,0,50])
					if scc gt 2*chunksz then plot2 = PLOT(xvals[scc-2*chunksz:scc-chunksz],wf[scc-2*chunksz:scc-chunksz,1],/overplot,color=[0,150,0])
					;shade current chunk
					if scc gt 2*chunksz then plot2 = PLOT(xvals[scc-chunksz:scc],wf[scc-chunksz:scc,0],/overplot,color=[100,0,240])
;					stop


					if wait ne 0 then wait,wait
				endif
			;endelse

			scc+=1
			stop_slow_plot=get_kbrd(0)
		endwhile
		stop_slow_plot=''
	endif else begin
		wset,1
		oplot,wf[s:e,0],wf[s:e,1] ;,color=chcolors[0]
	endelse

	if keyword_set(Bo) then begin
		wset,1
		oplot,[0,Bo[0]],[0,Bo[1]],color=0,thick=2
		oplot,[0],[0],color=0,thick=2,psym=7
	endif


	;---------------------
	;x vs z
	;---------------------
	hrange2 = max(wf[s:e,0],/nan) > max(wf[s:e,2],/nan)

	wset,1
	plot,[0,0],xtitle='X',ytitle='Z',$
	xstyle=1,ystyle=1,xrange=[-hrange2,hrange2],yrange=[-hrange2,hrange2],/NODATA,$
	position=[xl,y2+offset,xr,y3-offset],/noerase

	wset,2
	plot,wf[*,0],/nodata,title='wave'
	oplot,wf[*,0],color=50
	oplot,wf[*,2],color=250

	wset,3
	plot,sqrt(wf[*,0]^2 + wf[*,1]^2 + wf[*,2]^2)

	if keyword_set(slow) then begin
		stop_slow_plot = ''
		if keyword_set(Bo) then begin
			oplot,[0,Bo[0]],[0,Bo[2]],color=0,thick=2
			oplot,[0],[0],color=0,thick=2,psym=7
		endif
		wset,1
		plots,wf[s,0],wf[s,2] ;,color=chcolors[0]
		scc=s+1
		while scc le (nelem+s)-1 and stop_slow_plot ne 's' do begin
			if scc le chunksz+1 then plots,wf[scc,0],wf[scc,2],/continue else begin
				wset,1

				if scc mod plotnth eq 0 then begin
					plot,wf[scc-chunksz:scc,0],wf[scc-chunksz:scc,2],xtitle='X',ytitle='Z',$
					xstyle=1,ystyle=1,xrange=[-hrange2,hrange2],yrange=[-hrange2,hrange2],$
					position=[xl,y2+offset,xr,y3-offset]
					oplot,[wf[scc,0],0],[wf[scc,2],0],psym=2

					wset,2
					plot,wf[*,0],/nodata,title='wave: Blue=X, Red=Z'
					oplot,wf[*,0],color=50
					oplot,wf[*,2],color=250
					vline,scc-chunksz
					vline,scc
					wset,3
					plot,sqrt(wf[*,0]^2 + wf[*,1]^2 + wf[*,2]^2),title='modulus'
					vline,scc-chunksz
					vline,scc
					if wait ne 0 then wait,wait
				endif
			endelse

			scc+=1
			stop_slow_plot=get_kbrd(0)
		endwhile
		stop_slow_plot=''
	endif else begin
		wset,1
		oplot,wf[s:e,0],wf[s:e,2] ;,color=chcolors[0]
	endelse

	if keyword_set(Bo) then begin
		wset,1
		oplot,[0,Bo[0]],[0,Bo[2]],color=0,thick=2
		oplot,[0],[0],color=0,thick=2,psym=7
	endif


	;--------------------
	;y vs z
	;--------------------
	hrange3 = max(wf[s:e,1],/nan) > max(wf[s:e,2],/nan)

	wset,1
	plot,[0,0],xtitle='Y',ytitle='Z',$
	xstyle=1,ystyle=1,xrange=[-hrange3,hrange3],yrange=[-hrange3,hrange3],/NODATA,$
	position=[xl,y1+offset,xr,y2-offset],/noerase
	wset,2
	plot,wf[*,0],/nodata,title='wave'
	oplot,wf[*,1],color=150
	oplot,wf[*,2],color=250
	wset,3
	plot,sqrt(wf[*,0]^2 + wf[*,1]^2 + wf[*,2]^2)

	if keyword_set(slow) then begin
		stop_slow_plot = ''
		if keyword_set(Bo) then begin
			oplot,[0,Bo[1]],[0,Bo[2]],color=0,thick=2
			oplot,[0],[0],color=0,thick=2,psym=7
		endif
		wset,1
		plots,wf[s,1],wf[s,2] ;,color=chcolors[0]
		scc=s+1
		while scc le (nelem+s)-1 and stop_slow_plot ne 's' do begin
			if scc le chunksz+1 then plots,wf[scc,1],wf[scc,2],/continue else begin
				wset,1

				if scc mod plotnth eq 0 then begin
					plot,wf[scc-chunksz:scc,1],wf[scc-chunksz:scc,2],xtitle='Y',ytitle='Z',$
					xstyle=1,ystyle=1,xrange=[-hrange3,hrange3],yrange=[-hrange3,hrange3],$
					position=[xl,y1+offset,xr,y2-offset]
					oplot,[wf[scc,1],0],[wf[scc,2],0],psym=2

					wset,2
					plot,wf[*,0],/nodata,title='Green=Y, Red=Z'
					oplot,wf[*,1],color=150
					oplot,wf[*,2],color=250
					vline,scc-chunksz
					vline,scc
					wset,3
					plot,sqrt(wf[*,0]^2 + wf[*,1]^2 + wf[*,2]^2),title='modulus'
					vline,scc-chunksz
					vline,scc
					if wait ne 0 then wait,wait
				endif
			endelse

			scc+=1
			stop_slow_plot=get_kbrd(0)
		endwhile
		stop_slow_plot=''
	endif else begin
		wset,1
		oplot,wf[s:e,1],wf[s:e,2] ;,color=chcolors[0]
	endelse

	if keyword_set(Bo) then begin
		wset,1
		oplot,[0,Bo[1]],[0,Bo[2]],color=0,thick=2
		oplot,[0],[0],color=0,thick=2,psym=7
	endif



	;---------------------------------------------------------------------------------
	;Finish up
	;---------------------------------------------------------------------------------

	wset,1
	plot,wf[s:e,0],wf[s:e,1],xtitle='X',ytitle='Y',$
	xstyle=1,ystyle=1,xrange=[-hrange1,hrange1],yrange=[-hrange1,hrange1],$
	position=[xl,y3+offset,xr,y4-offset],/noerase

	if keyword_set(Bo) then begin
		oplot,[0,Bo[0]],[0,Bo[1]],color=0,thick=2
		oplot,[0],[0],color=0,thick=2,psym=7
	endif

	plot,wf[s:e,0],wf[s:e,2],xtitle='X',ytitle='Z',$
	xstyle=1,ystyle=1,xrange=[-hrange2,hrange2],yrange=[-hrange2,hrange2],$
	position=[xl,y2+offset,xr,y3-offset],/noerase

	if keyword_set(Bo) then begin
		oplot,[0,Bo[0]],[0,Bo[2]],color=0,thick=2
		oplot,[0],[0],color=0,thick=2,psym=7
	endif

	plot,wf[s:e,1],wf[s:e,2],xtitle='Y',ytitle='Z',$
	xstyle=1,ystyle=1,xrange=[-hrange3,hrange3],yrange=[-hrange3,hrange3],$
	position=[xl,y1+offset,xr,y2-offset],/noerase

	if keyword_set(Bo) then begin
		oplot,[0,Bo[1]],[0,Bo[2]],color=0,thick=2
		oplot,[0],[0],color=0,thick=2,psym=7
	endif



	wset,2
	plot,wf[*,0],/nodata,title='wave'
	oplot,wf[*,0],color=50
	oplot,wf[*,1],color=150
	oplot,wf[*,2],color=250

	wset,3
	plot,sqrt(wf[*,0]^2 + wf[*,1]^2 + wf[*,2]^2)




	if keyword_set(ps) then begin
		device,/close
		set_plot,'x'
	endif

end
