;+
;*****************************************************************************************
;
;  PROCEDURE : animate_hodogram.pro  
;  PURPOSE  :  Creates an MPEG movie of "slow" hodograms (also see hod_plot.pro)
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
;				start -> The initial array element of waveform to plot
;				final -> final array element of waveform to plot
;				chunksz -> the number of points to keep on the screen while
;					plotting the slow hodograms. Defaults to 10
;				framerate -> defaults to 60 frames/sec
;				plotnth -> plot every nth data point
;
;  EXAMPLES:    
;               
;
;
;   CHANGED:  1)  NA [MM/DD/YYYY   v1.0.0]
;
;   NOTES:      To find out if you're using good values of chunksz and plotnth variables
;				test with hod_plot.pro
;
;				For the non MPEG version see hod_plot.pro
;      
;				Forked from K. Kersten's animate_fft.pro
;               
;
;   CREATED:  08/02/2011
;   CREATED BY:  Aaron W. Breneman
;    LAST MODIFIED:  MM/DD/YYYY   v1.0.0
;    MODIFIED BY:  Aaron W. Breneman 
;						08/09/2011 -> added plotnth keyword
;
;*****************************************************************************************
;-

pro animate_hodogram,wf,Bo=Bo,start=s,final=e,chunksz=chunksz,framerate=framerate,plotnth=plotnth


if ~keyword_set(framerate) then framerate=60
start_time=systime(1)

; set up a randomly named image directory
imgdir='shot.'+strcompress(string(randomu(systime(1),/long)),/remove_all)
spawn,'mkdir '+imgdir

; figure out what to call the movie
if ~keyword_set(moviename) then begin
	moviename=''
	print,'Movie filename not set on input (set with moviename="example" keyword.)'
	print,'Please enter the base filename.' 
	read,moviename,prompt='(.mov extension will be added): '
endif

; set up the window or Z buffer
if keyword_set(fullseries) then begin
	!p.multi=[0,0,2]
	if ~keyword_set(xsize) then xsize=230
	if ~keyword_set(ysize) then ysize=550
endif else begin
	!p.multi=[0,0,2]
	if ~keyword_set(xsize) then xsize=230
	if ~keyword_set(ysize) then ysize=550
endelse

; plot to screen or z buffer (default)
if keyword_set(screen) then begin
	set_plot,'X'
	!p.charsize=0.8
	window,xsize=xsize,ysize=ysize
endif else begin
	set_plot,'Z'
	device,set_pixel_depth=24,set_resolution=[xsize,ysize]
	!p.charsize=0.65
endelse


device,decomposed=0
loadct,39,/silent

;plot position values
xl = 0.2
xr = 0.92
y1 = 0.2
y2 = 0.5
y3 = 0.6
y4 = 0.9


if ~keyword_set(plotnth) then plotnth=1
if ~keyword_set(chunksz) then chunksz=10.
if ~keyword_set(s) then s=0.
if ~keyword_set(e) then e = n_elements(wf[*,0])-1.
print,"****ENTER 's' TO STOP SLOW HODOGRAM****"

e = double(e)
s = double(s)
nelem = e-s

;-----------------------
;x vs y
;-----------------------

hrange1 = max(wf[s:e,0],/nan) > max(wf[s:e,1],/nan)

stop_slow_plot = ''
scc=s+1
counter = 0
while scc le (nelem+s)-1 and stop_slow_plot ne 's' do begin

  if scc mod 100 eq 0 then print,'...still working...',scc,' out of',nelem	
  if scc gt chunksz+1 then begin


	if scc mod plotnth eq 0 then begin
		plot,wf[scc-chunksz:scc,0],wf[scc-chunksz:scc,1],xtitle='Ex',ytitle='Ey',$
			xstyle=1,ystyle=1,xrange=[-hrange1,hrange1],yrange=[-hrange1,hrange1],$
			position=[xl,y3,xr,y4]
		oplot,[wf[scc,0],0],[wf[scc,1],0],psym=2
	
		if keyword_set(Bo) then begin
			oplot,[0,Bo[0]],[0,Bo[1]],color=0,thick=2
			oplot,[0],[0],color=0,thick=2,psym=7
		endif
	
		plot,wf[*,0],/nodata,title='wave: Blue=Ex, Green=Ey',position=[xl,y1,xr,y2],xstyle=1
		oplot,wf[*,0],color=50
		oplot,wf[*,1],color=150
		vline,scc-chunksz
		vline,scc	
		counter++
		sccstr = string(counter,format='(i04)')
		write_png,imgdir+'/'+sccstr+'.png',tvrd(true=1)
	endif


  endif

  scc+=1
  stop_slow_plot=get_kbrd(0)
endwhile
stop_slow_plot=''



print,'Encoding movie: '+moviename
encode_command='ffmpeg -r '+ strtrim(fix(framerate),2) + ' -i '+imgdir+'/%04d.png -qscale 2 '+moviename + '_Ex_Ey' + '.mov'

spawn,encode_command
print,''
print,'Movie complete: '+moviename + '_Ex_Ey' + '.mov'


print,'Cleaning up frame images...'
rmcommand='rm -rf '+imgdir
spawn,rmcommand
print,'Done!'

end_time=systime(1)

execution_time=end_time-start_time
print,''
print,'EXECUTION TIME: '+strcompress(string(execution_time),/remove_all)+' sec'
if execution_time gt 60. then print,'                      ( '+strcompress(string(long(execution_time)/60),/remove_all)+' min '+strcompress(string(long(execution_time) mod 60),/remove_all)+' sec )'
print,''






;-----------------------
;x vs z
;-----------------------

imgdir='shot.'+strcompress(string(randomu(systime(1),/long)),/remove_all)
spawn,'mkdir '+imgdir


hrange1 = max(wf[s:e,0],/nan) > max(wf[s:e,2],/nan)

stop_slow_plot = ''
scc=s+1
counter = 0
while scc le (nelem+s)-1 and stop_slow_plot ne 's' do begin

  if scc mod 100 eq 0 then print,'...still working...',scc,' out of',nelem
  if scc gt chunksz+1 then begin

	if scc mod plotnth eq 0 then begin

		plot,wf[scc-chunksz:scc,0],wf[scc-chunksz:scc,2],xtitle='Ex',ytitle='Ez',$
			xstyle=1,ystyle=1,xrange=[-hrange1,hrange1],yrange=[-hrange1,hrange1],$
			position=[xl,y3,xr,y4]
		oplot,[wf[scc,0],0],[wf[scc,1],2],psym=2
	
		if keyword_set(Bo) then begin
			oplot,[0,Bo[0]],[0,Bo[2]],color=0,thick=2
			oplot,[0],[0],color=0,thick=2,psym=7
		endif
	
		plot,wf[*,0],/nodata,title='wave: Blue=Ex, Green=Ez',position=[xl,y1,xr,y2],xstyle=1
		oplot,wf[*,0],color=50
		oplot,wf[*,2],color=150
		vline,scc-chunksz
		vline,scc	
		counter++
		sccstr = string(counter,format='(i04)')
		write_png,imgdir+'/'+sccstr+'.png',tvrd(true=1)
	endif

  endif

  scc+=1
  stop_slow_plot=get_kbrd(0)
endwhile
stop_slow_plot=''


print,'Encoding movie: '+moviename
encode_command='ffmpeg -r '+ strtrim(fix(framerate),2) + ' -i '+imgdir+'/%04d.png -qscale 2 '+moviename + '_Ex_Ez' + '.mov'

spawn,encode_command
print,''
print,'Movie complete: '+moviename + '_Ex_Ez' + '.mov'


print,'Cleaning up frame images...'
rmcommand='rm -rf '+imgdir
spawn,rmcommand
print,'Done!'

end_time=systime(1)

execution_time=end_time-start_time
print,''
print,'EXECUTION TIME: '+strcompress(string(execution_time),/remove_all)+' sec'
if execution_time gt 60. then print,'                      ( '+strcompress(string(long(execution_time)/60),/remove_all)+' min '+strcompress(string(long(execution_time) mod 60),/remove_all)+' sec )'
print,''




;-----------------------
;y vs z
;-----------------------

imgdir='shot.'+strcompress(string(randomu(systime(1),/long)),/remove_all)
spawn,'mkdir '+imgdir


hrange1 = max(wf[s:e,1],/nan) > max(wf[s:e,2],/nan)

stop_slow_plot = ''
scc=s+1
counter = 0
while scc le (nelem+s)-1 and stop_slow_plot ne 's' do begin

  if scc mod 100 eq 0 then print,'...still working...',scc,' out of',nelem
  if scc gt chunksz+1 then begin

	if scc mod plotnth eq 0 then begin

		plot,wf[scc-chunksz:scc,1],wf[scc-chunksz:scc,2],xtitle='Ey',ytitle='Ez',$
			xstyle=1,ystyle=1,xrange=[-hrange1,hrange1],yrange=[-hrange1,hrange1],$
			position=[xl,y3,xr,y4]
		oplot,[wf[scc,1],0],[wf[scc,2],0],psym=2
	
		if keyword_set(Bo) then begin
			oplot,[0,Bo[1]],[0,Bo[2]],color=0,thick=2
			oplot,[0],[0],color=0,thick=2,psym=7
		endif
	
		plot,wf[*,0],/nodata,title='wave: Blue=Ey, Green=Ez',position=[xl,y1,xr,y2],xstyle=1
		oplot,wf[*,1],color=50
		oplot,wf[*,2],color=150
		vline,scc-chunksz
		vline,scc	
		counter++
		sccstr = string(counter,format='(i04)')
		write_png,imgdir+'/'+sccstr+'.png',tvrd(true=1)
	endif

  endif

  scc+=1
  stop_slow_plot=get_kbrd(0)
endwhile
stop_slow_plot=''



print,'Encoding movie: '+moviename
encode_command='ffmpeg -r '+ strtrim(fix(framerate),2) + ' -i '+imgdir+'/%04d.png -qscale 2 '+moviename + '_Ey_Ez' + '.mov'
spawn,encode_command
print,''
print,'Movie complete: '+moviename + '_Ey_Ez' + '.mov'


print,'Cleaning up frame images...'
rmcommand='rm -rf '+imgdir
spawn,rmcommand
print,'Done!'

end_time=systime(1)

execution_time=end_time-start_time
print,''
print,'EXECUTION TIME: '+strcompress(string(execution_time),/remove_all)+' sec'
if execution_time gt 60. then print,'                      ( '+strcompress(string(long(execution_time)/60),/remove_all)+' min '+strcompress(string(long(execution_time) mod 60),/remove_all)+' sec )'
print,''



device,/close
set_plot,'x'

end
