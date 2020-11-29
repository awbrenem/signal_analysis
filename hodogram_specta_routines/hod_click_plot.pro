;plot hodograms by clicking on a start and endpoint on a waveform. 
;"waveform" is any [n,3] waveform. 

;Code will ask you to click on start and stop times for 
;each hodogram. Click on as many as you want. A hodogram will
;be created for each click pair (start/stop times). 



pro hod_click_plot,waveform 

    split_vec,waveform 

    get_data,waveform,t,d
    trangefull = [min(t),max(t)]

    ctime,t,v

    for j=0,n_elements(t) - 1 do begin 
        x = tsample(waveform,t[j:j+1],time=tt)
        store_data,waveform+'_tmp',tt,x
        plot_wavestuff,waveform+'_tmp',/hod,/nowindow

        tplot,waveform+'_'+['x','y','z'],trange=trangefull
        timebar,[t[j],t[j+1]],thick=2
        j++
        stop
    endfor

end
