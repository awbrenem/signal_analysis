"""
Functions

plot_scatter [NOT FINISHED] --> make scatter plot with colors for the dots
plot_spectrogram --> plot (pre-existing) spectrogram nicely
slice_spectrogram --> extract a slice (with option of time averaging)
"""


def plot_scatter():
    fig, axs = plt.subplots(2)
    logv = np.log(data['SpecBMax_lb'])

    p1 = axs[0].scatter(data['MLTrb'],data['Lrb'],c=data['EMFb'],vmin=0,vmax=np.max(data['EMFb']))
    p2 = axs[1].scatter(data['MLTrb'],data['Lrb'],c=logv,vmin=np.min(logv),vmax=np.max(logv))
    axs[0].set_ylabel('EMFISIS burst data\nsec\nL-shell')
    axs[1].set_ylabel('chorus lower band amp\nlog10(pT^2/Hz)\nL-shell')
    axs[0].set_xlabel('MLT')
    axs[1].set_xlabel('MLT')
    plt.colorbar(p1)
    plt.colorbar(p2)





def slice_spectrogram(tslice, tspec, spec, nsec=0):
    """
    Extract a slice out of a spectrogram with possible time averaging

    tslice = time of slice 
    tspec --> spectrogram time values 
    spec --> spectrogram power(amplitude) values
    nsec --> number of seconds to average over. Defaults to zero, which means no averaging over time

    Returns: (powavg) The power (for all freqs) at requested time averaged over "navg" time bins.
             (powarr) Also returns individual slices for each time bin used in average
             (tarr) The final times that lie b/t tslice and tslice + nsec

    NOTE: take the abs value of complex spectra before passing
    """

    import numpy as np

    #data cadence
    cadence = np.median(tspec - np.roll(tspec,1))
   
    #number of data points to average
    navg = int(nsec/cadence)
    if navg == 0: navg = 1

    #starting point
    goo = np.where(tspec >= tslice)

    #extract data from start to end point and average
    #powgoo = np.abs(spec)[:,goo[0][0]:goo[0][navg]]
    powgoo = spec[:,goo[0][0]:goo[0][navg]]

    powavg = [0] * np.shape(powgoo)[0]

    for i in range(np.shape(powgoo)[0]):
        powavg[i] = np.nanmean(powgoo[i,:])    
        #powavg[i] = np.average(powgoo[i,:])    

    #powarr = np.abs(spec)[:,goo[0][0:navg]]
    powarr = spec[:,goo[0][0:navg]]
    tarr = tspec[goo[0][0:navg]]

    return powavg, powarr, tarr




def plot_spectrogram(t,f,p,
                     vr=[0,1],xr=0, yr=0, 
                     minzval=0, maxzval=0,
                     yscale='linear', zscale='log', 
                     ax=False, show=True, invert=False, 
                     title='', xlabel='', ylabel='', 
                     plot_kwargs={'cmap':'turbo'},
                     plot_kwargs2={'origin':'lower','alpha':1,'interpolation':'nearest','aspect':'auto'},
                     colorbar_kwargs={},
                     colorbar=1,
                     xaxis2=0, bc=[0,0,0]):


    """
    Plot large-array spectograms quickly in Python using imshow. 
    This method is so much faster than pcolormesh, which can take forever. 

    NOTE: make sure the x-data are regularly gridded. This code will still work, but the FFT routine that 
    created the data will spit out incorrect timing if not. 

    E.g. 
    First take the FFT
        fs = 1/(np.mean(vlf12.tsec - vlf12.tsec.shift()))
        freq12, tspec12, power12 = signal.spectrogram(vlf12.amp, fs, nperseg=512, return_onesided=1)

    Now plot 
        plotv = {'vmin':0.5, 'vmax':0.8, 'xlim':[600,800], 'ylim':[6000,8000], 'yscale':'linear'}
        plot_spectrogram(tspec12, freq12, power12, plotv)

    Note: this can be a standalone plotting program, or it can plot alongside other plots
    by inputting the ax keyword from "fig, ax = plt.subplots()"

    Input:
        t -> array of x-axis values (e.g. times)
        f -> array of y-axis values (e.g. freqs)
        p -> [t,f] sized array of image values (e.g. sliding spectra from scipy.signal.spectrogram)
        yscale -> linear or log 
        zscale -> linear or log 
            linear = don't modify input spectral values
            log = take dB=10*log(p) of input spectral values
            NOTE: scipy.signal.spectrogram will return units**2/Hz (density) or units**2 (spectrum) for an input waveform of "units".
        xr -> xrange 
        yr -> yrange 
        ax -> from matplotlib's subplots (e.g. fig, ax = plt.subplots()) - only set if you want this plot to be a part of other plots.
        minzval -> any values less than this will be assigned NaN
        maxzval -> any values greater than this will be assigned NaN
        bc --> (bottom color). Set this as the lowest color value (RGB). Defaults to black
        
                        
            ---------------------------------------------
            Example using two subplots
                fig, ax = plt.subplots(2)
                ax[0].set_yscale('log')
                ax[0].set_ylim(100,12000)
                ps.plot_spectrogram(t,f,p1,ax=ax[0])
                ps.plot_spectrogram(t,f,p2,ax=ax[1])

            #NOTE: I've had trouble getting a y-log scale set using subplots when setting directly with this routine. 
            #It's best to set them outside of the plot_spectrogram call. Also note that Jupyter will sometimes not display 
            #multiple subplots inline. To display in a popup window use
            #   %matplotlib inline
            #   %matplotlib qt   

            --------------------------------------------- 

        invert -> flip the yscale upside down if set to True
    """


    import matplotlib.pyplot as plt
    import matplotlib
    import numpy as np
    from matplotlib import cm
    from matplotlib.colors import ListedColormap, LinearSegmentedColormap
    
    

    #Merge the two keyword dictionaries
    plot_kwargs.update(plot_kwargs2)


    #Modify the color bar so the lowest value is black
    cmapgoo = cm.get_cmap('turbo', 256)
    newcolors = cmapgoo(np.linspace(0, 1, 256))
    #redefine RGB colors
    bottom = np.array([bc[0],bc[1],bc[2], 1])
    #topval = np.array([255/256, 0, 125/256, 1]) #raspberry
    #topval = np.array([255/256, 0, 255/256, 1])
    newcolors[0, :] = bottom
    #newcolors[-1,:] = topval
    newcmp = ListedColormap(newcolors)
    plot_kwargs['cmap'] = newcmp



    #Imshow requires accurate specification of plot limits or else it will show x/y axes incorrectly. 
    #(1) if xr, yr are not set, then set them to the limits of the data
    #(2) if xr, yr are set, then:
    #   (a) if xr, yr are within the limits of the min/max xr and min/max yr, then reduce the arrays to be plotted
    #   (b) if xr, yr are outside the limits of the min/max xr and min/max yr, then plot the limits of the data. 


    if not xr: 
        xrplot = [np.min(t),np.max(t)]
    else: 
        xrplot = [0,0]
        xrplot[0] = max(np.min(t), xr[0])
        xrplot[1] = min(np.max(t), xr[1])
    if not yr: 
        yrplot = [np.min(f),np.max(f)]
    else: 
        yrplot = [0,0]
        yrplot[0] = max(np.min(f), yr[0])
        yrplot[1] = min(np.max(f), yr[1])




    #Subselect data within the desired timerange and frequency range
    #(required for plotting with imshow, which doesn't have a time base input)
    goodindt = np.where((t >= xrplot[0]) & (t <= xrplot[1]))[0]
    goodindf = np.where((f >= yrplot[0]) & (f <= yrplot[1]))[0]
    p2 = p[:,goodindt]
    p3 = p2[goodindf,:]




    if zscale == 'log':
        pn = 10.*np.log10(p3)
    else: 
        pn = p3.copy()


    if minzval != 0:
        pn[pn < minzval] = float("nan")
    if maxzval != 0:
        pn[pn > maxzval] = float("nan")


    if not ax:
        fig, ax = plt.subplots()


    ##invert image or not
    #origin='lower'
    #if invert:
    #    origin='upper'


    #Plot the spectrogram. Note that we need to force x,y axes to be defined exactly based on the range of data, otherwise image 
    #won't be displayed properly. 

    im = ax.imshow(pn,vmin=vr[0],vmax=vr[1],
                   extent=[xrplot[0],xrplot[1],yrplot[0],yrplot[1]],
                   **plot_kwargs)


    #Code to add top x-axis
    #if xaxis2:
    #    ax2 = ax.twiny()
    #    ax2.plot([11, 12, 31, 41, 15], [13, 51, 17, 11, 76], color='blue')


    #colorbar_kwargs = {'min':-180,'max':180}
    if colorbar == 1:
        plt.colorbar(im,ax=ax, **colorbar_kwargs)

#cbar.set_label('# of contacts', rotation=270)



    if title != '':
        ax.set_title(title)
    if xlabel != '':
        ax.set_xlabel(xlabel)
    if ylabel != '':
        ax.set_ylabel(ylabel)

    ax.set_yscale(yscale)    
    ax.set_xlim(xr)
    ax.set_ylim(yr)

    if show == True:
        plt.show()

    return pn
