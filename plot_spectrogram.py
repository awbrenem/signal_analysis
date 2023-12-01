"""
Functions

plot_spectrogram --> plot (pre-existing) spectrogram nicely
slice_spectrogram --> extract a slice (with option of time averaging)
"""



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
        powavg[i] = np.average(powgoo[i,:])    

    #powarr = np.abs(spec)[:,goo[0][0:navg]]
    powarr = spec[:,goo[0][0:navg]]
    tarr = tspec[goo[0][0:navg]]

    return powavg, powarr, tarr




def plot_spectrogram(t,f,p,
                     vr=[0,1], xr=0, yr=0, 
                     yscale='linear', zscale='log', 
                     ax=False, show=True, invert=False, 
                     title='', xlabel='', ylabel='', 
                     minzval=0, maxzval=0, cmap='turbo',
                     alpha=1):

    """
    Plot large-array spectograms quickly in Python using imshow. 
    This method is so much faster than pcolormesh, which can take forever. 

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
            NOTE: scipy.signal.spectrogram will return V**2/Hz (density) or V**2 (spectrum) for an input waveform of V.
        xr -> xrange 
        yr -> yrange 
        ax -> from matplotlib's subplots (e.g. fig, ax = plt.subplots()) - only set if you want this plot to be a part of other plots.
        minzval -> any values less than this will be assigned NaN
        maxzval -> any values greater than this will be assigned NaN
        
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


    if not xr: xr = [np.min(t),np.max(t)]
    if not yr: yr = [np.min(f),np.max(f)]


    if zscale == 'log':
        pn = 10.*np.log10(p)
    else: 
        pn = p.copy()


    if minzval != 0:
        pn[pn < minzval] = "nan"
    if maxzval != 0:
        pn[pn > maxzval] = "nan"


    if not ax:
        fig, ax = plt.subplots()


    #invert image or not
    origin='lower'
    if invert:
        origin='upper'


    #Plot the spectrogram. Note that we need to force x,y axes to be defined exactly based on the range of data, otherwise image 
    #won't be displayed properly. 

    #plt.set_cmap('turbo')
    #plt.set_cmap('terrain')
    #current_cmap = matplotlib.cm.get_cmap()
    #current_cmap.set_bad(color='red')

    #im = ax.imshow(pn,vmin=vr[0],vmax=vr[1],cmap='turbo',aspect='auto', extent=[np.min(t),np.max(t),np.min(f),np.max(f)], origin='lower', **plotkwargs)
    im = ax.imshow(pn,vmin=vr[0],vmax=vr[1],cmap=cmap,aspect='auto',
                   extent=[np.min(t),np.max(t),np.min(f),np.max(f)],
                   origin=origin,alpha=alpha)
    #im = ax.imshow(pn,vmin=vr[0],vmax=vr[1],aspect='auto', extent=[np.min(t),np.max(t),np.min(f),np.max(f)], origin='lower', **plotkwargs)

    if title != '':
        ax.set_title(title)
    if xlabel != '':
        ax.set_xlabel(xlabel)
    if ylabel != '':
        ax.set_ylabel(ylabel)

    #Now we can change to desired xrange and yrange
    plt.yscale(yscale)
    ax.set_ylim(yr[0],yr[1])
    ax.set_xlim(xr[0],xr[1])

    if show == True:
        plt.show()

    return pn
