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

        Example using two subplots
            fig, ax = plt.subplots(2)
            ps.plot_spectrogram(t,f,p1,ax=ax[0])
            ps.plot_spectrogram(t,f,p2,ax=ax[1])

    invert -> flip the yscale upside down if set to True
"""


def plot_spectrogram(t,f,p,vr=[0,1], yscale='linear', zscale='log', xr=0, yr=0, ax=False, show=True, invert=False, title='spec'):

    import matplotlib.pyplot as plt
    import matplotlib
    import numpy as np


    if not xr: xr = [np.min(t),np.max(t)]
    if not yr: yr = [np.min(f),np.max(f)]


    if zscale == 'log':
        pn = 10.*np.log10(p)
    else: 
        pn = p.copy()



    if not ax:
        fig, ax = plt.subplots()


    #invert image or not
    origin='lower'
    if invert:
        origin='upper'


    #Plot the spectrogram. Note that we need to force x,y axes to be defined exactly based on the range of data, otherwise image 
    #won't be displayed properly. 

    #plt.set_cmap('turbo')
    #current_cmap = matplotlib.cm.get_cmap()
    #current_cmap.set_bad(color='red')

    #im = ax.imshow(pn,vmin=vr[0],vmax=vr[1],cmap='turbo',aspect='auto', extent=[np.min(t),np.max(t),np.min(f),np.max(f)], origin='lower', **plotkwargs)
    im = ax.imshow(pn,vmin=vr[0],vmax=vr[1],cmap='turbo',aspect='auto', extent=[np.min(t),np.max(t),np.min(f),np.max(f)], origin=origin)
    #im = ax.imshow(pn,vmin=vr[0],vmax=vr[1],aspect='auto', extent=[np.min(t),np.max(t),np.min(f),np.max(f)], origin='lower', **plotkwargs)

    ax.set_title(title)


    #Now we can change to desired xrange and yrange
    plt.yscale(yscale)
    ax.set_ylim(yr[0],yr[1])
    ax.set_xlim(xr[0],xr[1])

    if show == True:
        plt.show()

    return pn
