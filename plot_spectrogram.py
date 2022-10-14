#Plot large-array spectograms quickly in Python using imshow. 
#This method is so much faster than pcolormesh, which can take forever. 

#E.g. 
#    fs = 1/(np.mean(vlf12.tsec - vlf12.tsec.shift()))

#    freq12, tspec12, power12 = signal.spectrogram(vlf12.amp, fs, nperseg=512, return_onesided=1)

#plotv = {'vmin':0.5, 'vmax':0.8, 'xlim':[600,800], 'ylim':[6000,8000], 'yscale':'linear'}
#    plot_spectrogram(tspec12, freq12, power12, plotv)




def plot_spectrogram(t,f,p,vr=[0,1], yscale='linear', pl=1, xr=0, yr=0):

    import matplotlib.pyplot as plt
    import numpy as np


    if not xr: xr = [np.min(t),np.max(t)]
    if not yr: yr = [np.min(f),np.max(f)]

    if pl:
        pn = 10.*np.log10(p)
    else: 
        #For linear scale, normalize data to 0-1
        #pn = (((p - np.min(p))) / (np.max(p) - np.min(p)))
        pn = p

    fig, ax = plt.subplots()


    
    im = ax.imshow(pn,vmin=vr[0],vmax=vr[1],cmap='turbo',aspect='auto', extent=[np.min(t),np.max(t),np.min(f),np.max(f)], origin='lower')
    #Force x,y axes to be defined exactly based on the range of data.
    ax.set_ylim(yr[0],yr[1])
    ax.set_xlim(xr[0],xr[1])
    #ax.set_ylim(np.min(f),np.max(f))
    #ax.set_xlim(np.min(t),np.max(t))

    plt.yscale(yscale)
    plt.xlim(xr)
    plt.ylim(yr)

    plt.show()

    return pn
