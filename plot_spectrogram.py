#Plot large-array spectograms quickly in Python using imshow. 
#This method is so much faster than pcolormesh, which can take forever. 

#E.g. 
#    fs = 1/(np.mean(vlf12.tsec - vlf12.tsec.shift()))

#    freq12, tspec12, power12 = signal.spectrogram(vlf12.amp, fs, nperseg=512, return_onesided=1)

#plotv = {'vmin':0.5, 'vmax':0.8, 'xlim':[600,800], 'ylim':[6000,8000], 'yscale':'linear'}
#    plot_spectrogram(tspec12, freq12, power12, plotv)




def plot_spectrogram(t,f,p, vmn=0, vmx=1, yscale='linear'):

    import matplotlib.pyplot as plt
    import numpy as np


    vmn, vmx = 0, 1
    xr = [np.min(t),np.max(t)]
    yr = [np.min(f),np.max(f)]

    plog = 10.*np.log10(p)



    #Normalize data to 0-1
    pn = (((plog - np.min(plog))) / (np.max(plog) - np.min(plog)))

    fig, ax = plt.subplots()
    im = ax.imshow(pn,vmin=vmn,vmax=vmx,cmap='turbo',aspect='auto', extent=[np.min(t),np.max(t),np.min(f),np.max(f)])
    #Force x,y axes to be defined exactly based on the range of data.
    ax.set_ylim(np.min(f),np.max(f))
    ax.set_xlim(np.min(t),np.max(t))

    plt.yscale('linear')
    plt.xlim(xr)
    plt.ylim(yr)


