"""
Functions for correlation analysis of signals

cross_spectral_density (gain and phase)
cross_spectral_density_spectrogram (gain only)
signal_coherence

"""


"""
Compute the coherence vs freq from scipy.signal.coherence
"""
def signal_coherence(wf1,wf2,fs,nperseg=1024,plot=False):
    from scipy import signal

    f, Cxy = signal.coherence(wf1, wf2, fs, nperseg=nperseg)

    if plot == True:
        import matplotlib.pyplot as plt
        plt.plot(f, Cxy)
        plt.xlabel('frequency [Hz]')
        plt.ylabel('Coherence')
        plt.yscale('linear')
        plt.show()

    return Cxy, f


"""
Calculate gain and phase from a cross spectral density analysis (matplotlib mlab csd)
https://stackoverflow.com/questions/21647120/how-to-use-the-cross-spectral-density-to-calculate-the-phase-shift-of-two-relate 
"""
def cross_spectral_density(wf1,wf2,fs,nperseg=256):
    from matplotlib import mlab
    import matplotlib.pyplot as plt
    import numpy as np

    # First create power spectral densities for normalization
    (ps1, f) = mlab.psd(wf1, Fs=fs, scale_by_freq=False, NFFT=nperseg)
    (ps2, f) = mlab.psd(wf2, Fs=fs, scale_by_freq=False, NFFT=nperseg)
    #plt.plot(f, ps1)
    #plt.plot(f, ps2)

    
    # Then calculate cross spectral density
    (csd, f) = mlab.csd(wf1, wf2,NFFT=nperseg, Fs=fs,sides='default', scale_by_freq=False)
    csd_norm = np.absolute(csd)**2 / (ps1 * ps2)

    
    fig = plt.figure()
    ax1 = fig.add_subplot(2, 1, 1)
    # Normalize cross spectral absolute values by auto power spectral density
    ax1.plot(f, csd_norm)
    ax2 = fig.add_subplot(2, 1, 2)
    angle = np.angle(csd, deg=True)  #from -180 to 180
    ax2.plot(f, angle)
    
    return csd_norm, angle, f
    


"""
Compute the sliding spectrogram of the cross spectral density (SciPy.signal.csd)

Units are amplitude**2/Hz
NOTE: there's no phase info with this method
"""
def cross_spectral_density_spectrogram(wf1,wf2,times,fs,timechunk,nperseg=1024,plot=False):

    import numpy as np 
    from scipy import signal

    nchunks = int(np.ceil((wf1.size/fs)/timechunk))
    tchunks = np.linspace(np.min(times), np.max(times), num=nchunks)

    #number of array elements for each timechunk
    ios = int(timechunk * fs)

    #Determine how large the FFT arrays will be
    freqs, Ptst = signal.csd(wf1[0:ios], wf2[0:ios], fs, nperseg=nperseg)
    Pxy = np.zeros((len(Ptst), nchunks))

    for i in range(nchunks): 
        freqs, Pxy[:,i] = signal.csd(wf1[i*ios:(i+1)*ios], wf2[i*ios:(i+1)*ios], fs, nperseg=nperseg)


    if plot == True: 
        import sys
        sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
        import plot_spectrogram as ps
        ps.plot_spectrogram(tchunks,freqs,Pxy,vr=[-100,-10], zscale='log')


    return Pxy, tchunks, freqs




"""
Auto correlation
***********************************************
WARNING: this will always decrease with lag due to the decreasing number of overlap values as window is slid along
PROBABLY FIND A BETTER WAY. 
"""
def auto_correlation(wf):
    import numpy as np 

    # Normalized data
    ndata = wf - np.mean(wf)

    acorr = np.correlate(ndata, ndata, 'full')[len(ndata)-1:] 
    acorr = acorr / np.var(wf) / len(ndata)

    return acorr


"""
******************************************
NOT IMPLEMENTED YET
******************************************
"""
def cross_correlation(wf):

    x = np.arange(128) / 128
    sig = np.sin(2 * np.pi * x)
    sig_noise = sig + rng.standard_normal(len(sig))
    
    
    corr = signal.correlate(sig_noise, sig)
    lags = signal.correlation_lags(len(sig), len(sig_noise))
    corr /= np.max(corr)

