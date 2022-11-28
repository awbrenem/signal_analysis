"""
Functions for correlation analysis of signals

cross_spectral_density_spectrogram

"""



"""
Auto correlation
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


"""
Compute the sliding spectrogram of the cross spectral density from SciPy.signal.csd

Units are amplitude**2/Hz
"""
def cross_spectral_density_spectrogram(wf1, wf2, times, fs, timechunk, nperseg=1024, plot=False):

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

