"""

Perform a piecewise FFT to produce a dynamic spectrum 
This is important for data that has sample rate variations. The Python FFT routines input the waveform and the sample rate, not the time series. 
Sample rate variations thus cause timing issues in the output. 

This code computes the FFT for small chunks and gloms them together for a final dynamic spectrum (sonogram)

fs_thres --> fraction of the median sample rate that current sample rate can be off by. 
            e.g. fs_thres = 0.02 (i.e. 2%) means that sample rates < 98% of the median sample rate or > 102% are not considered. 
            Spectral data at these times will be NaN'd.

full_freqs --> returns a freq spectrum for each time. Since the sample rate can vary, so do the freqs. 
                However, attempting to plot a spectrogram with these will essentially crash Python. 
                If not set, returns the median (over all times) frequencies 
            

Returns:     freqs_fin --> [f] (or [f,t] if 'full_freqs' set)  (final freqs)
             tcenter_fin --> [t] (final time values - center of bin)
             spec_fin --> [f,t]  (final sonogram values)
             fs_fin --> [t]    (sample rate at each time chunk)
             
"""

def fft_spectrum_piecewise(times, data, nfft=512, noverlap=8, fs_thres=0.02, full_freqs=0):

    import sys 
    sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/Endurance/')
    sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
    sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/plasma-physics-general/')
    import numpy as np 
    from matplotlib.mlab import psd
    import matplotlib.mlab as mlab



    #Find most common sample rate
    dt = (times - np.roll(times,1))[1:]
    fs_median = np.median(1/dt)


    nchunks = int(np.floor(len(times) / nfft))

    #create arrays of final data
    #one-sided spectrum
    spec_fin = np.empty((int(nfft/2 + 1), nchunks))
    freqs_fin = np.empty((int(nfft/2 + 1), nchunks))
    #two-sided spectrum
    #spec_fin = np.empty((int(nfft), nchunks))
    #freqs_fin = np.empty((int(nfft), nchunks))
    tcenter_fin = np.empty(nchunks)
    tleft_fin = np.empty(nchunks)
    fs_fin = np.empty(nchunks)


    for i in range(nchunks):

        dtmp = data[i*nfft:(i*nfft)+nfft]
        ttmp = times[i*nfft:(i*nfft)+nfft]
    

        #Test here for even sample rates
        #If sample rate at current time is close to median sample rate, then continue
        fs_fin[i] = 1/(ttmp[1]-ttmp[0])
        fs_delta = fs_fin[i] / fs_median
        if (fs_delta > (1-fs_thres)) & (fs_delta < (1+fs_thres)):
            spec_fin[:,i], freqs_fin[:,i] = psd(dtmp, NFFT=nfft, Fs=fs_fin[i], scale_by_freq=True, 
                                                window=mlab.window_hanning, noverlap=noverlap,sides='onesided')
            #spec_fin[:,i], freqs_fin[:,i] = psd(dtmp, NFFT=nfft, Fs=fs_fin[i], scale_by_freq=True, 
            #                                    window=mlab.window_hanning, noverlap=noverlap,sides='twosided')
        else:
            spec_fin[:,i] = np.nan


        #define time bins
        tcenter_fin[i] = (ttmp[-1] + ttmp[0])/2.
        tleft_fin[i] = ttmp[0]


    #if you don't need the freq arrays for each time, average the frequencies. 
    fmedian = np.zeros(int(nfft/2 + 1))

    if not full_freqs:
        for f in range(len(fmedian)):
            fmedian[f] = np.median(freqs_fin[f,:])
        freqs_fin = fmedian

    #Test the median freqs vs freqs at a specific time
    #plt.plot(freqs_fin[:,10])
    #plt.plot(fmedian[:])
    #df = freqs_fin[:,100] - fmedian
    #plt.plot(df)


    return freqs_fin, tcenter_fin, spec_fin, fs_fin

    """
    #-------------------------------------------------
    #turn the two-sided power array into a complex array
    indextmp = int(nfft/2)

    pcomplex = spec_fin[indextmp:,:] + spec_fin[0:indextmp,:] * 1j
    #pcomplex = complex(spec_fin[indextmp:,:],spec_fin[0:indextmp,:])

    #freqs_fin = freqs_fin[0:indextmp,:]
    freqs_fin = freqs_fin[indextmp:,:]


    return freqs_fin, tcenter_fin, pcomplex, fs_fin
    #return freqs_fin, tcenter_fin, spec_fin, fs_fin
    """