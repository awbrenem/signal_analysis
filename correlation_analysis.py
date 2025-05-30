"""
Functions for correlation analysis of signals.
NOTE: Best current bet for computing phase/coherence dynamic spectra is csd_spectrum_piecewise.py. 

psd - compute simple 1D power spectral density

phase_cc_timelag_analysis - calculate running phase vs time using a cc lag analysis. 

cross_spectral_density - compute 1D coherence and phase vs freq from calculation of CSD

interferometric_coherence_2D(Z1,Z2,N): [NOTE: phase calculation not working. Use csd_spectrum_piecewise.py] 
    essentially the same as cross_spectral_density_spectrogram but can return 
    gain/phase values at higher resolution. Better to use this one if you're able to input complex spectrograms.

csd_spectrum_piecewise - [****NOTE: current best routine to use] same as cross_spectral_density_spectrogram but is better b/c it 
    doesn't rely on an exactly constant sample rate (which can lead to timing discrepancies). It gets the phases, coherence, and timing correct.
    Note that the coherence array can have slightly different freqs than the phase array.

cross_spectral_density_spectrogram - [NOTE: use csd_spectrum_piecewise.py if you have complex input] compute coherence and phase in 2D spectrogram form from calculation of CSD.
    Maybe use this one if you don't have complex spectrograms to input. 

signal_coherence - Compute 1D coherence vs freq

auto_correlation 

cross_correlation 
"""


def psd(wf, tvals, fs, tr, nft=None, zlog=0, ScaleByFreq=True):
    """
    Compute simple power spectral density (real valued) from input waveform
    Returns units of units/sqrt(Hz). 

    To compare with waveform amplitudes at a set frequency you must then integrate over a freq range. 
        amp_f0-f1 = sum(i=0,nbins)[PSD(units/sqrt(Hz)) * sqrt(binsize)]

    E.g. for comparison to a wave with a particular amplitude that extends from 100-200 Hz.
    Assume freq bin sizes of 20 Hz. From 100-200 Hz there are five 20 Hz bins. 

        amp_100-200 = sum(i=0,20)[PSD(units/sqrt(Hz)) * np.sqrt(20)]   (units of "units")
    
    I've tested on some Endurance data and this seems to be working fine. 
    """

    #import scipy.signal
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.mlab import psd
    import matplotlib.mlab as mlab


    goot = np.where((tvals >= tr[0]) & (tvals <= tr[1]))
    wfz = wf[goot]


    #power spectral density (units**2 / Hz)
    S, f = psd(wfz, NFFT=nft, Fs=fs, scale_by_freq=ScaleByFreq, window=mlab.window_hanning)
    #S, f = plt.psd(wfz, Fs=fs, scale_by_freq=ScaleByFreq, NFFT=nft)



    if zlog == 1:
        S = 10*np.log10(S)
    else:
        S = np.sqrt(S)

    """
    plt.plot(f,S)
    if ScaleByFreq == False:
        if zlog == 1:
            plt.ylabel('PSD (dB of from power)')
        else:
            plt.ylabel('PSD (amplitude)')
    else:
        if zlog == 1:
            plt.ylabel('PSD (dB/Hz from power)')
        else:
            plt.ylabel('PSD (amplitude/sqrt(Hz))')


    plt.xlabel('frequency [Hz]')
    plt.show()
    """
    
    return S, f


"""
def psd(wf, tvals, fs, tr, nft=None, zlog=0):
    #Compute simple power spectral density (real valued) from input waveform
    #Returns dB of units**2/Hz

    #import scipy.signal
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.mlab as mlab

    goot = np.where((tvals >= tr[0]) & (tvals <= tr[1]))
    #window = np.hanning(len(goot[0]))
    #wfz = wf[goot]*window
    wfz = wf[goot]

    #units of units**2
    S, f = mlab.psd(wfz, Fs=fs, scale_by_freq=False, NFFT=nft)


    if zlog == 1:
        S = 20*np.log10(S)
    else:
        S = np.sqrt(S)


    plt.plot(f,S)
    if zlog:
        plt.ylabel('PSD (dB of units**2)')
    else:
        plt.ylabel('PSD (units**1)')
    plt.xlabel('frequency [Hz]')
    plt.show()

    return S, f
"""



def phase_cc_timelag_analysis(wf1,wf2,times,fs,ccstep=10):

    """
    Calculate running phase vs time using a cc lag analysis. 
    Uses the FFT to determine the frequency of peak power, so 
    works the best when the signal in question is fairly sinusoidal. 


    wf1, wf2 - input waveforms
    times - times of each waveform must be the same 
    fs - sample rate (#/sec)
    ccstep - helps to determine the stepsize

    Outputs the phases (deg) and times for the phase bins
    """

    import numpy as np
    from scipy import signal


    #normalize waveforms
    wf1 /= np.nanmax(wf1)
    wf2 /= np.nanmax(wf2)


    #determine the dominant wave period to turn lags into phases
    tmp = np.abs(np.fft.rfft(wf1))
    tmpf = np.fft.rfftfreq(len(times), d=1/fs)
    fmax = tmpf[np.argmax(tmp)]  #freq of max power 


    lagmax = (1/fmax) * fs      #number of data points in a single period at peak power


    stepsize = (1/fmax)*ccstep
    nsteppts = int(stepsize * fs)
    nsteps = int(np.floor((np.nanmax(times) - np.nanmin(times))/stepsize))


    bb = 0
    uu = nsteppts
    phase = np.zeros(nsteps)
    phaseT = np.zeros(nsteps)

    for qq in range(nsteps):

        v1t = wf1[bb:uu]
        v2t = wf2[bb:uu]
        v1tx = times[bb:uu]

        corr = signal.correlate(v1t,v2t)
        corr /= np.max(corr)
        lags = signal.correlation_lags(len(v1t), len(v2t)) / fs
        lagmax = lags[np.argmax(corr)]

        #delay time as fraction of signal period
        frac_period = 100 * lagmax / (1/fmax)
        #turn the delay time into a phase
        phase[qq] = 360*frac_period / 100

        #midpoint time for each chunk
        phaseT[qq] = (v1tx[len(v1tx)-1] + v1tx[0]) / 2

        bb += nsteppts
        uu += nsteppts


    #wrap phase to range -180 to 180 
    phase = np.degrees((np.radians(phase) + np.pi) % (2 * np.pi) - np.pi)

    return phase, phaseT




def signal_coherence(wf1,wf2,fs,nperseg=1024,plot=False):
    """
    Compute the coherence vs freq from scipy.signal.coherence
    """

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



def cross_spectral_density(wf1,wf2,fs,nperseg=256,plotshow=False):
    """
    Calculate 1D coherence and phase vs freq from a cross spectral density analysis (matplotlib mlab csd)
    Coherence defined as: 
        Cxy(f) = |Pxy(f)|^2 / (psd1 * psd2)
        with units of amplitude**2/Hz

    (for the spectrogram version use cross_spectral_density_spectrogram)

    https://stackoverflow.com/questions/21647120/how-to-use-the-cross-spectral-density-to-calculate-the-phase-shift-of-two-relate 

    Returns coherence, phase (rad), and frequency (Hz)

    """

    from matplotlib import mlab
    import matplotlib.pyplot as plt
    import numpy as np
    # First create power spectral densities for normalization
    ps1, f = mlab.psd(wf1, Fs=fs, scale_by_freq=False, NFFT=nperseg)
    ps2, f = mlab.psd(wf2, Fs=fs, scale_by_freq=False, NFFT=nperseg)
    
    # Then calculate cross spectral density
    csd, f = mlab.csd(wf1, wf2,NFFT=nperseg, Fs=fs,sides='default', scale_by_freq=False)
    coherence = np.absolute(csd)**2 / (ps1 * ps2)
    angle = np.angle(csd)  #-Pi to Pi

    if plotshow:
        fig = plt.figure()
        ax1 = fig.add_subplot(2, 1, 1)
        ax1.plot(f, csd_norm)
        ax2 = fig.add_subplot(2, 1, 2)
        ax2.plot(f, angle)

    return coherence, angle, f
    



def interferometric_coherence_2D(Z1,Z2,N, coh_min=0):
    """
    https://elisecolin.medium.com/why-2d-convolutions-everywhere-in-sar-imaging-583e046e4b1c

    Same result as cross_spectral_density_spectrogram but with better resolution. 
    That one will output phase and coherence arrays that are smaller than the power array taken from 
    FFT'ing the input waveforms. This routine, however, will output arrays of the same size. 

    ********************************
    NOTE: The phases are wrong. They're only going from -90 to 90, and are mostly centered on zero. 
    Conversely, the phases are correct for cross_spectral_density_spectrogram.
    I don't believe that this issue is being caused by the input as I've used two different programs to 
    input values here and get the same results. 
    ********************************

    (For use see interferometry_routines_call.py)

    Input:
        Z1, Z2 --> complex spectra ("images") of wave power
        N --> 2D blurring window. N must be large enough for the estimate to be robust, 
                but small enough not to blur the transitions between different backscattering zones.
                Try N~few 
    Returns:
        coherence 
        phase (radians)
                
    """
    import numpy as np 
    from scipy import signal

    win = np.ones((N,N))
    num = signal.convolve2d(Z1*np.conj(Z2), win, mode='same')
    den1 = signal.convolve2d(Z2*np.conj(Z2), win, mode='same')
    den2 = signal.convolve2d(Z1*np.conj(Z1), win, mode='same')

    gamma=num/np.sqrt(np.abs(den1)*np.abs(den2))  #this is the cross-spectral density CSD
    coherence=np.abs(gamma)
    #phase= -1 * np.angle(gamma)  #-1 to define phase sense as + in the direction of probe that measured Z1
    phase= np.angle(gamma)  
    phase2= np.arctan2(gamma.imag,gamma.real)  
    #phase3= np.arctan(gamma.imag/gamma.real)  


    #Remove values corresponding to low coherence, if desired
    if coh_min != 0:
        bad = np.where(coherence < coh_min)
        coherence[bad] = np.nan
        phase[bad] = np.nan


    return gamma,coherence,phase



"""

Perform a piecewise FFT to produce a dynamic cross spectral density spectrum.
Includes spectra of phase and coherence. 
(see fft_spectrum_piecewise.py)

This is important for data that has sample rate variations. The Python FFT routines input the waveform and the sample rate, not the time series. 
Sample rate variations thus cause timing issues in the output. 

This code computes the FFT for small chunks and gloms them together for a final dynamic spectrum (sonogram)

fs_thres --> fraction of the median sample rate that current sample rate can be off by. 
            e.g. fs_thres = 0.02 (i.e. 2%) means that sample rates < 98% of the median sample rate or > 102% are not considered. 
            Spectral data at these times will be NaN'd.

NOTE: the calculation of the phase requires 'twosided', and hence has both neg/pos freqs. 
I'm calculating the coherence separately using matplotlib's cohere function, which has only real freqs. 
To get this to work properly I have to calculate it at half the "nfft" value and pad it to the correct length.  
This is why the phase and coherence spectra have separate frequency arrays (they're slightly different).
I can't get the coherence calculation to work out any other way. 

NOTE: the frequency arrays have shape [nfft/2, ntimes]. This is b/c the sample rate can chance with each chunk
and thus so can the actual frequencies. For plotting, you'll want to choose a single freq array (else the plotting will take forever)

"""

def csd_spectrum_piecewise(times, data1, data2, nfft=512, noverlap=8, fs_thres=0.1):

    import sys 
    sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/Endurance/')
    sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
    sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/plasma-physics-general/')
    import numpy as np 
    import matplotlib.mlab as mlab



    #Find most common sample rate
    dt = (times - np.roll(times,1))[1:]
    fs_median = np.median(1/dt)

    #timerange = times[-1] - times[0]

    nchunks = int(np.floor(len(times) / nfft))
    #nchunks = int(np.floor(len(times) / (10*nfft)))
    #chunkpts = 10*nfft

    #create arrays of final data
    #spec_fin = np.empty((int(nfft/2 + 1), nchunks))
    #freqs_fin = np.empty((int(nfft/2 + 1), nchunks))
    spec_fin1 = np.empty((int(nfft), nchunks), dtype=np.complex128)
    spec_fin2 = np.empty((int(nfft), nchunks), dtype=np.complex128)
    csd_fin = np.empty((int(nfft), nchunks), dtype=np.complex128)
    phase_fin = np.empty((int(nfft), nchunks))
    coh_direct = np.empty((int(nfft/2), nchunks))
    freqs_coh = np.empty((int(nfft/2), nchunks))
    freqs_fin = np.empty((int(nfft), nchunks))
    tcenter_fin = np.empty(nchunks)
    tleft_fin = np.empty(nchunks)
    fs_fin = np.empty(nchunks)

    for i in range(nchunks):
        dtmp1 = data1[i*nfft:(i*nfft)+nfft]
        dtmp2 = data2[i*nfft:(i*nfft)+nfft]
        ttmp = times[i*nfft:(i*nfft)+nfft]

        #Test here for even sample rates
        #If sample rate at current time is close to median sample rate, then continue
        fs_fin[i] = 1/(ttmp[1]-ttmp[0])
        fs_delta = fs_fin[i] / fs_median
        if (fs_delta > (1-fs_thres)) & (fs_delta < (1+fs_thres)):
            spec_fin1[:,i], freqs_fin[:,i] = mlab.psd(dtmp1, NFFT=nfft, Fs=fs_fin[i], scale_by_freq=False, 
                                                window=mlab.window_hanning, noverlap=noverlap,sides='twosided')
            spec_fin2[:,i], freqs_fin[:,i] = mlab.psd(dtmp2, NFFT=nfft, Fs=fs_fin[i], scale_by_freq=False, 
                                                window=mlab.window_hanning, noverlap=noverlap,sides='twosided')

            # Then calculate cross spectral density
            csd_fin[:,i], freqs_fin[:,i] = mlab.csd(dtmp1, dtmp2, NFFT=nfft, Fs=fs_fin[i],sides='twosided',scale_by_freq=False)
            #tstx, fx = mlab.csd(dtmp1, dtmp2, NFFT=nfft, Fs=fs_fin[i],sides='twosided',scale_by_freq=False)
            #coh_fin[:,i] = np.abs(csd_fin[:,i])**2 / (spec_fin1[:,i] * spec_fin2[:,i])
            phase_fin[:,i] = np.angle(csd_fin[:,i])  #-Pi to Pi

            #Coherence arrays different sizes b/c only real freqs are considered
            cohtmp, cohf = np.asarray(mlab.cohere(dtmp1, dtmp2, int(nfft/2), pad_to=nfft,Fs=fs_fin[i]))
            coh_direct[:,i] = cohtmp[1:]
            freqs_coh[:,i] = cohf[1:]

        else:
            spec_fin1[:,i] = np.nan
            spec_fin2[:,i] = np.nan
            csd_fin[:,i] = np.nan
            #coh_fin[:,i] = np.nan
            phase_fin[:,i] = np.nan
            coh_direct[:,i] = np.nan

        #define time bins
        tcenter_fin[i] = (ttmp[-1] + ttmp[0])/2.
        tleft_fin[i] = ttmp[0]


    #-------------------------------------------------
    #turn the two-sided power array into a complex array
    #indextmp = int(nfft/2)

    #pcomplex1 = spec_fin1[indextmp:,:] + spec_fin1[0:indextmp,:] * 1j
    #pcomplex = complex(spec_fin[indextmp:,:],spec_fin[0:indextmp,:])

    #freqs_fin = freqs_fin[0:indextmp,:]
    #freqs_fin = freqs_fin[indextmp:,:]



    return freqs_fin[int(nfft/2):,:], freqs_coh, tcenter_fin, csd_fin[int(nfft/2):,:], coh_direct, phase_fin[int(nfft/2):,:], spec_fin1[int(nfft/2):,:], spec_fin2[int(nfft/2):,:], fs_fin







def cross_spectral_density_spectrogram(wf1,wf2,times,fs,timechunk,nperseg=1024,plot=False,coh_min=0):
    """
    NOTE: if you can input complex spectrograms, consider using interferometric_coherence_2D instead. 

        =***********************************************
    NOTE: THE TIMING VALUES ARE OFF IN THIS ROUTINE, WHEREAS THEY ARE FINE IN interferometric_coherence_2D.
    I'VE DISCOVERED THAT THIS IS DUE TO A VARYING SAMPLE FREQUENCY WITH TIME. 
    BEST TO SWITCH TO USING FFT_SPECTRUM_PIECEWISE.PY
        =***********************************************

    Compute the sliding spectrogram of the cross spectral density (Pxy(f):  see SciPy.signal.csd).
    Returns the coherence, phase (radians), and power spectrum (along with time and freq values). 
    Coherence defined as: 
        Cxy(f) = |Pxy(f)|^2 / (psd1 * psd2)
        with units of amplitude**2/Hz

        coh_min -> set to remove all coherence/phase data below this coherence threshold.
      
        timechunk --> (time, sec) determines how many datapoints to include in each CSD calculation. e.g. timechunk=0.2 sec
        nperseg --> FFT current "timechunk" using this size

    """

    from matplotlib import mlab
    import matplotlib.pyplot as plt
    import numpy as np 
    from scipy import signal

    nchunks = int(np.ceil((wf1.size/fs)/timechunk))
    tchunks = np.linspace(np.min(times), np.max(times), num=nchunks)

    #number of array elements for each timechunk
    ios = int(timechunk * fs)

    #Determine how large the FFT arrays will be
    Ptst, freqs = mlab.csd(wf1[0:ios], wf2[0:ios],NFFT=nperseg, Fs=fs,sides='default', scale_by_freq=True)

    #cross spectrum 
    Pxy = np.empty((len(Ptst), nchunks), dtype=complex)
    phase = np.zeros((len(Ptst), nchunks))
    coherence = np.transpose(phase)
    
    for i in range(nchunks): 
        Pxy[:,i], freqs = mlab.csd(wf1[i*ios:(i+1)*ios], wf2[i*ios:(i+1)*ios],NFFT=nperseg, Fs=fs,sides='default', scale_by_freq=True)
        ##individial power spectral densities
        #psd1, f = mlab.psd(wf1[i*ios:(i+1)*ios], Fs=fs, scale_by_freq=True, NFFT=nperseg)
        #psd2, f = mlab.psd(wf2[i*ios:(i+1)*ios], Fs=fs, scale_by_freq=True, NFFT=nperseg)
 
        phase[:,i] = np.angle(Pxy[:,i])
        #angle of complex PSD using arctan2 [see Eqn 2, Graham+16; doi:10.1002/2015JA021527]
        #phase[:,i] = np.degrees(np.arctan2(np.imag(Pxy[:,i]), np.real(Pxy[:,i])))

        coherence[i,:] = np.asarray(mlab.cohere(wf1[i*ios:(i+1)*ios],wf2[i*ios:(i+1)*ios],nperseg,fs))[0,:]

    coherence = np.transpose(coherence)

    #Remove values corresponding to low coherence, if desired
    if coh_min != 0:
        bad = np.where(coherence < coh_min)
        coherence[bad] = np.nan
        phase[bad] = np.nan


    if plot == True: 
        import sys
        sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
        import plot_spectrogram as ps
        import matplotlib.pyplot as plt
        ps.plot_spectrogram(tchunks,freqs,Pxy,vr=[-100,-10], zscale='log')

        fig,axs = plt.subplots(2)
        ps.plot_spectrogram(tchunks,freqs,coherence,vr=[0,1],zscale='linear',ax=axs[0])
        ps.plot_spectrogram(tchunks,freqs,np.abs(phase),vr=[0,180],zscale='linear',ax=axs[1])


    return coherence, phase, tchunks, freqs




def auto_correlation(wf):
    """
    Auto correlation
    ***********************************************
    WARNING: this will always decrease with lag due to the decreasing number of overlap values as window is slid along
    PROBABLY FIND A BETTER WAY. 
    """

    import numpy as np 

    # Normalized data
    ndata = wf - np.mean(wf)

    acorr = np.correlate(ndata, ndata, 'full')[len(ndata)-1:] 
    acorr = acorr / np.var(wf) / len(ndata)

    return acorr


def cross_correlation(wf):
    """
    ******************************************
    NOT IMPLEMENTED YET
    ******************************************
    """

    x = np.arange(128) / 128
    sig = np.sin(2 * np.pi * x)
    sig_noise = sig + rng.standard_normal(len(sig))
    
    
    corr = signal.correlate(sig_noise, sig)
    lags = signal.correlation_lags(len(sig), len(sig_noise))
    corr /= np.max(corr)

"""
import numpy as np
from scipy import signal
rng = np.random.default_rng()
x = rng.standard_normal(1000)
y = np.concatenate([rng.standard_normal(100), x])
correlation = signal.correlate(x, y, mode="full")
lags = signal.correlation_lags(x.size, y.size, mode="full")
lag = lags[np.argmax(correlation)]
"""

