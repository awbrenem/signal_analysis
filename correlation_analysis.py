"""
Functions for correlation analysis of signals

psd - compute simple 1D power spectral density
phase_cc_timelag_analysis - calculate running phase vs time using a cc lag analysis. 
cross_spectral_density - compute 1D coherence and phase vs freq from calculation of CSD
interferometric_coherence_2D(Z1,Z2,N): essentially the same as cross_spectral_density_spectrogram but can return 
    gain/phase values at higher resolution. Better to use this one.
cross_spectral_density_spectrogram - compute coherence and phase in 2D spectrogram form from calculation of CSD
signal_coherence - Compute 1D coherence vs freq
auto_correlation 
cross_correlation 
"""


def psd(wf, tvals, fs, tr):
    """
    Compute simple power spectral density
    """
        
    import scipy.signal
    import numpy as np
    import matplotlib.pyplot as plt

    goot = np.where((tvals >= tr[0]) & (tvals <= tr[1]))
    window = np.hanning(len(goot[0]))
    wfz = wf[goot]*window

    f, S = scipy.signal.periodogram(wfz, fs, scaling='density')


    plt.semilogy(f, S)
    #plt.ylim([1e-7, 1e2])
    #plt.xlim([0,100])
    plt.xlabel('frequency [Hz]')
    plt.ylabel('PSD [V**2/Hz]')
    plt.show()

    return S, f




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
    angle = np.angle(csd, deg=True)  #from -180 to 180

    if plotshow:
        fig = plt.figure()
        ax1 = fig.add_subplot(2, 1, 1)
        ax1.plot(f, csd_norm)
        ax2 = fig.add_subplot(2, 1, 2)
        ax2.plot(f, angle)

    return coherence, angle, f
    



def interferometric_coherence_2D(Z1,Z2,N):
    """
    https://elisecolin.medium.com/why-2d-convolutions-everywhere-in-sar-imaging-583e046e4b1c

    Same result as cross_spectral_density_spectrogram but with better resolution. 
    That one will output phase and coherence arrays that are smaller than the power array taken from 
    FFT'ing the input waveforms. This routine, however, will output arrays of the same size. 


    (For use see interferometry_routines_call.py)

    Input:
        Z1, Z2 --> complex spectra ("images") of wave power
        N --> 2D blurring window. N must be large enough for the estimate to be robust, 
                but small enough not to blur the transitions between different backscattering zones.
                Try N~few 

    """
    import numpy as np 
    from scipy import signal

    win = np.ones((N,N))
    num = signal.convolve2d(Z1*np.conj(Z2), win, mode='same')
    den1 = signal.convolve2d(Z2*np.conj(Z2), win, mode='same')
    den2 = signal.convolve2d(Z1*np.conj(Z1), win, mode='same')

    gamma=num/np.sqrt(np.abs(den1)*np.abs(den2))
    coherence=np.abs(gamma)
    phase=np.angle(gamma)

    return gamma,coherence,phase


def cross_spectral_density_spectrogram(wf1,wf2,times,fs,timechunk,nperseg=1024,plot=False,coh_min=0):
    """
    NOTE: consider using interferometric_coherence_2D instead. 

    Compute the sliding spectrogram of the cross spectral density (Pxy(f):  see SciPy.signal.csd).
    Returns the coherence, phase, and power spectrum (along with time and freq values). 
    Coherence defined as: 
        Cxy(f) = |Pxy(f)|^2 / (psd1 * psd2)
        with units of amplitude**2/Hz

        coh_min -> set to remove all coherence/phase data below this coherence threshold.
        
    """

    from matplotlib import mlab
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
    coherence = np.zeros((len(Ptst), nchunks))
    phase = np.zeros((len(Ptst), nchunks))

    for i in range(nchunks): 
        Pxy[:,i], freqs = mlab.csd(wf1[i*ios:(i+1)*ios], wf2[i*ios:(i+1)*ios],NFFT=nperseg, Fs=fs,sides='default', scale_by_freq=True)
        #individial power spectral densities
        psd1, f = mlab.psd(wf1[i*ios:(i+1)*ios], Fs=fs, scale_by_freq=True, NFFT=nperseg)
        psd2, f = mlab.psd(wf2[i*ios:(i+1)*ios], Fs=fs, scale_by_freq=True, NFFT=nperseg)
 
        #angle of complex PSD using arctan2 [see Eqn 2, Graham+16; doi:10.1002/2015JA021527]
        phase[:,i] = np.angle(Pxy[:,i],deg=True)
        coherence[:,i] = np.abs(Pxy[:,i])**2 / (psd1 * psd2)

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

