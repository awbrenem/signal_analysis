"""
Routines for frequency filtering waveform data

Import these routines as:
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
import filter_wave_frequency

Types of filters: 
-Butter bandpass filter (butter_bandpass_filter.py)
-Butter lowpass filter (butter_lowpass_filter.py)
-Butter highpass filter (butter_highpass_filter.py)
-FFT bandpass filter (fft_bandpass_filter.py)
"""



"""
-------------------------------------------------------------------------------------------
Butterworth bandpass filter routines 
-------------------------------------------------------------------------------------------
From https://stackoverflow.com/questions/12093594/how-to-implement-band-pass-butterworth-filter-with-scipy-signal-butter

NOTE: This method (using sosfiltfilt) preserves signal phase!! (tested and verified)
(see https://gist.github.com/junzis/e06eca03747fc194e322)

#Example: 
lowcut = 100
highcut = 1000
fs = 30000 #Samples/sec
wfbp = filter_wave_frequency.butter_bandpass_filter(wf, lowcut, highcut,fs,order=9)

#Test
f, t, p = signal.spectrogram(wfbp, fs, nperseg=16384,noverlap=16384/2,window='hann')
ps.plot_spectrogram(t,f,p,vr=[-100,-40],yr=yrspec, xr=trspec, yscale='linear')

"""

def butter_bandpass(lowcut, highcut, fs, order=5):
        from scipy.signal import butter
        nyq = 0.5 * fs
        low = lowcut / nyq
        high = highcut / nyq
        sos = butter(order, [low, high], analog=False, btype='band', output='sos')
        return sos

def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
        from scipy.signal import sosfiltfilt
        sos = butter_bandpass(lowcut, highcut, fs, order=order)
        #y = sosfilt(sos, data)   #does NOT preserve signal phase
        y = sosfiltfilt(sos, data) #preserves signal phase

        """
        # Plot the frequency response.
        from scipy.signal import sosfreqz
        w, h = sosfreqz(sos)
        import matplotlib.pyplot as plt
        import numpy as np
        plt.subplot(2, 1, 1)
        plt.plot(0.5*fs*w/np.pi, np.abs(h), 'b')
        plt.plot(lowcut, 0.5*np.sqrt(2), 'ko')
        plt.axvline(lowcut, color='k')
        plt.xlim(0.1, 0.5*fs)
        plt.title("Bandpass Filter Frequency Response")
        plt.xlabel('Frequency [Hz]')
        plt.xscale("log")
        plt.grid()
        plt.show()
        """

        return y


def butter_lowpass(lowcut, fs, order=5):
        from scipy.signal import butter
        nyq = 0.5 * fs
        low = lowcut / nyq
        sos = butter(order, low, analog=False, btype='low', output='sos')
        return sos

def butter_lowpass_filter(data, lowcut, fs, order=5):
        from scipy.signal import sosfiltfilt
        sos = butter_lowpass(lowcut, fs, order=order)
        y = sosfiltfilt(sos, data) #preserves signal phase
        return y

def butter_highpass(highcut, fs, order=5):
        from scipy.signal import butter
        nyq = 0.5 * fs
        high = highcut / nyq
        sos = butter(order, high, analog=False, btype='high', output='sos')
        return sos

def butter_highpass_filter(data, highcut, fs, order=5):
        from scipy.signal import sosfiltfilt
        sos = butter_lowpass(highcut, fs, order=order)
        y = sosfiltfilt(sos, data) #preserves signal phase
        return y



                
"""
-------------------------------------------------------------------------------------------
FFT bandpass (usable for purely real waveforms)

***NEEDS TO BE TESTED*********
***NEEDS TO BE TESTED*********
***NEEDS TO BE TESTED*********
-------------------------------------------------------------------------------------------

#Example: 
data --> 1D waveform data
tvals --> corresponding time values for waveform
lowcut = 100
highcut = 1000
fs = 30000 #Samples/sec

import numpy as np 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
import filter_wave_frequency
import plot_spectrogram as ps

data_filt = filter_wave_frequency.fft_bandpass_filter(data,lowcut,highcut,fs)

#Plot comparison of FFT'd data vis FFT'd and filtered data
window = np.hanning(len(data))
X = np.fft.rfft(data*window)
X_filtered = np.fft.rfft(data_filt)
plt.plot(freq,np.abs(X), '.',freq,np.abs(X_filtered2),'x')
plt.show()


#Plot comparison of original and filtered waveforms
plt.plot(tvals,data,tvals,data_filt)
plt.show()

#Plot dynamic spectral comparison of regular and filtered data
f, t, p = signal.spectrogram(data, fs, nperseg=16384,noverlap=16384/2,window='hann') #, return_onesided=1)
fz, tz, pz = signal.spectrogram(data_filt, fs, nperseg=16384,noverlap=16384/2,window='hann') #, return_onesided=1)
ps.plot_spectrogram(t,f,p,vr=[-200,0],xr=[10,200],yr=[0,8000], yscale='linear')
ps.plot_spectrogram(tz,fz,pz,vr=[-200,0],xr=[10,200],yr=[0,8000], yscale='linear')

"""

def fft_bandpass_filter(data,lowcut,highcut,fs):

        import numpy as np
        window = np.hanning(len(data))

        X = np.fft.rfft(data*window)
        X_filtered = X.copy()

        #Get frequencies (will be both + and -)
        freq = np.fft.fftfreq(len(X), d=1./fs)

        #Select desired frequency range and remove data outside of it
        cond = (np.abs(freq) > lowcut) & (np.abs(freq) < highcut)
        X_filtered = X_filtered * cond

        return np.fft.irfft(X_filtered)

        

