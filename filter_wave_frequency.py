"""
Routines for frequency filtering waveform data

Import these routines as:
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
import filter_wave_frequency

Types of filters: 
-Butter bandpass filter (butter_bandpass_filter.py)

"""



"""
-------------------------------------------------------------------------------------------
Butterworth bandpass filter routines 
-------------------------------------------------------------------------------------------
From https://stackoverflow.com/questions/12093594/how-to-implement-band-pass-butterworth-filter-with-scipy-signal-butter

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
        nyq = 0.5 * fs
        low = lowcut / nyq
        high = highcut / nyq
        sos = butter(order, [low, high], analog=False, btype='band', output='sos')
        return sos

def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
        from scipy.signal import butter, sosfilt, sosfreqz
        sos = butter_bandpass(lowcut, highcut, fs, order=order)
        y = sosfilt(sos, data)
        return y


"""
-------------------------------------------------------------------------------------------
FFT bandpass (usable for purely real waveforms)
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

        

