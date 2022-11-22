"""
Routines for bandpassing waveform data


bandpassing waveform data
From https://stackoverflow.com/questions/12093594/how-to-implement-band-pass-butterworth-filter-with-scipy-signal-butter



Example: 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
import bandpass_data_butter

lowcut = 100
highcut = 1000
fs = 30000 #Samples/sec

wfbp = bandpass_data_butter.butter_bandpass_filter(wf, lowcut, highcut,fs,order=9)


filter_wave_frequency

"""


from scipy.signal import butter, sosfilt, sosfreqz


def butter_bandpass(lowcut, highcut, fs, order=5):
        nyq = 0.5 * fs
        low = lowcut / nyq
        high = highcut / nyq
        sos = butter(order, [low, high], analog=False, btype='band', output='sos')
        return sos

def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
        sos = butter_bandpass(lowcut, highcut, fs, order=order)
        y = sosfilt(sos, data)
        return y


