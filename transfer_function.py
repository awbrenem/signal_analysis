"""
Calculate transfer function 

From Bode (gain/phase) plot:
(1) Find gain as |H(w)| = 10^B/10, where B is gain in dB from Bode plot
(2) H(w) = |H(w)| * exp(i*theta), where theta is the phase in radians

(see https://resources.pcb.cadence.com/blog/2021-understanding-a-circuit-transfer-function-from-a-bode-plot)
"""

import pandas as pd
import matplotlib.pyplot as plt
from scipy import signal
import numpy as np
import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/Endurance/')
import end_load_data as end
from scipy.interpolate import interp1d
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
import plot_spectrogram as ps




#Load Endurance data for testing
vlf = end.efield_vlf()
tvlf = vlf.tvlf
vlf12 = vlf.dvlf12_mvm

fs = vlf.samplerate


#FFT Endurance VLF data
vr = [-80,-60]
#vr = [-0.1,0.1]
xr = [100,500]
yr = [0,12000]

#FFT to get power spectral density (V^2/Hz). 
#Mode defaults to PSD (V^2/Hz). 
#Complex consists of amplitude (where magnitude = abs(complex))
fspec, tspec, powerc = signal.spectrogram(vlf12, fs, nperseg=16384,noverlap=16384/2,window='hann',return_onesided=True,mode='complex')
fspec, tspec, powerm = signal.spectrogram(vlf12, fs, nperseg=16384,noverlap=16384/2,window='hann',return_onesided=True,mode='magnitude')
fspec, tspec, powera = signal.spectrogram(vlf12, fs, nperseg=16384,noverlap=16384/2,window='hann',return_onesided=True,mode='angle')
#ft, tt, pt = signal.spectrogram(vlf12, fs, nperseg=16384,noverlap=16384/2,window='hann',return_onesided=True)

#Prove that abs(powerc) = powerm
pr = np.abs(powerc)
fig, axs3 = plt.subplots(2)
axs3[0].plot(tspec,pr[3000,:])
axs3[1].plot(tspec,powerm[3000,:])
plt.show()






path = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/gain_phase_files/'

#fn = "Endurance_Analog 1_HF12_1000-20000000-100.txt"
#fn = "Endurance_Analog 1_HF12 (1)_1000-20000000-100.txt"
#fn = "Endurance_Analog 1_HF12 (2)_1000-20000000-100.txt"
#fn = "Endurance_Analog 1_HF34 (3)_1000-20000000-100.txt"

#fn = "Endurance_Analog 1_V12D_10-10000-100.txt"
fn = "Endurance_Analog 1_VLF12D_6-30000-100.txt"

#fn = "Endurance_Analog 1_HF34 (3)_1000-20000000-100.txt"
with open(path + fn) as f:
    lines = f.readlines()

f = lines[0].split()  #freq in Hz
p = lines[1].split()  #phase in deg
g = lines[2].split()  #gain in dB
f = [float(i) for i in f]
p = [float(i) for i in p]
g = [float(i) for i in g]

#change to radians
prad = [np.deg2rad(i) for i in p]




#change gain from dB to linear scaler for calculation of transfer function
#From Steve Martin email on Nov 7, 2022: 
#Gain=10^(0.05 * (opchan+gainoffset))
offset = 0.
Hmag = [10**(0.05*i + offset) for i in g]

fig, axs = plt.subplots(3)
axs[0].plot(f,g)
axs[1].plot(f,Hmag)
axs[2].plot(f,p)
#axs[2].plot(f,prad)
axs[0].set_title('gain/phase; \n fn='+ fn)
axs[0].set_xscale('log')
axs[1].set_xscale('log')
axs[2].set_xscale('log')
axs[0].set_yscale('linear')
axs[1].set_ylim(0,20)
axs[0].set(ylabel='gain(dB)',xlabel='freq(kHz)')
axs[1].set(ylabel='gain',xlabel='freq(kHz)')
axs[2].set(ylabel='phase(deg)',xlabel='freq(kHz)')
#axs[:].set_xlim(-40,10)
plt.show()



#------------------------------------------------------------------
#Transfer function  H(w) = |H|*exp(i*theta)

#Interpolate transfer function to frequencies of FFT'd waveform data
interp = interp1d(f,Hmag,kind='cubic', bounds_error=False)
Hmag = interp(fspec)
interp2 = interp1d(f,prad,kind='cubic', bounds_error=False)
prad = interp2(fspec)

H = [Hmag[i] * np.exp(1j*prad[i]) for i in range(len(prad))]


fig, axs = plt.subplots(2)
axs[0].plot(f,p)
axs[1].plot(fspec, prad)
axs[0].set_xlim(0,15000)
axs[1].set_xlim(0,15000)
plt.show()


#Apply transfer function to the complex FFT data.
#(1) apply it to the complex FFT 
ttime = 2000

Htst = H
Htst = [1 + 0j for i in H]



#The transfer function correction applies to the POWER, not the amplitude

corrected = powerc[:,ttime]/H
#......
#corrected = powerc[:,ttime]/Htst

ratio = [np.abs(powerc[i,ttime])/np.abs(corrected[i]) for i in range(len(H))]


fig, axs = plt.subplots(2)
axs[0].plot(fspec,np.abs(powerc[:,ttime]),'.',fspec,np.abs(corrected),'x')
axs[1].plot(fspec,ratio)
axs[0].set_xscale('log')
axs[0].set_yscale('log')
axs[0].set_xlim(1,100000)
axs[1].set_xscale('log')
axs[1].set_yscale('log')
axs[1].set_xlim(1,100000)
axs[0].set_ylim(1e-7,1e-3)
#axs[1].set_ylim(1e-7,1e-3)
plt.show()















#Take Laplace transform:
https://stackoverflow.com/questions/38316225/numerical-laplace-transform-python
#Laplace from FFT:
https://www.tutorialspoint.com/relation-between-laplace-transform-and-fourier-transform
Therefore, the Laplace transform is just the complex Fourier transform of a signal. 
Hence, the Fourier transform is equivalent to the 
Laplace transform evaluated along the imaginary axis of the s-plane, i.e.,
L[x(t)]=F[x(t)e−σt]















"""
#HF channels 
Endurance_Analog 1_HF34 (3)_1000-20000000-100.txt
Endurance_Analog 1_HF12_1000-20000000-100.txt
Endurance_Analog 1_HF12 (2)_1000-20000000-100.txt
Endurance_Analog 1_HF12 (1)_1000-20000000-100.txt


#VLF digitial files
Endurance_Analog 1_VLF12D_6-30000-100.txt
Endurance_Analog 1_VLF13D_6-30000-100.txt
Endurance_Analog 1_VLF41D_6-30000-100.txt
Endurance_Analog 1_VLF24D_6-30000-100.txt
Endurance_Analog 1_VLF42D_6-30000-100.txt
Endurance_Analog 1_VLF32D_6-30000-100.txt
Endurance_Analog 1_VLF34D_6-30000-100.txt


Endurance_Analog 1_VLF12A_6-100000-100.txt

Endurance_Analog 1_V42D_10-10000-100.txt
Endurance_Analog 1_V41D_10-10000-100.txt
Endurance_Analog 1_V34D_10-10000-100.txt
Endurance_Analog 1_V34A_10-10000-100.txt
Endurance_Analog 1_V32D_10-10000-100.txt
Endurance_Analog 1_V24D_10-10000-100.txt
Endurance_Analog 1_V13D_10-10000-100.txt
Endurance_Analog 1_V12D_10-10000-100.txt
Endurance_Analog 1_V12A_10-10000-100.txt
Endurance_Analog 1_V4SD_10-10000-100.txt
Endurance_Analog 1_V4SA_10-10000-100.txt
Endurance_Analog 1_V3SD_10-10000-100.txt
Endurance_Analog 1_V3SA_10-10000-100.txt
Endurance_Analog 1_V2SD_10-10000-100.txt
Endurance_Analog 1_V2SA_10-10000-100.txt
Endurance_Analog 1_V1SD_10-10000-100.txt
Endurance_Analog 1_V1SA_10-10000-100.txt

"""




#Input frequency spectrum 
Xs = []

#Output frequency spectrum

