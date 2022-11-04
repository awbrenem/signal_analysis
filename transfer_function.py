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




path = '/Users/abrenema/Desktop/Research/Rocket_missions/Endurance/gain_phase_files/'

#fn = "Endurance_Analog 1_HF12_1000-20000000-100.txt"
#fn = "Endurance_Analog 1_HF12 (1)_1000-20000000-100.txt"
#fn = "Endurance_Analog 1_HF12 (2)_1000-20000000-100.txt"
#fn = "Endurance_Analog 1_HF34 (3)_1000-20000000-100.txt"

fn = "Endurance_Analog 1_V12D_10-10000-100.txt"


#fn = "Endurance_Analog 1_HF34 (3)_1000-20000000-100.txt"
with open(path + fn) as f:
    lines = f.readlines()

f = lines[0].split()  #freq in Hz
p = lines[1].split()  #phase in deg
g = lines[2].split()  #gain in dB
f = [float(i) for i in f]
p = [float(i) for i in p]
g = [float(i) for i in g]


fig, axs = plt.subplots(2)
axs[0].plot(f,p)
axs[1].plot(f,g)
for i in axs: i.set_xscale('log')

plt.show()


#change to radians and make imaginary
prad = [i*np.deg2rad for i in p]

#change back to linear scaler
Hmag = [10**(i/10.) for i in g]


#Transfer function  H(w) = |H|*exp(-1*i*theta)
H = [Hmag[i] * np.exp(-1j*prad[i]) for i in range(len(prad))]

fig2, axs2 = plt.subplots(2)
axs2[0].plot(f,Hmag)
axs2[0].set_yscale('log')
axs2[1].plot(f,np.imag(prad),'x',f,prad)
plt.show()



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
#Complex consists of amplitude (where magnitude = |amplitude|) and "phase" (not PSD)
fspec, tspec, powerc = signal.spectrogram(vlf12, fs, nperseg=16384,noverlap=16384/2,window='hann',return_onesided=True,mode='complex')
fspec, tspec, powerm = signal.spectrogram(vlf12, fs, nperseg=16384,noverlap=16384/2,window='hann',return_onesided=True,mode='magnitude')
fspec, tspec, powera = signal.spectrogram(vlf12, fs, nperseg=16384,noverlap=16384/2,window='hann',return_onesided=True,mode='angle')
#ft, tt, pt = signal.spectrogram(vlf12, fs, nperseg=16384,noverlap=16384/2,window='hann',return_onesided=True)


pr = np.real(powerc)
fig, axs3 = plt.subplots(2)
axs3[0].plot(tspec,pr[3000,:])
axs3[1].plot(tspec,powerm[3000,:])
axs3.set_ylim[0,0.1]
plt.show()





#Interpolate transfer function to frequencies of FFT'd waveform data
interp = interp1d(f,H,kind='cubic', bounds_error=False)
Hnew = interp(fspec)
fig, axs4 = plt.subplots(2)
axs4[0].plot(fspec, Hnew)
axs4[1].plot(f, np.real(H))
axs4[0].set_xlim(0,15000)

#axs4[0].plot(fspec, np.imag(Hnew))
#axs4[1].plot(f, np.imag(H))

axs4[0].plot(fspec, np.real(Hnew),'.',f,np.real(H),'-')
axs4[1].plot(fspec, np.imag(Hnew),'.',f,np.imag(H),'-')
plt.show()



#Apply transfer function to the complex FFT data.
#(1) apply it to the complex FFT 

pnewc = powerc[:,0]/Hnew

plt.plot(fspec,powerm[:,0],'.',fspec,pnewc)
plt.show()



#(2) apply it to the mag and phase FFT combined to create complex FFT

ttime = 1000

    #construct complex power
pnew = [powerm[i,ttime] * np.exp(-1j*powera[i,ttime]) for i in range(len(fspec))]


pnewCorr = pnew/Hnew

plt.plot(fspec,pnewCorr,'.',fspec,pnew)
plt.show()

for i in range(100):
    print(fspec[i], np.real(pnewCorr[i]),  np.real(pnew[i]))




num = np.abs(np.real(pnewCorr))
den = np.abs(np.real(pnew))
ratio = [num[i]/den[i] for i in range(len(num))]
plt.plot(fspec,ratio)
plt.show()




plt.plot(fspec,np.abs(np.real(pnew)),'x',fspec,powerm[:,ttime])
plt.xlim(6000,6200)
#plt.ylim(-0.01,0.01)
plt.yscale('log')
plt.show()



powerc2 = powerm[:,0] * np.exp(-1*powera[:,0])



pnewm = powerm[:,0]/np.real(Hnew)
pnewa = powera[:,0]/np.imag(Hnew)
pnewac = [complex(i) for i in pnewa]
pnew = [pnewm[i] * np.exp(-1*pnewac[i]) for i in range(len(pnewm))]


plt.plot(fspec,powerm[:,0],'.',fspec,np.real(pnew))
plt.ylim(-1,1)
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

