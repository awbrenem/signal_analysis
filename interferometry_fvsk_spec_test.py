"""
Example setup for calls to interferometry_routines.py 

NOTE: this technique only makes sense for parallel, spaced antennas. 
If the antennas are perpendicular, the ~90 deg phase shift in waves will be interpreted as a 
large k-value. 


TESTING: See tests from antenna_finite_wavelength_response.py 

NOTE: routine works great, with correct phase sense, with following two exceptions

1) phase identification can be difficult if the wave is propagating nearly along 
one of the interferometry axes (e.g. if along the y-hat' direction then V32 and V14 will
measure very little signal, meaning the waveform can be noise-dominated)

2) if wavelength is on the order or less than boom length than results can't be trusted.



                       (y-hat)
                         V3
                      /  |  \             
           (y-hat') z    |     z (x-hat')         
                  /      |       \        
               V2--------O--------V1 (x-hat)    
                  \      |       /
                    z    |     z
                      \  |  /
                         V4
          

z-points represent the centers of potential of the interferometry diagonals          

Coordinate system (system of input test wave)
x-hat --> E12 = V1 - V2 direction (positive to right)
y-hat --> E34 = V3 - V4 direction (positive upwards)

This code outputs the spectrum of k-values in the kx' and ky' directions, where
x'-hat --> Ex' = V1V3z - V4V2z (45 deg inclined from x-hat)
y'-hat --> Ey' = V3V2z - V1V4z (45 deg inclined from y-hat)

"""


import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/Endurance/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
from scipy import signal
import numpy as np 
import interferometry_routines as interf
import correlation_analysis
import plot_spectrogram as ps
import matplotlib.pyplot as plt
import pickle



#-------------------------------------------------------------
#Example: Load artificially created waveforms for testing 
#(uses pickle file created from antenna_finite_wavelength_response.py)
#-------------------------------------------------------------

#vAstr = 'VLF13D'
#vBstr = 'VLF42D'

#Load artificially generated test waveform
#goo = pickle.load(open('/Users/abrenema/Desktop/wavetst.pkl','rb'))
#goo = pickle.load(open('/Users/abrenema/Desktop/wavetst_f=50-200Hz.pkl','rb'))
goo = pickle.load(open('/Users/abrenema/Desktop/wavetst_f=20-60Hz.pkl','rb'))


#---Choose interferometry pair
#NOTE: positive sense of phase defined as pointing towards center of potential of "wfA"
#e.g. y-hat':  wfA = wf32; wfB = wf14
#x-hat'
#wfA = goo['wf13']
#wfB = goo['wf42']
#y-hat'
wfA = goo['wf32']
wfB = goo['wf14']

tdat = goo['tvals']
vals = goo['vals']
fs = 1/(tdat[1]-tdat[0])




wavelength = [vals['lambda_min'],vals['lambda_max']]
#Values in antenna (non-primed) coordinates (see header)
kx = vals['kx']
ky = vals['ky']


#plt.plot(tdat,wfA)
#plt.xlim(0,0.05)

#Plot FFT
#ft = np.fft.fft(wfA)
#timestep = tdat[1]-tdat[0]
#freqs = np.fft.fftfreq(ft.shape[-1], d=timestep)
#plt.plot(freqs,np.abs(ft))
#plt.xlim(0,1000)



#---------------------------------------------------------------------------


#Get complex power spectrum. This contains phase info that will be used to calculate phase differences
nps = 128
fspec, tspec, powercA = signal.spectrogram(wfA, fs, nperseg=nps,noverlap=nps/2,window='hann',return_onesided=True,mode='complex')
fspec, tspec, powercB = signal.spectrogram(wfB, fs, nperseg=nps,noverlap=nps/2,window='hann',return_onesided=True,mode='complex')


plt.plot(fspec,np.abs(powercA[:,10]))


cohmin = 0.9  #Best to limit bad coherence values at the onset. Otherwise get a lot of salt/pepper noise in final result


#NOTE: + sense of phase (rad) defined as pointing towards center of potential of "powercA"
g,coh,phase = correlation_analysis.interferometric_coherence_2D(powercA,powercB,3)



goo = coh < cohmin
coh[goo] = float("nan")
phase[goo] = float("nan")




#receiver_spacing = 2.26 #meters -- effective length of spaced receiver
d12 = 5 #probe spacing b/t V1 and V2
receiver_spacing = np.cos(np.radians(45))*d12 #meters -- effective length of spaced receiver


fkpowspec = 0 
kvals = 0
fvals = 0
#Plot full res results 

del(fkpowspec, kvals, fvals)
fkpowspec, kvals, fvals = interf.inter_fvsk(np.abs(powercA),tspec,fspec, 
                                         phase,tspec,fspec,
                                         receiver_spacing,
                                         mean_max='max',
                                         nkbins=200,
                                         klim=[-1,1])




#Find k-location of max power for every frequency
pmaxvals = np.zeros(np.shape(fkpowspec))
pmaxvals[:,:] = 1
for f in range(len(fspec)):
    if np.nanmax(fkpowspec[f,:]):
        cond = fkpowspec[f,:] == np.nanmax(fkpowspec[f,:])
    pmaxvals[f,:] = pmaxvals[f,:]*cond


#Plot the results (freq vs k)
vr = [-40,-10]
yr = [0,100]
kr = [-1,1]
titlegoo = 'Artificial wave test'
#xr = [tz-3*nsec,tz+3*nsec]


fig,axs = plt.subplots(6)
ps.plot_spectrogram(tspec,fspec,np.abs(powercA),vr=vr,yr=yr, yscale='linear',ax=axs[0],xlabel='time(s)',ylabel='power\nf(Hz)',title=titlegoo)
ps.plot_spectrogram(tspec,fspec,np.abs(powercB),vr=vr,yr=yr, yscale='linear',ax=axs[1],xlabel='time(s)',ylabel='power\nf(Hz)')
ps.plot_spectrogram(tspec,fspec,coh,vr=[0,1],zscale='linear',yr=yr,yscale='linear',ax=axs[2],xlabel='time(s)',ylabel='Coherence**2\nf(Hz)')
ps.plot_spectrogram(tspec,fspec,np.abs(np.degrees(phase)),vr=[0,5],zscale='linear',yr=yr,yscale='linear',ax=axs[3],xlabel='time(s)',ylabel='|Phase|(0-180deg)\nf(Hz)',cmap='Spectral')
ps.plot_spectrogram(kvals,fvals,fkpowspec,xr=kr,yr=yr,vr=vr,yscale='linear',ax=axs[4],minzval=0,maxzval=np.nanmax(fkpowspec),xlabel='k(rad/m)',ylabel='f(Hz)')
ps.plot_spectrogram(kvals,fvals,pmaxvals,vr=[0,2],xr=kr,yr=yr,zscale='linear',yscale='linear',ax=axs[5],minzval=-120,maxzval=10,xlabel='k(rad/m)',ylabel='f(Hz)',cmap='Greys')



#plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
#cax = plt.axes((0.85,0.1,0.075,0.8))
#plt.colorbar(cax=cax)
kvec_sim = np.sqrt(vals['kx']**2 + vals['ky']**2)
kvec_sim
#[print('|k| = ' + str(i)) for i in kvec_sim]


print('h')







#TEST COHERENCE CALCULATION 

ctst,ptst,ttst,ftst = cohtst(wfA,wfB,tdat,fs,nperseg=1024,coh_min=0)

print('h')














