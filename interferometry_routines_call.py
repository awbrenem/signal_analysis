"""
Example setup for calls to interferometry_routines.py 

(for testing see interferometry_fvsk_spec_test.py and 
antenna_finite_wavelength_response.py)

NOTE: this technique only makes sense for parallel, spaced antennas. 
If the antennas are perpendicular, the ~90 deg phase shift in waves will be interpreted as a 
large k-value. 

NOTE: phase identification can be difficult if:

1) the wave is propagating nearly along one of the interferometry axes (e.g. if along the y-hat' direction then V32 and V14 will
measure very little signal, meaning the waveform can be noise-dominated)

2) if wavelength is on the order or less than boom length than results can't be trusted.


                       (y-hat)
                         V3
                      /  |  \             
           (y-hat') x    |     x (x-hat')         
                  /      |       \        
               V2--------O--------V1 (x-hat)    
                  \      |       /
                    x    |     x
                      \  |  /
                         V4
          

x-points represent the centers of potential of the interferometry diagonals          

Coordinate system (system of input test wave)
x-hat --> E12 = V1 - V2 direction (positive to right)
y-hat --> E34 = V3 - V4 direction (positive upwards)

This code outputs the spectrum of k-values in the kx' and ky' directions, where
x'-hat --> Ex' = V1V3x - V4V2x (45 deg inclined from x-hat)
y'-hat --> Ey' = V3V2x - V1V4x (45 deg inclined from y-hat)



"""

import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/Endurance/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
from end_fields_loader import Endurance_Fields_Loader as EFL
import end_data_loader
from scipy import signal
import numpy as np 
import interferometry_routines as interf
import correlation_analysis
import plot_spectrogram as ps
import matplotlib.pyplot as plt



#-------------------------------------------------------------
#Example 1: Do a freq (yaxis) vs k-value (xaxis) vs |E|^2 (zaxis) interferometry analysis
#using Endurance data
#-------------------------------------------------------------

#---Choose interferometry pair
#NOTE: positive sense of phase defined as pointing towards center of potential of "wfA"
#e.g. y-hat':  wfA = wf32; wfB = wf14
#So, for Endurance the interferometry pairs are 
#   vAstr = 'VLF32D' and vBstr = 'VLF41D' (which will be flipped in sign)
#or vAstr = 'VLF13D' and vBstr = 'VLF24D' (which will be flipped in sign)

#x-hat'
#vAstr = 'VLF13D'
#vBstr = 'VLF24D'

#y-hat'
vAstr = 'VLF32D'
vBstr = 'VLF41D'

#main booms
#vAstr = 'VLF12D'
#vBstr = 'VLF34D'


#Load Endurance waveform
vA = EFL(vAstr)
vB = EFL(vBstr)
wfA, tdat = vA.load_data_gainphase_corrected()
wfB, tdat = vB.load_data_gainphase_corrected()



#Phase flip probe pair for "sense" consistency 
# 41 --> 14  and  24 --> 42
gooA = (vAstr == 'VLF24D' or vAstr == 'VLF41D')
if gooA:
    wfA *= -1
    vAstr = '-' + vAstr
gooB = (vBstr == 'VLF24D' or vBstr == 'VLF41D')
if gooB:
    wfB *= -1
    vBstr = '-' + vBstr


#sample rate
fs = vA.chnspecs['fs']



#---------------------------------------------------------------------------


#Get complex power spectrum. This contains phase info that will be used to calculate phase differences
nps = 1024
fspec, tspec, powercA = signal.spectrogram(wfA, fs, nperseg=nps,noverlap=nps/2,window='hann',return_onesided=True,mode='complex')
fspec, tspec, powercB = signal.spectrogram(wfB, fs, nperseg=nps,noverlap=nps/2,window='hann',return_onesided=True,mode='complex')



cohmin = 0.5  #Best to limit bad coherence values at the onset. Otherwise get a lot of salt/pepper noise in final result

#Reduce data to time range of interest [tz to tz+nsec]. 
#(e.g. select wave packet of interest)
#--Bernstein waves on upleg
#tz = 140 
#nsec = 20
#tz = 160 
#nsec = 1
#k-hat' = -15 @6300 Hz

#--unknown power just prior to two-stream
#tz = 882  
#nsec = 4
#--attitude adjustment maneuver (returns noise due to no coherence)
#tz = 595  
#nsec = 2
#--two-stream
#tz = 887
#nsec = 2
#--VLF hiss
tz = 472.3
nsec = 0.5



#NOTE: + sense of phase defined as pointing towards center of potential of "powercA"
#Nval = 3
Nval = 10
g,coh,phase = correlation_analysis.interferometric_coherence_2D(powercA,powercB,Nval)



goo = coh < cohmin
coh[goo] = float("nan")
phase[goo] = float("nan")

phase = np.degrees(phase)



#Reduce arrays to desired timerange
goo, powercAz, tspecz = ps.slice_spectrogram(tz,tspec,np.abs(powercA),nsec)
goo, powercBz, tpowercBz = ps.slice_spectrogram(tz,tspec,np.abs(powercA),nsec)
goo, cohz, ttmp = ps.slice_spectrogram(tz,tspec,coh,nsec)
goo, phasez, ttmp = ps.slice_spectrogram(tz,tspec,phase,nsec)





#Endurance has ~3m booms, so the 
receiver_spacing = 2.1 #meters -- effective length of spaced receiver (=3*cos(45))


fkpowspec = 0 
kvals = 0
fvals = 0
#Plot full res results 

del(fkpowspec, kvals, fvals)
fkpowspec, kvals, fvals = interf.inter_fvsk(np.abs(powercAz),tspecz,fspec, 
                                         phasez,tspecz,fspec,
                                         receiver_spacing,
                                         mean_max='max',
                                         nkbins=400,
                                         klim=[-50,50])


#Find k-location of max power for every frequency
pmaxvals = np.zeros(np.shape(fkpowspec))
pmaxvals[:,:] = 1
for f in range(len(fspec)):
    if np.nanmax(fkpowspec[f,:]):
        cond = fkpowspec[f,:] == np.nanmax(fkpowspec[f,:])
    pmaxvals[f,:] = pmaxvals[f,:]*cond



#Plot the results (freq vs k)
vr = [-45,-20]
yr = [5000,9000]
kr = [-20,20]
titlegoo = 'slice from '+ str(tz) + '-' + str(tz + nsec) + ' sec\n' + vAstr + ' and ' + vBstr
xr = [tz-3*nsec,tz+3*nsec]


fig,axs = plt.subplots(6)
ps.plot_spectrogram(tspec,fspec,np.abs(powercA),vr=vr,yr=yr,xr=xr, yscale='linear',ax=axs[0],xlabel='time(s)',ylabel='power\nf(Hz)',title=titlegoo)
ps.plot_spectrogram(tspec,fspec,np.abs(powercB),vr=vr,yr=yr,xr=xr, yscale='linear',ax=axs[1],xlabel='time(s)',ylabel='power\nf(Hz)')
ps.plot_spectrogram(tspec,fspec,coh,vr=[0,1],zscale='linear',xr=xr,yr=yr,yscale='linear',ax=axs[2],xlabel='time(s)',ylabel='Coherence**2\nf(Hz)')
ps.plot_spectrogram(tspec,fspec,np.abs(phase),vr=[0,180],zscale='linear',xr=xr,yr=yr,yscale='linear',ax=axs[3],xlabel='time(s)',ylabel='|Phase|(0-180deg)\nf(Hz)',cmap='Spectral')
ps.plot_spectrogram(kvals,fvals,fkpowspec,vr=vr,xr=kr,yr=yr,yscale='linear',ax=axs[4],minzval=-120,maxzval=10,xlabel='k(rad/m)',ylabel='f(Hz)')
ps.plot_spectrogram(kvals,fvals,pmaxvals,vr=[0,2],xr=kr,yr=yr,zscale='linear',yscale='linear',ax=axs[5],minzval=-120,maxzval=10,xlabel='k(rad/m)',ylabel='f(Hz)',cmap='Greys')

for i in range(3): axs[i].axvline(tz)
for i in range(3): axs[i].axvline(tz+nsec)


fig,axs = plt.subplots(2)
ps.plot_spectrogram(kvals,fvals,fkpowspec,vr=vr,xr=kr,yr=yr,yscale='linear',ax=axs[0],minzval=-120,maxzval=10,xlabel='k(rad/m)',ylabel='f(Hz)')
ps.plot_spectrogram(kvals,fvals,pmaxvals,vr=[0,2],xr=kr,yr=yr,zscale='linear',yscale='linear',ax=axs[1],minzval=-120,maxzval=10,xlabel='k(rad/m)',ylabel='f(Hz)',cmap='Greys')



#plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
#cax = plt.axes((0.85,0.1,0.075,0.8))
#plt.colorbar(cax=cax)

print('h')



print('h')




#TEST COHERENCE CALCULATION 

ctst,ptst,ttst,ftst = cohtst(wfA,wfB,tdat,fs,nperseg=1024,coh_min=0)

print('h')














