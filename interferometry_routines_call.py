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

#---interferometry pairs
#NOTE: positive sense of phase defined as pointing towards center of potential of "wfA"
#e.g. y-hat':  wfA = wf32; wfB = wf14
#So, for Endurance the interferometry pairs are 
#   vAstr = 'VLF32D' and vBstr = 'VLF41D' (which will be flipped in sign)
#or vAstr = 'VLF13D' and vBstr = 'VLF24D' (which will be flipped in sign)

#x-hat'
vAstrx = 'VLF13D'
vBstrx = 'VLF24D'

#y-hat'
vAstry = 'VLF32D'
vBstry = 'VLF41D'



#Load Endurance waveform
vAx = EFL(vAstrx)
vBx = EFL(vBstrx)
wfAx, tdatx = vAx.load_data_gainphase_corrected()
wfBx, tdatx = vBx.load_data_gainphase_corrected()
vAy = EFL(vAstry)
vBy = EFL(vBstry)
wfAy, tdaty = vAy.load_data_gainphase_corrected()
wfBy, tdaty = vBy.load_data_gainphase_corrected()



#Phase flip probe pair for "sense" consistency 
# 41 --> 14  and  24 --> 42
gooA = (vAstrx == 'VLF24D' or vAstrx == 'VLF41D')
if gooA:
    wfAx *= -1
    vAstrx = '-' + vAstrx
gooB = (vBstrx == 'VLF24D' or vBstrx == 'VLF41D')
if gooB:
    wfBx *= -1
    vBstrx = '-' + vBstrx

gooA = (vAstry == 'VLF24D' or vAstry == 'VLF41D')
if gooA:
    wfAy *= -1
    vAstry = '-' + vAstry
gooB = (vBstry == 'VLF24D' or vBstry == 'VLF41D')
if gooB:
    wfBy *= -1
    vBstry = '-' + vBstry



#sample rate
fs = vAx.chnspecs['fs']



#---------------------------------------------------------------------------


#Get complex power spectrum. This contains phase info that will be used to calculate phase differences
nps = 1024
fspecx, tspecx, powercAx = signal.spectrogram(wfAx, fs, nperseg=nps,noverlap=nps/2,window='hann',return_onesided=True,mode='complex')
fspecx, tspecx, powercBx = signal.spectrogram(wfBx, fs, nperseg=nps,noverlap=nps/2,window='hann',return_onesided=True,mode='complex')
fspecy, tspecy, powercAy = signal.spectrogram(wfAy, fs, nperseg=nps,noverlap=nps/2,window='hann',return_onesided=True,mode='complex')
fspecy, tspecy, powercBy = signal.spectrogram(wfBy, fs, nperseg=nps,noverlap=nps/2,window='hann',return_onesided=True,mode='complex')



cohmin = 0.5  #Best to limit bad coherence values at the onset. Otherwise get a lot of salt/pepper noise in final result

#Reduce data to time range of interest [tz to tz+nsec]. 
#(e.g. select wave packet of interest)
#--Bernstein waves on upleg
#vr = [-45,-20]
#ys = 'linear'
#yr = [5500,7500]
#kr = [-50,50]
#tz = 140 
#nsec = 20
#tz = 160 
#nsec = 1


#--unknown power just prior to two-stream
#vr = [-45,-20]
#ys = 'linear'
#yr = [0,2000]
#kr = [-50,50]
#tz = 882  
#nsec = 4
#--attitude adjustment maneuver (returns noise due to no coherence)
#tz = 595  
#nsec = 2
#--two-stream
#vr = [-45,-20]
#ys = 'linear'
#yr = [0,500]
#kr = [-50,50]
#tz = 888
#nsec = 2
#--VLF hiss
vr = [-45,-20]
ys = 'linear'
yr = [5000,9000]
kr = [-25,25]
tz = 472.
nsec = 0.5



#NOTE: + sense of phase defined as pointing towards center of potential of "powercA"
#Nval = 3
Nval = 10
gx,cohx,phasex = correlation_analysis.interferometric_coherence_2D(powercAx,powercBx,Nval)
gy,cohy,phasey = correlation_analysis.interferometric_coherence_2D(powercAy,powercBy,Nval)



goo = cohx < cohmin
cohx[goo] = float("nan")
phasex[goo] = float("nan")
phasex = np.degrees(phasex)

goo = cohy < cohmin
cohy[goo] = float("nan")
phasey[goo] = float("nan")
phasey = np.degrees(phasey)



#Reduce arrays to desired timerange
goo, powercAxz, tspecxz = ps.slice_spectrogram(tz,tspecx,np.abs(powercAx),nsec)
goo, powercBxz, tpowercBxz = ps.slice_spectrogram(tz,tspecx,np.abs(powercAx),nsec)
goo, cohxz, ttmp = ps.slice_spectrogram(tz,tspecx,cohx,nsec)
goo, phasexz, ttmp = ps.slice_spectrogram(tz,tspecx,phasex,nsec)

goo, powercAyz, tspecyz = ps.slice_spectrogram(tz,tspecy,np.abs(powercAy),nsec)
goo, powercByz, tpowercByz = ps.slice_spectrogram(tz,tspecy,np.abs(powercAy),nsec)
goo, cohyz, ttmp = ps.slice_spectrogram(tz,tspecy,cohy,nsec)
goo, phaseyz, ttmp = ps.slice_spectrogram(tz,tspecy,phasey,nsec)


#Endurance has ~3m booms, so the 
receiver_spacing = 2.1 #meters -- effective length of spaced receiver (=3*cos(45))


#Plot full res results 

fkpowspecx, kvalsx, fvalsx = interf.inter_fvsk(np.abs(powercAxz),tspecxz,fspecx, 
                                         phasexz,tspecxz,fspecx,
                                         receiver_spacing,
                                         mean_max='max',
                                         nkbins=200,
                                         klim=[-50,50])
fkpowspecy, kvalsy, fvalsy = interf.inter_fvsk(np.abs(powercAyz),tspecyz,fspecy, 
                                         phaseyz,tspecyz,fspecy,
                                         receiver_spacing,
                                         mean_max='max',
                                         nkbins=200,
                                         klim=[-40,40])



#Find k-location of max power for every frequency
pmaxvalsx = np.zeros(np.shape(fkpowspecx))
pmaxvalsx[:,:] = 1
for f in range(len(fspecx)):
    if np.nanmax(fkpowspecx[f,:]):
        cond = fkpowspecx[f,:] == np.nanmax(fkpowspecx[f,:])
    else:
        cond = 1
    pmaxvalsx[f,:] = pmaxvalsx[f,:]*cond

pmaxvalsy = np.zeros(np.shape(fkpowspecy))
pmaxvalsy[:,:] = 1
for f in range(len(fspecy)):
    if np.nanmax(fkpowspecy[f,:]):
        cond = fkpowspecy[f,:] == np.nanmax(fkpowspecy[f,:])
    else:
        cond = 1
    pmaxvalsy[f,:] = pmaxvalsy[f,:]*cond



#Plot the results (freq vs k)
titlegoo = 'slice from '+ str(tz) + '-' + str(tz + nsec) + ' sec\n' #+ vAstrx + ' and ' + vBstrx
xr = [tz-3*nsec,tz+3*nsec]

xptitle = 'V1V3-V4V2'
yptitle = 'V3V2-V1V4'


fig,axs = plt.subplots(6)
ps.plot_spectrogram(tspecx,fspecx,np.abs(powercAx),vr=vr,yr=yr,xr=xr, yscale=ys,ax=axs[0],xlabel='time(s)',ylabel='power\nf(Hz)',title=titlegoo)
ps.plot_spectrogram(tspecx,fspecx,np.abs(powercBx),vr=vr,yr=yr,xr=xr, yscale=ys,ax=axs[1],xlabel='time(s)',ylabel='power\nf(Hz)')
ps.plot_spectrogram(tspecx,fspecx,cohx,vr=[0,1],zscale='linear',xr=xr,yr=yr,yscale=ys,ax=axs[2],xlabel='time(s)',ylabel='Coherence**2\nf(Hz)')
ps.plot_spectrogram(tspecx,fspecx,np.abs(phasex),vr=[0,180],zscale='linear',xr=xr,yr=yr,yscale=ys,ax=axs[3],xlabel='time(s)',ylabel='|Phase|(0-180deg)\nf(Hz)',cmap='Spectral')
ps.plot_spectrogram(kvalsx,fvalsx,pmaxvalsx,vr=[0,1],xr=kr,yr=yr,zscale='linear',yscale=ys,ax=axs[4],minzval=0,maxzval=1,xlabel='k(rad/m)',ylabel="kx'\n"+xptitle+'\nf(Hz)',cmap='Greys')
ps.plot_spectrogram(kvalsx,fvalsx,fkpowspecx,vr=vr,xr=kr,yr=yr,yscale=ys,ax=axs[4],minzval=-120,maxzval=10,xlabel='k(rad/m)',ylabel="kx'\n"+xptitle+'\nf(Hz)',alpha=0.5)
ps.plot_spectrogram(kvalsy,fvalsy,pmaxvalsy,vr=[0,1],xr=kr,yr=yr,zscale='linear',yscale=ys,ax=axs[5],minzval=0,maxzval=1,xlabel='k(rad/m)',ylabel="ky'\n"+yptitle+'\nf(Hz)',cmap='Greys')
ps.plot_spectrogram(kvalsy,fvalsy,fkpowspecy,vr=vr,xr=kr,yr=yr,yscale=ys,ax=axs[5],minzval=-120,maxzval=10,xlabel='k(rad/m)',ylabel="ky'\n"+yptitle+'\nf(Hz)',alpha=0.5)


for i in range(4): axs[i].axvline(tz)
for i in range(4): axs[i].axvline(tz+nsec)

#Plot the limiting k-vector value based on receiver spacing. 
#i.e. location when wavelength equals interferometry receiver spacing
klim = 2*np.pi / receiver_spacing
for i in range(3): axs[i+3].axvline(klim)
for i in range(3): axs[i+3].axvline(-klim)


#----------------
#----------------
#----------------

fig,axs = plt.subplots(2)
ps.plot_spectrogram(kvalsx,fvalsx,pmaxvalsx,vr=[0,1],xr=kr,yr=yr,zscale='linear',yscale=ys,ax=axs[0],minzval=0,maxzval=1,xlabel='k(rad/m)',ylabel="kx'\n"+xptitle+'\nf(Hz)',cmap='Greys',title=titlegoo)
ps.plot_spectrogram(kvalsx,fvalsx,fkpowspecx,vr=vr,xr=kr,yr=yr,yscale=ys,ax=axs[0],minzval=-120,maxzval=10,xlabel='k(rad/m)',ylabel="kx'\n"+xptitle+'\nf(Hz)',alpha=0.5)
ps.plot_spectrogram(kvalsy,fvalsy,pmaxvalsy,vr=[0,1],xr=kr,yr=yr,zscale='linear',yscale=ys,ax=axs[1],minzval=0,maxzval=1,xlabel='k(rad/m)',ylabel="ky'\n"+yptitle+'\nf(Hz)',cmap='Greys')
ps.plot_spectrogram(kvalsy,fvalsy,fkpowspecy,vr=vr,xr=kr,yr=yr,yscale=ys,ax=axs[1],minzval=-120,maxzval=10,xlabel='k(rad/m)',ylabel="ky'\n"+yptitle+'\nf(Hz)',alpha=0.5)

axs[0].axvline(klim)
axs[0].axvline(-klim)
axs[1].axvline(klim)
axs[1].axvline(-klim)




print('h')
print('h')
print('h')
print('h')














