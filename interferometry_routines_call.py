"""
Example setup for calls to interferometry_routines.py 

NOTE: this technique only makes sense for parallel, spaced antennas. 
If the antennas are perpendicular, the ~90 deg phase shift in waves will be interpreted as a 
large k-value. 
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

#Choose boom pairs and load waveforms
vAstr = 'VLF32D'
vBstr = 'VLF41D'
#vAstr = 'VLF13D'
#vBstr = 'VLF24D'
#vAstr = 'VLF12D'
#vBstr = 'VLF34D'

#Load Endurance waveform
vA = EFL(vAstr)
vB = EFL(vBstr)
wfA, tdat = vA.load_data_gainphase_corrected()
wfB, tdat = vB.load_data_gainphase_corrected()


#Phase flip probe pair for "sense" consistency (e.g. 41 --> 14)
wfB = -wfB
vBstr = '-'+vBstr
#sample rate
fs = vA.chnspecs['fs']



#---------------------------------------------------------------------------


#Get complex power spectrum. This contains phase info that will be used to calculate phase differences
nps = 512
fspec, tspec, powercA = signal.spectrogram(wfA, fs, nperseg=nps,noverlap=nps/2,window='hann',return_onesided=True,mode='complex')
fspec, tspec, powercB = signal.spectrogram(wfB, fs, nperseg=nps,noverlap=nps/2,window='hann',return_onesided=True,mode='complex')



cohmin = 0.75  #Best to limit bad coherence values at the onset. Otherwise get a lot of salt/pepper noise in final result

#Reduce data to time range of interest [tz to tz+nsec]. 
#(e.g. select wave packet of interest)
#--Bernstein waves on upleg
#tz = 140 
#nsec = 20
tz = 160 
nsec = 1

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
#tz = 472.3
#nsec = 0.5



g,coh,phase = correlation_analysis.interferometric_coherence_2D(powercA,powercB,3)



goo = coh < 0.5
coh[goo] = float("nan")
phase[goo] = float("nan")

phase = np.degrees(phase)



#Reduce arrays to desired timerange
goo, powercAz, tspecz = ps.slice_spectrogram(tz,tspec,np.abs(powercA),nsec)
goo, powercBz, tpowercBz = ps.slice_spectrogram(tz,tspec,np.abs(powercA),nsec)
goo, cohz, ttmp = ps.slice_spectrogram(tz,tspec,coh,nsec)
goo, phasez, ttmp = ps.slice_spectrogram(tz,tspec,phase,nsec)






receiver_spacing = 2.26 #meters -- effective length of spaced receiver


fkpowspec = 0 
kvals = 0
fvals = 0
#Plot full res results 

del(fkpowspec, kvals, fvals)
fkpowspec, kvals, fvals = interf.inter_fvsk(np.abs(powercAz),tspecz,fspec, 
                                         phasez,tspecz,fspec,
                                         receiver_spacing,
                                         mean_max='max',
                                         nkbins=200,
                                         klim=[-2,2])


#Find k-location of max power for every frequency
#---QUITE SURE THIS LOOP IS WORKING
pmaxvals = np.zeros(np.shape(fkpowspec))
pmaxvals[:,:] = 1
for f in range(len(fspec)):
    cond = fkpowspec[f,:] == np.nanmax(fkpowspec[f,:])
    pmaxvals[f,:] = pmaxvals[f,:]*cond


#Plot the results (freq vs k)
vr = [-45,-20]
yr = [5000,7500]
kr = [-0.5,0.5]
titlegoo = 'slice from '+ str(tz) + '-' + str(tz + nsec) + ' sec\n' + vAstr + ' and ' + vBstr
xr = [tz-3*nsec,tz+3*nsec]


fig,axs = plt.subplots(6)
ps.plot_spectrogram(tspec,fspec,np.abs(powercA),vr=vr,yr=yr,xr=xr, yscale='linear',ax=axs[0],xlabel='time(s)',ylabel='power\nf(Hz)',title=titlegoo)
ps.plot_spectrogram(tspec,fspec,np.abs(powercB),vr=vr,yr=yr,xr=xr, yscale='linear',ax=axs[1],xlabel='time(s)',ylabel='power\nf(Hz)')
ps.plot_spectrogram(tspec,fspec,coh,vr=[0,1],zscale='linear',xr=xr,yr=yr,yscale='linear',ax=axs[2],xlabel='time(s)',ylabel='Coherence**2\nf(Hz)')
#ps.plot_spectrogram(tspec,fspec,np.abs(phase),vr=[0,180],zscale='linear',xr=xr,yr=yr,yscale='linear',ax=axs[3],xlabel='time(s)',ylabel='|Phase|(0-180deg)\nf(Hz)',cmap='PiYG')
ps.plot_spectrogram(tspec,fspec,np.abs(phase),vr=[0,180],zscale='linear',xr=xr,yr=yr,yscale='linear',ax=axs[3],xlabel='time(s)',ylabel='|Phase|(0-180deg)\nf(Hz)',cmap='Spectral')
ps.plot_spectrogram(kvals,fvals,fkpowspec,vr=vr,xr=kr,yr=yr,yscale='linear',ax=axs[4],minzval=-120,maxzval=10,xlabel='k(rad/m)',ylabel='f(Hz)')
ps.plot_spectrogram(kvals,fvals,pmaxvals,vr=[0,2],xr=kr,yr=yr,zscale='linear',yscale='linear',ax=axs[5],minzval=-120,maxzval=10,xlabel='k(rad/m)',ylabel='f(Hz)',cmap='Greys')

for i in range(3): axs[i].axvline(tz)
for i in range(3): axs[i].axvline(tz+nsec)

#plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
#cax = plt.axes((0.85,0.1,0.075,0.8))
#plt.colorbar(cax=cax)

print('h')



print('h')




#TEST COHERENCE CALCULATION 

ctst,ptst,ttst,ftst = cohtst(wfA,wfB,tdat,fs,nperseg=1024,coh_min=0)

print('h')














