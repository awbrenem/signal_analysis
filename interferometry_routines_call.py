"""
Example setup for calls to interferometry_routines.py 

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



"""
Example 1: Do a freq (yaxis) vs k-value (xaxis) vs |E|^2 (zaxis) interferometry analysis
using Endurance data
"""

#Choose boom pairs and load waveforms
vA = EFL('VLF32D')
vB = EFL('VLF41D')
#vA = EFL('VLF12D')
#vB = EFL('VLF34D')
wfA, tdat = vA.load_data_gainphase_corrected()
wfB, tdat = vB.load_data_gainphase_corrected()

#Phase flip probe pair for "sense" consistency (e.g. 41 --> 14)
wfB = -wfB

#sample rate
fs = vA.chnspecs['fs']

#Get complex power spectrum. This contains phase info that will be used to calculate phase differences
nps = 512
fspec, tspec, powercA = signal.spectrogram(wfA, fs, nperseg=nps,noverlap=nps/2,window='hann',return_onesided=True,mode='complex')
fspec, tspec, powercB = signal.spectrogram(wfB, fs, nperseg=nps,noverlap=nps/2,window='hann',return_onesided=True,mode='complex')

#Get cross spectrum 
#tchunk = 0.2  #sec
tchunk = 0.05  #sec

cohmin = 0.5  #Best to limit bad coherence values at the onset. Otherwise get a lot of salt/pepper noise in final result
coh, phase, tchunks, freqs = correlation_analysis.cross_spectral_density_spectrogram(wfA,wfB,
                                                                                     tdat,fs,tchunk,
                                                                                     coh_min=cohmin,
                                                                                     nperseg=512)


#Reduce data to time range of interest [tz to tz+nsec]. 
#(e.g. select wave packet of interest)
#--Bernstein waves on upleg
tz = 140 
nsec = 6
#--unknown power just prior to two-stream
#tz = 882  
#nsec = 4
#--attitude adjustment maneuver (returns noise due to no coherence)
#tz = 595  
#nsec = 2
#--two-stream
#tz = 887
#nsec = 2


pavg, phasearr, tphasearr = ps.slice_spectrogram(tz,tchunks,phase,nsec)
cavg, coharr, tgoo = ps.slice_spectrogram(tz,tchunks,coh,nsec)
powavg, powarr, tpowarr = ps.slice_spectrogram(tz,tspec,np.abs(powercA),nsec)

receiver_spacing = 2.26 #meters -- effective length of spaced receiver



#Get the f-k-power interferometry results
fkpowspec, kvals, fvals = interf.inter_fvsk(np.abs(powarr),tpowarr,fspec, 
                                         phasearr,tphasearr,freqs,
                                         receiver_spacing,
                                         mean_max='mean',
                                         klim=[-4,4])



#Plot the results (freq vs k)
vr = [-45,-20]
yr = [5000,8000]
kr = [-3,3]
titlegoo = 'slice from '+ str(tz) + '-' + str(tz + nsec) + ' sec'
xr = [tz-20,tz+20]

fig,axs = plt.subplots(4)
ps.plot_spectrogram(tspec,fspec,np.abs(powercA),vr=vr,yr=yr,xr=xr, yscale='linear',ax=axs[0],xlabel='time(s)',ylabel='power\nf(Hz)',title=titlegoo)
ps.plot_spectrogram(tchunks,freqs,coh,vr=[0,1],zscale='linear',xr=xr,yr=yr,yscale='linear',ax=axs[1],xlabel='time(s)',ylabel='Coherence**2\nf(Hz)')
ps.plot_spectrogram(tchunks,freqs,np.abs(phase),vr=[0,180],zscale='linear',xr=xr,yr=yr,yscale='linear',ax=axs[2],xlabel='time(s)',ylabel='|Phase|(0-180deg)\nf(Hz)',cmap='PiYG')
ps.plot_spectrogram(kvals,freqs,fkpowspec,vr=vr,xr=kr,yr=yr,yscale='linear',ax=axs[3],minzval=-120,maxzval=10,xlabel='k(rad/m)',ylabel='f(Hz)')
for i in range(3): axs[i].axvline(tz)
for i in range(3): axs[i].axvline(tz+nsec)




print('h')
