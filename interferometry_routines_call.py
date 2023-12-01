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



#-------------------------------------------------------------
#Method 1: Do a freq (yaxis) vs k-value (xaxis) vs |E|^2 (zaxis) interferometry analysis using Endurance data
#-------------------------------------------------------------
#Method 2: Do a freq (xaxis) vs k-value (yaxis) interferometry analysis using Endurance data
#-------------------------------------------------------------
--> results from both methods are plotted together at end. 

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


#---------------------------------------------------------
#Reduce data to time range of interest [tz to tz+nsec]  (e.g. select wave packet of interest)
#--VLF hiss
#tz = 475.
#nsec = 0.5

#--VLF hiss
tz = 145.
nsec = 0.5

vr = [-45,-20]
ys = 'linear'
fr = [5000,9000]
kr = [-2,2]




cohmin = 0.5  #Best to limit bad coherence values at the onset. Otherwise get a lot of salt/pepper noise in final result


#Endurance has ~3.2 m booms, so the 
receiver_spacing = 2.27 #meters -- effective length of spaced receiver





#-------------------------------------------------------------
#Method 1: Do a freq (yaxis) vs k-value (xaxis) vs |E|^2 (zaxis) interferometry analysis
#using Endurance data
#-------------------------------------------------------------

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



#Create f vs k values for plotting
fkpowspecx, kvalsx, fvalsx, pmaxvalsx = interf.inter_fvsk(np.abs(powercAxz),tspecxz,fspecx, 
                                         np.radians(phasexz),tspecxz,fspecx,
                                         receiver_spacing,
                                         mean_max='max',
                                         nkbins=200,
                                         klim=[-5,5])
fkpowspecy, kvalsy, fvalsy, pmaxvalsy = interf.inter_fvsk(np.abs(powercAyz),tspecyz,fspecy, 
                                         np.radians(phaseyz),tspecyz,fspecy,
                                         receiver_spacing,
                                         mean_max='max',
                                         nkbins=200,
                                         klim=[-5,5])

#Turn k-values into wavelength
wl1 = np.zeros(len(fvalsx))
for i in range(len(fvalsx)):
    tmp = pmaxvalsx[i,:]
    kxgoo = kvalsx[np.where(tmp == 1)]
    tmp = pmaxvalsy[i,:]
    kygoo = kvalsy[np.where(tmp == 1)]
    kmaggoo = np.sqrt(kxgoo**2 + kygoo*2)
    wl1[i] = 2*np.pi/kmaggoo



#-------------------------------------------------------------
#Method 2: Do a freq (xaxis) vs k-value (yaxis) interferometry analysis using Endurance data
#-------------------------------------------------------------


tchunk = 0.1  #sec
cohx2, phasex2, tchunks2, freqs2 = correlation_analysis.cross_spectral_density_spectrogram(wfAx,wfBx,tdatx,fs,tchunk,coh_min=cohmin,nperseg=512)
cohy2, phasey2, tchunks2, freqs2 = correlation_analysis.cross_spectral_density_spectrogram(wfAy,wfBy,tdaty,fs,tchunk,coh_min=cohmin,nperseg=512)


#Reduce arrays to desired timerange
pavgx2, phasearrx2, tarr_phasex2 = ps.slice_spectrogram(tz,tchunks2,phasex2,nsec)
cavgx2, coharrx2, tarr_cohx2 = ps.slice_spectrogram(tz,tchunks2,cohx2,nsec)
powavgx2, powarrx2, tarr_powx2 = ps.slice_spectrogram(tz,tspecx,np.abs(powercAx),nsec)

pavgy2, phasearry2, tarr_phasey2 = ps.slice_spectrogram(tz,tchunks2,phasey2,nsec)
cavgy2, coharry2, tarr_cohy2 = ps.slice_spectrogram(tz,tchunks2,cohy2,nsec)
powavgy2, powarry2, tarr_powy2 = ps.slice_spectrogram(tz,tspecy,np.abs(powercAy),nsec)


#Turn the phase values into k-values 
kx2 = [np.radians(i) / receiver_spacing for i in pavgx2]
ky2 = [np.radians(i) / receiver_spacing for i in pavgy2]
#and then into wavelength perp to Bo
kmag = [np.sqrt(kx2[i]**2 + ky2[i]**2) for i in range(len(kx2))]
wl2 = [2*np.pi/i for i in kmag]

yr = [5000,9000]
vr=[-60,-20]
xrspec = [tz,tz+nsec]




#Plot values from Method 2
fig,axs = plt.subplots(6,2, figsize=(10,10))  #,gridspec_kw={'height_ratios':[1,1,1,1,1,1,1,1,1]})
fig.subplots_adjust(bottom=0.1,right=0.8,left=0.2,top=0.9,hspace=0.4,wspace=0.4)
#fig,axs = plt.subplots(5,2,figsize=(12,7),gridspec_kw={'height_ratios':[1,1,1,1],'width_ratios':[1,0.5]})

ps.plot_spectrogram(tspecxz,fspecx,np.abs(powarrx2),vr=vr,xr=xrspec,yr=yr,yscale='linear',ax=axs[0,0],xlabel='time(sec)',ylabel='power spectrum\nfreq(Hz)')
ps.plot_spectrogram(tarr_cohx2,freqs2,coharrx2,vr=[0.9,1], zscale='linear',xr=xrspec,yr=yr,yscale='linear',ax=axs[1,0],xlabel='time(sec)',ylabel='coherence\nfreq(Hz)')
ps.plot_spectrogram(tarr_phasex2,freqs2,phasearrx2,vr=[-140,140], zscale='linear',xr=xrspec,yr=yr,yscale='linear',ax=axs[2,0],xlabel='time(sec)',ylabel='phase(deg)\nfreq(Hz)')
for i in range(3): axs[i,0].axvline(tz)
for i in range(3): axs[i,0].axvline(tz + nsec)
axs[0,1].plot(fspecx,powarrx2)
axs[0,1].plot(fspecx,powavgx2,'.',color='black')
axs[0,1].set_xlim(yr)
axs[0,1].set_ylim(0,np.nanmax(powarrx2))
axs[0,1].set_ylabel('powermax')
axs[0,1].set_xlabel('freq(Hz)')
axs[1,1].plot(freqs2,coharrx2)
axs[1,1].plot(freqs2,cavgx2,'.',color='black')
axs[1,1].set_xlim(yr)
axs[1,1].set_ylim(0,1)
axs[1,1].set_ylabel('coherence')
axs[1,1].set_xlabel('freq(Hz)')
axs[2,1].plot(freqs2,phasearrx2)
axs[2,1].plot(freqs2,pavgx2,'.',color='black')
axs[2,1].set_xlim(yr)
axs[2,1].set_ylim(-180,180)
axs[2,1].set_ylabel("phase (x')")
axs[2,1].set_xlabel('freq(Hz)')
axs[3,1].plot(freqs2,phasearry2)
axs[3,1].plot(freqs2,pavgy2,'.',color='black')
axs[3,1].set_xlim(yr)
axs[3,1].set_ylim(-180,180)
axs[3,1].set_ylabel("phase (y')")
axs[3,1].set_xlabel('freq(Hz)')
axs[4,1].plot(freqs2,kx2,'.',freqs2,ky2,'.')
axs[4,1].set_xlim(yr)
axs[4,1].set_ylim(min(np.nanmin(kx2),np.nanmin(ky2)),max([np.nanmax(kx2),np.nanmax(ky2)]))
axs[4,1].set_ylabel("kx'(blue)\nky'(orange)")
axs[4,1].set_xlabel('freq(Hz)')
axs[5,1].plot(freqs2,wl2,'.',color='black')
axs[5,1].set_xlim(yr)
axs[5,1].set_ylim(0,50)
axs[5,1].set_ylabel("wavelength(m)\nfrom |k|")
axs[5,1].set_xlabel('freq(Hz)')

fig.delaxes(axs[3,0])
fig.delaxes(axs[4,0])
fig.delaxes(axs[5,0])





#------------------------------------------------------------
#Plot final results
#------------------------------------------------------------

titlegoo = 'slice from '+ str(tz) + '-' + str(tz + nsec) + ' sec\n' #+ vAstrx + ' and ' + vBstrx
xr = [tz-3*nsec,tz+3*nsec]

xptitle = 'V1V3-V4V2'
yptitle = 'V3V2-V1V4'


fig,axs = plt.subplots(4,2,figsize=(12,7),gridspec_kw={'height_ratios':[1,1,1,1],'width_ratios':[1,0.5]})
fig.subplots_adjust(bottom=0.1,right=0.8,left=0.2,top=0.9,hspace=0.6,wspace=0.6)
ps.plot_spectrogram(tspecx,fspecx,np.abs(powercAx),vr=vr,yr=yr,xr=xr, yscale=ys,ax=axs[0,0],xlabel='time(s)',ylabel='power\nf(Hz)',title=titlegoo)
ps.plot_spectrogram(tspecx,fspecx,np.abs(powercBx),vr=vr,yr=yr,xr=xr, yscale=ys,ax=axs[1,0],xlabel='time(s)',ylabel='power\nf(Hz)')
ps.plot_spectrogram(tspecx,fspecx,cohx,vr=[0,1],zscale='linear',xr=xr,yr=yr,yscale=ys,ax=axs[2,0],xlabel='time(s)',ylabel='Coherence\nf(Hz)')
ps.plot_spectrogram(tspecx,fspecx,np.abs(phasex),vr=[0,180],zscale='linear',xr=xr,yr=yr,yscale=ys,ax=axs[3,0],xlabel='time(s)',ylabel='|Phase|(0-180deg)\nf(Hz)',cmap='Spectral')
for i in range(4): axs[i,0].axvline(tz)
for i in range(4): axs[i,0].axvline(tz+nsec)

ps.plot_spectrogram(kvalsx,fvalsx,pmaxvalsx,vr=[0,1],xr=krplot,yr=yr,zscale='linear',yscale=ys,ax=axs[0,1],minzval=0,maxzval=1,xlabel='k(rad/m)',ylabel="kx'\n"+xptitle+'\nf(Hz)',cmap='Greys')
ps.plot_spectrogram(kvalsx,fvalsx,fkpowspecx,vr=vr,xr=krplot,yr=yr,yscale=ys,ax=axs[0,1],minzval=-120,maxzval=10,xlabel='k(rad/m)',ylabel="kx'\n"+xptitle+'\nf(Hz)',alpha=0.5)
ps.plot_spectrogram(kvalsy,fvalsy,pmaxvalsy,vr=[0,1],xr=krplot,yr=yr,zscale='linear',yscale=ys,ax=axs[1,1],minzval=0,maxzval=1,xlabel='k(rad/m)',ylabel="ky'\n"+yptitle+'\nf(Hz)',cmap='Greys')
ps.plot_spectrogram(kvalsy,fvalsy,fkpowspecy,vr=vr,xr=krplot,yr=yr,yscale=ys,ax=axs[1,1],minzval=-120,maxzval=10,xlabel='k(rad/m)',ylabel="ky'\n"+yptitle+'\nf(Hz)',alpha=0.5)
#Plot the limiting k-vector value where short wavelength effects start to occur. 
#i.e. location when wavelength equals about twice the interferometry receiver spacing
klim = 2*np.pi / (2 * receiver_spacing)
for i in range(2): axs[i,1].axvline(klim)
for i in range(2): axs[i,1].axvline(-klim)

#Oplot the results from the initial (1d) analysis (Method 2)
axs[0,1].plot(kx2,freqs2,'*',color='black')
axs[1,1].plot(ky2,freqs2,'*',color='black')

axs[2,1].plot(wl1,fspecx,'.',color='black',markersize=2)
axs[2,1].plot(wl2,freqs2,'*',color='black')
axs[2,1].set_ylim(yr)
axs[2,1].set_xscale('log')
axs[2,1].set_xlim(0.1,2000)
axs[2,1].set_xlabel("wavelength(m)\n(from |k|)")


fig.delaxes(axs[3,1])



#--------------------------------
#Version focused on the k-spectra
#--------------------------------

fig,axs = plt.subplots(2)
ps.plot_spectrogram(kvalsx,fvalsx,pmaxvalsx,vr=[0,1],xr=kr,yr=fr,zscale='linear',yscale=ys,ax=axs[0],minzval=0,maxzval=1,xlabel='k(rad/m)',ylabel="kx'\n"+xptitle+'\nf(Hz)',cmap='Greys',title=titlegoo)
ps.plot_spectrogram(kvalsx,fvalsx,fkpowspecx,vr=vr,xr=kr,yr=fr,yscale=ys,ax=axs[0],minzval=-120,maxzval=10,xlabel='k(rad/m)',ylabel="kx'\n"+xptitle+'\nf(Hz)',alpha=0.5)
ps.plot_spectrogram(kvalsy,fvalsy,pmaxvalsy,vr=[0,1],xr=kr,yr=fr,zscale='linear',yscale=ys,ax=axs[1],minzval=0,maxzval=1,xlabel='k(rad/m)',ylabel="ky'\n"+yptitle+'\nf(Hz)',cmap='Greys')
ps.plot_spectrogram(kvalsy,fvalsy,fkpowspecy,vr=vr,xr=kr,yr=fr,yscale=ys,ax=axs[1],minzval=-120,maxzval=10,xlabel='k(rad/m)',ylabel="ky'\n"+yptitle+'\nf(Hz)',alpha=0.5)

axs[0].axvline(klim)
axs[0].axvline(-klim)
axs[1].axvline(klim)
axs[1].axvline(-klim)

#Oplot the results from the initial (1d) analysis
axs[0].plot(kx2,freqs2,'*',color='black')
axs[1].plot(ky2,freqs2,'*',color='black')



print('h')
print('h')
print('h')
print('h')

