"""
Routines for performing interferometry analysis 

"""


def inter_fvsk(powspec,tpowspec,fpowspec,
               phasespec,tphasespec,fphasespec,
               receiver_spacing,
               nkbins=100,klim=[-2,2],
               mean_max='max'):
    

    """
    Interferometry analysis f vs k plots, with wave power on z-axis. 
    (see interferometry_routines_call.py for example usage)

    Reference:
    Lalti, A., Khotyaintsev, Y. V., & Graham, D. B. (2023). 
    Short-Wavelength Electrostatic Wave Measurement Using MMS Spacecraft
    JGR: Space Phys https://doi. org/10.1029/2022JA031150

    Idea: use a coherence analysis to find the phase difference as a function of freq/time.
    This gives us the k-vector spectrum from k*d = delta-phase. 

    For a select time interval of interest, bin the k-values for each frequency to make 
    a freq vs k plot, with the z-values being wave power. 

    Can compare results to specific wave dispersion relations of f(k).

    NOTE: Generally the full resolution power array (what you get from FFTing the original waveform)
    will be much larger in size than the coherence array, which requires averaging over a number of 
    time steps. This difference is accounted for in this code.

    Inputs:
    powspec, tpowspec, fpowspec --> power spectrum (NOT complex) and assoc. time and freq values
    phasespec, tphasespec, fphasespec --> spectrum of phase values and assoc. time and freq values
    receiver_spacing --> separation of "centers of potential" from interferometry measurement

    nkbins --> number of bins in k-space
    klim --> range of k-values (rad/m)
    mean_max --> Since powspec is a larger array than phasespec (in general), I end up having to 
        reduce its dimensionality to size phasespec. Do this by taking the mean/max of relevant bins
    """


    import numpy as np

    #Define the range of k-values to populate
    kstep = (klim[1]-klim[0])/nkbins
    kvals = np.arange(klim[0],klim[1],kstep)


    #final values will be array of size [nfreqs,nkbins]
    nfreqs = np.shape(fphasespec)[0]
    powK = np.empty((nfreqs,nkbins)) #array of wave power sorted by freq and k-value 

     

    #time spacing
    tdelta = np.median(tphasespec - np.roll(tphasespec,1))

    #For each frequency
    for f in range(0,nfreqs-1,1):
        pslice = phasespec[f,:] #phase values for a particular freq and all times
        #kslice = [(3.1416/180)*i/receiver_spacing for i in pslice] #Change delta-phases into k-values (rad/m)
        kslice = [np.radians(i)/receiver_spacing for i in pslice] #Change delta-phases into k-values (rad/m)

        #continue if we have finite k values at current freq slice
        if np.sum(kslice != 0):
            #for current freq slice bin the k-values for all times
            for k in range(0,nkbins-1,1):

                #Extract the power values corresponding to finite k values.
                #--Note that the power array likely has more elements than phase array. So, we'll need to extract all power values
                #--in current frequency and k range
                gootimeIDX = np.where((kslice >= kvals[k]) & (kslice < kvals[k+1])) #time indices satisfying condition
                powgoo = np.empty(len(gootimeIDX[0])) #Can be multiple times for current freq that have k-values in current range (b/c we're considering all times)

                if len(gootimeIDX[0]) != 0:

                    #relevant freq range of (larger) power array to average over 
                    fgoo = [fpowspec[f],fpowspec[f+1]]
                    frange_powspecIDX = np.where((fpowspec >= fgoo[0]) & (fpowspec <= fgoo[1]))[0]
                    frange_powspecIDX = [np.min(frange_powspecIDX), np.max(frange_powspecIDX)]

                    #time index (t) defined from smaller phase array
                    for t in range(0,len(gootimeIDX[0]),1):
                        
                        #relevant time range and indices of (larger) power array to average over 
                        if t < len(gootimeIDX[0])-1:
                            tgoo = [tphasespec[t],tphasespec[t+1]]
                        else:
                            tgoo = [tphasespec[t],tphasespec[t]+tdelta] #populate last array element
                        trange_powspecIDX = np.where((tpowspec >= tgoo[0]) & (tpowspec <= tgoo[1]))[0]
                        trange_powspecIDX = [np.min(trange_powspecIDX), np.max(trange_powspecIDX)]

                        #For current time (defined from smaller phase array) can have 
                        #multiple freqs and times in larger array. Select these 
                        #and find their average/max value 
                        powarr_subset = powspec[frange_powspecIDX[0]:frange_powspecIDX[1],trange_powspecIDX[0]:trange_powspecIDX[1]]
                        if mean_max == 'max':
                            powgoo[t] = np.nanmax(powarr_subset)
                        else: 
                            powgoo[t] = np.mean(powarr_subset)


                    #now we need to sum over all times (defined from smaller phase array)
                    if np.sum(powgoo) > 1e-100:  #avoid absurdly low values
                        if mean_max == 'max':
                            powK[f,k] = np.nanmax(powgoo)            
                        else: 
                            powK[f,k] = np.mean(powgoo)            
                    

    return(powK, kvals, fphasespec)




