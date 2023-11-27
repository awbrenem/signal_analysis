"""
Routines for performing interferometry analysis 

NOTE: this technique only makes sense for parallel, spaced antennas. 
If the antennas are perpendicular, the ~90 deg phase shift in waves will be interpreted as a 
large k-value (k*d = delta-phase)


def inter_fvsk --> returns power spectrum of k-vectors vs freq


"""


def inter_fvsk(powspec,tpowspec,fpowspec,
               phasespec,tphasespec,fphasespec,
               receiver_spacing,
               nkbins=100,klim=[-2,2],
               mean_max='max'):
    

    """
    Interferometry analysis f vs k plots, with wave power on z-axis. 
    (see interferometry_routines_call.py for example usage)

    ------
    Reference:
    Lalti, A., Khotyaintsev, Y. V., & Graham, D. B. (2023). 
    Short-Wavelength Electrostatic Wave Measurement Using MMS Spacecraft
    JGR: Space Phys https://doi. org/10.1029/2022JA031150

    Graham, D. B., Y. V. Khotyaint-sev, A. Vaivads, and M. André (2016), 
    Electrostatic solitary wavesand electrostatic waves at the magnetopause,
    J. Geophys. Res.Space Physics,121, 3069–3092,doi:10.1002/2015JA021527.
    ------

    
    Idea: From an interferometry coherence analysis we can find the phase difference as a function of freq/time.
    This gives us the k-vector spectrum from k*d = delta-phase. 

    For a select time interval of interest, bin the k-values for each frequency to make 
    a freq vs k plot, with the z-values being wave power. 


    NOTE: Generally the full resolution power array (what you get from FFTing the original waveform)
    will be much larger in size than the coherence array, which requires averaging over a number of 
    time steps. This difference is accounted for in this code.

    NOTE: Can compare results to specific wave dispersion relations of f(k).
    Assuming a linear dispersion relation, can then calculate the phase speed with the slope.
        f = vph * (k/2*pi)

        
    Inputs:
    powspec, tpowspec, fpowspec --> power spectrum (NOT complex) and assoc. time and freq values
    phasespec, tphasespec, fphasespec --> spectrum of phase values (radians) and assoc. time (s) and freq (Hz) values
    receiver_spacing --> separation (m) of "centers of potential" from interferometry measurement

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
    powK = np.zeros((nfreqs,nkbins)) #array of wave power sorted by freq and k-value 

     

    #time spacing
    tdelta = np.median(tphasespec - np.roll(tphasespec,1))

    #------------------
    #SIMPLE VERSION ASSUMING THAT POWSPEC AND PHASESPEC HAVE SAME SIZE
    #---works (see interferometry_routines_call.py on how to make sure these are the same size)

    #For each frequency
    for f in range(0,nfreqs-1,1):
        pslice = phasespec[f,:] #phase values for a particular freq and all times
        #kslice = [(3.1416/180)*i/receiver_spacing for i in pslice] #Change delta-phases into k-values (rad/m)
        kslice = np.array([i/receiver_spacing for i in pslice]) #Change delta-phases into k-values (rad/m)

        #continue if we have finite k values at current freq slice
        if np.nansum(kslice) != 0:
            #for current freq slice bin the k-values for all times
            for k in range(0,nkbins-1,1):

                #Extract the power values corresponding to finite k values.

                gootimeIDX = np.where((kslice >= kvals[k]) & (kslice < kvals[k+1])) #time indices satisfying condition
                powgoo = np.zeros(len(gootimeIDX[0])) #Can be multiple times for current freq that have k-values in current range (b/c we're considering all times)


                if len(gootimeIDX[0]) != 0:

                    powgoo = powspec[f,gootimeIDX]
                
                    if mean_max == 'max':            
                        powK[f,k] = np.nanmean(powgoo)
                    else:
                        powK[f,k] = np.nanmax(powgoo)

                    
                    

                    """
                    #for each time that has a relevant k-value
                    for t in range(0,len(gootimeIDX[0]),1):
                        
                        ##relevant time range and indices of (larger) power array to average over 
                        #if t < len(gootimeIDX[0])-1:
                        #    tgoo = [tphasespec[t],tphasespec[t+1]]
                        #else:
                        #    tgoo = [tphasespec[t],tphasespec[t]+tdelta] #populate last array element
                        #trange_powspecIDX = np.where((tpowspec >= tgoo[0]) & (tpowspec <= tgoo[1]))[0]
                        #trange_powspecIDX = [np.min(trange_powspecIDX), np.max(trange_powspecIDX)]


                        powarr_subset = powspec[f,:]

                        #**********************************
                        #---ISSUE: EVERY TIME K INDEX IS UPDATED
                        #   POWARR_SUBSET GETS RESET TO 0.00010887 
                        #WHY IS THIS HAPPENING????
                        if mean_max == 'max':
                            if np.sum(powarr_subset) != 0:
                                #powgoo[t] = np.nanmax(powarr_subset)
                                powgoo[t] = np.median(powarr_subset)
                            else: 
                                powgoo[t] = float("nan")
                        else: 
                            powgoo[t] = np.mean(powarr_subset)
                        print('h')


                    #now we need to sum over all times (defined from smaller phase array)
                    if np.sum(powgoo) > 1e-100:  #avoid absurdly low values
                        if mean_max == 'max':
                            powK[f,k] = np.nanmax(powgoo)            
                        else: 
                            powK[f,k] = np.mean(powgoo)            

                    """


    """
    #For each frequency
    for f in range(0,nfreqs-1,1):
        pslice = phasespec[f,:] #phase values for a particular freq and all times
        #kslice = [(3.1416/180)*i/receiver_spacing for i in pslice] #Change delta-phases into k-values (rad/m)
        kslice = [i/receiver_spacing for i in pslice] #Change delta-phases into k-values (rad/m)

        #continue if we have finite k values at current freq slice
        if np.nansum(kslice) != 0:
            #for current freq slice bin the k-values for all times
            for k in range(0,nkbins-1,1):

                #Extract the power values corresponding to finite k values.


                #--Note that the power array can have more elements than phase array. 
                #--If so, we'll need to extract all power values in current frequency and k range.
                #--If not, then "else" loop is a simpler version of this "if" loop
                #*****I've tested this when the two have the same shape and it works fine. 

                gootimeIDX = np.where((kslice >= kvals[k]) & (kslice < kvals[k+1])) #time indices satisfying condition
                powgoo = np.zeros(len(gootimeIDX[0])) #Can be multiple times for current freq that have k-values in current range (b/c we're considering all times)

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
                        powarr_subset = powspec[frange_powspecIDX[0]:frange_powspecIDX[1],trange_powspecIDX[0]:trange_powspecIDX[1]][0]

                        #**********************************
                        #---ISSUE: EVERY TIME K INDEX IS UPDATED
                        #   POWARR_SUBSET GETS RESET TO 0.00010887 
                        #WHY IS THIS HAPPENING????
                        if mean_max == 'max':
                            if np.sum(powarr_subset) != 0:
                                #powgoo[t] = np.nanmax(powarr_subset)
                                powgoo[t] = np.median(powarr_subset)
                            else: 
                                powgoo[t] = float("nan")
                        else: 
                            powgoo[t] = np.mean(powarr_subset)
                        print('h')


                    #now we need to sum over all times (defined from smaller phase array)
                    if np.sum(powgoo) > 1e-100:  #avoid absurdly low values
                        if mean_max == 'max':
                            powK[f,k] = np.nanmax(powgoo)            
                        else: 
                            powK[f,k] = np.mean(powgoo)            
    """



    #Find k-location of max power for every frequency
    pmaxvals = np.zeros(np.shape(powK))
    pmaxvals[:,:] = 1
    for f in range(len(fphasespec)):
        if np.nanmax(powK[f,:]):
            cond = powK[f,:] == np.nanmax(powK[f,:])
        else:
            cond = 0
        pmaxvals[f,:] = pmaxvals[f,:]*cond


    return(powK, kvals, fphasespec, pmaxvals)




