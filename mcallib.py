# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 2017
@author: Eder Martioli
Description: Library based on MCAL V.3 by Vasco Neves
Laboratorio Nacional de Astrofisica, Brazil
"""

import os
import numpy as np
from pylab import *
from scipy import optimize
from scipy import integrate

######################
def get_fitsfilepaths(directory):
    
    """
    Generates a list of file names in a directory tree
    by walking the tree either top-down or bottom-up.
        
    Parameters
    ----------
    directory : directory path
        
    Returns
    -------
    file_paths: a list of file paths
    """
        
    file_paths = []  # List which will store all of the full filepaths.

    # Walk the tree:
    for root, directories, files in os.walk(directory):
        for filename in files:
			# Merge strings to form full filepaths
            filepath = os.path.join(root, filename)
            if filename.endswith(".fits") :
                file_paths.append(filepath)

    file_paths.sort()
    return file_paths
######################

############# get basename from FITS file path ###############
def getbasename(filepath) :
    
    basename = ''
    base = os.path.basename(filepath)
    
    dir = os.path.dirname(filepath)
    
    if ".fits.gz" in base :
        basename = os.path.splitext(os.path.splitext(base)[0])[0]
    elif ".fits" in base :
        basename = os.path.splitext(base)[0]
    elif ".txt" in base :
        basename = os.path.splitext(base)[0]
    elif ".spc.gz" in base :
        basename = os.path.splitext(os.path.splitext(base)[0])[0]
    else :
        print "Error: unknown extension in file: ", filepath
        exit()

    return basename, dir
###########################################

############# get Heliocentric Radial Velocity from header ####
def getEspaonsHelioRV(header) :
    comments = header['COMMENT']
    rv = 0.0
    for c in comments :
        if "Heliocentric velocity" in c :
            rv = float(c[c.find(": ")+1:c.find(" km")])
            break
    return rv
###########################################


############# FFT filtering ###############
def fft_filter(y, step=0.01, samp=20):
    """
        Perform Fast Fourier Transform filtering of data
        
        Parameters
        ----------
        y : data array
        step : step to define frequency sampling
        samp : number of sample points

        Returns
        -------
        yfftclean: filtered data array
    """

    fourier = np.fft.fft(y)
    n = len(fourier)
    freq = np.fft.fftfreq(n, d=step)
    iy = np.zeros(n, dtype=complex)
        
    for j in xrange(n):
        if -samp < freq[j] and freq[j] < samp:
            iy[j] = fourier[j]
        else:
            iy[j] = 0
    yfftclean = np.real(np.fft.ifft(iy, n))
    
    return yfftclean
######################

############# Measure Pseudo Equivalent Widths ###############
def measureEquivalentWidths(xx, yy, inputlinelist='lines.rdb', output = 'ew_out.npz') :
    
    """
        Measure pseudo equivalent widths in spectrum
        
        Parameters
        ----------
        xx : wavelength data array (angstrom)
        yy : flux data array (normalized)
        inputlinelist : input file containing list of line ranges [x0,xf]
        output : output file name
        
        Returns
        -------
        eqwidths: array of measured equivalent widths
    """
    
    xinif, xendf = np.loadtxt(inputlinelist, unpack=True, delimiter='\t',usecols=(0,1),dtype='float',skiprows=2)
    
    eqwidths = []
    
    for i in range(len(xinif)) :

        mask = np.where(np.logical_and(xx>=xinif[i], xx<=xendf[i]))
        
        maxlocf = mask[0][0]
        endmaxf = mask[0][-1] + 1
        
        maxloc_adjust = yy[maxlocf-2:maxlocf+3]
        endmax_adjust = yy[endmaxf-2:endmaxf+3]
        
        ind_temp = max(find(maxloc_adjust == max(maxloc_adjust)))
        maxlocfs = int(maxlocf + ind_temp - 2)
        ind_temp = min(find(endmax_adjust == max(endmax_adjust)))
        ind_endmax = int(endmaxf+ind_temp-2)
        
        endmaxfs = ind_endmax
        yinit = yy[maxlocfs]
        yendt = yy[endmaxfs]
        
        yit = yy[maxlocfs:endmaxfs+1]  # nao normalizado
        xint = xx[maxlocfs:endmaxfs+1]

        tam = yit.size

        if any(yinit < yit[0:5]) : minit = max(yit[0:5])
        else : minit = yinit
        
        if any(yendt < yit[-5:-1]) : mendt = max(yit[-5:-1])
        else : mendt = yendt

        yintt = yy[maxlocfs:endmaxfs+1]-max(minit,mendt)  #normalizado
        yinttnew = yy[maxlocfs:endmaxfs+1]

        indini = max(find (minit == yit))
        indend = min(find (mendt == yit))
        
        if (xint[indend]-xint[indini]) != 0.0:
            mt = (mendt-minit)/(xint[indend]-xint[indini])
            
            yeqt = minit + mt*(xint-xint[indini]) - max(minit,mendt) # condicao fronteira para o integral normalizada
            yeqplott = minit + mt*(xint-xint[indini])
        
            #### if yeqt < yintt
            eqw = 1000*np.trapz(1-yinttnew/yeqplott,dx=0.01)
        
            eqwidths.append(eqw)
        else :
            eqw = float('nan')
            eqwidths.append(eqw)

        #print i,"/",len(xinif), eqw

    
    np.save(output,eqwidths)
    return eqwidths
######################

############# Function to calculate Teff and [Fe/H]  ###############
def mcal(int_ew, calibmatrix) :
    
    """
        Calculate Effective Temperature (Teff) and Metallicity ([Fe/H])
        
        Parameters
        ----------
        int_ew : array of peak to peak equivalent widths
        calibmatrix : calibration matrix (n x N), where n=number of
                      coefficients and N=number of lines/features used 
                      in the fit
        
        Returns
        -------
        feh_fit,teff_fit: tuple containing the calcualted [Fe/H] and Teff
        """

    # removing nans from the EW file
    int_ew = nan_to_num(int_ew)
    
    ###loading the calibration matrix and errors
    file_cal = np.load(calibmatrix) ###calibration matrix
    coef = file_cal['coef'] #tellurics exclusion
    e_total = file_cal['e_total'] #tellurics exclusion
    
    ####2nd part of the calibration: refit with weights
    z2 = int_ew
    
    fun2 = lambda a,xx2,yy2,zz2 : a[0] + a[1]*xx2 + a[2]*yy2 + a[3]*zz2
    err2 = lambda a,xx2,yy2,zz2,z2,erro : (z2 - fun2(a,xx2,yy2,zz2))*((1/erro**2)/(sum(1/erro**2)))
    
    xx2fit = coef[:,0]
    yy2fit = coef[:,1]
    zz2fit = coef[:,2]
    
    a = array([0,0,0,3500])
    
    fit2 = optimize.leastsq (err2,a,args=(xx2fit,yy2fit,zz2fit,z2,e_total),full_output=1)
            
    rss = sum(fit2[2]['fvec']**2)/(z2.size-4)
    efitfeh = sqrt(diag(rss*fit2[1])[2])
    efitteff = sqrt(diag(rss*fit2[1])[3])
    coef2 = fit2[0]

    feh_fit = coef2[2]
    teff_fit = coef2[3]
    
    return (feh_fit,efitfeh),(teff_fit,efitteff)
######################

##### CALCULATE corrected Halpha to obtain a measure of activity ############
def calculateHalphaActivity(xx,yy) :

    aha =  1.8565921 #correction coefficients
    bha = -0.0870908 #correction coefficients
    
    mask = np.where(np.logical_and(xx>=6545., xx<=6585.))
    
    if len(mask[0]) == 0 :
        print "Error: spectrum does not cover H-alpha range! "
        exit()
    
    x = xx[mask]
    y = yy[mask]
    
    hlc = max(find (x <= 6562.01))
    hrc = min(find (x >= 6563.61))
    hl1 = max(find (x <= 6545.495))
    hl2 = min(find (x >= 6556.245))
    hr1 = max(find (x <= 6575.934))
    hr2 = min(find (x >= 6584.684))
    
    ha_core = sum(y[hlc:hrc])
    href1 = sum(y[hl1:hl2])
    href2 = sum(y[hr1:hr2])
    
    ha_rtest = ha_core/(href1+href2)
    
    ha = ha_rtest*aha + bha

    return ha
######################

##### Function to load Source RV from Coolsnap datasheet file
def getSourceRadialVelocity(odonumber="",targetName="",coolsnapfile="clichesfroids_log.dat") :
    sourceRV = 0.0
    try:
        if(os.path.exists(coolsnapfile)) :
            f = open(coolsnapfile, 'r')
            lines = f.readlines()
            for i in range(len(lines)) :
                data = (lines[i].rstrip('\n')).split(',')
                if str(data[2]).replace(" ","") == odonumber or str(data[3]).replace(" ","")==targetName:
                    sourceRV = float(data[19])
            f.close()
            #print odonumber, targetName, 'sourceRV:',sourceRV
    except :
        print 'Error: could not load sourceRV from Coolsnap data file:',coolsnapfile, ' odonumber=',odonumber, 'targetName=',targetName

    return sourceRV
######################


