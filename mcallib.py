# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 2017
@author: Eder Martioli
Description: Library based on MCAL V.3 by Vasco Neves
Laboratorio Nacional de Astrofisica, Brazil
"""

import os
import numpy as np
from matplotlib import pylab
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

######################
def get_spectrafilepaths(directory):
    
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
    
    for root, directories, files in os.walk(directory):
        for filename in files:
            # Merge strings to form full filepaths
            filepath = os.path.join(root, filename)
            if filename.endswith(".fits") or filename.endswith(".fits.gz") or \
                filename.endswith(".spc.gz") or filename.endswith(".txt") or \
                filename.endswith("iu.s.gz") or filename.endswith("iu.s") or \
                filename.endswith("in.s.gz") or filename.endswith("in.s") or \
                filename.endswith("pu.s.gz") or filename.endswith("pu.s") or \
                filename.endswith("pn.s.gz") or filename.endswith("pn.s"):

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
    elif "iu.s.gz" in base or "pu.s.gz" in base or "in.s.gz" in base or "pn.s.gz" in base:
        basename = os.path.splitext(os.path.splitext(base)[0])[0]
    elif "iu.s" in base or "pu.s" in base or "in.s" in base or "pn.s" in base:
        basename = os.path.splitext(base)[0]
    else :
        print "Error: unknown extension in file: ", filepath
        exit()

    return basename, dir
###########################################

############# get Heliocentric Radial Velocity from header ####
def getEspadonsHelioRV(header) :
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
        if -samp < freq[j] < samp:
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
        
        if (maxlocf-2 < 0) or (endmaxf+3 > len(xx)) :
            print "WARNING 1: index out of range! Setting EQ=NaN"
            eqwidths.append('nan')
            continue
        
        maxloc_adjust = yy[maxlocf-2:maxlocf+3]
        endmax_adjust = yy[endmaxf-2:endmaxf+3]
        
        xmaxloc = xx[maxlocf-2:maxlocf+3]
        xendmax = xx[endmaxf-2:endmaxf+3]
        
        """
        if 3760 < i < 3790:
            
            for j in range(len(xmaxloc)) :
                print i, j, xmaxloc[j], maxloc_adjust[j]

            for j in range(len(xendmax)) :
                print i, j, xendmax[j], endmax_adjust[j]
        """
        locmaxmask = np.where(maxloc_adjust == max(maxloc_adjust))
        ind_temp_max = locmaxmask[0][0]
        maxlocfs = int(maxlocf + ind_temp_max - 2)

        locendmask = np.where(endmax_adjust == max(endmax_adjust))
        ind_temp_end = locendmask[0][0]
        endmaxfs = int(endmaxf + ind_temp_end - 2)
        
        if (maxlocfs < 0) or (maxlocfs > len(xx)) or  (endmaxfs < 0) or (endmaxfs > len(xx)) :
            print "WARNING 2: index out of range! Setting EQ=NaN"
            eqwidths.append('nan')
            continue

        yinit = yy[maxlocfs]
        yendt = yy[endmaxfs]
        
        yit = yy[maxlocfs:endmaxfs+1]  # nao normalizado
        xint = xx[maxlocfs:endmaxfs+1]

        tam = len(yit)

        if tam < 5 :
            print "WARNING 3: index out of range! Setting EQ=NaN"
            eqwidths.append('nan')
            continue

        if any(yinit < yit[0:5]) : minit = max(yit[0:5])
        else : minit = yinit
        
        if any(yendt < yit[-5:-1]) : mendt = max(yit[-5:-1])
        else : mendt = yendt

        yintt = yy[maxlocfs:endmaxfs+1]-max(minit,mendt)  #normalizado
        yinttnew = yy[maxlocfs:endmaxfs+1]

        indini = max(pylab.find (minit == yit))
        indend = min(pylab.find (mendt == yit))
        
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

    eqwidths = np.array(eqwidths)
    np.save(output,eqwidths)
    return eqwidths
######################

############# Function to calculate Teff and [Fe/H]  ###############
def mcal(eqw, calibmatrix) :
    
    """
        Calculate Effective Temperature (Teff) and Metallicity ([Fe/H])
        
        Parameters
        ----------
        eqw : array of peak to peak equivalent widths
        calibmatrix : calibration matrix (n x N), where n=number of
                      coefficients and N=number of lines/features used 
                      in the fit
        
        Returns
        -------
        feh_fit,teff_fit: tuple containing the calcualted [Fe/H] and Teff
        """
    
    ### loading the calibration matrix and errors
    file_cal = np.load(calibmatrix) ###calibration matrix
    coef = file_cal['coef'] #tellurics exclusion
    e_total = file_cal['e_total'] #tellurics exclusion
    
    #### Refit with weights
    ncoef = 4

    fun = lambda a,x,y,z : a[0] + a[1]*x + a[2]*y + a[3]*z
    err = lambda a,x,y,z,w,erro : (w - fun(a,x,y,z))*((1/erro**2)/(sum(1/erro**2)))

    xfit = coef[:,0] #alpha
    yfit = coef[:,1] #beta
    zfit = coef[:,2] #gamma

    a = np.array([0,0,0,3500]) #initial guess

    fit = optimize.leastsq (err,a,args=(xfit,yfit,zfit,eqw,e_total),full_output=1)

    rss = sum(fit[2]['fvec']**2)/(eqw.size - ncoef)
    efitfeh = np.sqrt(np.diag(rss*fit[1])[2])
    efitteff = np.sqrt(np.diag(rss*fit[1])[3])

    coef = fit[0]
    feh_fit = coef[2]
    teff_fit = coef[3]

    return (feh_fit,efitfeh),(teff_fit,efitteff)
######################

############# Function to calculate the error from a fit  ###############
def fit_err(fit, z, ncoef) :
    """
        Function to calculate the error from a fit
        
        Parameters
        ----------
        fit: fit object returned from optimize.leastsq()
        z: data used in the fit
        ncoef: number of coefficients used in the fit
        
        Returns error
        -------
        """
    rss = sum(fit[2]['fvec']**2)/(z[0].size - ncoef)
    efit_coef_sum = 0.0
    for j in range(ncoef) :
        efit_coef = np.sqrt(np.diag(rss*fit[1])[j])
        efit_coef_sum += efit_coef**2
    return np.sqrt(efit_coef_sum)
######################

############# Function to calculate the calibration matrix for mcal  ###############
def calibrate_mcal(ewcal, Tcal, FeHcal, outputcalibmatrix) :
    
    """
        Function to calculate the calibration matrix for MCAL
        
        Parameters
        ----------
        ewcal: array (N x 1) of peak to peak equivalent widths measured on the calibration sample (without NaNs)
        Tcal: array (Ns x 1)of Effective Temperature (Teff) values for the calibration sample
        FeHcal: array (Ns x 1) of Metalliticity ([Fe/H]) values for the calibration sample
        outputcalibmatrix: output calibration matrix (dimension n x N), where n=number of coefficients and N=number of lines/features used in the fit
        
        Returns
        -------
        """
    
    ewcal = np.array(ewcal)
    
    fun = lambda a,x,y : a[0] + a[1]*x + a[2]*y
    err = lambda a,x,y,z : (z - fun(a,x,y))

    a = np.array([0,0,0])
    
    xx = np.array(FeHcal)
    yy = np.array(Tcal)
    z = ewcal.transpose()
    
    coeffs, e_total = [],[]
    ncoef = 3
    
    if z[0].size < ncoef :
        print "Number of calibration stars (", z[0].size,") must be greater than the number of coefficients (",ncoef,")"
        exit()

    for i in range(len(z)) :
        fit = optimize.leastsq(err,a,args=(xx,yy,z[i]),full_output=1)
        e_total.append(fit_err(fit,z,ncoef))
        coeffs.append(fit[0])

    coeffs = np.array(coeffs)
    e_total = np.array(e_total)
    np.savez(outputcalibmatrix,coef=coeffs,e_total=e_total)
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
    
    hlc = max(pylab.find (x <= 6562.01))
    hrc = min(pylab.find (x >= 6563.61))
    hl1 = max(pylab.find (x <= 6545.495))
    hl2 = min(pylab.find (x >= 6556.245))
    hr1 = max(pylab.find (x <= 6575.934))
    hr2 = min(pylab.find (x >= 6584.684))
    
    ha_core = sum(y[hlc:hrc])
    href1 = sum(y[hl1:hl2])
    href2 = sum(y[hr1:hr2])
    
    ha_rtest = ha_core/(href1+href2)
    
    ha = ha_rtest*aha + bha

    return ha
######################

##### Function to load Source RV from Coolsnap/Archive datasheet files
def getSourceRadialVelocity(odonumber="",targetName="",coolsnapfile="clichesfroids_log.dat",SSArchivefile="clichesfroids_log_archive.dat",PolarArchivefile="clichesfroids_log_archive_p.dat",CSPolFupfile="clichesfroids_log_coolsnap_foll.dat",CSSpSFupfile="clichesfroids_log_coolsnap_foll_i.dat") :
    sourceRV = 0.0
    
    try:
        if(os.path.exists(coolsnapfile)) :
            f = open(coolsnapfile, 'r')
            lines = f.readlines()
            for i in range(len(lines)) :
                data = (lines[i].rstrip('\n')).split(',')
                if str(data[2]).replace(" ","") == odonumber :
                    sourceRV = float(data[19])
                    sourceSNR = float(data[11])
            f.close()
        
        if(os.path.exists(SSArchivefile)) :
            f = open(SSArchivefile, 'r')
            lines = f.readlines()
            for i in range(len(lines)) :
                data = (lines[i].rstrip('\n')).split(',')
                if str(data[1]).replace(" ","") == odonumber :
                    sourceRV = float(data[18])
                    sourceSNR = float(data[11])
            f.close()
        
        if(os.path.exists(PolarArchivefile)) :
            f = open(PolarArchivefile, 'r')
            lines = f.readlines()
            for i in range(len(lines)) :
                data = (lines[i].rstrip('\n')).split(',')
                if str(data[2]).replace(" ","") == odonumber :
                    sourceRV = float(data[20])
                    sourceSNR = float(data[11])
            f.close()
        
        if(os.path.exists(CSPolFupfile)) :
            f = open(CSPolFupfile, 'r')
            lines = f.readlines()
            for i in range(len(lines)) :
                data = (lines[i].rstrip('\n')).split(',')
                if str(data[2]).replace(" ","") == odonumber :
                    sourceRV = float(data[20])
                    #sourceSNR = float(data[11]) for ESPaDOnS
                    sourceSNR = float(data[12])
            f.close()
        
        if(os.path.exists(CSSpSFupfile)) :
            f = open(CSSpSFupfile, 'r')
            lines = f.readlines()
            for i in range(len(lines)) :
                data = (lines[i].rstrip('\n')).split(',')
                if str(data[1]).replace(" ","") == odonumber :
                    sourceRV = float(data[18])
                    sourceSNR = float(data[11])
            f.close()

    except :
        print 'Error: could not load sourceRV for odonumber=',odonumber,' in Coolsnap/Archive data file:',coolsnapfile, '/', SSArchivefile, '/', PolarArchivefile, '/', CSPolFupfile, '/', CSSpSFupfile
    
    return sourceRV,sourceSNR
######################

##### Function to load input file containing sample of calibrators
def loadSampleOfCalibrators(input) :
    try:
        twomassid, id, teff, feh = [],[],[],[]
        f = open(input, 'r')
        lines = f.readlines()
        for i in range(len(lines)) :
            data = (lines[i].rstrip('\n')).split(' ')
            twomassid.append(data[0])
            id.append(data[1])
            teff.append(float(data[2]))
            feh.append(float(data[3]))
        f.close()
        return twomassid, id, teff, feh
    except :
        print 'Error: could not load file of sample of calibrators=',input
        exit()
######################
