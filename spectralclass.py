# -*- coding: utf-8 -*-

"""
Spectral Classes
---------------------------
Created on January 31 2017
@author: Eder Martioli
Laboratorio Nacional de Astrofisica, Brazil
"""
import os
import numpy as np
from astropy.io import fits
import mcallib
from scipy import constants

########## SPECTRUM CLASS ############
class Spectrum :
    'Common base class for a spectrum'
    
    def __init__(self, Filename, FFT=False):

        """
        Create a Spectrum object.
        
        Parameters
        ----------
        filename : string
            File to read the spectrum from.

        Examples
        --------
        >>> spc = Spectrum("spectrumfile.spc.gz")
        >>> spc = Spectrum("spectrumfile.fits")
        >>> spc = Spectrum("spectrumfile.txt")
        """
        self.heliorv = 0.0
        
        self.filepath = Filename
        
        basename = mcallib.getbasename(self.filepath)
        self.eqw_output = basename[1] + '/' + basename[0] + '_ew.npy'
        
        self.id = basename[0]
        self.filename = os.path.basename(self.filepath)

        try :
            if self.filepath.endswith(".fits") :
                self.wl,self.flux=self.loadSpectrumFromFITS(self.filepath)
            elif self.filepath.endswith(".spc.gz") :
                self.wl,self.flux=self.loadSpectrumFromSPC(self.filepath)
            elif self.filepath.endswith(".txt") :
                self.wl,self.flux=self.loadSpectrumFromTXT(self.filepath)
            else :
                print "Error: file type not supported for input spectrum: ",self.filepath
                exit()
        except :
            print "Error: could not open file: ",self.filepath
            exit()

        if FFT :
            self.flux = mcallib.fft_filter(self.flux)

    #--- Function to load spectrum from .spc.gz (OPERA) file
    def loadSpectrumFromSPC(self,Filename):
        wl,rvcorr,xcorr,rawflux,rawvar,normflux,calflux = np.loadtxt(Filename, unpack=True, comments='#', skiprows=11, usecols=(4,6,7,8,9,10,11), delimiter=' ')

        threshold = 10.0
        snr = rawflux / np.sqrt(rawvar)
        mask = np.where(snr > threshold)
                       
        indices = (wl[mask]).argsort()
        
        sortedwl = (wl[mask])[indices]
        #sortedrvcorrwl = (wl[mask] + rvcorr[mask])[indices]
        #sortednormflux = (normflux[mask])[indices]
        sortedcalflux = (calflux[mask])[indices]
        cleanwl = np.linspace(np.min(sortedwl), np.max(sortedwl), len(sortedwl))
        #cleannormflux = np.interp(cleanwl, sortedwl, sortednormflux)
        #cleanrvcorrnormflux = np.interp(cleanwl, sortedrvcorrwl, sortednormflux)
        cleancalflux = np.interp(cleanwl, sortedwl, sortedcalflux)
    
        cleanwl = cleanwl * 10.0
        return cleanwl, cleancalflux
    #------------

    #--- Function to load spectrum from .fits file
    def loadSpectrumFromFITS(self,Filename):
        
        wl,flux = [],[]
        hdu = fits.open(Filename)
        
        try :
            if hdu[0].header['INSTRUME'] == 'ESPaDOnS' :

                self.heliorv = mcallib.getEspaonsHelioRV(hdu[0].header)
                
                self.instrument = 'ESPaDOnS'
                self.object = hdu[0].header['OBJECT']
                
                sourceRV = mcallib.getSourceRadialVelocity(self.object)
                
                if hdu[0].header['INSTMODE'] == 'Polarimetry, R=65,000' :
                    # data[0],data[1] for normalized spectrum
                    # data[6],data[7] for unnormalized spectrum
                    wltmp = hdu[0].data[0]*(1.0 - sourceRV*1000.0/constants.c)
                    indices = wltmp.argsort()
                    wl = 10.0*wltmp[indices]
                    flux = (hdu[0].data[1])[indices]
                
                elif hdu[0].header['INSTMODE'] == 'Spectroscopy, star+sky, R=65,000':
                    # data[0],data[1] for normalized spectrum
                    # data[7],data[8] for unnormalized spectrum
                    wltmp = hdu[0].data[0]*(1.0 - sourceRV*1000.0/constants.c)
                    indices = wltmp.argsort()
                    wl = 10.0*wltmp[indices]
                    flux = (hdu[0].data[1])[indices]

        except :
            self.object = (self.id).rsplit('_')[0]
            self.instrument = 'HARPS'
            length = hdu[0].header['NAXIS1']
            start = hdu[0].header['CRVAL1']
            step = hdu[0].header['CDELT1']
            wl = np.linspace(start, start+length*step, length)
            flux = hdu[0].data

        return wl, flux
    #------------

    #--- Function to load spectrum from .txt file
    def loadSpectrumFromTXT(self,Filename):
        x,y = np.loadtxt(Filename, unpack=True, comments='#',usecols=(0,1),skiprows=0, delimiter=' ')
        return x,y
    #------------

    #--- Function to load/calculate Equivalent Widths
    def equivalentWidths(self, inputlinelist='lines.rdb', output=True, override=False, verbose=False) :
        try :
            if os.path.exists(self.eqw_output) and override == False :
                if verbose: print "Loading EWs from existing file: ",self.eqw_output
                self.eqwidths = np.load(self.eqw_output)
            else :
                if verbose: print "Calculating EWs from spectrum using list of lines: ",inputlinelist
                self.eqwidths = mcallib.measureEquivalentWidths(self.wl, self.flux, inputlinelist, self.eqw_output)
        except:
            print "Error: could not calculate Eq Widths. Input line list:",inputlinelist
            exit()
    #------------

    #--- Function to calculate Teff and [Fe/H]
    def calculateTeffAndFeH(self, calibmatrix, verbose=False) :
        try :
            if verbose: print "Calculating Teff and [Fe/H] using calibration: ", calibmatrix
            (self.FeH,self.eFeH),(self.Teff,self.eTeff) = mcallib.mcal(self.eqwidths, calibmatrix)
        except:
            print "Error: could not calculate Teff and [Fe/H]. Calibration matrix:",calibmatrix
            exit()
    #------------

    #--- Print spectrum data
    def resampling(self, wlsampling, wl0, wlf) :
        npoints = int((wlf-wl0)/wlsampling)
        wl_new = np.linspace(wl0, wlf, npoints)
        flux_new = np.interp(wl_new, self.wl, self.flux)
        self.wl = wl_new
        self.flux = flux_new
    #------------

    #--- Print spectrum information
    def info(self) :
        print "**************************"
        print "Info for spectrum: ",self.filename, " Object:", self.object
        print "Instrument:",self.instrument
        if self.instrument == 'ESPaDOnS' :
            print "RV=",self.heliorv,"km/s"
        print "wl0 =",self.wl[0],"A -- wlf =",self.wl[-1],"A"
        sampling = (self.wl[-1] - self.wl[0])/float(len(self.wl))
        print "sampling =",sampling," A/pixel"
        print "<F> =",self.flux.mean(),"+-",self.flux.std()
    
        eqwbad = np.count_nonzero(np.isnan(self.eqwidths))
        eqwgood = np.count_nonzero(~np.isnan(self.eqwidths))
        print "Using", np.count_nonzero(~np.isnan(self.eqwidths)), "of",len(self.eqwidths),"lines in list for EWs"
        print "**************************\n"
    #------------

    #--- Print spectrum data
    def printdata(self) :
        for i in range(len(self.wl)) :
            print self.wl[i],self.flux[i]
    #------------

    #--- Check activity
    def halphaActivity(self, verbose=False) :
        if verbose: print 'Calculating Halpha activity index...'

        self.halpha = mcallib.calculateHalphaActivity(self.wl,self.flux)
        
        if verbose :
            if self.halpha >= 0.25 :
                print 'WARNING: the star may be too active to use this calibration. Halpha =',self.halpha,'. Check the webpage or Neves (2014) for more details.'
            else :
                print 'No significant Halpha emission'
    #------------



