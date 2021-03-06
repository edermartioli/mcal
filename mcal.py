"""
    Created on Jan 31 2017
    
    Description: A wrapper to run MCAL V.3 by Vasco Neves
    
    @author: Eder Martioli <emartioli@lna.br>
    
    Laboratorio Nacional de Astrofisica, Brazil.
    
    Simple usage example:
    
    ./mcal.py --datadir=/Users/edermartioli/Desktop/CoolSnap/vasco-stuff/data/ --calibmatrix=./coef_calvb12.npz
    """

__version__ = "1.0"

__copyright__ = """
    Copyright (c) ...  All rights reserved.
    """

from optparse import OptionParser
import os,sys
import mcallib
from spectralclass import Spectrum

import time
start_time = time.time()

parser = OptionParser()
parser.add_option("-i", "--datadir", dest="datadir", help="Data directory",type='string',default="")
parser.add_option("-c", "--calibmatrix", dest="calibmatrix", help="Calibration matrix",type='string',default="coef_calvb12.npz")
parser.add_option("-f", action="store_true", dest="fft", help="activate FFT filtering",default=False)
parser.add_option("-v", action="store_true", dest="verbose", help="verbose",default=False)

try:
    options,args = parser.parse_args(sys.argv[1:])
except:
    print "Error: check usage with mcal.py -h ";sys.exit(1);

if options.verbose:
    print 'Data directory: ', options.datadir
    print 'Calibration matrix: ', options.calibmatrix

listOfStarSpectra = mcallib.get_spectrafilepaths(options.datadir)

print "--------------------------------------------------"
print "File","Object","[Fe/H]_corr","[Fe/H]","e[Fe/H]","Teff_corr","Teff","eTeff","H-alpha", "RV(km/s)"
print "--------------------------------------------------"

for file in listOfStarSpectra :

    if options.verbose :
        print "Processing spectrum: ", file

    spc = Spectrum(file, FFT=options.fft)

    if spc.instrument == 'ESPaDOnS' :
        spc.resampling(0.01,4000,10400)
    elif spc.instrument == 'NARVAL' :
        spc.resampling(0.01,4000,10400)

    spc.equivalentWidths(override=True,verbose=options.verbose)

    spc.calculateTeffAndFeH(options.calibmatrix, verbose=options.verbose)

    spc.halphaActivity(verbose=options.verbose)

    #spc.info()
    #spc.printdata()

    Teff_corr, Fe_H_corr = spc.TeffAndFeHCorr()

    print spc.filename, spc.object, round(Fe_H_corr,3), round(spc.FeH,3), round(spc.eFeH,3), int(round(Teff_corr,0)), int(round(spc.Teff,0)), int(round(spc.eTeff,0)), round(spc.halpha,2), round(spc.sourceRV,2)

print("--- Total time: %s seconds ---" % (time.time() - start_time))
