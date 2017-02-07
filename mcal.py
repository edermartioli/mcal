#!/opt/anaconda/bin/python
# -*- coding: iso-8859-1 -*-
"""
    Shebang options:
            #!/usr/bin/python
            #!/opt/anaconda/bin/python
            #!/Users/edermartioli/Local/Ureka/variants/common/bin/python
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

parser = OptionParser()
parser.add_option("-i", "--datadir", dest="datadir", help="Data directory",type='string',default="")
parser.add_option("-c", "--calibmatrix", dest="calibmatrix", help="Calibration matrix",type='string',default="")
parser.add_option("-f", action="store_true", dest="fft", help="activate FFT filtering",default=False)
parser.add_option("-v", action="store_true", dest="verbose", help="verbose",default=False)

try:
    options,args = parser.parse_args(sys.argv[1:])
except:
    print "Error: check usage with mcal.py -h ";sys.exit(1);

if options.verbose:
    print 'Data directory: ', options.datadir
    print 'Calibration matrix: ', options.calibmatrix

listOfStarSpectra = mcallib.get_fitsfilepaths(options.datadir)

print "---------------------------------------------"
print "File","Object","[Fe/H]","e[Fe/H]","Teff","eTeff","H-alpha"
print "---------------------------------------------"

for file in listOfStarSpectra :

    spc = Spectrum(file, FFT=options.fft)
    
    if spc.instrument == 'ESPaDOnS' :
        spc.resampling(0.01,4000,7000)

    #spc.printdata()

    spc.equivalentWidths(verbose=options.verbose)
    
    spc.calculateTeffAndFeH(options.calibmatrix, verbose=options.verbose)

    spc.halphaActivity(verbose=options.verbose)

    #spc.info()
    #spc.printdata()

    print spc.filename, spc.object, round(spc.FeH,3), round(spc.eFeH,3), int(round(spc.Teff,0)), int(round(spc.eTeff,0)), round(spc.halpha,2)
