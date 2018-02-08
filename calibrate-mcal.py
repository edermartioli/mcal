"""
    Created on Aug 25 2017
    
    Description: Tool to create calibration matrix from calibrators for MCAL to calculate Teff and [Fe/H]
    
    @author: Eder Martioli <emartioli@lna.br>
    
    Laboratorio Nacional de Astrofisica, Brazil.
    
    Simple usage example:
    
    python calibrate-mcal.py --inputsample=calibrators_db.dat --datadir=./data/ --outputcalibmatrix=./coefcal.npz
    """

__version__ = "1.0"

__copyright__ = """
    Copyright (c) ...  All rights reserved.
    """

from optparse import OptionParser
import os,sys
import mcallib
from spectralclass import Spectrum
import numpy as np
from scipy import optimize

#import matplotlib.pyplot as plt

import time
start_time = time.time()

parser = OptionParser()
parser.add_option("-i", "--inputsample", dest="inputsample", help="Input sample file",type='string',default="")
parser.add_option("-d", "--datadir", dest="datadir", help="Data directory",type='string',default="")
parser.add_option("-c", "--outputcalibmatrix", dest="outputcalibmatrix", help="Output calibration matrix",type='string',default="coef_calvb12.npz")
parser.add_option("-v", action="store_true", dest="verbose", help="verbose",default=False)

try:
    options,args = parser.parse_args(sys.argv[1:])
except:
    print "Error: check usage with mcal.py -h ";sys.exit(1);

if options.verbose:
    print 'Input sample file: ', options.inputsample
    print 'Data directory: ', options.datadir
    print 'Output calibration matrix: ', options.outputcalibmatrix

listOfStarSpectra = mcallib.get_spectrafilepaths(options.datadir)

twomassid, id, teff, feh = mcallib.loadSampleOfCalibrators(options.inputsample)

Tcal = []
FeHcal = []
ewcal = []

for file in listOfStarSpectra :

    spc = Spectrum(file)
    objname = spc.object.replace(" ","").replace("_","")

    if options.verbose :
        print "Processing spectrum: ", file, objname
    
    if spc.instrument == 'ESPaDOnS' :
        spc.resampling(0.01,4000,10400)

    for i in range(len(id)) :
        if objname in id[i] or objname in twomassid[i]:
            spc.equivalentWidths(override=False, verbose=options.verbose)
            ewcal_tmp = spc.eqwidths
            if not np.any(np.isnan(ewcal_tmp)) :
                ewcal.append(ewcal_tmp)
                Tcal.append(teff[i])
                FeHcal.append(feh[i])
                print "Selected Calibrator: ", file, objname, twomassid[i], id[i], teff[i], feh[i]
            else :
                print "Eq. Widths contain NaNs. Ignoring potential calibrator: ", file, objname

mcallib.calibrate_mcal(ewcal, Tcal, FeHcal, options.outputcalibmatrix)

print("--- Total time: %s seconds ---" % (time.time() - start_time))
