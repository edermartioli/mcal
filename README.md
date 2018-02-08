# mcal
MCAL is a software to perform measurements of the effective temperature (Teff) and metallicity ([Fe/H]) of cool M-dwarf stars. This version of MCAL is a wrapper to run an adapted version of the original MCAL V.3 by Vasco Neves.

Clone and run:

python mcal.py --datadir=./data/ --calibmatrix=coef_calvb12.npz

To recalibrate MCAL:

python calibrate-mcal.py --inputsample=calibrators_db.dat --datadir=./calibsample/ --outputcalibmatrix=coefcal.npz
