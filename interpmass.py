
from astropy.io import fits
from astropy.stats import sigma_clip
import matplotlib.patches as mpatches
import numpy as np
import matplotlib.pyplot as plt
import os, sys
import glob 
from scipy.interpolate import interp1d
from astropy.stats import bootstrap
from astropy.utils import NumpyRNGContext
from scipy import stats
import attenuationfunc as af
import chisq as cs
import stackspectra as ss
import emcee
from scipy import integrate
import pandas as pd
import scipy.optimize as op
from dynesty import plotting as dyplot
import dynesty
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

log_flux_Z001, wavelength = cs.loadspectrum('/home/s1535092/summerproject_updated/S99_models/S99-v00-Z001-IMF2.3_100myr.spectrum')

flux_Z001 = 10**(log_flux_Z001)
	
log_flux_Z002, wavelength = cs.loadspectrum('/home/s1535092/summerproject_updated/S99_models/S99-v00-Z002-IMF2.3_100myr.spectrum')

flux_Z002 = 10**(log_flux_Z002)

log_flux_Z008, wavelength = cs.loadspectrum('/home/s1535092/summerproject_updated/S99_models/S99-v00-Z008-IMF2.3_100myr.spectrum')

flux_Z008 = 10**(log_flux_Z008)

log_flux_Z014, wavelength = cs.loadspectrum('/home/s1535092/summerproject_updated/S99_models/S99-v00-Z014-IMF2.3_100myr.spectrum')

flux_Z014 = 10**(log_flux_Z014)

log_flux_Z040, wavelength = cs.loadspectrum('/home/s1535092/summerproject_updated/S99_models/S99-v00-Z040-IMF2.3_100myr.spectrum')

flux_Z040 = 10**(log_flux_Z040)


plt.plot(wavelength, flux_Z001, color='black', lw=1.5)
plt.plot(wavelength, flux_Z002, color='red', lw=1.5)
plt.plot(wavelength, flux_Z008, color='blue', lw=1.5)
plt.plot(wavelength, flux_Z014, color='yellow', lw=1.5)
plt.plot(wavelength, flux_Z040, color='pink', lw=1.5)
plt.tick_params(direction='in', which='both', length=6, width=0.5, grid_alpha=0.5, grid_linewidth=0.5)
plt.xlabel("$\lambda$ [$\mathrm{\AA}$]", size=16)
plt.ylabel("f$_\lambda$ [$erg$ $s^{-1}cm^{-2}\mathrm{\AA}^{-1}$]", size=16)
plt.show()
