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

def main():
	
	#sort out intrinsic flux values and calculate observed values
	#at different values of Av
	'''
	obs_data = np.loadtxt('observed_data')

	master_grid = obs_data[:,0]
	flux_obs = obs_data[:,1]
	error_obs = obs_data[:,2]
	'''
	
	obs_data = np.loadtxt('vandels_high-mass_data')
	
	master_grid = obs_data[:,0]
	flux_obs = obs_data[:,1]
	error_obs = obs_data[:,2]
	
	Z001_data = np.loadtxt('vandels_high-mass_data')
	
	masked_wav = Z001_data[:,0]
	Z001_model = Z001_data[:,1]
	
	Z002_data = np.loadtxt('vandels_high-mass_data')

	Z002_model = Z002_data[:,1]
	
	Z008_data = np.loadtxt('vandels_high-mass_data')

	Z008_model = Z008_data[:,1]

	Z014_data = np.loadtxt('vandels_high-mass_data')

	Z014_model = Z014_data[:,1]
	
	Z040_data = np.loadtxt('vandels_high-mass_data')

	Z040_model = Z040_data[:,1]
	
	master_grid_bump = np.arange(1000, 5500, 1)
	A_cullen, A_calzetti = af.attenuation(master_grid_bump/1.e4)
	A_mod = af.attenuation_mod(master_grid_bump/1.e4, d, B)
	
	plt.plot(master_grid_bump, A_calzetti, color='black', lw=1.5)
	plt.plot(master_grid_bump, A_mod, color='red', lw=1.5)
	plt.tick_params(direction='in', which='both', length=6, width=0.5, grid_alpha=0.5, grid_linewidth=0.5)
	plt.title("Absolute Attentuation vs Wavelength")
	plt.xlabel("$\lambda$ [$\mathrm{\AA}$]", size=16)
	plt.ylabel("$A_\lambda / A_v$", size=16)
	plt.xlim(1000, 5500)
	plt.show()

	plt.plot(master_grid, flux_obs, color='black', lw=1.)
	for i in range(len(unmasked_wav)):
		plt.plot(unmasked_wav[i],f_obs_unmask[i], color='deepskyblue', lw=1.5)
	plt.plot(masked_wav, Z001_model, color='limegreen', lw=1.5)
	plt.plot(masked_wav, Z002_model, color='limegreen', lw=1.5)
	plt.plot(masked_wav, Z008_model, color='limegreen', lw=1.5)
	plt.plot(masked_wav, yM, color='limegreen', lw=1.5)
	plt.plot(masked_wav, yM, color='limegreen', lw=1.5)
	plt.plot(master_grid, error_obs, color='red', lw=1.)
	plt.tick_params(direction='in', which='both', length=6, width=0.5, grid_alpha=0.5, grid_linewidth=0.5)
	plt.xlabel("$\lambda$ [$\mathrm{\AA}$]", size=16)
	plt.ylabel("f$_\lambda$ [$erg$ $s^{-1}cm^{-2}\mathrm{\AA}^{-1}$]", size=16)
	plt.xlim(1200, 2000)
	plt.show()

if __name__ == "__main__": 
    main()

