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
	
	obs_data = np.loadtxt('observed_data')

	master_grid = obs_data[:,0]
	flux_obs = obs_data[:,1]
	error_obs = obs_data[:,2]
	
	log_flux_int, wavelength = cs.loadspectrum('/home/s1535092/summerproject_updated/S99_models/S99-v00-Z002-IMF2.3_100myr.spectrum')

	flux_int = 10**(log_flux_int)

	masked_wav, unmasked_wav, f_int_mask, f_obs_mask, f_obs_unmask, sigma_mask = cs.maskwavelength(wavelength, flux_int, flux_obs, error_obs)
	
	#define free parameters 
	A_v = np.array(np.linspace(0, 6.0, 31))
	d = np.array(np.linspace(-0.5, 0.5, 11))
	B = np.array(np.linspace(0, 1.0, 31))

	theta = [A_v, d, B]
	flux_mod, A_v_list , B_list , d_list = cs.fluxmod(f_int_mask, masked_wav/1.e4)

	beta = cs.normalisingfactor(flux_mod, f_obs_mask, sigma_mask)

	chi_squared = cs.chisquared(flux_mod, f_obs_mask, sigma_mask, beta)

	y = []
	for i in range(len(beta)):
		_b = beta[i]
		_fm = np.array(flux_mod[i])
		y_data = _b * _fm
		y.append(y_data)

	flux_best_fit = y[np.argmin(chi_squared)]
		
	print(''+str(min(chi_squared))+','+str(np.array(A_v_list[np.argmin(chi_squared)]))+','+str(np.array(B_list[np.argmin(chi_squared)]))+','+str(np.array(d_list[np.argmin(chi_squared)]))+','+str(cs.chi2_likelihood(flux_mod, f_obs_mask, sigma_mask, beta)))

	plt.plot(master_grid, flux_obs, color='black', lw=1.)
	for i in range(len(unmasked_wav)):
		plt.plot(unmasked_wav[i],f_obs_unmask[i], color='deepskyblue', lw=1.5)
	plt.plot(masked_wav, flux_best_fit, color='limegreen', lw=1.5)
	plt.plot(master_grid, error_obs, color='red', lw=1.)
	plt.tick_params(direction='in', which='both', length=6, width=0.5, grid_alpha=0.5, grid_linewidth=0.5)
	plt.xlabel("$\lambda$ [$\mathrm{\AA}$]", size=16)
	plt.ylabel("f$_\lambda$ [$erg$ $s^{-1}cm^{-2}\mathrm{\AA}^{-1}$]", size=16)
	plt.xlim(1200, 2000)
	plt.show()

if __name__ == "__main__": 
    main()


