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
#from dynesty import plotting as dyplot
#import dynesty
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

def main():
	
	#sort out intrinsic flux values and calculate observed values
	#at different values of Av
	'''
	obs_data1 = np.loadtxt('observed_data')

	master_grid1 = obs_data1[:,0]
	flux_obs1 = obs_data1[:,1]
	error_obs1 = obs_data1[:,2]
	'''

	obs_data = np.loadtxt('vandelsz45_low-mass_data')

	master_grid = obs_data[:,0]
	flux_obs = obs_data[:,1]
	error_obs = obs_data[:,2]

	
	Z001_data = cs.getdataz45('/Users/carolynmichael/Documents/summerproject_updated 2/S99_models/S99-v00-Z001-IMF2.3_100myr.spectrum', obs_data)
	
	Z002_data = cs.getdataz45('/Users/carolynmichael/Documents/summerproject_updated 2/S99_models/S99-v00-Z002-IMF2.3_100myr.spectrum', obs_data)
	
	Z008_data = cs.getdataz45('/Users/carolynmichael/Documents/summerproject_updated 2/S99_models/S99-v00-Z008-IMF2.3_100myr.spectrum', obs_data)
	
	Z014_data = cs.getdataz45('/Users/carolynmichael/Documents/summerproject_updated 2/S99_models/S99-v00-Z014-IMF2.3_100myr.spectrum', obs_data)
	
	model_data = [Z001_data, Z002_data, Z008_data, Z014_data]

	data = [Z001_data[0], Z001_data[3], Z001_data[5]]
	unmasked_data = [Z001_data[1], Z001_data[4]]


	#define free parameters 
	A_v = float(input('A_v value: '))
	d = float(input('delta value: '))
	B = float(input('B value: '))
	z = float(input('z value: '))

	theta = [A_v, d, B, z]

	flux_mod = cs.fluxmod2(model_data, data[0]/1.e4, theta)

	beta = cs.normalisingfactor2(flux_mod, data[1], data[2])
	
	yM = flux_mod*beta

	print('the chi-squared value is: ' + str(cs.chisquared2(flux_mod, data[1], data[2], beta)))

	np.savetxt('bestfit_9.0', X=np.column_stack((data[0], yM)))
		
	'''
	master_grid_bump = np.arange(1000, 5500, 1)
	A_cullen, A_calzetti = af.attenuation(master_grid_bump/1.e4)
	A_mod_low = af.attenuation_mod(master_grid_bump/1.e4, 0.09, 0.97)
	A_mod_med = af.attenuation_mod(master_grid_bump/1.e4, 0.08, 0.94)
	A_mod_high = af.attenuation_mod(master_grid_bump/1.e4, 0.04, 0.76)
	Az45_mod_low = af.attenuation_mod(master_grid_bump/1.e4, -0.08, 0.51)
	Az45_mod_high = af.attenuation_mod(master_grid_bump/1.e4, -0.06, 0.51)


	fig, ax = plt.subplots()
	ax.plot(master_grid_bump, A_calzetti, color='black', lw=1)
	ax.plot(master_grid_bump, A_mod_low, color='green', lw=1, label='low stellar masses')
	ax.plot(master_grid_bump, A_mod_med, color='orange', lw=1, label='medium stellar masses')
	ax.plot(master_grid_bump, A_mod_high, color='red', lw=1, label='high stellar masses')
	ax.plot(master_grid_bump, Az45_mod_low, color='blue', lw=1, label='low stellar masses z45')
	ax.plot(master_grid_bump, Az45_mod_high, color='purple', lw=1, label='high stellar masses z45')
	ax.set_xlabel("Rest Wavelength ($\mathrm{\AA}$)", size=16)
	ax.set_ylabel("$A_\lambda / A_v$", size=16)
	ax.set_xticks(np.arange(1000, 5501, 250), minor=True)
	ax.set_yticks(np.arange(0, 4.1, 0.25), minor=True)
	ax.set_xticks(np.arange(1000, 5501, 1000))
	ax.set_yticks(np.arange(0, 4.1, 1))
	ax.tick_params(direction='in', which='major', length=6, width=0.5, grid_alpha=0.5, grid_linewidth=0.5, bottom=True, top=True, left=True, right=True)
	ax.tick_params(direction='in', which='minor', length=4, width=0.5, grid_alpha=0.5, grid_linewidth=0.5, bottom=True, top=True, left=True, right=True)
	ax.set_ylim(0.75, 4)
	ax.set_xlim(1000, 5500)
	ax.legend()
	plt.show()

	'''


	fig, ax = plt.subplots()
	ax.plot(master_grid, flux_obs, color='black', lw=1.)
	for i in range(len(unmasked_data[1])):
		ax.plot(unmasked_data[0][i],unmasked_data[1][i], color='deepskyblue', lw=1)
	ax.plot(data[0], yM, color='limegreen', lw=1.5)
	ax.plot(master_grid, error_obs, color='red', lw=1.)
	ax.set_xlabel("Rest Wavelength ($\mathrm{\AA}$)", size=16)
	ax.set_ylabel("f$_\lambda$ ($erg$ $s^{-1}cm^{-2}\mathrm{\AA}^{-1}$)", size=16)
	ax.set_xticks(np.arange(1200, 2001, 50), minor=True)
	ax.set_xticks(np.arange(1200, 2001, 200))
	ax.set_yticks(np.arange(0, 8*10**(-19), 0.125*10**(-19)), minor=True)
	ax.tick_params(direction='in', which='major', length=6, width=0.5, grid_alpha=0.5, grid_linewidth=0.5, bottom=True, top=True, left=True, right=True)
	ax.tick_params(direction='in', which='minor', length=4, width=0.5, grid_alpha=0.5, grid_linewidth=0.5, bottom=True, top=True, left=True, right=True)
	#plt.xlim(1200, 2000)
	ax.set_ylim(0, 3.5*10**(-19))
	ax.set_xlim(1200, 2000)
	plt.show()
	'''

	fig, ax = plt.subplots()
	ax.plot(master_grid1, flux_obs1, color='black', lw=1.)
	for i in range(len(unmasked_data[1])):
		ax.plot(unmasked_data[0][i],unmasked_data[1][i], color='deepskyblue', lw=1)
	ax.plot(master_grid1, error_obs1, color='red', lw=1.)
	ax.set_xlabel("Rest Wavelength ($\mathrm{\AA}$)", size=16)
	ax.set_ylabel("f$_\lambda$ ($erg$ $s^{-1}cm^{-2}\mathrm{\AA}^{-1}$)", size=16)
	ax.set_xticks(np.arange(1200, 2001, 50), minor=True)
	ax.set_xticks(np.arange(1200, 2001, 200))
	ax.set_yticks(np.arange(0, 8*10**(-19), 0.125*10**(-19)), minor=True)
	ax.tick_params(direction='in', which='major', length=6, width=0.5, grid_alpha=0.5, grid_linewidth=0.5, bottom=True, top=True, left=True, right=True)
	ax.tick_params(direction='in', which='minor', length=4, width=0.5, grid_alpha=0.5, grid_linewidth=0.5, bottom=True, top=True, left=True, right=True)
	#plt.xlim(1200, 2000)
	ax.set_ylim(0, 3.5*10**(-19))
	ax.set_xlim(1200, 2000)
	plt.show()
	'''


if __name__ == "__main__": 
    main()

