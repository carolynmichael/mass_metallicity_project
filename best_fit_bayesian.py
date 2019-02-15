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
	
	obs_data = np.loadtxt('vandelsz45_high-mass_data')

	master_grid = obs_data[:,0]
	flux_obs = obs_data[:,1]
	error_obs = obs_data[:,2]
	
	Z001_data = cs.getdataz45('/Users/carolynmichael/Documents/summerproject_updated 2/S99_models/S99-v00-Z001-IMF2.3_100myr.spectrum', obs_data)
	
	Z002_data = cs.getdataz45('/Users/carolynmichael/Documents/summerproject_updated 2/S99_models/S99-v00-Z002-IMF2.3_100myr.spectrum', obs_data)
	
	Z008_data = cs.getdataz45('/Users/carolynmichael/Documents/summerproject_updated 2/S99_models/S99-v00-Z008-IMF2.3_100myr.spectrum', obs_data)
	
	Z014_data = cs.getdataz45('/Users/carolynmichael/Documents/summerproject_updated 2/S99_models/S99-v00-Z014-IMF2.3_100myr.spectrum', obs_data)
	
	model_data = [Z001_data, Z002_data, Z008_data, Z014_data]
	
	#now do bayesian analysis on model fitting
	data = [Z001_data[0], Z001_data[3], Z001_data[5]]
	
	#cs.interpolate2(model_data)
	
	#nested sampling approach to bayesian model selection
	
	def log_likelihood(theta, data=data, model_data=model_data):
		x, y, sigma_y = data
		y_M = cs.fluxmod2(model_data, x/1.e4, theta)
		B = cs.normalisingfactor2(y_M, y, sigma_y)
		yM = y_M*B
		
		return -0.5 * np.sum(np.log(2 * np.pi * sigma_y ** 2) + (y - yM) ** 2 / sigma_y ** 2)

	def prior_transform(utheta):
		uAv, ud, uB, uz = utheta
		Av = 6.0 * uAv
		d = 2*(ud - 0.5)
		B = 1 * uB
		z = 0.013*uz + 0.001

		return Av, d, B, z

	def compute_ns(log_likelihood, prior_transform, data):
		dsampler = dynesty.DynamicNestedSampler(log_likelihood, prior_transform, ndim=4, bound='multi', sample='rwalk', update_interval=3.)
		#dsampler = dynesty.DynamicNestedSampler(log_likelihood, prior_transform, ndim=3)
		dsampler.run_nested()
		dres = dsampler.results
		return dres
	
	dres = compute_ns(log_likelihood, prior_transform, data)
	
	titles = ['$A_v$', '$\delta$', 'B', '$Z_*$']
	
	fig, axes = dyplot.cornerplot(dres, show_titles=True, labels=titles, fig=plt.subplots(4, 4, figsize=(15, 15)))
	#axes = axes.flatten()
	#axes.tick_params(direction='in', which='major', length=6, width=0.5, grid_alpha=0.5, grid_linewidth=0.5, bottom=True, top=True, left=True, right=True)
	#axes.tick_params(direction='in', which='minor', length=4, width=0.5, grid_alpha=0.5, grid_linewidth=0.5, bottom=True, top=True, left=True, right=True)
	plt.show()
	
if __name__ == "__main__": 
    main()
