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
from astropy.io import ascii
from astropy.table import Table, Column, MaskedColumn
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

def main():
	
	#split spectra lists by mass values at third quartiles 
	
	spectra_id_lists, mass_lists = ss.pickoutspectraID2('/Users/carolynmichael/Documents/summerproject_updated 2/vandels_spectra_z45/', 'vandels_masses')
	
	data_low = Table([spectra_id_lists[0], mass_lists[0]], names=['Spectrum_ID', 'Masses'])
	ascii.write(data_low, 'lowvaluesz45.dat', overwrite=True)
	data_low = Table([spectra_id_lists[1], mass_lists[1]], names=['Spectrum_ID', 'Masses'])
	ascii.write(data_low, 'highvaluesz45.dat',overwrite=True)
	'''
	data_low = Table([spectra_id_lists[2], mass_lists[2]], names=['Spectrum_ID', 'Masses'])
	ascii.write(data_low, 'highvalues.dat')
	'''

	#extract data from filtered files
	wl_rest_low, flux_low, error_low = ss.readpickedspectra(spectra_id_lists[0], '/Users/carolynmichael/Documents/summerproject_updated 2/vandels_spectra_z45/')

	wl_rest_med, flux_med, error_med = ss.readpickedspectra(spectra_id_lists[1], '/Users/carolynmichael/Documents/summerproject_updated 2/vandels_spectra_z45/')


	wl_rest_list_low = []
	flux_list_low = []
	wl_rest_list_med = []
	flux_list_med = []
	error_list_low = []
	error_list_med = []

	for i in range(len(wl_rest_low)):
		mask = (1200.0<wl_rest_low[i]) & (wl_rest_low[i]<1601.0)
		wl_rest_list_low.append(wl_rest_low[i][mask])
		flux_list_low.append(flux_low[i][mask])
		error_list_low.append(error_low[i][mask])

	for i in range(len(wl_rest_med)):
		mask = (1200.0<wl_rest_med[i]) & (wl_rest_med[i]<1601.0)
		wl_rest_list_med.append(wl_rest_med[i][mask])	
		flux_list_med.append(flux_med[i][mask])		
		error_list_med.append(error_med[i][mask])
	'''
	for i in range(len(wl_rest_list_low)):
		print(min(wl_rest_list_low[i]))
		print(max(wl_rest_list_low[i]))
	
	for i in range(len(wl_rest_list_med)):
		print(min(wl_rest_list_med[i]))
		print(max(wl_rest_list_med[i]))
	'''

	#create master spectra:
	master_grid = np.arange(1200, 1600, 1)
	master_grid_extended = np.arange(950, 2500, 1)
	master_grid_bump = np.arange(1000, 5500, 1)

	fluxes_low_list, error_low_list = ss.interpolatespectra(wl_rest_list_low, flux_list_low, error_list_low, master_grid)
	fluxes_med_list, error_med_list = ss.interpolatespectra(wl_rest_list_med, flux_list_med, error_list_med, master_grid)
	#fluxes_high_list, error_high_list = ss.interpolatespectra(wl_rest_list_high, flux_list_high, error_list_high, master_grid)

	#find median values, discounting any nan values
	flux_low = np.nanmedian(fluxes_low_list, axis=0)
	flux_error_low = np.nanmedian(error_low_list, axis=0)
	
	flux_med = np.nanmedian(fluxes_med_list, axis=0)
	flux_error_med = np.nanmedian(error_med_list, axis=0)
	
	#flux_high = np.nanmedian(fluxes_high_list, axis=0)
	#flux_error_high = np.nanmedian(error_high_list, axis=0)

	#correct error values
	booterr_low_list = ss.bootstrap_err_correction(fluxes_low_list)
	booterr_med_list = ss.bootstrap_err_correction(fluxes_med_list)
	#booterr_high_list = ss.bootstrap_err_correction(fluxes_high_list)
	
	
	np.savetxt('vandelsz45_low-mass_data', X=np.column_stack((master_grid, flux_low, booterr_low_list)))
	np.savetxt('vandelsz45_high-mass_data', X=np.column_stack((master_grid, flux_med, booterr_med_list)))
	#np.savetxt('vandels_high-mass_data', X=np.column_stack((master_grid, flux_high, booterr_high_list)))

if __name__ == "__main__": 
    main()


