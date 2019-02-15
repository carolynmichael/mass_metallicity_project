from astropy.io import fits
from astropy.stats import sigma_clip
import numpy as np
import matplotlib.pyplot as plt
import os, sys
import glob 
from scipy.interpolate import interp1d
from astropy.stats import bootstrap
from astropy.utils import NumpyRNGContext
import pandas as pd
import numpy.ma as ma
import csv
from astropy.io import ascii
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

def readmultiplespectra(path):

	#create empty lists for finding averages
	wl_rest_list = []
	flux_list = []
	flux_error_list = []

	#open all spectra files:
	for filename in glob.glob(os.path.join(path, '*.fits')):
		hdu = fits.open(filename)
		data = hdu[1].data
		hdr = hdu[0].header
		# get the redshift:
		redshift = float(hdr['HIERARCH PND Z'])

		# extract the flux and error (units are ergs/s/cm2/A):
		flux = hdu['EXR1D'].data
		flux_error = hdu['NOISE'].data

		# extract the wavelegnth information:
		wl0 = hdr['CRVAL1'] # intial wavelegnth
		dwl = hdr['CDELT1'] # wavelegnth interval
		naxis = flux.shape[0] # number of data points

		wlf = wl0 + (dwl * naxis) # work out the final wavelgnth

		wl_obs = np.arange(wl0, wlf, dwl) # setup the wavelegnth grid (observed frame)
		wl_rest = wl_obs / (1. + redshift) # shift to rest frame
		
		wl_rest_list.append(wl_rest)
		flux_list.append(flux)
		flux_error_list.append(flux_error)
	
	return wl_rest_list, flux_list, flux_error_list
	
def pickoutspectraID(path_vandels, path_masses):
	filenames = []
	#open all spectra files:
	filenames = []
	for filename in glob.glob(os.path.join(path_vandels, '*.fits')):
		filenames.append(os.path.basename(filename))
	
	data = ascii.read(path_masses)  
	spectra_id = np.array(data['col1'])
	z = data['col2']
	lmass_fpp = data['col3']
	lmass_rjm = data['col4']
	
	mask = np.in1d(spectra_id, filenames)
	mask = np.array(mask)
	spectra = spectra_id[mask]
	masses = lmass_rjm[mask]
	mass_ranges = pd.qcut(masses, 3)
	masses_lab = pd.qcut(masses, 3).codes

	mask_low = masses_lab == 0
	mask_med = masses_lab == 1
	mask_high = masses_lab == 2
	'''
	spectra_m_low = [spectra[np.array(mask_low)], masses[np.array(mask_low)]]
	spectra_m_med = [spectra[np.array(mask_med)], masses[np.array(mask_med)]]
	spectra_m_high = [spectra[np.array(mask_high)], masses[np.array(mask_high)]]
	'''
	
	spectra_id_lists = [spectra[np.array(mask_low)], spectra[np.array(mask_med)], spectra[np.array(mask_high)]]
	mass_lists = [masses[np.array(mask_low)], masses[np.array(mask_med)], masses[np.array(mask_high)]]
	
	return spectra_id_lists, mass_lists
	
def pickoutspectraID2(path_vandels, path_masses):
	filenames = []
	#open all spectra files:
	filenames = []
	for filename in glob.glob(os.path.join(path_vandels, '*.fits')):
		filenames.append(os.path.basename(filename))
	
	data = ascii.read(path_masses)  
	spectra_id = np.array(data['col1'])
	z = data['col2']
	lmass_fpp = data['col3']
	lmass_rjm = data['col4']
	
	mask = np.in1d(spectra_id, filenames)
	mask = np.array(mask)
	spectra = spectra_id[mask]
	masses = lmass_rjm[mask]
	mass_ranges = pd.qcut(masses, 2)
	masses_lab = pd.qcut(masses, 2).codes

	mask_low = masses_lab == 0
	mask_high = masses_lab == 1
	'''
	spectra_m_low = [spectra[np.array(mask_low)], masses[np.array(mask_low)]]
	spectra_m_med = [spectra[np.array(mask_med)], masses[np.array(mask_med)]]
	spectra_m_high = [spectra[np.array(mask_high)], masses[np.array(mask_high)]]
	'''
	
	spectra_id_lists = [spectra[np.array(mask_low)], spectra[np.array(mask_high)]]
	mass_lists = [masses[np.array(mask_low)], masses[np.array(mask_high)]]
	
	return spectra_id_lists, mass_lists
	
	

def readpickedspectra(file_list, path):

	#create empty lists for finding averages
	wl_rest_list = []
	flux_list = []
	flux_error_list = []

	#open all spectra files:
	for filename in file_list:
		filename = os.path.join(path, filename)
		hdu = fits.open(filename)
		data = hdu[1].data
		hdr = hdu[0].header
		# get the redshift:
		redshift = float(hdr['HIERARCH PND Z'])

		# extract the flux and error (units are ergs/s/cm2/A):
		flux = hdu['EXR1D'].data
		flux_error = hdu['NOISE'].data

		# extract the wavelegnth information:
		wl0 = hdr['CRVAL1'] # intial wavelegnth
		dwl = hdr['CDELT1'] # wavelegnth interval
		naxis = flux.shape[0] # number of data points

		wlf = wl0 + (dwl * naxis) # work out the final wavelgnth

		wl_obs = np.arange(wl0, wlf, dwl) # setup the wavelegnth grid (observed frame)
		wl_rest = wl_obs / (1. + redshift) # shift to rest frame
		
		wl_rest_list.append(wl_rest)
		flux_list.append(flux)
		flux_error_list.append(flux_error)
	
	return wl_rest_list, flux_list, flux_error_list
	

def interpolatespectra(wl_rest_list, flux_list, flux_error_list, master_grid):
	final_fluxes_list = []
	final_flux_error_list = []

	#interpolate spectra:
	for i in range(len(wl_rest_list)):
		#interpolate each x,y onto master grid
		f = interp1d(wl_rest_list[i], flux_list[i], fill_value="extrapolate")
		flux_final = f(master_grid)
		#interpolate y errors onto master grid
		f_error = interp1d(wl_rest_list[i], flux_error_list[i], fill_value="extrapolate")
		flux_error_final = f_error(master_grid)
		#append to empty lists
		final_fluxes_list.append(flux_final)
		final_flux_error_list.append(flux_error_final)
	
	return final_fluxes_list, final_flux_error_list
	
def bootstrap_err_correction(final_fluxes_list):
	#rearrange initial flux lists to hold all flux values 
	#at the Nth wavelength pixel in the ith element of the list
	final_fluxes_list = np.array(final_fluxes_list)

	new_flux_list = [0.0]*len(final_fluxes_list[0])

	for i in range(len(final_fluxes_list[0])):
		new_flux_list[i] = final_fluxes_list[:,i]


	#for all wavelengths at each N wavelength pixel, resample 
	#N times and after finding mean of each resample, append
	#stdev to list

	boot_error_list = []
	for i in range(len(new_flux_list)):
		with NumpyRNGContext(1):
			flux_bootresult = bootstrap(new_flux_list[i], 1000, bootfunc=np.median)
			flux_error_bootresult = np.nanstd(flux_bootresult)
			boot_error_list.append(flux_error_bootresult)
	
	return boot_error_list
