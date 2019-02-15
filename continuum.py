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

def loadspectrum(filename):
	
	data = np.loadtxt(filename)

	wavelength = data[:,0]
	log1_flux_int = data[:,1]
	norm = data[:,2]

	flux_int = 10**(log1_flux_int)

	log_flux_int = np.log10(flux_int)
	
	f_cont = flux_int/norm
	
	data_list = [wavelength, np.log10(flux_int), np.log10(f_cont)]

	return data_list



def main():

	Z001_data = loadspectrum('/Users/carolynmichael/Documents/summerproject_updated 2/S99_models/S99-v00-Z001-IMF2.3_100myr.spectrum')
	
	Z002_data = loadspectrum('/Users/carolynmichael/Documents/summerproject_updated 2/S99_models/S99-v00-Z002-IMF2.3_100myr.spectrum')
	
	Z008_data = loadspectrum('/Users/carolynmichael/Documents/summerproject_updated 2/S99_models/S99-v00-Z008-IMF2.3_100myr.spectrum')
	
	Z014_data = loadspectrum('/Users/carolynmichael/Documents/summerproject_updated 2/S99_models/S99-v00-Z014-IMF2.3_100myr.spectrum')
	
	Z040_data = loadspectrum('/Users/carolynmichael/Documents/summerproject_updated 2/S99_models/S99-v00-Z040-IMF2.3_100myr.spectrum')
	


	fig = plt.figure()

	ax = fig.add_subplot(111) 
	# Turn off axis lines and ticks of the big subplot
	ax.spines['top'].set_color('none')
	ax.spines['bottom'].set_color('none')
	ax.spines['left'].set_color('none')
	ax.spines['right'].set_color('none')
	ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
	ax.set_ylabel('log$f_\lambda$', size = 16)
	ax.set_xlabel('Rest Wavelength ($\mathrm{\AA}$)', size = 16)

	ax1 = fig.add_subplot(321)
	ax1.plot(Z001_data[0], Z001_data[1], linewidth = 0.5, color = 'blue', label = '$Z_*$ = 0.001')
	ax1.plot(Z001_data[0], Z001_data[2], color = 'orange', linewidth = 0.5)
	ax1.set_xticks(np.arange(1000, 2500, 125), minor=True)
	ax1.set_yticks(np.arange(39.5, 42, 0.125), minor=True)
	ax1.set_xticks(np.arange(1000, 2501, 500))
	ax1.set_yticks(np.arange(39.75, 42, 0.5))
	ax1.tick_params(direction='in', which='major', length=6, width=0.5, grid_alpha=0.5, grid_linewidth=0.5, bottom=True, top=True, left=True, right=True)
	ax1.tick_params(direction='in', which='minor', length=4, width=0.5, grid_alpha=0.5, grid_linewidth=0.5, bottom=True, top=True, left=True, right=True)
	ax1.set_ylim(39.5, 41)
	#ax1.set_xticklabels(['1000', '1500', '2000', '2500'])
	ax1.set_xlim(950, 2501)
	ax1.text(1900, 40.5, '$Z_*$ = 0.001')

	ax2 = fig.add_subplot(322)
	ax2.plot(Z002_data[0], Z002_data[1], linewidth = 0.5, label = '$Z_*$ = 0.002')
	ax2.plot(Z002_data[0], Z002_data[2], linewidth = 0.5)
	ax2.set_xticks(np.arange(1000, 2500, 125), minor=True)
	ax2.set_yticks(np.arange(39.5, 42, 0.125), minor=True)
	ax2.set_xticks(np.arange(1000, 2501, 500))
	ax2.set_yticks(np.arange(39.75, 42, 0.5))
	ax2.tick_params(direction='in', which='major', length=6, width=0.5, grid_alpha=0.5, grid_linewidth=0.5, bottom=True, top=True, left=True, right=True)
	ax2.tick_params(direction='in', which='minor', length=4, width=0.5, grid_alpha=0.5, grid_linewidth=0.5, bottom=True, top=True, left=True, right=True)
	ax2.set_ylim(39.5, 41)
	ax2.set_xlim(950, 2501)
	ax2.text(1900, 40.5, '$Z_*$ = 0.002')

	ax3 = fig.add_subplot(323)	
	ax3.plot(Z008_data[0], Z008_data[1], linewidth = 0.5, label = '$Z_*$ = 0.008')
	ax3.plot(Z008_data[0], Z008_data[2], linewidth = 0.5)
	ax3.set_xticks(np.arange(1000, 2500, 125), minor=True)
	ax3.set_yticks(np.arange(39.5, 42, 0.125), minor=True)
	ax3.set_xticks(np.arange(1000, 2501, 500))
	ax3.set_yticks(np.arange(39.75, 42, 0.5))
	ax3.tick_params(direction='in', which='major', length=6, width=0.5, grid_alpha=0.5, grid_linewidth=0.5, bottom=True, top=True, left=True, right=True)
	ax3.tick_params(direction='in', which='minor', length=4, width=0.5, grid_alpha=0.5, grid_linewidth=0.5, bottom=True, top=True, left=True, right=True)
	ax3.set_ylim(39.5, 41)
	ax3.set_xlim(950, 2501)
	ax3.text(1900, 40.5, '$Z_*$ = 0.008')

	ax4 = fig.add_subplot(324)
	ax4.plot(Z014_data[0], Z014_data[1], linewidth = 0.5, label = '$Z_*$ = 0.014')
	ax4.plot(Z014_data[0], Z014_data[2], linewidth = 0.5)	
	ax4.set_xticks(np.arange(1000, 2500, 125), minor=True)
	ax4.set_yticks(np.arange(39.5, 42, 0.125), minor=True)
	ax4.set_xticks(np.arange(1000, 2501, 500))
	ax4.set_yticks(np.arange(39.75, 42, 0.5))
	ax4.tick_params(direction='in', which='major', length=6, width=0.5, grid_alpha=0.5, grid_linewidth=0.5, bottom=True, top=True, left=True, right=True)
	ax4.tick_params(direction='in', which='minor', length=4, width=0.5, grid_alpha=0.5, grid_linewidth=0.5, bottom=True, top=True, left=True, right=True)
	ax4.set_ylim(39.5, 41)
	ax4.set_xlim(950, 2501)
	ax4.text(1900, 40.5, '$Z_*$ = 0.014')

	ax5 = fig.add_subplot(325)
	ax5.plot(Z040_data[0], Z040_data[1], linewidth = 0.5, label = '$Z_*$ = 0.040')
	ax5.plot(Z040_data[0], Z040_data[2], linewidth = 0.5)
	ax5.set_xticks(np.arange(1000, 2500, 125), minor=True)
	ax5.set_yticks(np.arange(39.5, 42, 0.125), minor=True)
	ax5.set_xticks(np.arange(1000, 2501, 500))
	ax5.set_yticks(np.arange(39.75, 41.75, 0.5))
	ax5.tick_params(direction='in', which='major', length=6, width=0.5, grid_alpha=0.5, grid_linewidth=0.5, bottom=True, top=True, left=True, right=True)
	ax5.tick_params(direction='in', which='minor', length=4, width=0.5, grid_alpha=0.5, grid_linewidth=0.5, bottom=True, top=True, left=True, right=True)
	ax5.set_ylim(39.5, 41)
	ax5.set_xlim(950, 2501)	
	ax5.text(1900, 40.5, '$Z_*$ = 0.040')

	plt.show()


if __name__ == "__main__": 
    main()

