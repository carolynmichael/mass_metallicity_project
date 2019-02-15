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

    def loadspectrum(filename):

        data = np.loadtxt(filename)

        master_grid = data[:,0]
        flux_obs = data[:,1]
        error_obs = data[:,2]

        data_list = [master_grid, flux_obs, error_obs]

        return data_list

    def loadspectrum2(filename):

        data = np.loadtxt(filename)

        master_grid = data[:,0]
        flux_obs = data[:,1]

        data_list = [master_grid, flux_obs]

        return data_list

    obs_data = np.loadtxt('vandels_low-mass_data')

    obs_data1 = np.loadtxt('vandels_med-mass_data')

    obs_data2 = np.loadtxt('vandels_high-mass_data')

    obs_data3 = np.loadtxt('vandelsz45_low-mass_data')

    obs_data4 = np.loadtxt('vandelsz45_high-mass_data')

    
    #observed data for each mass bin
    data1 = loadspectrum('vandels_low-mass_data')
    data2 = loadspectrum('vandels_med-mass_data')
    data3 = loadspectrum('vandels_high-mass_data')
    data4 = loadspectrum('vandelsz45_low-mass_data')
    data5 = loadspectrum('vandelsz45_high-mass_data')

    #model data for each mass bin
    data1_mod = loadspectrum2('bestfit_8.91')
    data2_mod = loadspectrum2('bestfit_9.47')
    data3_mod = loadspectrum2('bestfit_9.98')
    data4_mod = loadspectrum2('bestfit_9.0')
    data5_mod = loadspectrum2('bestfit_9.68')



    Z001_data = cs.getdata('/Users/carolynmichael/Documents/summerproject_updated 2/S99_models/S99-v00-Z001-IMF2.3_100myr.spectrum', obs_data)

    Z002_data = cs.getdata('/Users/carolynmichael/Documents/summerproject_updated 2/S99_models/S99-v00-Z002-IMF2.3_100myr.spectrum', obs_data)

    Z008_data = cs.getdata('/Users/carolynmichael/Documents/summerproject_updated 2/S99_models/S99-v00-Z008-IMF2.3_100myr.spectrum', obs_data)

    Z014_data = cs.getdata('/Users/carolynmichael/Documents/summerproject_updated 2/S99_models/S99-v00-Z014-IMF2.3_100myr.spectrum', obs_data)

    model_data = [Z001_data, Z002_data, Z008_data, Z014_data]

    data = [Z001_data[0], Z001_data[3], Z001_data[5]]
    unmasked_data = [Z001_data[1], Z001_data[4]]


    Z001_data_med = cs.getdata('/Users/carolynmichael/Documents/summerproject_updated 2/S99_models/S99-v00-Z001-IMF2.3_100myr.spectrum', obs_data1)

    Z002_data_med = cs.getdata('/Users/carolynmichael/Documents/summerproject_updated 2/S99_models/S99-v00-Z002-IMF2.3_100myr.spectrum', obs_data1)

    Z008_data_med = cs.getdata('/Users/carolynmichael/Documents/summerproject_updated 2/S99_models/S99-v00-Z008-IMF2.3_100myr.spectrum', obs_data1)

    Z014_data_med = cs.getdata('/Users/carolynmichael/Documents/summerproject_updated 2/S99_models/S99-v00-Z014-IMF2.3_100myr.spectrum', obs_data1)

    model_data = [Z001_data_med, Z002_data_med, Z008_data_med, Z014_data_med]

    data_med = [Z001_data_med[0], Z001_data_med[3], Z001_data_med[5]]
    unmasked_data_med = [Z001_data_med[1], Z001_data_med[4]]

    Z001_data_high = cs.getdata('/Users/carolynmichael/Documents/summerproject_updated 2/S99_models/S99-v00-Z001-IMF2.3_100myr.spectrum', obs_data2)

    Z002_data_high = cs.getdata('/Users/carolynmichael/Documents/summerproject_updated 2/S99_models/S99-v00-Z002-IMF2.3_100myr.spectrum', obs_data2)

    Z008_data_high = cs.getdata('/Users/carolynmichael/Documents/summerproject_updated 2/S99_models/S99-v00-Z008-IMF2.3_100myr.spectrum', obs_data2)

    Z014_data_high = cs.getdata('/Users/carolynmichael/Documents/summerproject_updated 2/S99_models/S99-v00-Z014-IMF2.3_100myr.spectrum', obs_data2)

    model_data = [Z001_data_high, Z002_data_high, Z008_data_high, Z014_data_high]

    data_high = [Z001_data_high[0], Z001_data_high[3], Z001_data_high[5]]
    unmasked_data_high = [Z001_data_high[1], Z001_data_high[4]]

    Z001_data_z45 = cs.getdataz45('/Users/carolynmichael/Documents/summerproject_updated 2/S99_models/S99-v00-Z001-IMF2.3_100myr.spectrum', obs_data3)

    Z002_data_z45 = cs.getdataz45('/Users/carolynmichael/Documents/summerproject_updated 2/S99_models/S99-v00-Z002-IMF2.3_100myr.spectrum', obs_data3)

    Z008_data_z45 = cs.getdataz45('/Users/carolynmichael/Documents/summerproject_updated 2/S99_models/S99-v00-Z008-IMF2.3_100myr.spectrum', obs_data3)

    Z014_data_z45 = cs.getdataz45('/Users/carolynmichael/Documents/summerproject_updated 2/S99_models/S99-v00-Z014-IMF2.3_100myr.spectrum', obs_data3)

    model_data = [Z001_data_z45, Z002_data_z45, Z008_data_z45, Z014_data_z45]

    data_z45 = [Z001_data_z45[0], Z001_data_z45[3], Z001_data_z45[5]]
    unmasked_data_z45 = [Z001_data_z45[1], Z001_data_z45[4]]

    Z001_data_zh = cs.getdataz45('/Users/carolynmichael/Documents/summerproject_updated 2/S99_models/S99-v00-Z001-IMF2.3_100myr.spectrum', obs_data4)

    Z002_data_zh = cs.getdataz45('/Users/carolynmichael/Documents/summerproject_updated 2/S99_models/S99-v00-Z002-IMF2.3_100myr.spectrum', obs_data4)

    Z008_data_zh = cs.getdataz45('/Users/carolynmichael/Documents/summerproject_updated 2/S99_models/S99-v00-Z008-IMF2.3_100myr.spectrum', obs_data4)

    Z014_data_zh = cs.getdataz45('/Users/carolynmichael/Documents/summerproject_updated 2/S99_models/S99-v00-Z014-IMF2.3_100myr.spectrum', obs_data4)

    model_data = [Z001_data_zh, Z002_data_zh, Z008_data_zh, Z014_data_zh]

    data_z45 = [Z001_data_zh[0], Z001_data_zh[3], Z001_data_zh[5]]
    unmasked_data_zh = [Z001_data_zh[1], Z001_data_zh[4]]


    fig = plt.figure()

    ax = fig.add_subplot(111) 
    # Turn off axis lines and ticks of the big subplot
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
    ax.set_xlabel("Rest Wavelength ($\mathrm{\AA}$)", size=16)
    ax.set_ylabel("f$_\lambda$ ($erg$ $s^{-1}cm^{-2}\mathrm{\AA}^{-1}$)", size=16)

    ax1 = fig.add_subplot(321)
    ax1.plot(data1[0], data1[1], color='black', lw=0.5)
    for i in range(len(unmasked_data[1])):
        ax1.plot(unmasked_data[0][i],unmasked_data[1][i], color='deepskyblue', lw=0.5)
    ax1.plot(data1_mod[0], data1_mod[1], color='limegreen', lw=0.5)
    ax1.plot(data1[0], data1[2], color='red', lw=0.5)
    ax1.set_xticks(np.arange(1200, 2001, 50), minor=True)
    ax1.set_xticks(np.arange(1200, 2001, 200))
    ax1.set_yticks(np.arange(0, 8*10**(-19), 0.25*10**(-19)), minor=True)
    ax1.set_yticks(np.arange(0, 8*10**(-19), 1*10**(-19)))
    ax1.tick_params(direction='in', which='major', length=6, width=0.5, grid_alpha=0.5, grid_linewidth=0.5, bottom=True, top=True, left=True, right=True)
    ax1.tick_params(direction='in', which='minor', length=4, width=0.5, grid_alpha=0.5, grid_linewidth=0.5, bottom=True, top=True, left=True, right=True)
    ax1.set_ylim(0, 3.5*10**(-19))
    ax1.set_xlim(1200, 2000)
    ax1.text(1500, 2.5*10**(-19), '$log_{10}(M_*/M_{\odot})$ = 8.91')

    ax2 = fig.add_subplot(322)
    ax2.plot(data2[0], data2[1], color='black', lw=0.5)
    for i in range(len(unmasked_data_med[1])):
        ax2.plot(unmasked_data_med[0][i],unmasked_data_med[1][i], color='deepskyblue', lw=0.5)
    ax2.plot(data2_mod[0], data2_mod[1], color='limegreen', lw=0.75)
    ax2.plot(data2[0], data2[2], color='red', lw=0.5)
    ax2.set_xticks(np.arange(1200, 2001, 50), minor=True)
    ax2.set_xticks(np.arange(1200, 2001, 200))
    ax2.set_yticks(np.arange(0, 8*10**(-19), 0.25*10**(-19)), minor=True)
    ax2.set_yticks(np.arange(0, 8*10**(-19), 1*10**(-19)))
    ax2.tick_params(direction='in', which='major', length=6, width=0.5, grid_alpha=0.5, grid_linewidth=0.5, bottom=True, top=True, left=True, right=True)
    ax2.tick_params(direction='in', which='minor', length=4, width=0.5, grid_alpha=0.5, grid_linewidth=0.5, bottom=True, top=True, left=True, right=True)
    #plt.xlim(1200, 2000)
    ax2.set_ylim(0, 3.5*10**(-19))
    ax2.set_xlim(1200, 2000)
    ax2.text(1500, 2.5*10**(-19), '$log_{10}(M_*/M_{\odot})$ = 9.47')

    ax3 = fig.add_subplot(323)	
    ax3.plot(data3[0], data3[1], color='black', lw=0.5)
    for i in range(len(unmasked_data_high[1])):
        ax3.plot(unmasked_data_high[0][i],unmasked_data_high[1][i], color='deepskyblue', lw=0.5)
    ax3.plot(data3_mod[0], data3_mod[1], color='limegreen', lw=0.75)
    ax3.plot(data3[0], data3[2],  color='red', lw=0.5)
    ax3.set_xticks(np.arange(1200, 2001, 50), minor=True)
    ax3.set_xticks(np.arange(1200, 2001, 200))
    ax3.set_yticks(np.arange(0, 8*10**(-19), 0.25*10**(-19)), minor=True)
    ax3.set_yticks(np.arange(0, 8*10**(-19), 1*10**(-19)))
    ax3.tick_params(direction='in', which='major', length=6, width=0.5, grid_alpha=0.5, grid_linewidth=0.5, bottom=True, top=True, left=True, right=True)
    ax3.tick_params(direction='in', which='minor', length=4, width=0.5, grid_alpha=0.5, grid_linewidth=0.5, bottom=True, top=True, left=True, right=True)
    #plt.xlim(1200, 2000)
    ax3.set_ylim(0, 3.5*10**(-19))
    ax3.set_xlim(1200, 2000)
    ax3.text(1500, 2.5*10**(-19), '$log_{10}(M_*/M_{\odot})$ = 9.98')

    ax4 = fig.add_subplot(324)
    ax4.plot(data4[0], data4[1], color='black', lw=0.5)
    for i in range(len(unmasked_data_z45[1])):
        ax4.plot(unmasked_data_z45[0][i],unmasked_data_z45[1][i], color='deepskyblue', lw=0.5)
    ax4.plot(data4_mod[0], data4_mod[1], color='limegreen', lw=0.75)
    ax4.plot(data4[0], data4[2], color='red', lw=0.5)
    ax4.set_xticks(np.arange(1200, 2001, 25), minor=True)
    ax4.set_xticks(np.arange(1200, 2001, 100))
    ax4.set_yticks(np.arange(0, 8*10**(-19), 0.25*10**(-19)), minor=True)
    ax4.set_yticks(np.arange(0, 8*10**(-19), 1*10**(-19)))
    ax4.tick_params(direction='in', which='major', length=6, width=0.5, grid_alpha=0.5, grid_linewidth=0.5, bottom=True, top=True, left=True, right=True)
    ax4.tick_params(direction='in', which='minor', length=4, width=0.5, grid_alpha=0.5, grid_linewidth=0.5, bottom=True, top=True, left=True, right=True)
    #plt.xlim(1200, 2000)
    ax4.set_ylim(0, 3.5*10**(-19))
    ax4.set_xlim(1200, 1601)
    ax4.text(1350, 2.5*10**(-19), '$log_{10}(M_*/M_{\odot})$ = 9.00')

    ax5 = fig.add_subplot(325)
    ax5.plot(data5[0], data5[1], color='black', lw=0.5)
    for i in range(len(unmasked_data_zh[1])):
        ax5.plot(unmasked_data_zh[0][i],unmasked_data_zh[1][i], color='deepskyblue', lw=0.5)
    ax5.plot(data5_mod[0], data5_mod[1], color='limegreen', lw=0.75)
    ax5.plot(data5[0], data5[2], color='red', lw=0.5)
    ax5.set_xticks(np.arange(1200, 2001, 25), minor=True)
    ax5.set_xticks(np.arange(1200, 2001, 100))
    ax5.set_yticks(np.arange(0, 8*10**(-19), 0.25*10**(-19)), minor=True)
    ax5.set_yticks(np.arange(0, 8*10**(-19), 1*10**(-19)))
    ax5.tick_params(direction='in', which='major', length=6, width=0.5, grid_alpha=0.5, grid_linewidth=0.5, bottom=True, top=True, left=True, right=True)
    ax5.tick_params(direction='in', which='minor', length=4, width=0.5, grid_alpha=0.5, grid_linewidth=0.5, bottom=True, top=True, left=True, right=True)
    #plt.xlim(1200, 2000)
    ax5.set_ylim(0, 3.5*10**(-19))
    ax5.set_xlim(1200, 1601)
    ax5.text(1350, 2.5*10**(-19), '$log_{10}(M_*/M_{\odot})$ = 9.68')

    plt.show()


if __name__ == "__main__": 
    main()

