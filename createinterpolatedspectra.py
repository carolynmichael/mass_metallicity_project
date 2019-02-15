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
import bayesian as b
import emcee
from scipy import integrate
import pandas as pd
import scipy.optimize as op
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)


# Read name of output file from command line

'''
if len(sys.argv)!=4:
	print("Wrong number of arguments.")
	print("Usage: " + sys.argv[0] + " <output file>")
	quit()
else:
	compositespectrum = sys.argv[1]
	attenutationcurve = sys.argv[2]
	attenuationbump = sys.argv[3]


outfile1 = open(compositespectrum, "w")
outfile2 = open(attenutationcurve, "w")
outfile3 = open(attenuationbump, "w")
'''

wl_rest_list, flux_list, flux_error_list = ss.readmultiplespectra('/Users/carolynmichael/Documents/summerproject_updated 2/vandels_spectra/')

#create master spectra:
master_grid = np.arange(1200, 2000, 1)
master_grid_extended = np.arange(950, 2500, 1)
master_grid_bump = np.arange(1000, 5500, 1)

y = np.mean(flux_list[10])
y1 = np.mean(flux_list[150])
y2 = np.mean(flux_list[250])
y3 = np.mean(flux_list[80])

print(len(flux_list))

ylim = 7*y
y1lim = 7*y1
y2lim = y2
y3lim = y3

'''
fig, ax = plt.subplots()
ax.plot(wl_rest_list[150], flux_list[150], color='black', lw=1.)
ax.plot(wl_rest_list[150], flux_error_list[150], color='red', lw=1.)
ax.set_xlabel("Rest Wavelength ($\mathrm{\AA}$)", size=16)
ax.set_ylabel("$f_\lambda$ ($erg$ $s^{-1}cm^{-2}\mathrm{\AA}^{-1}$)", size=16)
ax.set_xticks(np.arange(1200, 2001, 50), minor=True)
ax.set_xticks(np.arange(1200, 2001, 200))
ax.set_yticks(np.arange(-8*10**-19, 8*10**-19, 0.5*10**-19), minor=True)
ax.tick_params(direction='in', which='major', length=6, width=0.5, grid_alpha=0.5, grid_linewidth=0.5, bottom=True, top=True, left=True, right=True)
ax.tick_params(direction='in', which='minor', length=4, width=0.5, grid_alpha=0.5, grid_linewidth=0.5, bottom=True, top=True, left=True, right=True)
ax.set_ylim(-ylim, ylim)
ax.set_xlim(1200, 2000)
plt.show()
'''
fig = plt.figure()

ax = fig.add_subplot(111) 
# Turn off axis lines and ticks of the big subplot
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
ax.set_xlabel("Rest Wavelength ($\mathrm{\AA}$)", size=16)
ax.set_ylabel("$f_\lambda$ ($erg$ $s^{-1}cm^{-2}\mathrm{\AA}^{-1}$)", size=16)

ax1 = fig.add_subplot(221)
ax1.plot(wl_rest_list[150], flux_list[150], color='black', lw=0.5)
ax1.plot(wl_rest_list[150], flux_error_list[150], color='red', lw=0.5)
ax1.set_xticks(np.arange(1200, 2001, 50), minor=True)
ax1.set_xticks(np.arange(1200, 2001, 200))
ax1.set_yticks(np.arange(-10*10**-19, 10*10**-19,0.625*10**-19), minor=True)
ax1.tick_params(direction='in', which='major', length=6, width=0.5, grid_alpha=0.5, grid_linewidth=0.5, bottom=True, top=True, left=True, right=True)
ax1.tick_params(direction='in', which='minor', length=4, width=0.5, grid_alpha=0.5, grid_linewidth=0.5, bottom=True, top=True, left=True, right=True)
ax1.set_ylim(-y1lim, y1lim)
ax1.set_xlim(1200, 2000)

ax2 = fig.add_subplot(222)
ax2.plot(wl_rest_list[80], flux_list[80], color='black', lw=0.5)
ax2.plot(wl_rest_list[80], flux_error_list[80], color='red', lw=0.5)
ax2.set_xticks(np.arange(1200, 2001, 50), minor=True)
ax2.set_xticks(np.arange(1200, 2001, 200))
ax2.set_yticks(np.arange(-15*10**-19, 15*10**-19,1.25*10**-19), minor=True)
ax2.tick_params(direction='in', which='major', length=6, width=0.5, grid_alpha=0.5, grid_linewidth=0.5, bottom=True, top=True, left=True, right=True)
ax2.tick_params(direction='in', which='minor', length=4, width=0.5, grid_alpha=0.5, grid_linewidth=0.5, bottom=True, top=True, left=True, right=True)
ax2.set_ylim(-1.5*ylim, 1.5*ylim)
ax2.set_xlim(1200, 2000)

ax3 = fig.add_subplot(223)	
ax3.plot(wl_rest_list[10], flux_list[10], color='black', lw=0.5)
ax3.plot(wl_rest_list[10], flux_error_list[10], color='red', lw=0.5)
ax3.set_xticks(np.arange(1200, 2001, 50), minor=True)
ax3.set_xticks(np.arange(1200, 2001, 200))
ax3.set_yticks(np.arange(-15*10**-19, 15*10**-19, 1.25*10**-19), minor=True)
ax3.tick_params(direction='in', which='major', length=6, width=0.5, grid_alpha=0.5, grid_linewidth=0.5, bottom=True, top=True, left=True, right=True)
ax3.tick_params(direction='in', which='minor', length=4, width=0.5, grid_alpha=0.5, grid_linewidth=0.5, bottom=True, top=True, left=True, right=True)
ax3.set_ylim(-1.5*ylim, 1.5*ylim)
ax3.set_xlim(1200, 2000)

ax4 = fig.add_subplot(224)
ax4.plot(wl_rest_list[250], flux_list[250], color='black', lw=0.5)
ax4.plot(wl_rest_list[250], flux_error_list[250], color='red', lw=0.5)
ax4.set_xticks(np.arange(1200, 2001, 50), minor=True)
ax4.set_xticks(np.arange(1200, 2001, 200))
ax4.set_yticks(np.arange(-15*10**-19, 15*10**-19, 0.5*10**-19), minor=True)
ax4.tick_params(direction='in', which='major', length=6, width=0.5, grid_alpha=0.5, grid_linewidth=0.5, bottom=True, top=True, left=True, right=True)
ax4.tick_params(direction='in', which='minor', length=4, width=0.5, grid_alpha=0.5, grid_linewidth=0.5, bottom=True, top=True, left=True, right=True)
ax4.set_ylim(-0.5*ylim, 0.5*ylim)
ax4.set_xlim(1200, 2000)

plt.show()

'''
final_fluxes_list, final_flux_error_list = ss.interpolatespectra(wl_rest_list, flux_list, flux_error_list, master_grid)

#find median values, discounting any nan values
flux = np.nanmedian(final_fluxes_list, axis=0)
flux_error = np.nanmedian(final_flux_error_list, axis=0)


#check mean calculation by sigma clipping
flux_mean = np.nanmean(final_fluxes_list, axis=0)
flux_error_mean = np.nanmean(final_flux_error_list, axis=0)
#flux_check = sigma_clip(flux_mean, sigma=3)
#flux_error_check = sigma_clip(flux_error_mean, sigma=3)

#correct error values
boot_error_list = ss.bootstrap_err_correction(final_fluxes_list)


#calculation absolute attenuation for each model		
A_cullen, A_calzetti = af.attenuation(master_grid_bump/1.e4)
A_mod_bump = af.attenuation_mod(master_grid_bump/1.e4, 0.2, 0.5)
A_mod = af.attenuation_mod(master_grid_bump/1.e4, 0.2, 0)
A_mod_bump2 = af.attenuation_mod(master_grid_bump/1.e4, -0.5, 3)
A_mod2 = af.attenuation_mod(master_grid_bump/1.e4, -0.5, 0)

np.savetxt('observed_data', X=np.column_stack((master_grid, flux, boot_error_list)))


outfile2.write('{} {} {}'.format(master_grid_extended, A_cullen, A_calzetti))

outfile3.write('{} {} {} {} {}'.format(master_grid_bump, A_mod, A_mod_bump, A_mod2, A_mod_bump2))


fig, ax = plt.subplots()
ax.plot(master_grid, np.array(flux), color='black', lw=1.)
ax.plot(master_grid, np.array(boot_error_list), color='red', lw=1.)
ax.set_xlabel("Rest Wavelength ($\mathrm{\AA}$)", size=16)
ax.set_ylabel("$f_\lambda$ ($erg$ $s^{-1}cm^{-2}\mathrm{\AA}^{-1}$)", size=16)
ax.set_xticks(np.arange(1200, 2001, 50), minor=True)
ax.set_xticks(np.arange(1200, 2001, 200))
ax.set_yticks(np.arange(0, 3.5*10**(-19), 0.125*10**(-19)), minor=True)
ax.tick_params(direction='in', which='major', length=6, width=0.5, grid_alpha=0.5, grid_linewidth=0.5, bottom=True, top=True, left=True, right=True)
ax.tick_params(direction='in', which='minor', length=4, width=0.5, grid_alpha=0.5, grid_linewidth=0.5, bottom=True, top=True, left=True, right=True)
ax.set_ylim(0, 3.5*10**(-19))
ax.set_xlim(1200, 2000)
plt.show()


fig, ax = plt.subplots()
ax.plot(master_grid_bump, A_cullen, color='black', lw=1.)
ax.plot(master_grid_bump, A_calzetti, color='black', lw=1.)
ax.set_xlabel("Rest Wavelength ($\mathrm{\AA}$)", size=16)
ax.set_ylabel("$A_\lambda / A_v$", size=16)
ax.set_xticks(np.arange(1200, 5501, 25), minor=True)
ax.set_yticks(np.arange(0, 4, 0.25), minor=True)
plt.show()

#outfile2.write('{} {} {}'.format(master_grid_extended, A_cullen, A_calzetti))


fig, ax = plt.subplots()
ax.plot(master_grid_bump, A_mod, color='black', lw=1.5)
ax.plot(master_grid_bump, A_mod_bump, color='black', lw=1.5, dashes=[6, 2])
ax.plot(master_grid_bump, A_mod2, color='black', lw=1.5)
ax.plot(master_grid_bump, A_mod_bump2, color='black', lw=1.5, dashes=[6, 2])
ax.set_xlabel("Rest Wavelength ($\mathrm{\AA}$)", size=16)
ax.set_ylabel("$A_\lambda / A_v$", size=16)
ax.set_xticks(np.arange(1200, 5501, 25), minor=True)
ax.set_yticks(np.arange(0, 4, 0.25), minor=True)
plt.show()


#outfile3.write('{} {} {} {} {}'.format(master_grid_bump, A_mod, A_mod_bump, A_mod2, A_mod_bump2))

#outfile4.write('{} {}'.format(masked_wav, flux_best_fit))
fig, ax = plt.subplots()
ax.plot(master_grid, flux, color='black', lw=1.)
for i in range(len(unmasked_wav)):
	ax.plot(unmasked_wav[i],f_obs_unmask[i], color='deepskyblue', lw=1.5)
ax.plot(masked_wav, flux_best_fit, color='limegreen', lw=1.5)
ax.plot(master_grid, boot_error_list, color='red', lw=1.)
ax.set_xlabel("Rest Wavelength ($\mathrm{\AA}$)", size=16)
ax.set_ylabel("f$_\lambda$ ($erg$ $s^{-1}cm^{-2}\mathrm{\AA}^{-1}$)", size=16)
ax.set_xticks(np.arange(1200, 2001, 25), minor=True)
ax.set_yticks(np.arange(0, 3*10**(-19), 0.25*10**(-19)), minor=True)
#ax.tick_params(direction='in', which='both', length=6, width=0.5, grid_alpha=0.5, grid_linewidth=0.5)
#plt.xlim(1200, 2000)
plt.show()
'''

