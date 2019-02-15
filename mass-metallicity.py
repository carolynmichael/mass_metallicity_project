from astropy.io import fits
import matplotlib.patches as mpatches
import numpy as np
import matplotlib.pyplot as plt
import os, sys
from scipy.interpolate import interp1d
from astropy.stats import bootstrap
from astropy.utils import NumpyRNGContext
import attenuationfunc as af
import chisq as cs
import stackspectra as ss
import glob
from astropy.io import ascii
import math as m
import itertools
import matplotlib.cm as cm
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

def readtables(path):
	filenames = []
	masses = []
	for filename in glob.glob(os.path.join(path, '*.dat')):
		filenames.append(os.path.basename(filename))
		data = ascii.read(filename)  
		lmass_rjm = np.array(data['Masses'])
		masses.append(np.mean(lmass_rjm))
	
	return masses, filenames

mean_masses, filenames = readtables('/Users/carolynmichael/Documents/summerproject_updated 2/pythonfits/filteredmass_data')

filenames_sorted = [filenames[2], filenames[0], filenames[3], filenames[4], filenames[1]]

mean_masses.sort()

masses_rat = 10**(np.array(mean_masses))/(2*10**30)

masses = np.array(masses_rat, dtype=float)


z = [0.0014, 0.0015, 0.0018, 0.0023, 0.0021]

z_error = [[0.00017, 0.00025], [0.00023, 0.00035], [0.00032, 0.00056], [0.0006, 0.0010], [0.0003, 0.0005]]

y_err = []

for i in range(5):
	data_err = []
	for j in range(2):
		data_err.append(z_error[i][j]/z[i])
	y_err.append(data_err)

z_zsun = np.log10(np.array(z)/0.02)

t = np.arange(len(z_zsun))

y1 = 0.4*(np.array(mean_masses) - 10) + 0.67*m.exp(0.0*(-0.50)) - 1.04

y2 = 0.4*(np.array(mean_masses) - 10) + 0.67*m.exp(3.9*(-0.50)) - 1.04

mean_masses_z34 = [mean_masses[0], mean_masses[2], mean_masses[4]]
mean_masses_z45 = [mean_masses[1], mean_masses[3]]

z_zsun_z34 = [z_zsun[0], z_zsun[2], z_zsun[4]]
z_zsun_z45 = [z_zsun[1], z_zsun[3]]

print(mean_masses_z45)

from scipy.stats import linregress
print(linregress(mean_masses, y2))
print(linregress(mean_masses, y1))

y2_z34 = [y2[0], y2[2], y2[4]]
y2_z45 = [y2[1], y2[3]]

residuals = np.array(z_zsun) - np.array(y2)

residuals_z34 = np.array(z_zsun_z34) - np.array(y2_z34)
residuals_z45 = np.array(z_zsun_z45) - np.array(y2_z45)

np.savetxt('y_obs_exp', X=np.column_stack((mean_masses, y2, z_zsun)))

y = 0.3999999999999999*np.arange(0, 15, 1) + -4.908069177613189
y_2 = 0.3999999999999999*np.arange(0, 15, 1) + -4.37

fig, ax = plt.subplots()
for i in range(5):
	ax.errorbar(mean_masses[i], z_zsun[i], yerr = np.array([[y_err[i]]]).T, color = 'black', lw=0.5, linestyle="None", zorder=1)
#ax.plot(np.unique(mean_masses), np.poly1d(np.polyfit(mean_masses, z_zsun, 1))(np.unique(mean_masses)), color='black', lw=1.5, dashes=[6, 2], zorder=2)
plt.plot(np.arange(0,15,1), y_2, lw=1, zorder=2, color = 'red', label = '$ z = 0.0$')
#ax.plot(mean_masses, y2, color='blue', lw=1, zorder=2, label = 'X. Ma et al.')
ax.plot(np.arange(0, 15, 1), y, lw=1, zorder=2, color = 'blue', label = '$ z = 3.9$')
ax.scatter(mean_masses_z34, z_zsun_z34, color = 'green', zorder=3)
ax.scatter(mean_masses_z45, z_zsun_z45, color = 'orange', zorder=3)
ax.set_xlabel("$log_{10} (M_*/M_{\odot})$", size=16)
ax.set_ylabel("$log_{10} (Z_*/Z_{\odot})$", size=16)
ax.set_xticks(np.arange(8.8, 10.2, 0.1), minor=True)
ax.set_yticks(np.arange(-1.5, -0.2, 0.05), minor=True)
ax.set_xticks(np.arange(8.8, 10.2, 0.4))
ax.set_yticks(np.arange(-1.5, -0.2, 0.2))
ax.tick_params(direction='in', which='major', length=6, width=0.5, grid_alpha=0.5, grid_linewidth=0.5, bottom=True, top=True, left=True, right=True)
ax.tick_params(direction='in', which='minor', length=4, width=0.5, grid_alpha=0.5, grid_linewidth=0.5, bottom=True, top=True, left=True, right=True)
ax.set_ylim(-1.5, -0.2)
ax.set_xlim(8.7, 10.2)
ax.legend()	
plt.show()

fig, ax = plt.subplots()
for i in range(5):
	ax.errorbar(mean_masses[i], residuals[i], yerr = np.array([[y_err[i]]]).T, color = 'black', lw=0.5, linestyle="None", zorder=1)
ax.scatter(mean_masses_z34, residuals_z34, color = 'green', zorder=3)
ax.scatter(mean_masses_z45, residuals_z45, color = 'orange', zorder=3)
ax.axhline(linewidth=1, color='black')
ax.set_xlabel("$log_{10} (M_*/M_{\odot})$", size=16)
ax.set_ylabel("Residuals", size=16)
ax.set_xticks(np.arange(8.8, 10.2, 0.1), minor=True)
ax.set_yticks(np.arange(-1.0, 1.0, 0.05), minor=True)
ax.set_xticks(np.arange(8.8, 10.2, 0.4))
ax.set_yticks(np.arange(-1.0, 1.0, 0.2))
ax.tick_params(direction='in', which='major', length=6, width=0.5, grid_alpha=0.5, grid_linewidth=0.5, bottom=True, top=True, left=True, right=True)
ax.tick_params(direction='in', which='minor', length=4, width=0.5, grid_alpha=0.5, grid_linewidth=0.5, bottom=True, top=True, left=True, right=True)
ax.set_ylim(-0.6, 0.6)
ax.set_xlim(8.7, 10.2)
ax.legend()	
plt.show()


'''
marker = itertools.cycle(('+', 'o', '*')) 

fig = plt.figure()
ax = fig.add_subplot(111)

for i in range(3):
	for q,p in zip(mean_masses, z_zsun):
		plt.errorbar(q, p, yerr = np.array([[z_error_scaled[i]]]).T, lw=0.5, linestyle = '', marker=marker.next())
		#ax.plot(q,p, linestyle = '', marker=marker.next())

for i in range(3):
	plt.errorbar(mean_masses[i], z_zsun[i], yerr = np.array([[z_error_scaled[i]]]).T, fmt='ko', lw=0.5)

plt.plot(np.unique(mean_masses), np.poly1d(np.polyfit(mean_masses, z_zsun, 1))(np.unique(mean_masses)), color='black', lw=1.5, dashes=[6, 2])
plt.xlabel("log$(M_*)$", size=16)
plt.ylabel("log$(Z_*/Z_{\odot})$", size=16)
plt.show()



for marker in ['o', '>', 'd']:
    plt.plot(mean_masses, z_zsun, marker, label="marker='{0}'".format(marker))
plt.legend(numpoints=1)
plt.plot(mean_masses, z_zsun, color='black', lw=1.5, dashes=[6, 2])
#plt.xlabel("$M_*$/$M_{\odot}$", size=16)
plt.xlabel("log$(M_*)$", size=16)
plt.ylabel("log$(Z_*/Z_{\odot})$", size=16)
plt.show()
'''
