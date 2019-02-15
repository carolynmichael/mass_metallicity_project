import numpy as np
import numpy.ma as ma
from scipy import stats
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

def loadspectrum(filename):
	
	data = np.loadtxt(filename)

	wavelength = data[250:1050,0]
	flux_density_wav = data[250:1050,1]

	return flux_density_wav, wavelength

	
def maskwavelength(wav, flux_int, f_obs, boot_error_list):
	
	mask7 = (wav>1220) & (wav<1255)
	mask8 = (wav>1269) & (wav<1292)
	mask9 = (wav>1311) & (wav<1329)
	mask10 = (wav>1339) & (wav<1364)
	mask11 = (wav>1372) & (wav<1390)
	mask12 = (wav>1395) & (wav<1399)
	mask13 = (wav>1403) & (wav<1522)
	mask14 = (wav>1527) & (wav<1532)
	mask15 = (wav>1535) & (wav<1542)
	mask16 = (wav>1551) & (wav<1607)
	mask17 = (wav>1609) & (wav<1658)
	mask18 = (wav>1674) & (wav<1709)
	mask19 = (wav>1710) & (wav<1741)
	mask20 = (wav>1742) & (wav<1752)
	mask21 = (wav>1753) & (wav<1846)
	mask22 = (wav>1863) & (wav<1879)
	mask23 = (wav>1884) & (wav<1904)
	mask24 = (wav>1919) & (wav<2000)

	mask = mask7 | mask8 | mask9 | mask10 | mask11 | mask12 | mask13 | mask14 | mask15 | mask16 | mask17 | mask18 | mask19 | mask20 | mask21 | mask22 | mask23 | mask24

	masked_wav = wav[mask]
	f_int_mask = flux_int[mask]
	f_obs_mask = f_obs[mask]
	sigma_mask = np.array(boot_error_list)[mask]

	f_obs_unmask = [f_obs[0:21],f_obs[55:70],f_obs[92:112],f_obs[129:140],
	f_obs[164:173],f_obs[190:196],f_obs[199:204],f_obs[322:328],f_obs[332:336],
	f_obs[342:352],f_obs[407:410],f_obs[458:475],f_obs[509:511],f_obs[541:543],
	f_obs[552:554],f_obs[646:664],f_obs[679:685],f_obs[704:720]]

	unmasked_wav = [wav[0:21],wav[55:70],wav[92:112],wav[129:140],
	wav[164:173],wav[190:196],wav[199:204],wav[322:328],wav[332:336],
	wav[342:352],wav[407:410],wav[458:475],wav[509:511],wav[541:543],
	wav[552:554],wav[646:664],wav[679:685],wav[704:720]]

	return masked_wav, unmasked_wav, f_int_mask, f_obs_mask, f_obs_unmask, sigma_mask

def maskwavelengthz45(wav, flux_int, f_obs, boot_error_list):
	
	mask7 = (wav>1220) & (wav<1255)
	mask8 = (wav>1269) & (wav<1292)
	mask9 = (wav>1311) & (wav<1329)
	mask10 = (wav>1339) & (wav<1364)
	mask11 = (wav>1372) & (wav<1390)
	mask12 = (wav>1395) & (wav<1399)
	mask13 = (wav>1403) & (wav<1522)
	mask14 = (wav>1527) & (wav<1532)
	mask15 = (wav>1535) & (wav<1542)
	mask16 = (wav>1551) & (wav<1601)

	maskz45 = mask7 | mask8 | mask9 | mask10 | mask11 | mask12 | mask13 | mask14 | mask15 | mask16

	masked_wav = wav[maskz45]
	f_int_mask = flux_int[maskz45]
	f_obs_mask = f_obs[maskz45]
	sigma_mask = np.array(boot_error_list)[maskz45]

	f_obs_unmask = [f_obs[0:21],f_obs[55:70],f_obs[92:112],f_obs[129:140],
	f_obs[164:173],f_obs[190:196],f_obs[199:204],f_obs[322:328],f_obs[332:336],
	f_obs[342:352]]

	unmasked_wav = [wav[0:21],wav[55:70],wav[92:112],wav[129:140],
	wav[164:173],wav[190:196],wav[199:204],wav[322:328],wav[332:336],
	wav[342:352]]

	return masked_wav, unmasked_wav, f_int_mask, f_obs_mask, f_obs_unmask, sigma_mask
	
	
def getdata(filename, obs_data):
		
	log_flux_int, wavelength = loadspectrum(filename)

	flux_int = 10**(log_flux_int)
	
	master_grid = obs_data[:,0]
	flux_obs = obs_data[:,1]
	error_obs = obs_data[:,2]

	masked_wav, unmasked_wav, f_int_mask, f_obs_mask, f_obs_unmask, sigma_mask = maskwavelength(master_grid, flux_int, flux_obs, error_obs)
	
	data = masked_wav, unmasked_wav, f_int_mask, f_obs_mask, f_obs_unmask, sigma_mask
	
	return data

def getdataz45(filename, obs_data):
		
	log_flux_int, wavelength = loadspectrum(filename)

	mask_int = (1200<wavelength) & (wavelength<1601)

	flux_int = 10**(log_flux_int)
	
	master_grid = obs_data[:,0]
	flux_obs = obs_data[:,1]
	error_obs = obs_data[:,2]

	masked_wav, unmasked_wav, f_int_mask, f_obs_mask, f_obs_unmask, sigma_mask = maskwavelengthz45(master_grid, flux_int[mask_int], flux_obs, error_obs)
	
	data = masked_wav, unmasked_wav, f_int_mask, f_obs_mask, f_obs_unmask, sigma_mask
	
	return data
	
def fluxmod(f_int_mask, rest_wl):
	
	rest_wl_m = rest_wl*1.e4

	A_v = np.array(np.linspace(0, 3.0, 31))
	d = np.array(np.linspace(-1.0, 1.0, 11))
	B = np.array(np.linspace(0, 3.0, 31))

	flux_mod = []
	A_v_list = []
	B_list = []
	d_list = []

	for i in range(len(A_v)):	
		for j in range(len(B)):
			for k in range(len(d)):
				R_mod = 4.05/(5.05*(0.8**d[k]) - 4.05)
				D = (B[j]*(rest_wl_m**2)*(350**2))/((rest_wl_m**2 - (2175**2))**2 + (rest_wl_m**2)*(350**2))
				k_calzetti = 2.659*(-2.156 + 1.509*rest_wl**(-1) - 0.198*rest_wl**(-2) + 0.011*rest_wl**(-3)) + 4.05
				A_mod = (k_calzetti*(R_mod/4.05)*((rest_wl_m/5500)**d[k]) + D)/R_mod
				A_lambda = A_mod*A_v[i]
				flux_ob = []
				A_v_list.append(A_v[i])
				B_list.append(B[j])
				d_list.append(d[k])
				for l in range(len(f_int_mask)):
					flux_obs = f_int_mask[l]*10**(-0.4*A_lambda[l])
					flux_ob.append(flux_obs)
				flux_mod.append(flux_ob)
	

	return flux_mod, A_v_list, B_list, d_list
	
def interpolate(model_data, z):
	
	f_models = []
	for datalist in model_data:
		f_models.append(datalist[2])
	
	z_models = [0.001, 0.002, 0.008, 0.014]
	
	'''
	interp_func = []
	for i in range(len(z_models)):
		if z == z_models[i]:
			return f_models[i]
		else:	
			for j in range(i+1):
				if z_models[i] < z < z_models[j]:
					perc_z = (z - z_models[i])/(z_models[j] - z_models[i])
					interp_func = f_models[j]*perc_z
					print(interp_func)
					return interp_func

	'''
	for i in range(len(z_models)):
		if z == z_models[i]:
			#print('z is exact' + str(z))
			return f_models[i]
			
		elif 0.001 <= z <= 0.002:
			#print('z is between 0.001 and 0.002: z=' + str(z))
			perc_z = (z - 0.001)/(0.002 - 0.001)
			interp_func1 = (f_models[1] - f_models[0])*perc_z + f_models[0]
			return interp_func1
			
		elif 0.002 <= z <= 0.008:
			#print('z is between 0.002 and 0.008: z=' + str(z))
			perc_z = (z - 0.002)/(0.008 - 0.002)
			interp_func2 = (f_models[2] - f_models[1])*perc_z + f_models[1]
			return interp_func2
		
		elif 0.008 <= z <= 0.014:
			#print('z is between 0.008 and 0.014: z=' + str(z))
			perc_z = (z - 0.008)/(0.014 - 0.008)
			interp_func3 = (f_models[3] - f_models[2])*perc_z + f_models[2]
			return interp_func3
			
def interpolate2(model_data):
	
	f_models = []
	for datalist in model_data:
		f_models.append(datalist[2])
	
	z_models = [0.001, 0.002, 0.008, 0.014]
	
	'''
	interp_func = []
	for i in range(len(z_models)):
		if z == z_models[i]:
			return f_models[i]
		else:	
			for j in range(i+1):
				if z_models[i] < z < z_models[j]:
					perc_z = (z - z_models[i])/(z_models[j] - z_models[i])
					interp_func = f_models[j]*perc_z
					print(interp_func)
					return interp_func

	'''
	for i in range(len(z_models)):
		if z == z_models[i]:
			#print('z is exact' + str(z))
			return f_models[i]
			
		elif 0.001 <= z <= 0.002:
			#print('z is between 0.001 and 0.002: z=' + str(z))
			perc_z = (z - 0.001)/(0.002 - 0.001)
			interp_func1 = (f_models[1] - f_models[0])*perc_z + f_models[0]
			plt.plot(model_data[0][0], f_models[0])
			plt.plot(model_data[0][0], f_models[1])
			plt.plot(model_data[0][0], interp_func1)
			plt.show()
			return interp_func1
			
		elif 0.002 <= z <= 0.008:
			#print('z is between 0.002 and 0.008: z=' + str(z))
			perc_z = (z - 0.002)/(0.008 - 0.002)
			interp_func2 = (f_models[2] - f_models[1])*perc_z + f_models[1]
			plt.plot(model_data[0][0], f_models[1])
			plt.plot(model_data[0][0], f_models[2])
			plt.plot(model_data[0][0], interp_func2)
			plt.show()
			return interp_func2
		
		elif 0.008 <= z <= 0.014:
			#print('z is between 0.008 and 0.014: z=' + str(z))
			perc_z = (z - 0.008)/(0.014 - 0.008)
			interp_func3 = (f_models[3] - f_models[2])*perc_z + f_models[2]
			plt.plot(model_data[0][0], f_models[2])
			plt.plot(model_data[0][0], f_models[3])
			plt.plot(model_data[0][0], interp_func3)
			plt.show()
			return interp_func3

	
def fluxmod2(model_data, rest_wl, theta):
	
	rest_wl_m = rest_wl*1.e4
	A_v, d, B, z = theta
	R_mod = 4.05/(5.05*(0.8**d) - 4.05)
	D = (B*(rest_wl_m**2)*(350**2))/((rest_wl_m**2 - (2175**2))**2 + (rest_wl_m**2)*(350**2))
	k_calzetti = 2.659*(-2.156 + 1.509*rest_wl**(-1) - 0.198*rest_wl**(-2) + 0.011*rest_wl**(-3)) + 4.05
	A_mod = (k_calzetti*(R_mod/4.05)*((rest_wl_m/5500)**d) + D)/R_mod
	A_lambda = A_mod*A_v

	f_z = interpolate(model_data, z)

	flux_mod = f_z*10**(-0.4*A_lambda)


	return flux_mod
	
	
def normalisingfactor2(flux_mod, flux, sigma):
	
	a_tot = []
	b_tot = []
	

	for j in range(len(flux)):
		a = (flux_mod[j]*flux[j])/(sigma[j])**2
		b = (flux_mod[j]/sigma[j])**2
		a_tot.append(a)
		b_tot.append(b)
	
	a_total = np.sum(a_tot)
	b_total = np.sum(b_tot)

	B = a_total/b_total

	return B

def normalisingfactor(flux_mod, flux, sigma):
	a_total = [0.0]*len(flux_mod)
	b_total = [0.0]*len(flux_mod)
	
	for i in range(len(flux_mod)):
		for j in range(len(flux)):
			a = (flux_mod[i][j]*flux[j])/(sigma[j])**2
			b = (flux_mod[i][j]/sigma[j])**2
			a_total[i] += a
			b_total[i] += b
	
	B = [0.0]*len(flux_mod)
	
	for i in range(len(a_total)):
		B[i] = a_total[i]/b_total[i] 
	
	return B

def chisquared(flux_mod, f, sigma, B):
		
	chi = [0.0]*len(flux_mod)

	for i in range(len(flux_mod)):
		for j in range(len(f)):
			chi_sq = ((B[i]*flux_mod[i][j] - f[j])/sigma[j])**2
			chi[i] += chi_sq
	
	return chi
	

def chisquared2(flux_mod, f, sigma, B):
		
	chi = 0.0

	for i in range(len(flux_mod)):
		chi_sq = ((B*flux_mod[i] - f[i])/sigma[i])**2
		chi += chi_sq
	
	return chi


def compute_dof(theta, maskedwavelength):
    return len(maskedwavelength) - len(theta)

def chi2_likelihood(flux_mod, f, sigma, B):
    chi2 = min(chisquared(flux_mod, f, sigma, B))
    #dof = compute_dof(theta, maskedwavelength)
    dof = int(624)
    return stats.chi2(dof).pdf(chi2)
