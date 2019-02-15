import numpy as np

def attenuation(rest_wl):
	A_cullen = (0.587*rest_wl**(-1) - 0.020*rest_wl**(-2))
	A_calzetti = (2.659*(-2.156 + 1.509*rest_wl**(-1) - 0.198*rest_wl**(-2) + 0.011*rest_wl**(-3)) + 4.05)/4.05
	return A_cullen, A_calzetti


def attenuation_mod(rest_wl, d, B):
	

	R_mod = 4.05/(5.05*(0.8**d) - 4.05)
	#convert list of wavelengths into microns
	rest_wl_m = rest_wl*1.e4
	D = (B*(rest_wl_m**2)*(350**2))/((rest_wl_m**2 - (2175**2))**2 + (rest_wl_m**2)*(350**2))
	k_calzetti = 2.659*(-2.156 + 1.509*rest_wl**(-1) - 0.198*rest_wl**(-2) + 0.011*rest_wl**(-3)) + 4.05
	A_mod = (k_calzetti*(R_mod/4.05)*((rest_wl_m/5500)**d) + D)/R_mod


	return A_mod


