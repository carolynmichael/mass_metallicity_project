import numpy as np
import emcee
from scipy import integrate
import chisq as cs
import scipy.optimize as op
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

def log_prior(theta):
	# size of theta determines the model.
	# flat prior over a large range
	A_v, d, B = theta
	if 0.0 < A_v < 6.0 and -0.5 < d < 0.5 and 0.0 < B < 1.0:
		return 0.0
	else:
		return -np.inf


def log_likelihood(theta, data, f_int_mask):
	x, y, sigma_y = data
	y_M = cs.fluxmod2(f_int_mask, x/1.e4, theta)
	B = cs.normalisingfactor2(y_M, y, sigma_y)
	yM = y_M*B
	return -0.5 * np.sum(np.log(2 * np.pi * sigma_y ** 2) + (y - yM) ** 2 / sigma_y ** 2)

def log_posterior(theta, data, f_int_mask):
	theta = np.asarray(theta)
	'''
	if not np.isfinite(log_prior(theta)):
		return -np.inf
	'''
	return log_prior(theta) + log_likelihood(theta, data, f_int_mask)
	
def maxlikelihood(true_vals, data, theta, f_int_mask):
	A_v_true, d_true, B_true = true_vals
	x, y, yerr = data
	nll = lambda *args: -log_likelihood(*args)
	result = op.minimize(nll, [A_v_true, d_true, B_true], args=(x, y, yerr))
	Av_ml, d_ml, B_ml = result["x"]
	return result["x"]

def compute_mcmc(data, f_int_mask, true_vals, log_posterior=log_posterior, nwalkers=50, nburn=1000, nsteps=2000):
	ndim = 3.0 # this determines the model
	rng = np.random.RandomState(0) 
	starting_guesses = []
	for i in range(nwalkers):
		starting_guess = [np.random.uniform(low=0.0, high=6.0), np.random.uniform(low=-0.5, high=0.5), np.random.uniform(low=0.0, high=1.0)]
		starting_guesses.append(starting_guess)
	#result["x"] = maxlikelihood(true_vals, data, theta, f_int_mask)
	#pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
	sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=[data, f_int_mask])
	sampler.run_mcmc(np.array(starting_guesses), nsteps)
	trace = sampler.chain[:, nburn:, :].reshape(-1, ndim)
	return trace

def integrate_posterior_3D(log_posterior, xlim, ylim, zlim, data, f_int_mask, theta):
	func = lambda theta2, theta1, theta0: np.exp(log_posterior([theta0, theta1, theta2], data, f_int_mask))
	print(func(theta[2], theta[1], theta[0]))
	print(log_posterior([theta[0], theta[1], theta[2]], data, f_int_mask))
	return integrate.tplquad(func, xlim[0], xlim[1], lambda x: ylim[0], lambda x: ylim[1], lambda x, y: zlim[0], lambda x, y: zlim[1])
