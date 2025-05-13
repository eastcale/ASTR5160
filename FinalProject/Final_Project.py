#!/usr/bin/env python3

import numpy as np
from Tasks.week10.classif import line
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import emcee
import corner
from astropy.io import fits

def quadratic(x, a0, a1, a2):
	"""
	Given x-data and 2nd order parameters, returns
	y-values for a quadratic curve.

	Parameters
	----------
	x : array of floats
		The x-data to get values for

	a0 : float
		The 0th order term in the quadratic fit

	a1 : float
		The 1st order term in the quadratic fit

	a2 : float
		The 2nd order term in the quadratic fit

	Returns
	------
	y : array of floats
		y-values for each x-value given based on a quadratic
		function.
	"""

	y = a2 * x **2 + a1 * x + a0

	return y

def ln(params, mod, x, y, yerr):
    """
    Finds the likelihood of either a linear or quadratic model with 
    given parameters relative to given data.

    Parameters
    ----------
    mod : str
    	A string specifying the desired model to use.
    	Options are 'lin' and 'quad' (for linear and quadratic).

    params : array or list of floats
    	If linear: The slope (m) and y-intercept (b) for the proposed model.
        Must be in the form [m, b].

        If quadratic: The 3 fit paramters a0, a1, and a2
        (for y = a2*x**2 + a1*x + a_0). Must be in the form
        [a0, a1, a2].

    x : array or list of floats
        The x-data to compare the model to.

    y : array or list of floats
        The y-data to compare the model to.

    yerr : array or list of floats
        The error on the given y-data.

    Returns
    -------
    ln : float
        The likelihood of the model relative to the data.
    """

    mod = str(mod)

    if mod == 'lin':
    	m, b = params
    	model = m * x + b
    elif mod == 'quad':
    	a0, a1, a2 = params
    	model = a2 * x**2 + a1 * x + a0
    else:
    	raise ValueError("Invalid mod input. Please use either 'lin' or 'quad'.")

    ln = -0.5 * np.sum((y - model) ** 2 / (yerr ** 2) + np.log(yerr ** 2))

    return ln

def prior(mod, params, param_range):
    """
    Finds the flat (log) prior for a model given its expected range
    of parameters.

    Parameters
    ----------
    mod : str
        A string specifying the desired model to use.
        Options are 'lin' or 'quad'.

    params : array or list of floats
        If linear:
            The slope (m) and y-intercept (b) for the proposed model.
            Must be in the form [m, b].
        If quadratic:
            The fit parameters a0, a1, and a2 for the function
            a2 * x**2 + a1 * x + a0.
            Must be in the form [a0, a1, a2]

    param_range : list of lists of floats
        The allowed range for the parameters.
        If linear:
            Must be in the form:
            [[m_lower_limit, m_upper_limit], [b_lower_limit, b_upper_limit]]
        If quadratic:
            Must be in the form
            [[a0_lower, a0_upper], [a1_lower, a1_upper], [a2_lower, a2_upper]]

    Returns
    -------
    prior : float
        The flat (log) prior for the model.
    """

    mod = str(mod)

    if mod == 'lin':
        m, b = params
        m_range, b_range = param_range

        if m_range[0] < m < m_range[1] and b_range[0] < b < b_range[1]:
            prior = 0
        else:
            prior = -np.inf
    elif mod == 'quad':
        a0, a1, a2 = params
        a0_range, a1_range, a2_range = param_range

        if (a0_range[0] < a0 < a0_range[1] and
            a1_range[0] < a1 < a1_range[1] and
            a2_range[0] < a2 < a2_range[1]):
            prior = 0
        else:
            prior = -np.inf
    else:
        raise ValueError("Invalid mod input. Please use either 'lin' or 'quad'.")

    return prior

def posterior(params, mod, param_range, x, y, yerr):
    """
    Finds the posterior probability for a given linear model
    relative to given data.

    Parameters
    ----------
    mod : str
    	A string specifying the desired model to fit.
    	Options are 'lin' or 'quad'.

    params : array or list of floats
    	If linear:
	        The slope (m) and y-intercept (b) for the proposed model.
	        Must be in the form [m, b].
	    If quadratic:
	    	The fit parameters for the function a2*x**2+a1*x+a0.
	    	Must be in the form [a0, a1, 2]

    param_range : list of lists of floats
    	The allowed range for the parameters.
	    If linea	r:
	        Must be in the form:
	        [[m_lower_limit, m_upper_limit], [b_lower_limit, b_upper_limit]]
	    If quadratic:
	    	Must be in the form
	    	[[a0_lower, a0_upper], [a1_lower, a1_upper], [a2_lower, a2_upper]]

    x : array or list of floats
        The x-data to compare the model to.

    y : array or list of floats
        The y-data to compare the model to.

    yerr : array or list of floats
        The error on the given y-data.

    Returns
    -------
    P : float
        The posterior probability for the given model
        relative to the given data.
    """

    lp = prior(mod, params, param_range)

    if not np.isfinite(lp):
    	P = -np.inf
    else:
    	P = lp + ln(params, mod, x, y, yerr)

    return P

def max_likelihood(mod, x, y, yerr, guess, seed = 42):
    """
    Applies maximum likelihood methods to find a 
    give model fit to given data.

    Parameters
    ----------
    x : array or list of floats
        The x-data to fit the model to.

    y : array or list of floats
        The y-data to fit the model to.

    yerr : array or list of floats
        The error on the given y-data.

    guess : list or array of floats
        A guess on the fit parameters to start.
       	If linear, must be in the form [m, b].
       	If quadratic, must be in the form [a0, a1, a2].

    seed : int
        Sets the seed for random number generation. Default is 42.

    Returns
    -------
    sol : array
    	An array corresponding to the given parameters.
    	If linear, it has the form [m, b]. If quadratic, it has the form
    	[a0, a1, a2].
    """

    #CTE Setting seed to allow for reproducable results
    np.random.seed(seed)

    #CTE lambda is a way to define an "anonymous"
    #CTE function that isn't previously defined
    #CTE *args sets a variable number of arguments that can be passed
    #CTE "nll" stands for Negative Log Likelihood
    nll = lambda *args: -ln(*args)

    #CTE Finding the minimum -ln, which should
    #CTE correspond to maximum |ln|
    soln = minimize(nll, guess, args=(mod, x, y, yerr))

    sol = soln.x

    return sol

def plot_res(flat_samp, labels, mod, fit, data):
	"""
	Given a flattened sample as well as best fit
	parameters to a set of data, creates a corner
	plot of the sample and plots the data with a
	best fit model (and error).

	Parameters
	----------
	flat_samp : NumPy array
		The flattened sample from an MCMC chain,
		given by the function emcee_mcmc()

	labels : list of str
		A list of strings corresponding to the parameters
		(in the same order as returned by emcee_mcmc())
		that are being fit (e.g., if mod = 'lin', then labels = ['m', 'b'])

	mod : str
		The model that was fit to the data. Options are
		'lin' or 'quad'.

	fit : list of floats
		The best fit parameters in the form [[value, lower_erorr, upper_error]].
		e.g., for a linear fit:
		[[m, m_lower_error, m_upper_error], [b, b_lower_error, b_upper_error]] 
	
	data : list of NumPy arrays
		The data and data error that was fit.
		Must be in the form [x, y, yerr]

	Returns
	------
	None
	"""

	#CTE Making the corner plot
	c_fig = corner.corner(flat_samp, labels=labels)
	c_fig.savefig(mod + '-corner.png')
	print(f"{mod + '-corner.png'} saved.")

	#CTE MAking the fit plot
	f_fig, f_ax = plt.subplots()

	f_ax.scatter(x, y, color='navy')
	f_ax.errorbar(x, y, yerr=yerr, linestyle='', color='navy')

	x_plot = np.linspace(f_ax.get_xlim()[0], f_ax.get_xlim()[1], 1000)

	if mod == 'lin':
		m, b = fit[0][0], fit[1][0]
		f_ax.plot(x_plot, line(x_plot, m, b), color='red')
		f_ax.fill_between(x_plot, 
			line(x_plot, m-fit[0][1], b-fit[1][1]), line(x_plot, m+fit[0][2], b+fit[1][2]), color='red', alpha=0.2)
	elif mod == 'quad':
		a0, a1, a2 = fit[0][0], fit[1][0], fit[2][0]
		f_ax.plot(x_plot, quadratic(x_plot, a0, a1, a2), color='red')
		f_ax.fill_between(x_plot,
			quadratic(x_plot, a0-fit[0][1], a1-fit[1][1], a2-fit[2][1]),
			quadratic(x_plot, a0+fit[0][2], a1+fit[1][2], a2+fit[2][1]), color='red', alpha=0.2)
	f_fig.savefig(mod + '-fit.png')
	print(f"{mod + '-fit.png'} saved.")
	plt.show()

def emcee_mcmc(mod, iters, x, y, yerr, guess, param_range):
	"""
	Fits a specifed model to a given data based on flat
	priors via Markov Chain Monte Carlo methods. 
	Must be provided a fit guess and an expected
	range for each parameter.

	Parameters
	----------
	mod : str
		The model to fit. Options are 'lin' or 'quad'.

	iters : int
		The number of iterations in the MCMC chain. Suggested is 5000.

	x : NumPy array
		The x-data to fit a model to.

	y : NumPy array
		The y-data to fit a model to.

	yerr : NumPy array
		The error on the y-data to fit.

	guess : list of floats
		A list of guesses for each parameter in the form
		[param1_guess, param2_guess, ...]

	param_range : list of lists of floats
		A list of possible ranges for each parameter in the form
		[[param1_min, param1_max], [param2_min, param2_max], ...]

	Returns
	------
	flat_samples : NumPy array
		A flattened array of parameters sampled in the MCMC chain.

	fit : list of floats
		The best fit parameters and their errors in the form:
		[[param1_val, param1_lower_error, param1_upper_error], ...]
	"""

	#CTE Finding max likelihood parameters
	sol = max_likelihood(mod, x, y, yerr, guess)

	if mod == 'lin':
		m_ml, b_ml = sol[0], sol[1]

		#CTE Initialize random walk in a Gaussian ball
		#CTE around the max-likelihood result
		pos = [m_ml, b_ml] + 1e-4 * np.random.randn(32, 2)

		labels = ['m', 'b']
	elif mod == 'quad':
		a0_ml, a1_ml, a2_ml = sol[0], sol[1], sol[2]
		post = [a0_ml, a1_ml, a2_ml] + 1e-4 * np.random.randn(32, 3)
		labels = ['a0', 'a1', 'a2']

	#CTE Setting the number of walkers to the
	#CTE number of random points made
	#CTE and the numper of dimensions to the number
	#CTE of fitting parameters
	nwalkers, ndim = pos.shape

	#CTE Sampling the distribution according to nwalkers and ndim
	#CTE Given the function for the posterior probabilities and the
	#CTE parameters needed for it (excluding params, as this is what is being sampled)
	sampler = emcee.EnsembleSampler(
		nwalkers, ndim, posterior, args=(mod, param_range, x, y, yerr))
	sampler.run_mcmc(pos, iters, progress = True)

	#CTE Autocorrelation time is essientially
	#CTE how many steps it takes for the chain to "forget"
	#CTE where it started
	tau = sampler.get_autocorr_time()
	flat_samples = sampler.get_chain(discard=100, thin=int(tau[0] / 2), flat = True)

	#CTE Getting fit parameters based on teh median,
	#CTE and errors based on the 16th and 84th percentiles
	if mod == 'lin':
		m_emcee = ([np.percentile(flat_samples[:,0], [16, 50, 84])[1],
			np.diff(np.percentile(flat_samples[:,0], [16, 50, 84]))[0],
			np.diff(np.percentile(flat_samples[:,0], [16, 50, 84]))[1]])

		b_emcee = ([np.percentile(flat_samples[:,1], [16, 50, 84])[1],
			np.diff(np.percentile(flat_samples[:,1], [16, 50, 84]))[0],
			np.diff(np.percentile(flat_samples[:,1], [16, 50, 84]))[1]])

		fit = [m_emcee, b_emcee]
	elif	mod == 'quad':
		a0_emcee = ([np.percentile(flat_samples[:,0], [16, 50, 84])[1],
			np.diff(np.percentile(flat_samples[:,0], [16, 50, 84]))[0],
			np.diff(np.percentile(flat_samples[:,0]), [16, 50, 84])[1]])

		a1_emcee = ([np.percentile(flat_samples[:,1], [16, 50, 84])[1],
			np.diff(np.percentile(flat_samples[:,1], [16, 50, 84]))[0],
			np.diff(np.percentile(flat_samples[:,1]), [16, 50, 84])[1]])

		a2_emcee = ([np.percentile(flat_samples[:,2], [16, 50, 84])[1],
			np.diff(np.percentile(flat_samples[:,2], [16, 50, 84]))[0],
			np.diff(np.percentile(flat_samples[:,2]), [16, 50, 84])[1]])

		fit = [a0_emcee, a1_emcee, a2_emcee]

	#CTE Plotting corner plot and fit plot
	plot_res(flat_samples, labels, mod, fit, [x, y, yerr])

	return flat_samples, fit

if __name__ == '__main__':
	import argparse

	parser = argparse.ArgumentParser(description = "This code uses the emcee package to fit\
		both a linear and quadratic model to a pre-determined set of data. No user input is needed.")

	data = fits.open('/d/scratch/ASTR5160/final/dataxy.fits')

	x = np.array(data[1].data['x'])
	y = np.array(data[1].data['y'])
	yerr = np.array(data[1].data['yerr'])

	lin_guess = [-0.75, 2.5]
	quad_guess = [7.5, -2, 0.05]

	print("\nStarting linear fit...")
	lin_samp, lin_fit = emcee_mcmc('lin', 5000, x, y, yerr, lin_guess, [[-5, -0.01], [0, 10]])
	print("\nFor a general parameter q, the best fit is reported:\nq = value [lower_error, upper_error]...")
	print(f"Best linear fit [m, b]:\nm = {lin_fit[0][0]} [{lin_fit[0][1]}, {lin_fit[0][2]}]\nb = {lin_fit[1][0]} [{lin_fit[1][1]}, {lin_fit[1][2]}]")

	print("\nStarting quadratic fit...")
	quad_samp, quad_fit = emcee_mcmc('quad', 5000, x, y, yerr, quad_guess, [[0, 10], [-5, 5], [0, 5]])
	print("\nFor a general parameter q, the best fit is reported:\nq = value [lower_error, upper_error]...")
	print(f"Best quadratic fit [a0, a1, a2]:\na0 = {quad_fit[0][0]} [{quad_fit[0][1]}, {quad_fit[0][2]}]\na1 = {quad_fit[1][0]} [{quad_fit[1][1]}, {quad_fit[1][2]}]\na2 = {quad_fit[2][0]} [{quad_fit[2][1]}, {quad_fit[2][2]}]")

	print("\nFor starters, the value of a2 (the term that makes the fit quadratic rather than linear)\n\
	is already very small (nearly zero), which already supports not using a quadratic fit.\n\
	The corner plot for the quadratic fit supports this idea. Looking at the distribution for a2\n\
	in the corner plot, the value is rather tightly constrained (ranging from approximately 0 to ~0.10)\n\
	near zero, with the most common value around 0.06. In this case, a2 makes a rather negligible\n\
	contribution to the fit and could be dropped in favor of a less complex fit.")