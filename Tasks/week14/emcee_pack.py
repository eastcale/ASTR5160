#!/usr/bin/env python3

import numpy as np
from Tasks.week10.classif import line
from Modules.plots import save_show
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import emcee
import corner

def ln_linear(params, x, y, yerr):
    """
    Finds the likelihood of a linear model with given parameters
    relative to given data.

    Parameters
    ----------
    params : array or list of floats
        The slope (m) and y-intercept (b) for the proposed model.
        Must be in the form [m, b].

    x : array or list of floats
        The x-data to compare the model to.

    y : array or list of floats
        The y-data to compare the model to.

    yerr : array or list of floats
        The error on the given y-data.

    Returns
    -------
    ln : float
        The likelihood of the linear model (m, b) relative to the data.
    """
    m, b = params
    model = m * x + b

    ln = -0.5 * np.sum((y - model) ** 2 / (yerr ** 2) + np.log(yerr ** 2))

    return ln

def prior_linear(params, param_range):
    """
    Finds the flat (log) prior for a model given its expected range
    of parameters.

    Parameters
    ----------
    params : array or list of floats
        The slope (m) and y-intercept (b) for the proposed model.
        Must be in the form [m, b].

    param_range : list of lists of floats
        The allowed range for the parameters. Must be in the form
        [[m_lower_limit, m_upper_limit], [b_lower_limit, b_upper_limit]]

    Returns
    -------
    prior : float
        The flat (log) prior for the linear model (m, b).
	"""

    m, b = params
    m_range, b_range = param_range

    if m_range[0] < m < m_range[1] and b_range[0] < b < b_range[1]:
        prior = 0.0
    else:
    	prior = -np.inf

    return prior

def posterior_linear(params, param_range, x, y, yerr):
    """
    Finds the posterior probability for a given linear model
    relative to given data.

    Parameters
    ----------
    params : array or list of floats
        The slope (m) and y-intercept (b) for the proposed model.
        Must be in the form [m, b].

    param_range : list of lists of floats
        The allowed range for the parameters. Must be in the form
        [[m_lower_limit, m_upper_limit], [b_lower_limit, b_upper_limit]]

    x : array or list of floats
        The x-data to compare the model to.

    y : array or list of floats
        The y-data to compare the model to.

    yerr : array or list of floats
        The error on the given y-data.

    Returns
    -------
    P : float
        The posterior probability for the given linear model (m, b)
        relative to the given data.
    """

    lp = prior_linear(params, param_range)

    if not np.isfinite(lp):
        P =  -np.inf
    else:
    	P = lp + ln_linear(params, x, y, yerr) 

    return P

def least_squares_linear(x, y, yerr):
    """
    Applies least squares methods to find a linear
    fit to given data.

    Parameters
    ----------
    x : array or list of floats
        The x-data to fit the model to.

    y : array or list of floats
        The y-data to fit the model to.

    yerr : array or list of floats
        The error on the given y-data.

    Returns
    -------
    mb : list of floats
        The best fit parameters (slope: m and y-intercept: b)
        in the form [m, b].

    mb_err : list of floats
        The error on the best fit parameters in the form
        [m_err, b_err].
    """

    #CTE Creating a Vandermonde Matrix
    #CTE Vandermonde Matrix:
    #      [1    x_0    x_0^2    ...    x_0^n]
    #      [1    x_1    x_1^2    ...    x_1^n]
    # V = [.                                     ]
    #      [.                                     ]
    #      [1    x_m    x_m^2    ...    x_m^n]
    #CTE In this case, m = len(x) and n = 1
    #CTE That way we can retrieve a linear fit
    A = np.vander(x, 2)

    #CTE Setting up the problem to solve:
    #CTE (A.T dot A)(x) = (y)
    #CTE Where x contains the best fit params
    #CTE Not sure why we are dividing by variance here,
    ATA = np.dot(A.T, A / (yerr**2)[:, None])

    #CTE Finding the covariance matrix
    cov = np.linalg.inv(ATA)

    #CTE Solving for the best fit parameters
    params = np.linalg.solve(ATA, np.dot(A.T, y / yerr**2))

    #CTE Defining parameters explicitly
    mb = [params[0], params[1]]
    mb_err = [np.sqrt(cov[0, 0]), np.sqrt(cov[1, 1])]

    return mb, mb_err

def max_likelihood_linear(x, y, yerr, guess, seed = 42):
    """
    Applies maximum likelihood methods to find a linear
    fit to given data.

    Parameters
    ----------
    x : array or list of floats
        The x-data to fit the model to.

    y : array or list of floats
        The y-data to fit the model to.

    yerr : array or list of floats
        The error on the given y-data.

    guess : list or array of floats
        A guess on the linear fit parameters to start.
        Must be in the form [m_guess, b_guess].

    seed : int
        Sets the seed for random number generation. Default is 42.

    Returns
    -------
    m : float
        The slope of the linear best fit.

    b : float
        The y-intercept of the linear best fit.
    """

    #CTE Setting seed to allow for reproducable results
    np.random.seed(seed)

    #CTE lambda is a way to define an "anonymous"
    #CTE function that isn't previously defined
    #CTE *args sets a variable number of arguments that can be passed
    #CTE not sure why this variable is called nll though
    nll = lambda *args: -ln_linear(*args)

    #CTE Finding the minimum -ln_linear, which should
    #CTE correspond to maximim |ln_linear|
    soln = minimize(nll, guess, args=(x, y, yerr))

    m, b = soln.x

    return m, b

def emcee_mcmc(iters, x, y, yerr, guess, param_range):
    """
    Applies Markov Chain Monte Carlo methods via the emcee
    package to find a linear fit to given data. Plots the sampled
    parameters, a corner plot of the samples, and the best fit line
    in the data space.

    Parameters
    ----------
    iters : int
        The number of iterations for the chain. Recommended is 5000.

    x : array or list of floats
        The x-data to fit the model to.

    y : array or list of floats
        The y-data to fit the model to.

    yerr : array or list of floats
        The error on the given y-data.

    guess : list or array of floats
        A guess on the linear fit parameters to start.
        Must be in the form [m_guess, b_guess].

    param_range : list of lists of floats
        The allowed range for the parameters. Must be in the form
        [[m_lower_limit, m_upper_limit], [b_lower_limit, b_upper_limit]].

    Returns
    -------
    m_emcee : list of floats
        The slope of the linear best fit and the errors based
        on the 16th, 50th, and 84th percentiles in the form:
        [m_value, m_lower_error, m_upper_error].

    b_emcee : list of floats
        The y-intercept of the linear best fit and the errors based
        on the 16th, 50th, and 84th percentiles in the form:
        [b_value, b_lower_error, b_upper_error].
    """

    #CTE Finding the maximum likelihood parameters
    m_ml, b_ml = max_likelihood_linear(x, y, yerr, guess)

    #CTE Initialize random walk in a Gaussian ball
    #CTE around the max-likelihood result
    pos = [m_ml, b_ml] + 1e-4 * np.random.randn(32, 2)

    #CTE Setting the number of walkers to the 
    #CTE number of random points made
    #CTE and the number of dimensions to the number
    #CTE of fitting parameters
    nwalkers, ndim = pos.shape

    #CTE Sampling the distribution according to nwalkers and ndim
    #CTE Given the function for the posterior probabilities and the
    #CTE parameters needed for it (excluding params, as this is what is being sampled)
    sampler = emcee.EnsembleSampler(nwalkers, ndim, posterior_linear,
                                    args=(param_range, x, y, yerr))
    sampler.run_mcmc(pos, iters)

    #CTE Showing the sample of values
    fig, axes = plt.subplots(2, figsize=(10, 7), sharex=True)
    samples = sampler.get_chain()
    labels = ["m", "b"]
    for i in range(ndim):
        ax = axes[i]
        ax.plot(samples[:, :, i], "k", alpha=0.3)
        ax.set_xlim(0, len(samples))
        ax.set_ylabel(labels[i])
        ax.yaxis.set_label_coords(-0.1, 0.5)

    axes[-1].set_xlabel("step number")

    save_show(fig, 'sample.png')

    #CTE Autocorrelation time is essentially how many steps it
    #CTE takes for the chain to "forget" where it started.
    #CTE Apparently (???) it is not unreasonable to discard
    #CTE the first 100 steps and thin by half the autocorrelation time
    tau = sampler.get_autocorr_time()
    #CTE "Flatten" the sample so every iteration has the same weight (??)
    flat_samples = sampler.get_chain(discard=100, thin=int(tau[0] / 2), flat=True)

    #CTE Making a corner plot of the results
    fig = corner.corner(flat_samples, labels=labels)
    save_show(fig, 'corner.png')

    #CTE Reporting m as the value at the 50th percentile for the sample
    m_val = np.percentile(flat_samples[:, 0], [16, 50, 84])[1]
    #CTE Reporting m_down as the difference between the 50th and 16th percentile values
    m_down = np.diff(np.percentile(flat_samples[:, 0], [16, 50, 84]))[0]
    # CTE Reporting m_up as the difference between the 84th and 50th percentile values
    m_up = np.diff(np.percentile(flat_samples[:, 0], [16, 50, 84]))[1]
    #CTE Reporting m with error in the form [m, m_down, m_up]
    m_emcee = [m_val, m_down, m_up]

    #CTE Reporting b as the value at the 50th percentile for the sample
    b_val = np.percentile(flat_samples[:, 1], [16, 50, 84])[1]
    #CTE Reporting b_down as the difference between the 50th and 16th percentile values
    b_down = np.diff(np.percentile(flat_samples[:, 1], [16, 50, 84]))[0]
    # CTE Reporting b_up as the difference between the 84th and 50th percentile values
    b_up = np.diff(np.percentile(flat_samples[:, 1], [16, 50, 84]))[1]
    #CTE Reporting b with error in the form [b, b_down, b_up]
    b_emcee = [b_val, b_down, b_up]

    #CTE Showing the best fit in the data space.
    fig, ax = plt.subplots()

    ax.scatter(x, y, color='navy')
    ax.errorbar(x, y, yerr, color='navy', linestyle='')

    ax.plot(np.linspace(min(x), max(x), 1000), line(np.linspace(min(x), max(x), 1000),
                                                     m_emcee[0], b_emcee[0]), color='red', linestyle='--')

    ax.fill_between(np.linspace(min(x), max(x), 1000), 
                    line(np.linspace(min(x), max(x), 1000), m_emcee[0] - m_emcee[1], b_emcee[0] - b_emcee[1]),
                    line(np.linspace(min(x), max(x), 1000), m_emcee[0] + m_emcee[2], b_emcee[0] + b_emcee[2]),
                    color='red', alpha=0.25)

    ax.set_xlabel('x', fontsize=15)
    ax.set_ylabel('y', fontsize=15)

    save_show(fig, 'mcmc-bestfit.png')

    return m_emcee, b_emcee

if __name__ == '__main__':
	print('This code is adapted from the emcee package documentation')

	#CTE Reading in data and getting relevant statistics
	dat = np.loadtxt('/d/scratch/ASTR5160/week13/line.data')

	means = np.array([np.mean(i) for i in [dat[:, j] for j in range(len(dat[0]))]])
	variances = np.array([np.var(i, ddof = 1) for i in [dat[:,j] for j in range(len(dat[0]))]])

	#CTE Setting data for algorithms to work with
	x = np.arange(0.5, 10.5, 1)
	y = means
	yerr = np.sqrt(variances)

	#CTE Defining model x-values for plotting
	x_model = np.linspace(0.5, 9.5, 1000)

	#CTE Least Squares (LS) Fitting
	mb_ls, mb_ls_err = least_squares_linear(x, y, yerr)

	#CTE MCMC Fitting
	m_mcmc, b_mcmc = emcee_mcmc(5000, x, y, yerr, [4, 5], [[1,5],[0,10]])

	#CTE Plotting the difference between LS and MCMC.
	#CTE Thought they would differ more but they are almost exactly the same
	fig, ax = plt.subplots()

	ax.scatter(x, y, color='navy', zorder=10)
	ax.errorbar(x, y, yerr, color='navy', zorder=10, linestyle='')

	ax.plot(x_model, line(x_model, mb_ls[0], mb_ls[1]), color='orange', label = 'LS')
	ax.fill_between(x_model, 
	                line(x_model, mb_ls[0]+mb_ls_err[0], mb_ls[1]+mb_ls_err[1]), 
	                line(x_model, mb_ls[0]-mb_ls_err[0], mb_ls[1]-mb_ls_err[1]), 
	                color = 'orange', alpha = 0.25)

	ax.plot(x_model, line(x_model,
			m_mcmc[0], b_mcmc[0]), color = 'green', linestyle = '--', label = 'MCMC')

	ax.fill_between (x_model, line(np.linspace(min(x), max(x), 1000),
			m_mcmc[0]-m_mcmc[1], b_mcmc[0]-b_mcmc[1]), line(np.linspace(min(x), max(x), 1000),
			m_mcmc[0]+m_mcmc[2], b_mcmc[0]+b_mcmc[2]), color = 'red', alpha = 0.25)

	ax.legend(loc = 'best')

	ax.set_xlabel('x', fontsize = 15)
	ax.set_ylabel('y', fontsize = 15)

	save_show(fig, 'fits.png')