#!/usr/bin/env python3

import numpy as np
from Tasks.week10.classif import line
from Modules.plots import save_show
import matplotlib.pyplot as plt

def ln(x_dat, y_dat, variance, y_model):
    """
    Finds the likelihood of a linear model with given parameters
    relative to given data.

    Parameters
    ----------
    x_dat : array or list of floats
        The x-data to compare the model to.

    y_dat : array or list of floats
        The y-data to compare the model to.

    variance : array or list of floats
        The variance of the given y-data.

    y_model : array or list of floats
        The y-data of the linear model to compare.

    Returns
    -------
    ln : float
        The likelihood of the linear model (m, b) relative to the data.
    """

    ln = -1/2 * np.sum([((i - j)**2)/(k**2) + np.log(2*np.pi*(k**2)) 
                        for i, j, k in zip(y_dat, y_model, variance)])

    return ln

def post(m, b, x_dat, y_dat, variance, y_model):
    """
    Finds the posterior probability given a flat prior
    for a linear model with given parameters relative 
    to given data.

    Parameters
    ----------
    m : list of floats
    	The value and upper/lower limits for the proposed slope.
    	Must be in the form: [m_value, m_lower_limit, m_upper_limit]

    b : float
    	The value and upper/lower limits for the proposed y-intercept.
    	Must be in the form: [b_value, b_lower_limit, b_upper_limit]

    x_dat : array or list of floats
        The x-data to compare the model to.

    y_dat : array or list of floats
        The y-data to compare the model to.

    variance : array or list of floats
        The variance of the given y-data.

    y_model : array or list of floats
        The y-data of the linear model to compare.

    Returns
    -------
    posterior : float
        The posterior probability of the proposed linear model.
    """

    ln_L = ln(x_dat, y_dat, variance, y_model)

    if m[1] <= m[0] <= m[-1] and b[1] <= b[0] <= b[-1]:
        posterior = ln_L
    else:
        posterior = - np.inf

    return posterior

def mcmc_run(x_data, y_data, variance, param_range, guess, prop_step, n_iters = 5000):
    """
    Runs a Metropolis-Hastings Algorithm to find a linear best fit to given data.

    Parameters
    ----------
    x_data : array or list of floats
        The x-data to compare the model to.

    y_data : array or list of floats
        The y-data to compare the model to.

    variance : array or list of floats
        The variance of the given y-data.

    param_range : list of floats
        The expected range for each parameter. Must be in the form:
        [[m_lower_limit, m_upper_limit], [b_lower_limit, b_upper_limit]]

    guess : list of floats
        An inital guess for each parameter. Must be in the form:
        [m_guess, b_guess]

    prop_step : float
        The standard deviation of the Gaussian proposal function

    n_iters : int
        The number of iterations in the chain. Default is 5000 (recommended)

    Returns
    -------
    params : list of lists of floats
        A collection of all accepted parameters in the form:
        [[m_accepted], [b_accepted]]

    accept_rate : float
        The acceptance rate (in decimal form) of the algorithm.
        Ideally, this rate should be ~0.30

    likelihood : list of floats
        The likelihood for each accepted model.
    """

    # CTE Setting box for flat prior
    m_range = param_range[0]
    b_range = param_range[1]

    # CTE Initializing parameter lists, starting with initial guesses
    ms = [guess[0]]
    bs = [guess[1]]

    # CTE Initializing posterior probability list, starting with the probability
    # CTE of the initial guess
    posteriors = [post([ms[0], m_range[0], m_range[1]], [bs[0], b_range[0], b_range[1]],
                       x_data, y_data, variance, line(x_data, ms[0], bs[0]))]

    # CTE Initializing likelihood list to append to
    likelihood = []

    # CTE Starting count of accepted parameters
    accept = 0

    # CTE Iterating n_iters number of times
    iters = 0
    while iters < n_iters:
        # CTE Random walk in the parameter space
        m = np.random.normal(ms[-1], prop_step)
        b = np.random.normal(bs[-1], prop_step)

        # CTE Finding the posterior probability of the random guess
        posterior = post([m, m_range[0], m_range[1]], [b, b_range[0], b_range[1]],
                         x_data, y_data, variance, line(x_data, m, b))

        # CTE If the posterior probability is np.nan or +/- np.inf, then
        # CTE then R = 0. Otherwise, R is the np.exp(difference of P_new and P_old)
        if not np.isfinite(posterior):
            R = 0
        else:
            R = np.exp(posterior - posteriors[-1])

        # CTE Rejection/Acceptance
        # CTE If R is greater than or equal to 1, then always accept parameters
        # CTE Otherwise, generate a random number. If R is greater than the random
        # CTE number, accept the proposal. Otherwise, reject and move on.
        if R >= 1:
            accept += 1
            ms.append(m)
            bs.append(b)
            posteriors.append(posterior)
            likelihood.append(np.exp(ln(x_data, y_data, variance, line(x_data, m, b))))
        else:
            rand = np.random.rand()

            if R > rand:
                accept += 1
                ms.append(m)
                bs.append(b)
                posteriors.append(posterior)
                likelihood.append(np.exp(ln(x_data, y_data, variance, line(x_data, m, b))))
            else:
                pass

        iters += 1

    params = [ms, bs]
    accept_rate = (accept / n_iters)

    return params, accept_rate, likelihood

if __name__ == '__main__':
	#CTE Reading in data and getting relevant statistics
	dat = np.loadtxt('/d/scratch/ASTR5160/week13/line.data')

	means = [np.mean(i) for i in [dat[:, j] for j in range(len(dat[0]))]]
	variances = [np.var(i, ddof = 1) for i in [dat[:,j] for j in range(len(dat[0]))]]

	#CTE Setting x_dat based on the bins in dat
	x_dat = np.arange(0.5, 10.5, 1)
	#CTE Setting y_dat based on the average y values in each bin
	y_dat = means

	#CTE Making the MCMC chain
	mbs, acceptance_rate, likelihoods = mcmc_run(
		x_dat, y_dat, variances, [[1, 5], [0, 10]], [3, 4.5], 0.5)

	print(f'Acceptance rate is {(acceptance_rate) * 100}%')
	print(f'Highest likelihood pair is (m, b)={mbs[0][np.argmax(likelihoods)], mbs[1][np.argmax(likelihoods)]}')