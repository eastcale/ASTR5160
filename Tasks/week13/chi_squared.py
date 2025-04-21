#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from Modules.plots import save_show, ticks2
from Tasks.week10.classif import line
import itertools
from scipy.stats import chi2

def chisq(o, sig, e):
	"""
    Calculates chi-squared for a given array of model values,
	based on observed values and variances.

    Parameters
    ----------
    o : array or list of floats
    	The average dependent (y-axis) values for
    	every independent (x-axis) value in a set of data.

    sig : array or list of floats
    	The variance of each dependent (y-axis) value for
    	every independent (x-axis) value in a set of data.

    e : array or list of lists of floats
		Predicted dependent values (y-axis) for every 
		independent (x-axis) value based on a model. If multiple
		models are used, separate the values into individual lists.
		Each list should be the same length (corresponding to the number
		of independent values in a data set). i.e., if there are three models then:
		e = [[model1 values], [model2 values], [model3 values]]

    Returns
    -------
    chisq : list of floats
        A chi-squared value for every model given.
    """    

	chisq = [np.sum([((j - k) ** 2)/(l ** 2) for j, k, l in zip(o, i, sig)]) for i in e]

	return chisq

if __name__ == '__main__':
	#CTE Reading in data
	dat = np.loadtxt('/d/scratch/ASTR5160/week13/line.data')

	#CTE Getting means and variances for each bin
	means = [np.mean(i) for i in [dat[:, j] for j in range(len(dat[0]))]]
	variances = [np.var(i, ddof = 1) for i in [dat[:,j] for j in range(len(dat[0]))]]

	#CTE Making arrays for my grid of m and b guesses
	ms = np.arange(1, 4, 0.5)
	bs = np.arange(4, 7, 0.5)

	#CTE itertools.product() is going to find every combination of
	#CTE ms and bs to build a grid (36 combinations in my case)
	mbs = list(itertools.product(ms, bs))

	#CTE Making an x array for plotting (placed at the center of each bin)
	x_dat = np.arange(0.5, 10.5, 1)

	#CTE Generating a list of 36 lists that have the y-values for every (m, b) combination
	y_expect = np.array([line(x_dat, i[0], i[1]) for i in mbs])

	#CTE Calculating chisq for each of my (m, b) combinations
	chisq = chisq(means, variances, y_expect)

	#CTE Finding the best fit (m, b) parameters
	best_fit = mbs[np.argmin(chisq)]

	#CTE Plotting chi-squared for each (m, b) combination
	fig, ax = plt.subplots(figsize = (10, 9))

	ax.scatter(np.arange(0, 36, 1), chisq)

	ax.set_xticks(np.arange(0, 36, 1))
	ax.set_xticklabels([str(i) for i in mbs], fontsize = 10)

	ax.set_xlabel('(m, b)', fontsize = 15)
	ax.set_ylabel('$\chi^2$', fontsize = 15)

	ax.tick_params(labelrotation = 270)

	save_show(fig, 'chisq.png')

	#CTE Plotting all of my lines
	#CTE Highlighitng the line with the lowest chi-squared
	fig, ax = plt.subplots()

	[ax.scatter(len(dat[:,i])*[x_dat[i]], dat[:,i], color='navy') for i in range(len(dat[0]))]
	[ax.plot(x_dat, i, alpha=0.25) for i in y_expect]
	ax.plot(x_dat, y_expect[np.argmin(chisq)], label=f'Best Fit : (m, b) = {best_fit}', color = 'red')

	ticks2(ax, [0, 10, 1, 0.1], [0, 40, 5, 1])

	ax.legend(loc = 'best')

	save_show(fig, 'lines.png')