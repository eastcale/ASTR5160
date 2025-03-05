import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit as cf
import os
from PIL import Image as im
import math

def line(x, m, b):
	"""
	Given x, m, and b, return y according to y = m*x+b (a line)
	
	Parameters
	----------
	x : float or array of floats
		A value or values for x (the independent variable) 
		to generate the result, y (the dependent variable)
	
	m : float
		A value for the slope of the line
	
	b : float
		A value for the y-intercept of the line
	
	Returns
    	-------
    	y : float or array of floats
    		A value or values of y generated according to y=m*x+b
	"""
	y = m*x+b
	
	return y

def scatter(m, b):
	"""
	Given m, and b, return y (determined by line()) scattered in the +/- y axis by a magnitude given by
	a random value pull from a Gaussian centered at 0 and a standard deviation of 0.5
	
	Parameters
	----------
	m : float
		A value for the slope of the line given to the function line()
	
	b : float
		A value for the y-intercept of the line given to the function line()
	
	Returns
    	-------
    	x : array of floats
    		The independent variable in the function y = m*x+b
    	
    	y : array of floats
    		The dependent variable in the function y = m*x+b
    	
    	err : array of floats
    		The amount by which the original y values from line() were scattered
	"""
	x = np.random.uniform(0, 10, 10)
	
	y = line(x, m, b)
	
	err = np.random.normal(0, 0.5, np.size(y))
	
	return x, y+err, err
	
def min_max(x, err):
	"""
	Given an array, finds ideal minimum and maximum integer 
	values that include all points in the array
	
	Parameters
	----------
	x : array of floats
		The array this function will operate on.
	
	Returns
    	-------
    	xmin : int
    		The ideal minimum value that includes np.min(x)
    	
    	xmax : int
    		The ideal maximum value that includes np.max(x)
	"""
	xmin = math.floor(np.min(x)-np.abs(err[0]))
	xmax = math.ceil(np.max(x)+np.abs(err[-1]))
	
	return xmin, xmax


def tick_array(x, err):
	"""
	Given an array and the error on its points, generates arrays to determine ideal
	major and minor tick labels for the axis the array will be plotted on.
	
	Parameters
	----------
	x : array of floats
		The array this function will operate on.
		
	err : array of floats
		The error, if any, on the points in x.
	
	Returns
    	-------
    	maj_ticks : array of floats
    		An array of <= 11 evenly spaced tick marks for plotting 
    		(intended to be major tick marks)
    		
    	min_ticks : array of floats
    		An array of <= 51 evenly spaced tick marks for plotting
    		(intended to be minor tick marks)
	"""
	xmin, xmax = min_max(x, err) #CTE Using min_max() to find ideal upper/lower bounds
		
	rang = xmax - xmin #CTE Finding the range
	
	t_size = math.ceil(rang / 10) #CTE Dividing the range by 10 and rounding to the next integer 
				      #CTE to set first instance of tick spacing
	
	#CTE Comparing first instance of tick spacing to typical tick sizes and determines the nearest 
	#CTE and smallest typical tick size.
	if t_size <= 0.25:
		tick_size = 0.25
	elif t_size <= 0.5:
		tick_size = 0.5
	elif t_size <= 1:
		tick_size = 1
	elif t_size <= 2:
		tick_size = 2
	elif t_size <= 5:
		tick_size = 5
	elif t_size <= 10:
		tick_size = 10
	elif t_size <= 20:
		tick_size = 20
	elif t_size <= 50:
		tick_size = 50

	xmin_new = tick_size * math.floor(xmin / tick_size) #CTE Adjusting lower bounds to include new tick size
	xmax_new = tick_size * math.ceil(xmax / tick_size) #CTE Adjusting upper bounds to include new tick size
	
	maj_ticks = np.arange(xmin_new, xmax_new+tick_size, tick_size) #CTE Generating an array of tick locations
								       #CTE determined by the tick size
	
	min_t_size = (maj_ticks[1] - maj_ticks[0]) / 5 #CTE Dividing the major tick spacing into 5 smaller increments
	
	min_ticks = np.arange(xmin_new, xmax_new+min_t_size, min_t_size) #CTE Generating an array of tick locations
									 #CTE determined by min_t_size
									 
	return maj_ticks, min_ticks
	
def ticks(ax, x, y, xerr, yerr):
	"""
	Given an axis, adjust x and y limits and ticks to fit the min/max of each parameter
	and customizes them according to a specific style [the one I, CTE, prefer :)].
	
	Parameters
	----------
	ax : type?
		The axes on which to act.
		
	x : array of floats
		The array plotted on the x axis.
		
	y : array of floats
		The array plotted on the y axis.
		
	xerr : array of floats
		The error, if any, on the x axis data.
		
	yerr : array of floats
		The error, if any, on the y axis data.
 
	Returns
    	-------
    	None
	"""
	xticksmaj, xticksmin = tick_array(x, xerr) #CTE Uses tick_array() to determine major/minor ticks for x axis
	yticksmaj, yticksmin = tick_array(y, yerr) #CTE Uses tick_array() to determine major/minor ticks for y axis
	
	ax.set_xlim(min(xticksmaj), max(xticksmaj))

	ax.set_xticks(xticksmaj)
	ax.set_xticks(xticksmin, minor=True)

	ax.set_ylim(min(yticksmaj), max(yticksmaj))

	ax.set_yticks(yticksmaj)
	ax.set_yticks(yticksmin, minor=True)

	ax.tick_params(axis='both', which='major', direction='in', length=5)
	ax.tick_params(axis='both', which='minor', direction='in', length=3)

	axy=ax.twinx()
	axx=ax.twiny()

	axx.set_xlim(min(xticksmaj), max(xticksmaj))

	axx.set_xticks(xticksmaj)
	axx.set_xticks(xticksmin, minor=True)

	axx.set_xticklabels('')

	axy.set_ylim(min(yticksmaj), max(yticksmaj))
		
	axy.set_yticks(yticksmaj)
	axy.set_yticks(yticksmin, minor=True)

	axy.set_yticklabels('')

	axx.tick_params(axis='both', which='major', direction='in', length=5)
	axx.tick_params(axis='both', which='minor', direction='in', length=3)

	axy.tick_params(axis='both', which='major', direction='in', length=5)
	axy.tick_params(axis='both', which='minor', direction='in', length=3)
	
def save_show(fig, file_name):
	"""
	Intended to be used after a figure has been made. Given a file name, 
	checks if the file exists and deletes it if so. Saves a new one and
	displays it.
	
	Parameters
	----------
	fig : matplotlib.pyplot figure
		The figure to save

	file_name : string
		The desired name of the file to be saved.
		
	Returns
    	-------
    	None
	"""
	if os.path.exists(file_name): #CTE If the resulting image exists, delete it and save a new one in the pwd
		os.remove(file_name)
		fig.savefig(file_name)
	else:
		fig.savefig(file_name) #CTE If the image does not exist, save a new one in the pwd
		
	im.open(file_name).show() #CTE Show the saved image using any viewer available

if __name__ == '__main__':

	try:
		m = float(input("Slope: ")) #CTE Asks user for slope input
		b = float(input("y-intercept: ")) #CTE Asks user for y-intercept input
	except ValueError:
		raise ValueError('Expected a number. Please input a number.')

	x, y, err = scatter(m, b) #CTE Runs scatter() with the m/b user values
	yerr = np.array(np.size(err)*[0.5]) #CTE Explicitly setting the y error to be 0.5 for each point (as specified in HW0 doc)
	xerr = np.array(np.size(err)*[0]) #CTE Explicitly setting x error to be 0 for each point 

	params, cov=cf(line, x, y) #CTE Generating best fit parameters from scipy.optimize.curve_fit()

	fig, ax=plt.subplots()

	ax.plot(x, line(x, m, b), color='black', label='Original: m = {:.3f}; b = {:.3f}'.format(m, b), zorder=1)
	ax.plot(x, line(x, *params), color='red', label='Model: m = {:.3f}; b = {:.3f}'.format(*params), zorder=2)
	ax.scatter(x, y, color='navy', label='Data', zorder=3)
	ax.errorbar(x, y, yerr=np.abs(yerr), linestyle='', color='navy', zorder=4)

	ticks(ax, x, y, xerr, yerr) #CTE Adjusting the axes using ticks()

	ax.legend(loc='best')

	save_show(fig, 'line-plot.png') #CTE Saving and displaying the result of most recent plot.
