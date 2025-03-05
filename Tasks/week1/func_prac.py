#CTE Importing necessary packages
import numpy as np
import os
import matplotlib.pyplot as plt
from PIL import Image as im

def quad(x):
	"""
	Given x, return y according to y=x**2+3*x+8
	
	Parameters
	----------
	x : float or array of floats
		A value or values for x to generate the result, y
	
	Returns
    	-------
    	float
    		A value of values of y generated according to y=x**2+3*x+8
	"""
	y = x**2+3*x+8
	
	return y
	
def quad_plot(low_lim, up_lim):
	"""
	Given upper and lower limits, generate a plot of y(x) according to the func quad(x)
	
	Parameters
	----------
	low_lim : float
		The lower limit on the x axis for the plot
	
	up_lim : float
		The upper limit on the x axis for the plot
	
	Returns
    	-------
    	Image
    		Creates and displays an image named 'quad-plot.png' in the pwd
	"""
	
	x = np.linspace(low_lim, up_lim, 1000) #CTE Generating x array from lower/upper limits
	
	fig, ax = plt.subplots()
	
	ax.plot(x, quad(x), color='navy')
	
	ax.set_xlim(low_lim, up_lim)
	
	if os.path.exists('quad-plot.png'): #CTE If the resulting image exists, delete it and save a new one
		os.remove('quad-plot.png')
		fig.savefig('quad-plot.png')
	else:
		fig.savefig('quad-plot.png') #CTE If the image does not exist, save a new one
		
	im.open('quad-plot.png').show() #CTE Use whatever image viewer is available to view the plot

lower_lim=float(input("Lower Limit: ")) #CTE Requests user input for the lower limit of the graph
upper_lim=float(input("Upper Limit: ")) #CTE Requests user input for the upper limit of the graph

quad_plot(lower_lim, upper_lim)
