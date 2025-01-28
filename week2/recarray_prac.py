from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image as im 
import os
import math

def save_show(file_name):
	"""
	Intended to be used after a figure has been made. Given a file name, 
	checks if the file exists and deletes it if so. Saves a new one and
	displays it.
	
	Parameters
	----------
	file_name : string
		The desired name of the file to be saved.
		
	Returns
    	-------
    	None
	"""
	if os.path.exists(file_name):
		yn = input("Overwrite existing version of {} (y/n)? ".format(file_name))
		if yn == "y":
			os.remove(file_name)
			fig.savefig(file_name)
		elif yn == "n":
			print("Did not save new version of {}.".format(file_name))
	else:
		fig.savefig(file_name) #CTE If the image does not exist, save a new one in the pwd
	
	dis = input("Display {} (y/n)? ".format(file_name))
	 
	if dis == "y":
		im.open(file_name).show() #CTE Show the saved image using any viewer available
	else:
		pass
	
def fits_write(file_name, tab):
	"""
	Used to save a recarray (as a table) to a fits file. If the file already exists,
	asks user if they would like to overwrite. Does not write anything if user says no.
	
	Parameters
	----------
	file_name : string
		The desired name of the fits file to be saved.
		
	tab : astropy.table.Table
		The table to be written to a fits file
	"""
	if os.path.exists(file_name):
		yn = input("Overwrite existing version of {} (y/n)? ".format(file_name))
		if yn == "y":
			tab.write(file_name, overwrite=True)
		elif yn == "n":
			print("Did not save new version of {}.".format(file_name))
	else:
		tab.write(file_name)
		
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

objs = Table.read("/d/scratch/ASTR5160/week2/struc.fits")

ext_more = objs["EXTINCTION"][:,0] >= 0.22

ra, dec = objs["RA"], objs["DEC"]

ra_ext_more, dec_ext_more = ra[ext_more], dec[ext_more]

fig, ax=plt.subplots()

ax.plot(ra, dec, "bx")

ax.plot(ra_ext_more, dec_ext_more, "rx", label="Ext. $\geq$ 0.22")

ticks(ax, ra, dec, np.array(np.size(ra)*[0]), np.array(np.size(ra)*[0])) #CTE Using ticks() to set tick parameters;
									 #CTE the two arrays are to set the errors for each
									 #CTE axis to zero

ax.set_xlabel("Right Ascension")
ax.set_ylabel("Declination")

ax.invert_xaxis() #CTE Flipping the x axis to fit with standard RA convention

ax.legend(loc = 'best')

save_show('recarray-prac.png')

randnum1 = np.random.randint(1, 100, 100)
randnum2 = np.random.randint(1, 100, 100)
randnum3 = np.random.randint(1, 100, 100)

randnum = np.reshape(np.concatenate((randnum1, randnum2, randnum3)), (100, 3)) #CTE Stitching and reshaping the three arrays

ra_cop = np.copy(ra)
dec_cop = np.copy(dec)

tab = Table([ra_cop, dec_cop, randnum], names=["RA", "DEC", "RANDNUM"])

fits_write("recarray-prac.fits", tab)
