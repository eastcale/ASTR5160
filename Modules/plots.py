import os
import matplotlib.pyplot as plt
from PIL import Image as im
import numpy as np
import math

if __name__ == "__main__":
    print("This is meant to be used in conjunction with other code to customize plots. Please import this module from a different file.")

def min_max(x, err):
    """
    Given an array, finds ideal minimum and maximum integer 
    values that include all points in the array
	
    Parameters
    ----------
    x : array of floats
        The array this function will operate on.
        
    err : str or array of floats
        The error, of any, on the array. If no error, pass 'none'
	
    Returns
    -------
    xmin : int
        The ideal minimum value that includes np.min(x)
    	
    xmax : int
        The ideal maximum value that includes np.max(x)
    """
    
    if err == 'none':
        errs = np.size(x)*[0]
    else:
    	errs = err
    	
    xmin = math.floor(np.min(x)-np.abs(errs[0]))
    xmax = math.ceil(np.max(x)+np.abs(errs[-1]))

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
        The error, if any, on the points in x. If no error, pass 'none'.
	
    Returns
    -------
    maj_ticks : array of floats
        An array of <= 11 evenly spaced tick marks for plotting 
        (intended to be major tick marks)
    		
    min_ticks : array of floats
        An array of <= 51 evenly spaced tick marks for plotting
        (intended to be minor tick marks)
    """
    
    if err == 'none':
        errs = np.size(x)*[0]
    else:
        errs = err
    
    xmin, xmax = min_max(x, errs) #CTE Using min_max() to find ideal upper/lower bounds
    	
    rang = xmax - xmin #CTE Finding the range
	
    t_size = (rang / 10) #CTE Dividing the range by 10 
	                 #CTE to set first instance of tick spacing       
	
    #CTE Comparing first instance of tick spacing to typical tick sizes and determines the nearest 
    #CTE and smallest typical tick size.
    if t_size <= 0.01:
        tick_size = 0.01
    elif t_size <= 0.1:
        tick_size = 0.1
    elif t_size <= 0.25:
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
	
def ticks(ax, x, y, xerr='none', yerr='none'):
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
	The error, if any, on the x axis data. If no error, pass 'none'.
		
    yerr : array of floats
	The error, if any, on the y axis data. If no error, pass 'none'.
 
    Returns
    -------
    None
    """
    if xerr == 'none':
        xerrs = np.size(x)*[0]
    else:
        xerrs = xerr
   
    if yerr == 'none':
        yerrs = np.size(x)*[0]
    else:
        yerrs = yerr
   
    xticksmaj, xticksmin = tick_array(x, xerrs) #CTE Uses tick_array() to determine major/minor ticks for x axis
    yticksmaj, yticksmin = tick_array(y, yerrs) #CTE Uses tick_array() to determine major/minor ticks for y axis
	
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
	
def ticks2(ax, x_ticks, y_ticks):
    """
    Given an axis, set ticks in my (CTE) prefered style with given limits
    and increments
    
    Parameters
    ----------
    ax : type?
        The axis on which to act.
        
    xticks : list of floats or list of ints
        The desired x axis parameters in the format:
        [minimum, maximum, major_ticksize, minor_ticksize]

    yticks : list of floats or list of ints
        The desired y axis parameters in the format:
        [minimum, maximum, major_ticksize, minor_ticksize]
 
    Returns
    -------
    None
    """
    ax.set_xlim(x_ticks[0], x_ticks[1])
    ax.set_ylim(y_ticks[0], y_ticks[1])

    xmaj = np.arange(x_ticks[0], x_ticks[1]+x_ticks[2], x_ticks[2])
    xmin = np.arange(x_ticks[0], x_ticks[1]+x_ticks[3], x_ticks[3])

    ymaj = np.arange(y_ticks[0], y_ticks[1]+y_ticks[2], y_ticks[2])
    ymin = np.arange(y_ticks[0], y_ticks[1]+y_ticks[3], y_ticks[3])

    ax.set_xticks(xmaj)
    ax.set_xticks(xmin, minor=True)

    ax.set_yticks(ymaj)
    ax.set_yticks(ymin, minor=True)

    ax.tick_params(axis='both', which='major', direction='in', length=5)
    ax.tick_params(axis='both', which='minor', direction='in', length=3)

    axy=ax.twinx()
    axx=ax.twiny()

    axx.set_xlim(x_ticks[0], x_ticks[1])
    axy.set_ylim(y_ticks[0], y_ticks[1])
    
    axx.set_xticks(xmaj)
    axx.set_xticks(xmin, minor=True)

    axx.set_xticklabels('')
        
    axy.set_yticks(ymaj)
    axy.set_yticks(ymin, minor=True)

    axy.set_yticklabels('')

    axx.tick_params(axis='both', which='major', direction='in', length=5)
    axx.tick_params(axis='both', which='minor', direction='in', length=3)

    axy.tick_params(axis='both', which='major', direction='in', length=5)
    axy.tick_params(axis='both', which='minor', direction='in', length=3)
    

def save_show(fig, file_name, both = 'n'):
    """
    Intended to be used after a figure has been made. Given a file name, 
    checks if the file exists and deletes it if so. Saves a new one and
    displays it.
    
    Parameters
    ----------
    file_name : string
        The desired name of the file to be saved.
        The only accepted formats are currently .png or .svg files.
        
    Returns
    -------
    None
    """
    
    if file_name.find('.png') != -1:
        file = file_name
        fig.savefig(file)
        im.open(file).show()
    elif file_name.find('.png') == -1:
        file = file_name.replace('.svg', '.png')
        fig.savefig(file_name)
        fig.savefig(file)
        im.open(file).show()
        if both == 'n':
            os.remove(file)
