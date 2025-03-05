from astropy import units as u
from astropy.coordinates import SkyCoord as sc
from astropy.coordinates import EarthLocation as el
from astropy.coordinates import Angle as ang
from astropy.coordinates import AltAz as AltAz
from astropy.time import Time as ti
import os
import matplotlib.pyplot as plt
from PIL import Image as im
import numpy as np
import math
import dustmaps
from dustmaps.config import config
from dustmaps.sfd import SFDQuery

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
	
    t_size = (rang / 10) #CTE Dividing the range by 10 
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
	
dustdir = "/d/scratch/ASTR5160/data/dust/v0_1/maps"
config["data_dir"] = dustdir
sfd = SFDQuery()

ugriz = np.array([4.239, 3.303, 2.285, 1.698, 1.263]) #CTE SDSS ugriz R_v values from Schlafly & Finkbeiner 2011

q1 = sc(ra=246.933*u.deg, dec=40.795*u.deg) #CTE Coordinates for the first quasar
q2 = sc(ra=236.562*u.deg, dec=2.440*u.deg)

q1_gal = q1.galactic #CTE Converting ICRS coords to Galactic
q2_gal = q2.galactic 


q1_ugriz = np.array([18.82, 18.81, 18.73, 18.82, 18.90]) #CTE u,g,r,i,z magnitudes for the first quasar
q2_ugriz = np.array([19.37, 19.10, 18.79, 18.73, 18.63])

q1_ebv = sfd(q1_gal) #CTE E(B-V) for the first quasar
q2_ebv = sfd(q2_gal) #CTE E(B-V) for the second quasar

q1_A = q1_ebv * ugriz #CTE Extinction correction
q2_A = q2_ebv * ugriz

q1_ugriz_cor = np.array([i-j for i, j in zip(q1_ugriz, q1_A)]) #CTE Extinction corrected ugriz
q2_ugriz_cor = np.array([i-j for i, j in zip(q2_ugriz, q2_A)])

q1_gr = q1_ugriz[1] - q1_ugriz[2] #CTE g-r color without ext. correction
q1_gr_cor = q1_ugriz_cor[1] - q1_ugriz_cor[2] #CTE g-r color with ext. correction

q2_gr = q2_ugriz[1] - q2_ugriz[2]
q2_gr_cor = q2_ugriz_cor[1] - q2_ugriz_cor[2]

q1_ri = q1_ugriz[2] - q1_ugriz[3] #CTE r-i color without ext. correction
q1_ri_cor = q1_ugriz_cor[2] - q1_ugriz_cor[3] #CTE r-i color with ext. correction

q2_ri = q2_ugriz[2] - q2_ugriz[3]
q2_ri_cor = q2_ugriz_cor[2] - q2_ugriz_cor[3]

fig, ax = plt.subplots()

ax.scatter(q1_gr, q1_ri, color='r', label='Q1', alpha=0.3)
ax.scatter(q1_gr_cor, q1_ri_cor, color='r', label='Q1: Ext. Corrected', marker='*')

ax.scatter(q2_gr, q2_ri, color='orange', label='Q2', alpha=0.3)
ax.scatter(q2_gr_cor, q2_ri_cor, color='orange', label='Q2: Ext. Corrected', marker='*')

ticks(ax, [q1_gr, q1_gr_cor, q2_gr, q2_gr_cor], [q1_ri, q1_ri_cor, q2_ri, q2_ri_cor], [0, 0, 0, 0], [0, 0, 0, 0])

ax.set_xlabel('$(r-i)$', fontsize=15)
ax.set_ylabel('$(g-r)$', fontsize=15)

ax.legend(loc='best')

save_show(fig, 'quasar_cols.png')

ra_grid_q2 = np.arange(236.6-(0.1*50), 236.6+(0.1*50), 0.1) #CTE 100 RA coords. centered on 236.6 in steps of 0.1 deg
dec_grid_q2 = np.arange(2.4-(0.1*50), 2.4+(0.1*50), 0.1) #CTE 100 DEC coords. centered on 2.4 in steps of 0.1 deg

x_q2, y_q2 = np.meshgrid(ra_grid_q2, dec_grid_q2) #CTE 100x100 grid centered on ra=236.6, dec=2.4

coos_q2 = sc(ra = x_q2*u.deg, dec = y_q2*u.deg) #CTE Creating coord. objects for each point in the 100x100 grid

ebvs_q2 = sfd(coos_q2) #CTE E(B-V) for each of the coords. in the grid

ra_grid_q1 = np.arange(246.9-(0.13*50), 246.9+(0.13*50), 0.13) #CTE 100 RA coords. centered on 246.9 in steps of 0.13 deg 
dec_grid_q1 = np.arange(40.8-(0.1*50), 40.8+(0.1*50), 0.1) #CTE 100 DEC coords. centered on 40.8 in steps of 0.1 deg

x_q1, y_q1 = np.meshgrid(ra_grid_q1, dec_grid_q1) #CTE 100x100 grid centered on ra=246.9, dec=40.8

coos_q1 = sc(ra = x_q1*u.deg, dec = y_q1*u.deg) #CTE Creating coord. objects for each point in the 100x100 grid

ebvs_q1 = sfd(coos_q1) #CTE E(B-V) for each of the coords. in the grid

fig, ax = plt.subplots()

conts_q2 = ax.contourf(x_q2, y_q2, ebvs_q2, cmap='ocean')
fig.colorbar(conts_q2)

ax.scatter(q2.ra.deg, q2.dec.deg, color='orange', marker='*', s=50)

ax.invert_xaxis()

ax.set_xlabel("Right Ascension (deg)", fontsize=15)
ax.set_ylabel("Declination (deg)", fontsize=15)

save_show(fig, 'q2_contour.png')

fig, ax = plt.subplots()

conts_q1 = ax.contourf(x_q1, y_q1, ebvs_q1, cmap='ocean')
fig.colorbar(conts_q1)

ax.scatter(q1.ra.deg, q1.dec.deg, color='red', marker='*', s=50)

ax.invert_xaxis()

ax.set_xlabel("Right Ascension (deg)", fontsize=15)
ax.set_ylabel("Declination (deg)", fontsize=15)

save_show(fig, 'q1_contour.png')
