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

def to_cart(ra, dec):
    """
    Converts given Equatorial Coordinates, converts to Cartesian Coordinates
    """
    x = np.cos(ra)*np.cos(dec)
    y = np.sin(ra)*np.cos(dec)
    z = np.sin(dec)
    
    return [x, y, z]
    
def to_decdeg(d, m, s):
    """
    Given an angle in degrees, arcminutes,
    and arcseconds, converts to decimal degrees.

    Parameters
    ----------
    d : float
        The degree value for the given angle.
	
    m : float
        The arcminutes value for the given angle.

    s : float
        The arcseconds value for the given angle.

    Returns
    -------
    ddeg : float
        The resulting angle in decimal degrees.
    """
	
    ddeg = (d)+(m/60)+(s/3600)

    return ddeg
    
def time_to(current_time, time):
    """
    Given the current time,
    finds the number of hours until a desired time
    
    Parameters
    ----------
    current_time : astropy.time.Time
        The current time given by 
        astropy.time.Time.now()

    time : list of floats
        The desired time in the format 
        [hr, min, sec]

    Returns
    -------
    tot : float
        The number of hours until the desired time
	
    """
    hr_to = (time[0] - current_time.ymdhms[3])
    min_to = (time[1] - current_time.ymdhms[4])
    sec_to = (time[2] - current_time.ymdhms[5])

    tot = ((hr_to) + (min_to/60) + (sec_to/3600))*u.hour
	
    return tot

def zenith_to_icrs(o, l):
    """
    Given an observation time and a location, finds the Equatorial Coordinates
    for zenith at that time and place.
    """
   
    altaz_icrs = sc(alt = 90*u.deg, az = 0*u.deg, obstime = o, location = l, frame = 'altaz').icrs
    
    return altaz_icrs
    
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

bg_coo = sc('05 55 10.3054 +07 24 25.4304', unit=(u.hourangle, u.deg)) #CTE Betelgeuse coordinates from Simbad

cart_hand = to_cart(bg_coo.ra.rad, bg_coo.dec.rad) #CTE Important to note that the angles are in radians, not degrees

bg_coo_cart = np.copy(bg_coo) #CTE Making a copy of bg_coo so I can keep both versions.

bg_coo_cart.representation_type = "cartesian" #CTE The SkyCoord result agrees with by hand conversions.

gal_cent = sc('00 00 00 00 00 00', unit=(u.deg, u.deg), frame='galactic')

gal_cent_icrs = gal_cent.icrs #CTE The Galactic center is on the very edge of Sagittarius

wiro = el(lat = to_decdeg(41, 5, 49)*u.deg, lon = -1*to_decdeg(105, 58, 33)*u.deg, height = 2943) #CTE Geographical coordinates of WIRO
time = ti.now()

utc_offset = -7 * u.hour #CTE UTC is 7 hours ahead of MST

current_time = ti.now() + utc_offset #CTE The current time in MST

time_to_midnight = time_to(current_time, [24, 0, 0]) #CTE Finding the number hours from now to midnight

tonight_midnight = current_time + time_to_midnight #CTE Setting a time to tonight at midnight

obs_time = tonight_midnight #CTE Setting the first observation time to tonight at midnight

year_midnight = obs_time + 365*u.day #CTE Finding the date/time a year from now at midnight

ls = [] #CTE Empty list for l galactic coordinates to populate
bs = [] #CTE Empty list for b galactic coordinates to populate

while obs_time != year_midnight: #CTE Converting Equatorial zenith coordinates to Galactic Coordinates for each night this year
    coos = zenith_to_icrs(obs_time, wiro)
    coos_gal = coos.galactic
    ls.append(coos_gal.l.deg)
    bs.append(coos_gal.b.deg)
    obs_time += 1*u.day #CTE Stepping to the next day
    
#CTE Question: Is there a better way to do this without a loop? ^^^^
    
plot_times = np.arange(math.ceil(tonight_midnight.mjd), math.ceil(year_midnight.mjd), 1)

l_fig, l_ax=plt.subplots()

l_ax.scatter(plot_times, ls, color='navy', s=10)

l_ax.set_ylabel('Galactic Longitude (l)', fontsize=15)
l_ax.set_xlabel('MJD', fontsize=15)

ticks(l_ax, plot_times, ls, np.size(plot_times)*[0], np.size(ls)*[0])

save_show(l_fig, 'lvt.png')

b_fig, b_ax=plt.subplots()
b_ax.scatter(plot_times, bs, color='navy', s=10)

b_ax.set_ylabel('Galactic Latitude (l)', fontsize=15)
b_ax.set_xlabel('MJD', fontsize=15)

ticks(b_ax, plot_times, bs, np.size(plot_times)*[0], np.size(bs)*[0])

save_show(b_fig, 'bvt.png')

sol = (sc("21 11 18 -16 13 36", unit = (u.hourangle, u.deg), distance = 0.9858*u.au)).heliocentrictrueecliptic
moon = (sc("01 52 19 +14 12 15", unit= (u.hourangle, u.deg), distance = (368912*u.km).to(u.au))).heliocentrictrueecliptic
merc = (sc("20 58 35.5 -19 17 51", unit= (u.hourangle, u.deg), distance = 1.410*u.au)).heliocentrictrueecliptic
venus = (sc("23 54 44.6 02 05 07", unit= (u.hourangle, u.deg), distance = 0.499*u.au)).heliocentrictrueecliptic
mars = (sc("07 25 45.8 +26 13 17", unit= (u.hourangle, u.deg), distance = 0.699*u.au)).heliocentrictrueecliptic
jup = (sc("04 37 41.2 +21 36 13", unit= (u.hourangle, u.deg), distance = 4.590*u.au)).heliocentrictrueecliptic
sat = (sc("23 16 47.9 -06 43 54", unit= (u.hourangle, u.deg), distance = 10.441*u.au)).heliocentrictrueecliptic
uran = (sc("03 22 20.4 +18 15 59", unit= (u.hourangle, u.deg), distance = 19.393*u.au)).heliocentrictrueecliptic
nep = (sc("23 53 38.4 -02 04 18", unit= (u.hourangle, u.deg), distance = 30.614*u.au)).heliocentrictrueecliptic
pluto = (sc("20 19 27.6 -22 59 11", unit= (u.hourangle, u.deg), distance = 36.151*u.au)).heliocentrictrueecliptic

plans = [sol, moon, merc, venus, mars, jup, sat, uran, nep, pluto]
p_names = ['Sun', 'Moon', 'Mercury', 'Venus', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto']

lons = np.array([i.lon.deg for i in plans])
lats = np.array([i.lat.deg for i in plans])

p_fig, p_ax=plt.subplots()

for i, j in zip(plans, p_names):
    p_ax.scatter(i.lon.deg, i.lat.deg, label=j)
    
p_ax.legend(bbox_to_anchor=(1, 0.75))

#p_ax.scatter(lons, lats, color='navy', s=10)

p_ax.set_ylabel('Ecliptic Latitude (Degrees)', fontsize=15)
p_ax.set_xlabel('Ecliptic Latitude (Degrees)', fontsize=15)

p_ax.hlines(0, 400, 1, color='grey', linestyle='--', alpha = 0.5)

ticks(p_ax, lons, lats, np.size(lons)*[0], np.size(lats)*[0])

p_fig.tight_layout()

save_show(p_fig, 'planets.png')
