import sys
sys.path.insert(0, "/d/users/caleb/Py_Modules")
from plots import min_max, tick_array, ticks, save_show
import matplotlib.pyplot as plt
import numpy as np
from numpy.random import random
from astropy import units as u
from astropy.coordinates import SkyCoord as sc
from astropy import wcs
import math

def ang_dist(v1, v2):
	"""
    Given two three-dimensional Cartesian points,
    finds the angular distance in degrees between them
	
    Parameters
    ----------
    v1 : astropy.coordinates.SkyCoord array
        The Cartesian coordinate array for the first point
        
 	v2 : astropy.coordinates.SkyCoord array
 		The Cartesian coordinate array for the second point

    Returns
    -------
    ang : float
    	The angle, in degrees, between v1 and v2
    """

	v1_vec = np.array([v1.x, v1.y, v1.z]) #CTE Creating a workable array with the x/y/z
										  #CTE components of v1
	v2_vec = np.array([v2.x, v2.y, v2.z])

	v1_mag = np.linalg.norm(v1_vec) #CTE Using np.linalg.norm to find the magnitude of v1
	v2_mag = np.linalg.norm(v2_vec)

	mag_prod = v1_mag*v2_mag

	d_prod = np.dot(v1_vec, v2_vec) #CTE Using np.dot to take the dot product of v1 and v2

	ang = math.degrees(np.arccos(d_prod/mag_prod))

	return ang

def rand_bounds(up, low, s, seed):
	"""
	Given upper and lower bounds, generates 's' random
	floats between them. Setting a seed keeps the same
	array every time the script is run.
	
    Parameters
    ----------
    up : float
    	The upper bound to generate a random value; the result is
    	less than but never equal to the upper bound

    low : float
    	The lower bound to generate a random value; the result can
    	be greater than or equal to the lower bound

    s : int
    	The number of random values between low and up to be generated

    seed : int or str
    	Initializes np.random.random to a fixed state and guarantees
    	the same result every time. If passed a string, does not set a seed

    Returns
    -------
    nums : array of floats
    	An array of randomly generated floats that are in the bounds [low, up)
    	and has 's' terms.  	
	"""

	if isinstance(seed, str): #CTE If passed a string, don't set a seed
		pass
	else:
		np.random.seed(seed)

	nums = (up - low) * np.random.random(size=s) + low

	return nums

coo1 = sc(ra=263.75*u.deg, dec=-17.9*u.deg)
coo2 = sc('20 24 59.9 +10 06 00', unit=(u.hourangle, u.deg))

apy_dist = (coo2.separation(coo1)).deg #CTE Finding the angular distance between coo1 
									   #CTE and coo2 using astropy's separation() function

coo1.representation_type = 'cartesian' #CTE Converting coo1 to Cartesian coordinates
coo2.representation_type = 'cartesian'

manual_dist = ang_dist(coo2, coo1) #CTE Using my ang_dist() function to manually
								   #CTE find the angular distance between coo1 and coo2

ras1 = rand_bounds(3.000000001, 2, 100, 'None') #CTE Generating an array of 100 RA
												#CTE coordinates from 2-3 hrs

decs1 = rand_bounds(2.000000001, -2, 100, 'None') #CTE Generating an array of 100 Dec
												  # coordinates from -2 to 2 deg

ras2 = rand_bounds(3.000000001, 2, 100, 'None')
decs2 = rand_bounds(2.000000001, -2, 100, 'None')

ran_coo1 = sc(ra=ras1*u.hourangle, dec=decs1*u.deg)
ran_coo2 = sc(ra=ras2*u.hourangle, dec=decs2*u.deg)

#CTE Finding pairs of coordinates in ran_coo1 and ran_coo2 that are within 10' of eachother
id1, id2, d2, d3 = ran_coo2.search_around_sky(ran_coo1, (10/60)*u.deg)

ran_coo1_near = ran_coo1[id1] #CTE Indexing ran_coo1 to give only the pairs
ran_coo2_near = ran_coo2[id2]

fig, ax = plt.subplots()

ax.scatter(ran_coo1.ra, ran_coo1.dec, marker='o', color='navy')
ax.scatter(ran_coo2.ra, ran_coo2.dec, marker='s', color='orange')

ax.scatter(ran_coo1_near.ra, ran_coo1_near.dec, marker='o', color='red')
ax.scatter(ran_coo2_near.ra, ran_coo2_near.dec, marker='s', color='red')

ax.set_aspect('equal')

ax.invert_xaxis()

ax.set_xlabel('Right Ascension (deg)', fontsize=15)
ax.set_ylabel('Declinaiton (deg)', fontsize=15)

save_show(fig, 'randpts.png')

all_ra = np.concatenate([ras1, ras2])
all_dec = np.concatenate([decs1, decs2])

all_coo = sc(ra = all_ra*u.hourangle, dec = all_dec*u.deg)

obs_coo = sc(ra='02h20m05s', dec='-00d06m12s', unit=(u.hourangle, u.deg))

ii = obs_coo.separation(all_coo) < 1.8*u.deg #CTE Finding all points within 1.8 deg of obs_coo

in_obs = all_coo[ii]

fig, ax = plt.subplots()

ax.scatter(obs_coo.ra.deg, obs_coo.dec.deg, marker='x', color='black', s=50)

ax.scatter(all_coo.ra.deg, all_coo.dec.deg, marker = 'o', color='navy')
ax.scatter(in_obs.ra.deg, in_obs.dec.deg, marker = 'o', color='red')

ax.set_aspect('equal')

ax.invert_xaxis()

ax.set_xlabel('Right Ascension (deg)', fontsize=15)
ax.set_ylabel('Declinaiton (deg)', fontsize=15)

save_show(fig, 'spec-obs.png')