import sys
sys.path.insert(0, "/d/users/caleb/Py_Modules")
from plots import min_max, tick_array, ticks, save_show
import matplotlib.pyplot as plt
import numpy as np
from numpy.random import random
from astropy import units as u
from astropy.coordinates import SkyCoord as sc
import math
import healpy as hp

def pixels(index, ras, decs, pix_locations):
	"""
    Given an index corresponding to a pixel in a HEALPix
    grid, finds the Equatorial coordinates that are in that pixel
    from a pregenerated set of coordinates.
	
    Parameters
    ----------
    index : int
    	The HEALPix pixel on which to run this function
	
	ras ; array of floats
		Array of pregenerated right ascension coordinates

	decs : array of floats
		Array of pregenerated declination coordinates

	pix_locations : array of ints
		List of pixel locations for each point in 'ras'/'decs'

    Returns
    -------
    ii : Boolean array
		Array of True/False for each point in 'ras'/'decs'. True corresponds to points where
		'pix_locations == index'

	ra : array of floats
		'ras' evaluated at all 'True' indices in 'ii'

	dec ; array of floats
		'decs' evaluated at all 'True' indices in 'ii'
    """

	ii = pix_locations == index
	ra = ras[ii]
	dec = decs[ii]

	return ii, ra, dec

#CTE Generating 1,000,000 random points that populate a sphere equally in area
ras = 360. * (random(1000000))
decs = (180 / np.pi) * np.arcsin(1. - random(1000000) * 2.)

nside = 1 #CTE subdivides sphere into Npix = 12 * Nside**2 = 12 pixels

nside_area = hp.nside2pixarea(nside, degrees = True) #CTE area of a pixel with 12 * nside**2 equal area pixels

pix_locs = hp.ang2pix(nside, ras, decs, lonlat = True) #CTE returns collection of pixel indices for each point

vals, cts = np.unique(pix_locs, return_counts = True) #CTE vals is unique values (pixel indices in this case)
																 #CTE cts is number of points corresponding to each val

pt_perc = (cts / 1000000) * 100 #CTE Each of the 12 pixels contains ~ 8.33% of points (consistent w/ equal area)

ii_2, ra_2, dec_2 = pixels(2, ras, decs, pix_locs) #CTE All ra/dec points inside nside=1 pixel 2
ii_5, ra_5, dec_5 = pixels(5, ras, decs, pix_locs)
ii_8, ra_8, dec_8 = pixels(8, ras, decs, pix_locs)

nside2 = 2

pix_locs2 = hp.ang2pix(nside2, ras, decs, lonlat = True)
vals2, cts2 = np.unique(pix_locs2, return_counts = True)
pt_perc2 = (cts2 / 1000000) * 100

ii_22, ra_22, dec_22 = pixels(2, ras, decs, pix_locs2)
ii_52, ra_52, dec_52 = pixels(5, ras, decs, pix_locs2)
ii_82, ra_82, dec_82 = pixels(8, ras, decs, pix_locs2)

in_5 = np.unique(pix_locs2[ii_5]) #CTE Finding the pixels from nside=2 that are in pixel 5 of nside=1
										 #CTE by evaluating pix_locs2 at the indices that are in pixel 5 of pix_locs
										 #CTE result of this code is [14, 21, 22, 30]

fig, ax = plt.subplots()

ax.scatter(ras, decs, color='grey', alpha=0.5)
ax.scatter(ra_2, dec_2, color='red')
ax.scatter(ra_5, dec_5, color='green')
ax.scatter(ra_8, dec_8, color='navy')

ax.set_xlim(0, 360)
ax.set_ylim(-90, 90)

save_show(fig, 'nside1-258.png')

fig, ax = plt.subplots()

ax.scatter(ras, decs, color='grey', alpha=0.5)
ax.scatter(ra_22, dec_22, color='red')
ax.scatter(ra_52, dec_52, color='green')
ax.scatter(ra_82, dec_82, color='navy')

ax.set_xlim(0, 360)
ax.set_ylim(-90, 90)

save_show(fig, 'nside2-258.png')
