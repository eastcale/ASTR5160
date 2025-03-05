import sys
import os
import numpy as np
import pymangle
lap_path = '/Users/calebeastlund/Desktop/'
dep_path = '/d/dor1/caleb/Classes/AstrTechII/Notes/week5/week5-code/'
if os.path.exists(lap_path):
	sys.path.insert(0, lap_path)
elif os.path.exists(dep_path):
	sys.path.insert(0, dep_path)
	sys.path.insert(0, "/d/users/caleb/Py_Modules")
from plots import min_max, tick_array, ticks, save_show
from sph_caps import ra_cap, dec_cap, sph_cap, cap_to_str, to_poly, to_cart
import matplotlib.pyplot as plt
from numpy.random import random

#CTE FYI FOR ADAM: I made signifcant changes to my sph_caps.py
#CTE file from week5 in light of this assignment. Many of the functions I
#CTE call have been altered since they were last evaluated

def mk_rect(x, y, ax, col):
	"""
	Given the corners of a rectangle, plots a rectangle on
	the given axis in a given color.

	Parameters
    ----------
    x : list or array of floats
    	The x coordinates for the vertices of the rectangle from left to right
    	with no duplicates. e.g. if the corners are 
    	([x1, y1], [x2, y1], [x1, y2], [x2, y2]), then 'x' is just [x1, x2]

    y : list or array of floats
    	The y coordinates for the vertices of the rectangle from bottom to top
    	with no duplicates. e.g. if the corners are 
    	([x1, y1], [x2, y1], [x1, y2], [x2, y2]), then 'y' is just [y1, y2]

    ax : matplotlib.pyplot axs
		The axis to plot the rectangle on.

    col : str
    	The color to plot the rectangle in. Must be a matplotlib color.

    Returns
    -------
	None
	"""

	ax.vlines(x[0], y[0], y[1], color=col)
	ax.vlines(x[1], y[0], y[1], color=col)

	ax.hlines(y[0], x[0], x[1], color=col)
	ax.hlines(y[1], x[0], x[1], color=col)

def rand_radec(num):
	"""
	Generates a given number of random points distributed
	across a sphere.

	Parameters
    ----------
    num : int
    	The number of points to generate.

    Returns
    -------
    ras : array of floats
    	The right ascension for each random point.

    decs : array of floats
    	The declination for each random point.
	"""

	ras = 360. * (random(num))
	decs = (180 / np.pi) * np.arcsin(1. - random(num) * 2.)

	return ras, decs

if __name__ == '__main__':

	#CTE Defining the caps for the first rectangle
	p1c1 = ra_cap(5, 1)
	p1c2 = ra_cap(6, -1)
	p1c3 = dec_cap(30, 1) 
	p1c4 = dec_cap(40, -1)

	#CTE Making the first rectangle a polygon 
	p1c1_str = cap_to_str(p1c1, 8)
	p1c2_str = cap_to_str(p1c2, 8)
	p1c3_str = cap_to_str(p1c3, 8)
	p1c4_str = cap_to_str(p1c4, 8)

	p1c_list = [[p1c1_str, p1c2_str, p1c3_str, p1c4_str]]

	#CTE Making a polgyon file INCLUDING the area of the rectangle
	to_poly('rect1.ply', 1, [4], p1c_list, [5, 30], [6, 40], 'y')

	#CTE Making a mask for the first rectangle
	r1mask = pymangle.Mangle('rect1.ply')

	#CTE Following the same process as above for the second rectangle
	p2c1 = ra_cap(11, 1)
	p2c2 = ra_cap(12, -1)
	p2c3 = dec_cap(60, 1)
	p2c4 = dec_cap(70, -1)

	p2c1_str = cap_to_str(p2c1, 8)
	p2c2_str = cap_to_str(p2c2, 8)
	p2c3_str = cap_to_str(p2c3, 8)
	p2c4_str = cap_to_str(p2c4, 8)

	p2c_list = [[p2c1_str, p2c2_str, p2c3_str, p2c4_str]]

	#CTE Making a polgyon file INCLUDING the area of the rectangle
	to_poly('rect2.ply', 1, [4], p2c_list, [11, 60], [12, 40], 'y')

	r2mask = pymangle.Mangle('rect2.ply')

	#CTE Generating 1,000,000 random points that populate a sphere equally in area
	ras, decs = rand_radec(1000000)

	#CTE Generating 10,000 random points in each of the rectangles
	p1_ra, p1_dec = r1mask.genrand(10000)
	p2_ra, p2_dec = r2mask.genrand(10000)

	#CTE Finding all points in (ras, decs) that are in each rectangle
	inr1_in = r1mask.contains(ras, decs)
	inr2_in = r2mask.contains(ras, decs)

	#CTE Indexing (ras, decs) based on the information found above
	inr1_ra, inr1_dec = ras[inr1_in], decs[inr1_in]
	inr2_ra, inr2_dec = ras[inr2_in], decs[inr2_in]

	#CTE Plotting points in each rectangle based on the mask.contains method
	fig, ax = plt.subplots()

	mk_rect([5*15, 6*15], [30, 40], ax, 'navy')
	mk_rect([11*15, 12*15], [60, 70], ax, 'red')

	ax.scatter(ras, decs, color='green', alpha=0.1, s=0.01)

	ax.scatter(inr1_ra, inr1_dec, color = 'navy', alpha=0.5, s=0.01)
	ax.scatter(inr2_ra, inr2_dec, color = 'red', alpha=0.5, s=0.01)

	ax.set_xlim(0, 360)
	ax.set_ylim(-90, 90)

	save_show(fig, 'rects-contains.png')

	#CTE Plotting points in each rectangle based on the mask.genrand() method
	fig, ax = plt.subplots()

	mk_rect([5*15, 6*15], [30, 40], ax, 'navy')
	mk_rect([11*15, 12*15], [60, 70], ax, 'red')

	ax.scatter(ras, decs, color='green', alpha=0.1, s=0.01)

	ax.scatter(p1_ra, p1_dec, color = 'navy', alpha=0.5, s=0.01)
	ax.scatter(p2_ra, p2_dec, color = 'red', alpha=0.5, s=0.01)

	ax.set_xlim(0, 360)
	ax.set_ylim(-90, 90)

	save_show(fig, 'rects-genrand.png')