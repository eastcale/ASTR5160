#!/usr/bin/env python3

import sys
import matplotlib.pyplot as plt
sys.path.insert(0, "/d/users/caleb/Py_Modules")
from plots import min_max, tick_array, ticks, save_show, ticks2
from matplotlib_inline.backend_inline import set_matplotlib_formats
set_matplotlib_formats('svg')
import numpy as np
import pandas as pd
from astropy.table import QTable, Table, Column
from astropy.io import fits
import os
import glob

def decode_sweep_name(sweepname, nside=None, inclusive=True, fact=4):
    """Retrieve RA/Dec edges from a full directory path to a sweep file

    Parameters
    ----------
    sweepname : :class:`str`
        Full path to a sweep file, e.g., /a/b/c/sweep-350m005-360p005.fits
    nside : :class:`int`, optional, defaults to None
        (NESTED) HEALPixel nside
    inclusive : :class:`book`, optional, defaults to ``True``
        see documentation for `healpy.query_polygon()`
    fact : :class:`int`, optional defaults to 4
        see documentation for `healpy.query_polygon()`

    Returns
    -------
    :class:`list` (if nside is None)
        A 4-entry list of the edges of the region covered by the sweeps file
        in the form [RAmin, RAmax, DECmin, DECmax]
        For the above example this would be [350., 360., -5., 5.]
    :class:`list` (if nside is not None)
        A list of HEALPixels that touch the  files at the passed `nside`
        For the above example this would be [16, 17, 18, 19]
    """
    # ADM extract just the file part of the name.
    sweepname = os.path.basename(sweepname)

    # ADM the RA/Dec edges.
    ramin, ramax = float(sweepname[6:9]), float(sweepname[14:17])
    decmin, decmax = float(sweepname[10:13]), float(sweepname[18:21])

    # ADM flip the signs on the DECs, if needed.
    if sweepname[9] == 'm':
        decmin *= -1
    if sweepname[17] == 'm':
        decmax *= -1

    if nside is None:
        return [ramin, ramax, decmin, decmax]

    pixnum = hp_in_box(nside, [ramin, ramax, decmin, decmax],
                       inclusive=inclusive, fact=fact)

    return pixnum

def is_in_box(objs, radecbox):
    """Determine which of an array of objects are inside an RA, Dec box.

    Parameters
    ----------
    objs : :class:`~numpy.ndarray`
        An array of objects. Must include at least the columns "RA" and "DEC".
    radecbox : :class:`list`
        4-entry list of coordinates [ramin, ramax, decmin, decmax] forming the
        edges of a box in RA/Dec (degrees).

    Returns
    -------
    :class:`~numpy.ndarray`
        ``True`` for objects in the box, ``False`` for objects outside of the box.

    Notes
    -----
        - Tests the minimum RA/Dec with >= and the maximum with <
    """
    ramin, ramax, decmin, decmax = radecbox

    # ADM check for some common mistakes.
    if decmin < -90. or decmax > 90. or decmax <= decmin or ramax <= ramin:
        msg = "Strange input: [ramin, ramax, decmin, decmax] = {}".format(radecbox)
        log.critical(msg)
        raise ValueError(msg)

    ii = ((objs["RA"] >= ramin) & (objs["RA"] < ramax)
          & (objs["DEC"] >= decmin) & (objs["DEC"] < decmax))

    return ii

def sdss_cmatch(ra, dec, file):
	"""
    Given an array of RAs and DECs, 
    queries SDSS using sdssDR9query.py for objects near each point 
    and writes results to a file.
	
    Parameters
    ----------
    ra : array of floats
    	Right ascension coordinates in degrees,
	
	dec ; array of floats
		Declination coordinates in degrees

	file : str
		The file name to write to. If existing, will append

    Returns
    -------
    None
    """

	[os.system(f'python sdssDR9query.py {i} {j} >> {file}.txt') for i, j in zip(ra, dec)]

def in_sweeps(directory, table):
	"""
    Given a directory containing sweep files and a table with 
    atleast RA and DEC columns, returns the path to any sweep file
    that contains atleast one point from the given table.
	
    Parameters
    ----------
    directory : str
    	The full path to the directory containing the sweep files.
	
	table ; astropy.table Table
		Table of data; must contain atleast RA and DEC columns

    Returns
    -------
    in_sweeps ; list of str
    	A list of paths to any sweep file that contains atleast one of the given points.
    """

    #CTE Setting the directory to search through
	sweeps_dir = directory
	sweep_list = np.array(glob.glob(sweeps_dir + "/*.fits"))

	#CTE Creating a list of boxes from the sweep names
	boxes = [decode_sweep_name(i) for i in sweep_list]

	#CTE Checking every point in every box,
	#CTE and creating a list that is len(boxes) long, where
	#CTE each term is len(table) long containing True or False
	#CTE for if a point is in the box.
	in_box = [is_in_box(table, i) for i in boxes]

	#CTE Creating an array of True/False that is len(boxes) long
	#CTE that assigns True to any box that has atleast one point from table
	true_false = np.array([np.any(i) for i in in_box])

	#CTE Finding the box indices for all Trues
	indexes = np.where(true_false == True)[0]

	#CTE Evaluating the sweep list at all Trues from above
	in_sweeps = sweep_list[indexes]

	return in_sweeps

if __name__ == '__main__':
	#CTE Reading in the FIRST radio data
	first_tab = Table.read('/d/scratch/ASTR5160/data/first/first_08jul16.fits')

	#CTE Giving columns for first_tab
	first_cols = first_tab.columns

	#CTE Getting RA and DEC coordinates from first_tab
	ra, dec = first_tab['RA'], first_tab['DEC']

	#CTE Getting the first 100 rows of first_tab
	first_tab100 = first_tab[0:100]

	#CTE Getting the first 100 points of first_tab
	ra_100, dec_100 = ra['RA'], dec['DEC']

	#CTE Cross matching the first 100 points of first_tab with SDSS
	#CTE This is commented out because it takes so long!!
	# sdss_cmatch(ra_100, dec_100, 'sdss_cmatch')

	first_sweeps = in_sweeps('/d/scratch/ASTR5160/data/legacysurvey/dr9/north/sweep/9.0', first_tab100)

	#CTE Plotting all points in first_tab
	fig, ax = plt.subplots()

	ax.scatter(ra, dec, color='navy', s = .0001)

	ticks(ax, ra, dec)

	ax.set_xlim(-10, 370)
	ax.set_xlabel('Right Ascension (deg)', fontsize = 15)
	ax.set_ylabel('Declination (deg)', fontsize = 15)

	save_show(fig, 'first_08jul16.png')

