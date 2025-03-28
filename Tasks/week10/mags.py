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
from astropy.coordinates import SkyCoord as sc
import astropy.units as u

def to_ugriz(UBVRI):
    """
    Given UBVRI magnitudes, converts to ugriz.

    Parameters
    ----------
    UBVRI : list of floats
        A list of magnitudes in the format [U, B, V, R, I]

    Returns
    -------
    ugriz : list of floats
        A list of converted ugriz magnitudes.
    
    """

	ug = 1.28 * (UBVRI[0] - UBVRI[1]) + 1.13
	gr = 1.02 * (UBVRI[1] - UBVRI[2]) - 0.22
	ri = 0.91 * (UBVRI[3] - UBVRI[4]) - 0.20
	rz = 1.72 * (UBVRI[3] - UBVRI[4]) - 0.41
	g = UBVRI[2] + 0.60 * (UBVRI[1] - UBVRI[2]) - 0.12
	r = UBVRI[2] - 0.42 * (UBVRI[1] - UBVRI[2]) + 0.11

	u = ug + g
	i = r - ri
	z = r - rz

    ugriz = [u, g, r, i, z]

	return ugriz

def perc_diff(a, b):
    """
    Calculates the percent difference between to values.

    Parameters
    ----------
    a : flaot
        The first value.

    b : float
        The second value.

    Returns
    -------
    p_diff : float
        The percent difference between a and b
    """

    p_diff = (abs(a-b) / ((a + b) / 2)) * 100

	return p_diff

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

def f_to_m(flux):
    """
    Converts a given flux in nanomaggies to a magnitude.
    
    Parameters
    ----------
    flux : float
        The flux in nanomaggies.

    Returns
    -------
    m ; float or np.nan
        The converted magnitude. If flux <= 0, then a target is not detected
        and the magnitude is undefined (np.nan).
    """

    if flux <= 0:
        m = np.nan
    else:
        m = 22.5 - 2.5 * np.log10(flux)

    return m


V = 15.256
BV = 0.873
B = BV + V
UB = 0.320
U = UB + B
VR = 0.505
R = V - VR
RI = 0.511
I = R - RI

UBVRI = [U, B, V, R, I]

ugriz_nav = [17.30, 15.70, 15.19, 14.71, 14.55]

ugriz = to_ugriz(UBVRI)

#CTE Values converted from the table are close to the SDSS Nav values
tab_diffs = [perc_diff(i, j) for i, j in zip(ugriz, ugriz_nav)]

radec = {"RA":248.8583, "DEC":9.79806}
coo = sc(ra=radec['RA']*u.deg, dec=radec['DEC']*u.deg)

sweeps_dir = '/d/scratch/ASTR5160/data/legacysurvey/dr9/south/sweep/9.0'

sweep_file = in_sweeps(sweeps_dir, radec)

tab = Table.read(sweep_file[0])

sw_coo = sc(ra = tab['RA'], dec = tab['DEC'])

targ = (coo.separation(sw_coo)).argmin()

sw_g = f_to_m(tab['FLUX_G'][targ])
sw_r = f_to_m(tab['FLUX_R'][targ])
sw_z = f_to_m(tab['FLUX_Z'][targ])

sw_w1 = f_to_m(tab['FLUX_W1'][targ])
sw_w2 = f_to_m(tab['FLUX_W2'][targ])
sw_w3 = f_to_m(tab['FLUX_W3'][targ])
sw_w4 = f_to_m(tab['FLUX_W4'][targ])

sw_diffs = [perc_diff(i, j) for i, j in zip([sw_g, sw_r, sw_z], [ugriz_nav[1], ugriz_nav[2], ugriz_nav[-1]])]