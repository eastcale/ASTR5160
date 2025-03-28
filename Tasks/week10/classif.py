#!/usr/bin/env python3

import sys
import matplotlib.pyplot as plt
sys.path.insert(0, "/d/users/caleb/Py_Modules")
from plots import min_max, tick_array, ticks, save_show, ticks2
from matplotlib_inline.backend_inline import set_matplotlib_formats
set_matplotlib_formats('svg')
import numpy as np
import pandas as pd
from astropy.table import QTable, Table, Column, vstack
from astropy.io import fits
import os
import glob
from astropy.coordinates import SkyCoord as sc
import astropy.units as u

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

def star_or_qso(r_W1, g_z):
    """
    Given r-W1 and g-z colors, determines if a target is likely
    to be a star or a QSO.
    
    Parameters
    ----------
    r_W1 : float
        The r-W1 color of a target.
    
    g_z ; float
        The g-z color of a target

    Returns
    -------
    classif ; str
        The likely classification of the target. Either 'Star' or 'QSO'
    """
    
    c1 = line(g_z, 1, -1)
    c2 = line(g_z, -2.5, -0.5)
    
    if r_W1 > c1 and r_W1 > c2:
        classif = 'QSO'
    else:
        classif = 'Star'
    
    return classif

def line(x, m, b):
    """
    Given x points, a slope, and a y-intercept, generates
    a set of y points corresponding to a line.
    
    Parameters
    ----------
    x : list or array of floats
        The given x-points to find corresponding y-points

    m : float
        The slope of the line.

    b : float
        The y-intercept of the line.

    Returns
    -------
    y ; list or array of floats
        The corresponding y-points determined by a line equation.
    """

    y = m * x + b

    return y

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

def mag_append(tab, match, magnitude, flux, transmission):
    """
    Goes through the necessary conversions to append magnitudes to
    a given table.

    Parameters
    ----------
    tab : type ?
        The table to append to

    match : type ?
        The indices to append magnitudes at.

    magnitude : str
        The magnitude column to append to.

    flux : float
        The flux to convert and append.

    transmission : str
        The transmission column to assist in magnitude conversions.

    Returns
    ------
    None
    """
    tab[magnitude][match[1]] = [f_to_m(i) for i in ((sw_stack[match[0]][flux])
                                                   /(sw_stack[match[0]][transmission]))]

#CTE Reading in the star and qso files
direct = '/d/scratch/ASTR5160/week10/'
star_file = direct + 'stars-ra180-dec30-rad3.fits'
qso_file = direct + 'qsos-ra180-dec30-rad3.fits'

#CTE Setting the directory for sweep files
sweeps_dir = '/d/scratch/ASTR5160/data/legacysurvey/dr9/south/sweep/9.0'

#CTE Initializing star table and empty magnitude columns to be populated later.
stars = Table.read(star_file)
stars['g_mag'] = np.nan
stars['r_mag'] = np.nan
stars['z_mag'] = np.nan
stars['W1_mag'] = np.nan
stars['W2_mag'] = np.nan

#CTE Initializing qso table and empty magnitude columns to be populated later.
qsos = Table.read(qso_file)
qsos['g_mag'] = np.nan
qsos['r_mag'] = np.nan
qsos['z_mag'] = np.nan
qsos['W1_mag'] = np.nan
qsos['W2_mag'] = np.nan

#CTE Stars/QSOs in same 4 sweep files
star_insw = in_sweeps(sweeps_dir, stars)
qso_insw = in_sweeps(sweeps_dir, qsos)

#CTE Creating coordinate objects for stars and QSOs
st_coos = sc(ra = stars['RA']*u.deg, dec = stars['DEC']*u.deg)
qso_coos = sc(ra = qsos['RA']*u.deg, dec = qsos['DEC']*u.deg)

sw_stack = vstack([Table(fits.open(i, memmap=True)[1].data) for i in star_insw])

sw_coos = sc(sw_stack['RA']*u.deg, sw_stack['DEC']*u.deg)

stars_sw = st_coos.search_around_sky(sw_coos, 0.5*u.arcsec)
qsos_sw = qso_coos.search_around_sky(sw_coos, 0.5*u.arcsec)

mag_append(stars, stars_sw, 'g_mag', 'FLUX_G', 'MW_TRANSMISSION_G')
mag_append(stars, stars_sw, 'r_mag', 'FLUX_R', 'MW_TRANSMISSION_R')
mag_append(stars, stars_sw, 'z_mag', 'FLUX_Z', 'MW_TRANSMISSION_Z')
mag_append(stars, stars_sw, 'W1_mag', 'FLUX_W1', 'MW_TRANSMISSION_W1')
mag_append(stars, stars_sw, 'W2_mag', 'FLUX_W2', 'MW_TRANSMISSION_W2')

mag_append(qsos, qsos_sw, 'g_mag', 'FLUX_G', 'MW_TRANSMISSION_G')
mag_append(qsos, qsos_sw, 'r_mag', 'FLUX_R', 'MW_TRANSMISSION_R')
mag_append(qsos, qsos_sw, 'z_mag', 'FLUX_Z', 'MW_TRANSMISSION_Z')
mag_append(qsos, qsos_sw, 'W1_mag', 'FLUX_W1', 'MW_TRANSMISSION_W1')
mag_append(qsos, qsos_sw, 'W2_mag', 'FLUX_W2', 'MW_TRANSMISSION_W2')

x_plot = np.linspace(-2, 7, 1000)

st_rw1 = stars['r_mag'] - stars['W1_mag']
st_gz = stars['g_mag'] - stars['z_mag']

qso_rw1 = qsos['r_mag'] - qsos['W1_mag']
qso_gz = qsos['g_mag'] - qsos['z_mag']

fig, ax = plt.subplots()

ax.scatter(st_gz, st_rw1, marker = '*', color = 'navy', s = 5, label = 'Stars')
ax.scatter(qso_gz, qso_rw1, color = 'red',  s = 5, label = 'QSOs')

ax.plot(x_plot, line(x_plot, 1, -1), linestyle = '--', alpha = 0.5, color = 'grey')
ax.plot(x_plot, line(x_plot, -2.5, -0.5), linestyle = '--', alpha = 0.5, color = 'grey')

ticks2(ax, [-2, 7, 1, 0.25], [-5, 5, 1, 0.25])

ax.set_xlabel('g-z (mag)', fontsize = 15)
ax.set_ylabel('r-W1 (mag)', fontsize = 15)

# save_show(fig, 'star_qso_col.svg')