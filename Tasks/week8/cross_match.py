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

def sdss_cmatch(ra, dec, file):
	[os.system(f'python sdssDR9query.py {i} {j} >> {file}.txt') for i, j in zip(ra_100, dec_100)]

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
	sweeps_dir = directory
	sweep_list = np.array(glob.glob(sweeps_dir + "/*.fits"))

	boxes = [decode_sweep_name(i) for i in sweep_list]

	in_box = [is_in_box(table, i) for i in boxes]

	true_false = np.array([np.any(i) for i in in_box])

	indexes = np.where(true_false == True)[0]

	in_sweeps = sweep_list[indexes]

	return in_sweeps

if __name__ == '__main__':
	first_tab = Table.read('/d/scratch/ASTR5160/data/first/first_08jul16.fits')
	first_cols = first_tab.columns
	ra, dec = first_tab['RA'], first_tab['DEC']

	first_tab100 = first_tab[0:100]

	ra_100, dec_100 = ra[0:100], dec[0:100]

	fig, ax = plt.subplots()

	ax.scatter(ra, dec, color='navy', s = .0001)

	ticks(ax, ra, dec)

	ax.set_xlim(-10, 370)
	ax.set_xlabel('Right Ascension (deg)', fontsize = 15)
	ax.set_ylabel('Declination (deg)', fontsize = 15)

	save_show(fig, 'first_08jul16.png')

	#CTE This is commented out because it takes so long!!
	# sdss_cmatch(ra_100, dec_100, 'sdss_cmatch')

	first_sweeps = in_sweeps('/d/scratch/ASTR5160/data/legacysurvey/dr9/north/sweep/9.0', first_tab100)