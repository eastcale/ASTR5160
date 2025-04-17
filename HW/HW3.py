#!/usr/bin/env python3

import matplotlib.pyplot as plt
from Modules.plots import ticks2
import numpy as np
from astropy.table import Table
from astropy.io import fits
import astropy.units as u
from astropy.coordinates import SkyCoord as sc
import os
import glob
from time import sleep
import time
from PIL import Image as im
import argparse

#CTE Until otherwise said, the following code is adapted
#CTE from code written by Adam Myers (ADM)

class sdssQuery:
    """
    NAME: sdssQuery
 
    PURPOSE: class that can be initialized using Python's urllib tools to
    send an SQL command to SDSS web services
 
    CALLING SEQUENCE: from the UNIX command line:
      
      python sdssDR9query.py ra dec

    INPUTS: ra and dec shoud be sent from the command line:

      ra - Right Ascension of position to query around in SDSS DR9
      dec - declination of position to query around in SDSS DR9

    OUTPUTS: the result of the SQL command called "query" in the
    code, below, is executed by the SDSS DR9 SQL API and
    printed at the command line.
    
    COMMENTS: This is adapted from an example provide by the SDSS 
    in an early tutorial on web services.

    Note that the SQL command passed as "query" can be changed to
    any valid SDSS string.
 
    EXAMPLES: At the Unix command line:

      python sdssDR9query.py 145.285854 34.741254
      
    should return:

      145.28585385,34.74125418,21.132605,20.059248,19.613457,19.379345,19.423323,7.7489487E-4
    """
    # ADM this is the URL of the SDSS "web services" API.
    url='http://skyserver.sdss3.org/dr9/en/tools/search/x_sql.asp'
    # ADM always return the output in .csv format
    format = 'csv'

    # ADM initialize the class with a null query.
    def __init__(self):
        self.query = ''
        self.cleanQuery = ''

    # ADM use Python's urllib module to initialize a query string.
    def executeQuery(self):
        from urllib.parse import urlencode
        from urllib.request import urlopen
        self.filterQuery()
        params = urlencode({'cmd': self.cleanQuery, 'format':self.format})
        return urlopen(self.url + '?%s' % params)

    # ADM this cleans up the syntax in the query string.
    def filterQuery(self):
        from os import linesep
        self.cleanQuery = ''
        tempQuery = self.query.lstrip()
        for line in tempQuery.split('\n'):
            self.cleanQuery += line.split('--')[0] + ' ' + linesep;

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

def run_q(query):
    """
    Given a prompt, queries SDSS DR9 for the requested data
    
    Parameters
    ----------
    query : str
        An SDSS DR9 query in the required format.
        
    Returns
    -------
    res : str
        The result of the SDSS DR9 query in the format 'ra,dec,u,i'.
        If no object is detected, the result is 'No objects have been found'
    """

    # ADM execute the query.
    qry.query = query
    for line in qry.executeQuery():
        result = line.strip()

    # ADM NEVER remove this line! It won't speed up your code, it will
    # ADM merely overwhelm the SDSS server (a denial-of-service attack)!
    sleep(1)

    # ADM the server returns a byte-type string. Convert it to a string.
    res = result.decode()

    return res

def sdss_ui_query(ra, dec):
    """
    Given a list of right ascension and declination coordinates, 
    initializes an SDSS DR9 query for u and i magnitudes for each
    coordinate.
    
    Parameters
    ----------
    ra : list of floats
        A list of right ascension coordinates to query
        SDSS DR9 for.

    dec : list of floats
        A list of declination coordinates to query
        SDSS DR9 for.
        
    Returns
    -------
    qry : class sdssQuery()
        A variable to initialize an sdssQuery() object. MUST be set
        to qry.

    query : list of strings
        A list of SDSS DR9 query prompts
    """

    # ADM initialize the query.
    qry = sdssQuery()

    # ADM the query to be executed. You can substitute any query, here!
    #CTE Querying by coordinates for ra, dec, u, and i data.
    query = ["""SELECT top 1 ra,dec,u,i FROM PhotoObj as PT
    JOIN dbo.fGetNearbyObjEq(""" + i + """,""" + j + """,0.02) as GNOE
    on PT.objID = GNOE.objID ORDER BY GNOE.distance""" for i, j in zip(ra, dec)]
    
    return qry, query

#CTE End of Adapted Code

def decode_res(result):
    """
    Intended to be used in conjunction with decode_res2. Given an SDSS
    query result in the form 'ra,dec,u,i', begins the process of decoding
    the result into workable arrays.
    
    Parameters
    ----------
    result : str
        An SDSS query result in the form 'ra,dec,u,i'
        
    Returns
    -------
    info : list of floats
        The SDSS query result separated into floats for each parameter
        in the form [ra, dec, u, i].
    """

    #CTE Splitting the result at commas
    split = result.split(',')

    #CTE Indexing the split for each parameter
    ra = float(split[0])
    dec = float(split[1])
    u = float(split[2])
    i = float(split[3])

    #CTE Generating a list of data
    info = [ra, dec, u, i]
    
    return info

def decode_res2(result):
    """
    Intendec to be used in conjunction with decode_res. Given an SDSS
    query result in the form 'ra,dec,u,i', decodes the data into workable
    lists.
    
    Parameters
    ----------
    result : str
        A list of SDSS query result in the form ['ra,dec,u,i', 'ra,dec,u,i', ...]
        
    Returns
    -------
    ra : list of floats
        A list of right ascension coordinates from the SDSS query results.

    dec : list of floats
        A list of declination coordinates from the SDSS query.

    u_mag : list of floats
        A list of u-magnitudes from the SDSS query results.

    i_mag : list of floats
        A list of i-magnitudes from the SDSS query results.
    """
    
    #CTE Finding all results that did NOT come up with no objects found
    no_detection = result != 'No objects have been found'
    
    #CTE Using decode_res() to get a list of results for every other target
    res = [decode_res(i) for i in result[no_detection]]

    #CTE Separating the lists into individual parameters
    ra = [i[0] for i in res]
    dec = [i[1] for i in res]
    u_mag = [i[2] for i in res]
    i_mag = [i[3] for i in res]
    
    #CTE Finding indices of objects that had no result
    no_detections_ins = np.where(no_detection == False)[0]
    
    #CTE Inserting np.nans at the indices where objects were not detected
    ra_ins = [ra.insert(i, np.nan) for i in no_detections_ins]
    dec_ins = [dec.insert(i, np.nan) for i in no_detections_ins]
    u_ins = [u_mag.insert(i, np.nan) for i in no_detections_ins]
    i_ins = [i_mag.insert(i, np.nan) for i in no_detections_ins]
    
    return ra, dec, u_mag, i_mag

def save_show(fig, file_name):
    """
    Intended to be used after a figure has been made. Saves a .png
    of the given figure; overwrites previous saves if the same name is given!
    
    Parameters
    ----------
	fig : matplotlib.pyplot figure
		The figure to save.

    file_name : string
        The desired name of the file to be saved.
        The only accepted formats are currently .png or .svg files.
        
    Returns
    -------
    None
    """
    
    fig.savefig(file_name)
    im.open(file_name).show()


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

def match_sweep(path1, sweep_dir, cen):
    """
    Given file paths to a data set and a directory containing
    Legacy Survey sweep files, begins the process of cross matching
    both data sets that have data centered at a given point
    with a given radius. Prints the time taken to match files
    
    Parameters
    ----------
    path1 : str
        The full path to a data file to be cross matched
        to Legacy Survey sweep files. Must have atleast
        'RA' and 'DEC' columns.
    
    sweep_dir ; str
        The full path to the directory containing the Legacy Survey
        sweep files.

    cen : list or array of floats
        The center coordinate and radius to query around. Must be in the form
        [RA, DEC, radius] all in degrees

    Returns
    -------
    path1_coos : array of astropy.SkyCoord coordinate objects
        ALL coordinates contained within the first data file.

    path1_dict : dictionary
        A dictionary containing the coordinates of only the path1 objects
        that are within the range specified by 'cen'

    sweep_match : list of str
        A list of paths to all sweep files that contain atleast one object from path1_dict

    sweep_stack : array of Legacy Survey sweep file data
        ALL data from the sweep files specified by sweep_match; if multiple sweep
        files are read, they are stacked into one array (still indexable by columns)
    """

    print('Matching sweep files...')
    time1 = time.time()
    
    #CTE Reading in the first data set
    tab = Table.read(path1)

    #CTE Getting coordinate objects 
    path1_coos = sc(tab['RA'], tab['DEC'])
    
    #CTE Getting indices for objects in first data set that are within the specified range
    in_range = sc(cen[0]*u.deg, cen[1]*u.deg).separation(path1_coos) < cen[2]*u.deg

    #CTE Getting coordinates for the objects within the specified range
    path1_picks = path1_coos[in_range]
    path1_dict = {'RA':path1_picks.ra/u.deg, 'DEC':path1_picks.dec/u.deg}
    
    #CTE Finding any sweep file that contains at least one object from path1_dict
    sweep_match = in_sweeps(sweep_dir, path1_dict)
    
    #CTE Collecting all data from sweep_match files into one array
    sweep_stack = np.hstack([fits.open(i, memmap=True)[1].data for i in sweep_match])

    print(f'---> Matching sweep files took {time.time() - time1} seconds')
    
    return path1_coos, path1_dict, sweep_match, sweep_stack

def sweep_pick(pick_dict, sweep_stack, col_min, r_max, dist = 1*u.arcsec):
    """
    Meant to be used in conjunction with match_sweep. Given a set of coordinates
    and an array of Legacy Survey sweep file data, cross matches the two and 
    finds all common objects that are both above a minimum W1-W2 color
    and below a given r-magnitude.
    
    Parameters
    ----------
    pick_dict : dictionary
        A dictionary of coordinates to cross match with Legacy Survey
        sweep files

    sweep_stack ; array of Legacy Survey sweep file data
        Array of Legacy Survey sweep file data to cross match with
        pick_dict

    col_min : float
        A minimum W1-W2 color to limit data. Will exclude objects that
        have a color below this limit AND an r-magnitude above the following limit.

    r_max : float
        A maximum r-magnitude to limit data. Will exclude objects that
        have a magnitude above this limit AND a color below the previous limit.

    dist : float
        A cross matching distance; finds coordinates that are within this distance from eachother.
        Must specify units. Default is 1 arcsecond (recommended)

    Returns
    -------
    sweep_picks_fin : array of ints
        Indexes corresponding to the sweep file objects that meet the
        specified criteria

    picked : array of Legacy Survey sweep file data
        Data for all sweep file objects that meet the specified criteries

    picked_ra : list of floats
        List of right ascension coordinates for all sweep file objects that
        meet the specified criteria.

    picked_dec : list of floats
        List of declination coordinates for all sweep file objects that meet
        the specified criteria.
    """
    print('Selecting targets...')
    time1 = time.time()

    #CTE Cross matching the given coordinates to the given Legacy Survey sweep files.
    #CTE Finds all coordinates that are within 'dist' apart from eachother.
    pick_crit_1 = sc(pick_dict['RA']*u.deg, pick_dict['DEC']*u.deg).search_around_sky(
    sc(sweep_stack['RA']*u.deg, sweep_stack['DEC']*u.deg), dist)
    
    #CTE Getting all the cross matched sweep data.
    sweep_picks1 = sweep_stack[pick_crit_1[0]]
    
    #CTE Calculating colors and r-magnitudes for all the cross matched data
    w1 = np.array([f_to_m(i) for i in sweep_picks1['FLUX_W1']])
    w2 = np.array([f_to_m(i) for i in sweep_picks1['FLUX_W2']])
    col = w1 - w2
    r = np.array([f_to_m(i) for i in sweep_picks1['FLUX_R']])
    
    #CTE Finding where cross matched data meets given criteria
    sweep_picks_fin = np.where((col > col_min) & (r < r_max))[0]
    
    #CTE Getting only targets from the cross matched data
    #CTE that meet the given criteria.
    picked = sweep_picks1[sweep_picks_fin]
    picked_ra = [str(i) for i in picked['RA']]
    picked_dec = [str(i) for i in picked['DEC']]

    print(f'{len(sweep_picks_fin)} sources have been found with r < 22 and W1 - W2 > 0.5.')
    
    print(f'---> Selecting targets took {time.time() - time1} seconds.')

    return sweep_picks_fin, picked, picked_ra, picked_dec

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

    #CTE Making a copy of flux to do operations on
    m = np.copy(flux)
    m = m.astype(float)

    #CTE Setting magnitudes for all objects
    #CTE with negative fluxes to np.nan
    #CTE as a negative flux means no detection
    no_det = flux <= 0
    m[no_det] = np.nan
    
    #CTE Converting all other fluxes to magnitudes
    m[~no_det] = [22.5 - 2.5*np.log10(i) for i in m[~no_det]]

    return m

def m_to_f(mag):
    """
    Converts a given magnitude to a flux in nanomaggies.
    
    Parameters
    ----------
    mag : float
        The magnitude to convert.

    Returns
    -------
    f : float or -1
        The converted flux. If mag == np.nan, then a target 
        is not detected and the flux is negative (by convention).
    """

    #CTE Making a copy of mag to do operations on
    f = np.copy(mag)

    #CTE Setting flux for all targets with magnitude np.nan
    #CTE to -1, as negative flux means no detection.
    no_det = np.isnan(mag)
    f[no_det] = -1
    
    #CTE Converting all other magnitudes to fluxes.
    f[~no_det] = [10**(-1*(i-22.5)/2.5) for i in f[~no_det]]

    return f

def get_fluxes(index, tab, u_mag, i_mag):
    """
    Given a table, retrieves SDSS ugriz and WISE W1W2W3W4
    fluxes in nanonmaggies from a Legacy Survey sweep file
    at a given index. u and i magnitudes must be provided.
    
    Parameters
    ----------
    index : int
        The index of the target for which fluxes are needed.

    tab : array of Legacy Survey sweep file data
        The table from which to pull fluxes.

    u_mag : list of floats
        A list of u-magnitudes. Must be the same length as the columns
        in tab and have matching indices (i.e., u_mag[0] should provide the 
        u-magnitude for a target, and tab['FLUX_R'][0] should provide the r flux for
        the SAME target)

    i_mag : list of floats
        A list of i-magnitudes. Must be the same length as the columns
        in tab and have matching indices (i.e., i_mag[0] should provide the 
        i-magnitude for a target, and tab['FLUX_R'][0] should provide the r flux for
        the SAME target)

    Returns
    -------
    coos : astropy.SkyCoord coordinate object
        The coordinates for the object in question

    sdss_fluxes : array of floats
        An aray of SDSS ugriz fluxes for the object in the form
        np.array([u, g, r, i, z]).

    sdss_lams : array of ints
        An array of wavelengths in angstroms for the SDSS ugriz
        passbands. Will always be np.array([3543, 4770, 6231, 7625, 9134])

    wise_fluxes : array of floats
        An array of WISE W1W2W3W4 fluxes for the object in the form
        np.array([W1, W2, W3, W4]).

    wise_lams : array of ints
        An array of wavelengths in angstroms for the WISE
        W1W2W3W4 passbands. Will always be 
        np.array([33680, 46180, 120820, 221940])
    """

    #CTE Getting the coordinates for the object at index
    coos = sc(tab['RA'][index]*u.deg, tab['DEC'][index]*u.deg)
    
    #CTE Getting the u and i magnitudes for the object at index
    u_mag = u_mag[index]
    i_mag = i_mag[index]
    
    #CTE Indexing G, R, and Z fluxes, and converting u and i magnitudes
    #CTE to fluxes. Writing into array in the form np.array([u, g, r, i, z])
    sdss_fluxes = np.array([
        m_to_f(u_mag), tab['FLUX_G'][index], 
        tab['FLUX_R'][index], m_to_f(i_mag), tab['FLUX_Z'][index]])

    #CTE Setting the wavelenghts for ugriz passbands
    sdss_lams = np.array([3543, 4770, 6231, 7625, 9134])
    
    #CTE Indexing all W1W2W3W4 fluxes and writing
    #CTE into array in the form np.array([W1, W2, W3, W4])
    wise_fluxes = np.array([
        tab['FLUX_W1'][index], tab['FLUX_W2'][index], 
        tab['FLUX_W3'][index], tab['FLUX_W4'][index]])
    
    #CTE Setting the wavlength for the W1W2W3W4 passpbands
    wise_lams = np.array([33680, 46180, 120820, 221940])
    
    return coos, sdss_fluxes, sdss_lams, wise_fluxes, wise_lams

if __name__ == '__main__':
    #CTE Initializing ArgumentParser description
    parser = argparse.ArgumentParser(description = 'This module cross-matches\
    FIRST data collected on July 8, 2016 with data from Legacy Survey sweep files\
    and SDSS DR9 and collects objects in a circle centered on (RA, DEC) = (163, 50) with a\
    radius of 3 degrees. The module selects targets with SDSS r-magnitudes less than 22 and\
    WISE W1 - W2 colors greater than 0.5. The module also plots a photometric SED for a target\
    dubbed "ubrite1", that is the brightest of the selected targets in the SDSS u-band. This\
    code makes use of code written by Dr. Adam Myers at the University of Wyoming')
    args = parser.parse_args()

    #CTE Starting a 'timer'
    start = time.time()

    #CTE Setting paths to the FIRST data file and the directory containing relevant sweep files
    first_path = '/d/scratch/ASTR5160/data/first/first_08jul16.fits'
    sweep_path = '/d/scratch/ASTR5160/data/legacysurvey/dr9/north/sweep/9.0/'

    #CTE Finding sweep files that have objects in a circle with radius 3 degrees
    #CTE centered at (RA, DEC) = (163, 50) degrees
    all_first_coos, first_pick_dict, first_sweeps, sweep_stack = match_sweep(
        first_path, sweep_path, [163, 50, 3])

    #CTE Finding objects from the sweep files that are also in the FIRST data
    #CTE and have colors above 0.5 and r magnitudes below 22
    pick_ins, sweep_picks, sweep_picks_ra, sweep_picks_dec = sweep_pick(
        first_pick_dict, sweep_stack, 0.5, 22)

    #CTE Querying SDSS to get u and i magnitudes for the 
    #CTE objects found in the previous line
    print('Querying SDSS...')
    q_start = time.time()

    qry, query = sdss_ui_query(sweep_picks_ra, sweep_picks_dec)

    quer = np.array([run_q(i) for i in query])

    query_time = time.time()

    #CTE Decoding the result of the SDSS query to get ra, dec, u_mag
    #CTE and i_mag arrays.
    ra, dec, u_mag, i_mag = decode_res2(quer)

    print(f'{len(np.where(np.isnan(ra) == False)[0])} objects have been retrieved from SDSS DR9')

    print(f'---> SDSS Query took {query_time-q_start} seconds.')

    #CTE Finding the index of the brightest object in the u-band
    bright = np.nanargmin(u_mag)

    #CTE Getting the SDSS and WISE fluxes (in nanomaggies)
    #CTE for the brightest object in the u-band.
    ubrite1, ubrite1_sdss, sdss_lams, ubrite1_wise, wise_lams = get_fluxes(bright, sweep_picks, u_mag, i_mag)

    #CTE Plotting SDSS and WISE fluxes (in nanomaggies)
    #CTE as a function of wavelength (in angstroms)
    fig, ax = plt.subplots()

    ax.scatter(sdss_lams, ubrite1_sdss, color='navy')
    ax.scatter(wise_lams, ubrite1_wise, color='navy')

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlabel('$\log(\lambda \ / \ \mathrm{\AA})$', fontsize = 10)
    ax.set_ylabel('$\log(f_{\lambda} \ / \ \mathrm{nanomaggy})$', fontsize = 10)

    save_show(fig, 'query-res.png')

    #CTE Stopping the timer and printing total run time
    end = time.time()

    print('\nWith a redshift of ~1.04, the spectrum of ubrite1 shows prominent emission')
    print('at ~5700 angstroms, corresponding to the [MgII] emission line, which is typically')
    print('at ~2800 angstroms. The redshift puts this emission in the r-band, which caused')
    print('ubrite1 to be bright in the r-band (and thus, be detected in our survey')
    print('Similarly, the [CIII] 1908 emission line has been redshifted into')
    print('the u-band and has caused ubrite1 to be bright in the u-band')
    print('(and thus, be detected in our survey). These traits, paired the strong emssion in WISE')
    print('filters, makes ubrite1 very likely to be a QSO')

    print(f'\nWhole script took {end-start} seconds.')