#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

def mk_rect(x, y, ax, col, lab):
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

	ax.vlines(x[0], y[0], y[1], color=col, label = lab)
	ax.vlines(x[1], y[0], y[1], color=col)

	ax.hlines(y[0], x[0], x[1], color=col)
	ax.hlines(y[1], x[0], x[1], color=col)

def sph_area(ra, dec, no_dot = None):
    """
    Given rectangular bounds, calculates the number of square degrees
    subtended by the rectangle. Optionally, if no_dot != None, also
    calculates the expected number of dots in the rectangle if there are
    no_dot total random points evenly distributed across the sphere.

    Parameters
    ----------
    ra : list or array of floats
        The ra coordinates for the vertices of the rectangle from left to right
        with no duplicates. e.g. if the corners are 
        ([x1, y1], [x2, y1], [x1, y2], [x2, y2]), then 'dec' is just [x1, x2]

    dec : list or array of floats
        The dec coordinates for the vertices of the rectangle from bottom to top
        with no duplicates. e.g. if the corners are 
        ([x1, y1], [x2, y1], [x1, y2], [x2, y2]), then 'dec' is just [y1, y2]

    n_dot : int
        The total number, if any, of evenly distributed across the sphere.

    Returns
    -------
    area : float
        The area, in square degrees, subtended by the rectangle.

    dots_area : float
        The expected number of dots, if any, that lie within the given rectangle.
    """

    #CTE Find the area of the given box
    area = (180 / np.pi ) ** 2 * (ra[1] - ra[0]) * (np.sin(dec[1]) - np.sin(dec[0]))

    #CTE Calculate the total area of a sphere
    sphere = 4 * np.pi * (180 / np.pi) ** 2

    #CTE Find the number of dots expected to be within the box
    dots_area = ( area / sphere ) * no_dot

    return area, dots_area

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

def ran_sphere(n_dots):
    """
    Generates n_dots randomly generated points that are evenly distributed in area
    across a sphere.

    Parameters
    ----------
    n_dots : int
        The number of dots to generate.

    Returns
    -------
    ra : array of floats
        The randomly generated right ascension coordinates.

    dec : array of floats
        The randomly generated declination coordinates.

    coos_dict : dictionary
        ra and dec combined into a single dictionary with keys "RA" and "DEC".
    """
    from numpy.random import random

    #CTE Randomly generate RA and DEC evenly distributed in area
    ra = 2*np.pi*(random(n_dots)-0.5)
    dec = np.arcsin(1.-random(n_dots)*2.)

    #CTE Stitch the coordinates into an array
    coos_dict = {"RA":np.degrees(ra), "DEC":np.degrees(dec)}

    return ra, dec, coos_dict

def pts_in_box(coos_dict, boxes):
    """
    Given a dictionary of coordinates with atleast the keys
    "RA" and "DEC", finds what points, if any, are within given boxes.

    Parameters
    ----------
    coos_dict : dictionary
        A dictionary containing atleast the keys "RA" and "DEC"
        with right ascension and declination coordinates, respectively.
    
    boxes : list of lists or list of arrays
        A list with each term containing a 4-entry list or array of coordinates 
        [ramin, ramax, decmin, decmax] forming the
        edges of a box in RA/Dec (degrees).

    Returns
    -------
    in_boxes : list of lists
        A list containing True/False for every point that lies within each of the boxes
        (True corresponding to a point that is within the box)

    ra_box : list of lists
        A list with each term corresponding to the right ascension coordinates that lie within
        each of the boxes.

    dec_box : list of lists
        A list with each term corresponding to the declination coordinates that lie within
        each of the boxes.
    """

    #CTE Find which objects are in each of the given boxes
    in_boxes = [is_in_box(coos_dict, i) for i in boxes]

    #CTE Finding the coordinates for each object in each box
    ra_box = [ra[i] for i in in_boxes]
    dec_box = [dec[i] for i in in_boxes]

    return in_boxes, ra_box, dec_box

def aitoff(ra = None, dec = None):
    """
    Initializes an Aitoff Projection plot in the preferred style. Optionally, plots
    RA and DEC if any are given.

    Parameters
    ----------
    ra : list or array of floats or None
        The right ascension coordinates to plot. If none are given, none are plotted

    dec : list or array of floats or None
        The declination coordinates to plot. If none are given, none are plotted

    Returns
    -------
    fig : type ?
        The figure created by this function.

    ax : type ? 
        The axis on which plotting is done.
    """

    x_lab = [str(i)+'h' for i in [14, 16, 18, 20, 22, 0, 2, 4, 6, 8, 10]]
    y_lab = [str(i)+'$^{\circ}$' for i in np.arange(-75, 90, 15)]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection = 'aitoff')

    ax.scatter(ra, dec, s = 1, color = 'tab:gray')

    ax.set_xticklabels(x_lab, weight=800)
    ax.set_yticklabels(y_lab, weight=800)

    ax.grid(color='black', alpha=0.5)

    return fig, ax

if __name__ == '__main__':
    from PIL import Image as im
    import argparse

    #CTE Setting arguments for argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('path', 
        help='Where do you want to save the file? Include the full path, including file name') 

    args = parser.parse_args()

    #C.TE Setting boxes
    a_min = np.radians(4*[-20])
    a_max = np.radians(4*[20])

    d_min = np.radians([-55, -25, 20, 50])
    d_max = np.radians([-30, 0, 45, 75])

    a_bound = [[i, j] for i, j in zip(a_min, a_max)]
    d_bound = [[i, j] for i, j in zip(d_min, d_max)]

    boxes = [np.degrees([i+j][0]) for i, j in zip(a_bound, d_bound)]

    #CTE Finding area of boxes and expected no of dots in each box
    n_dots_tot = 10000
    areas = [sph_area(i, j, n_dots_tot) for i, j in zip(a_bound, d_bound)]
    exp_dot_areas = [i[1] for i in areas]

    #CTE Randomly populating a sphere
    ra, dec, coos_dict = ran_sphere(n_dots_tot)

    #CTE Finding points in each box
    in_boxes, ra_box, dec_box = pts_in_box(coos_dict, boxes)

    #CTE Actual number of dots in each box
    dots_in_box = [np.size(i) for i in ra_box]

    #CTE Plotting
    cols = ['red', 'orange', 'blue', 'green']
    a_labs = [f'{round(i[0], 2)}'+' sq. deg.'+f' Expected = {round(i[1], 1)}'+f' Actual = {j}'
            for i, j in zip(areas, dots_in_box)]

    fig, ax = aitoff(ra, dec)

    [mk_rect([a_min[i], a_max[i]], [d_min[i], d_max[i]], ax, j, k) for i, j, k in zip(range(4), cols, a_labs)]

    [ax.scatter(ra[i], dec[i], color=j, s = 1) for i, j in zip(in_boxes, cols)]

    ax.legend(bbox_to_anchor = (1, 0), fontsize = 8)
    
    file = args.path
    fig.savefig(file)
    print(f'File saved at {file}')
    im.open(file).show()