#!/usr/bin/env python3

import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord as sc
import astropy.units as u
from Tasks.week8.cross_match import *
from HW.HW3 import f_to_m

def decode_type(obj):
    """
    Given an object pulled from a Legacy Sky Survey sweep file,
    determines what type of object it is.

    Parameters
    ----------
    obj : fits.data object or astropy.Table
        The object to classify. Must be in the format of a Legacy Sky Survey
        sweep file object.

    Returns
    -------
    obj_type : str
        The given object's type
    """    

    #CTE Type dictionary with definitions pulled from sweep file documentation
    type_dict = {'PSF':'stellar', 'REX':'round exp galaxy with variable radius', 
                 'EXP':'exponential galaxy', 'DEV':'deVauc', 
                 'SER':'Sersic', 'DUP':'Gaia source fit by different model'}
    
    obj_type = type_dict[obj['TYPE']]
    
    return obj_type


def circle(x_cent, y_cent, r):
    """
    Given a circle's center and radius, generates 20
    [x.y] coordinates around a circle.

    Parameters
    ----------
    x_cent : float
        The x-coordinate of the center of the circle.

    y_cent : float
        The y-coordinate of the center of the circle.

    r : The radius of the circle

    Returns
    -------
    coos : list of arrays of floats
        A list of x and y coordinate arrays
    """    

    theta = np.linspace(0, 2*np.pi, 20)
    x = x_cent + r * np.cos(theta)
    y = y_cent + r * np.sin(theta)
        
    coos = [x, y]
            
    return coos

if __name__ == '__main__':
    #CTE Making coordinates in the two formats I will need to work with
    coos_dict = {'RA' : 188.53667, 'DEC' : 21.04572}
    coos = sc(coos_dict['RA']*u.deg, coos_dict['DEC']*u.deg)

    #CTE Finding which sweep file the given coordinate is in
    sweep = in_sweeps('/d/scratch/ASTR5160/data/legacysurvey/dr9/south/sweep/9.0', coos_dict)[0]

    #CTE Reading in sweep file data
    sweep_dat = fits.open(sweep)[1].data
    sweep_coos = sc(sweep_dat['RA']*u.deg, sweep_dat['DEC']*u.deg)

    #CTE Finding the closes object to the given coordinates
    targ = sweep_dat[(coos.separation(sweep_coos)).argmin()]
    targ_type = decode_type(targ)

    #CTE Saturated in all 3 bands
    #CTE Definitely saturated in Legacy Viewer
    #CTE According to SIMBAD, it is a blazar candidate
    targ_bitmasks = [targ['ALLMASK_G'], targ['ALLMASK_R'], targ['ALLMASK_Z']]

    #CTE Making a circle with a radius of 3 degrees centered
    #CTE on RA, DEC = 180, 30
    circle_coos = {'RA':circle(180, 30, 3)[0], 'DEC':circle(180, 30, 3)[1]}
    cen_coo = sc(180*u.deg, 30*u.deg)

    #CTE Finding which sweep files contain circle_coos and reading in their data
    circle_sweep = in_sweeps('/d/scratch/ASTR5160/data/legacysurvey/dr9/south/sweep/9.0', circle_coos)
    sweep_stack = np.hstack([fits.open(i, memmap=True)[1].data for i in circle_sweep])
    sweep_stack_coos = sc(sweep_stack['RA']*u.deg, sweep_stack['DEC']*u.deg)

    #CTE Limiting sweep data to targets within 3 degrees of cen_coo
    sweep_stack_coo_lim = sweep_stack[(cen_coo.separation(sweep_stack_coos)) < 3*u.deg]

    #CTE Further limiting data to targets with r-magnitude < 20
    sweep_r_lim = sweep_stack_coo_lim[np.where(f_to_m(sweep_stack_coo_lim['FLUX_R']) < 20)[0]]

    #CTE FURTHER limiting to targets with 'TYPE'='PSF}'
    psf_objs = sweep_r_lim[np.where(sweep_r_lim['TYPE'] == b'PSF')[0]]
    psf_coos = sc(psf_objs['RA']*u.deg, psf_objs['DEC']*u.deg)

    #CTE Reading in QSO data
    qso_file = '/d/scratch/ASTR5160/week10/qsos-ra180-dec30-rad3.fits'
    qso_dat = Table.read(qso_file)
    qso_coos = sc(qso_dat['RA']*u.deg, qso_dat['DEC']*u.deg)

    #CTE Coordinate matching QSO targets to PSF targets
    qsos = qso_dat[qso_coos.search_around_sky(psf_coos, 0.5*u.arcsec)[1]]

    print(np.size(qsos))