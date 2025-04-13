#!/usr/bin/env python3

import numpy as np
from sklearn.datasets import load_iris
import matplotlib.pyplot as plt
from sklearn import neighbors
import sys
sys.path.insert(0, "/d/users/caleb/Py_Modules")
from plots import min_max, tick_array, ticks, save_show, ticks2
from astropy.table import Table, Column, hstack
from astropy.io import fits
import os
import glob
from astropy.coordinates import SkyCoord as sc
import astropy.units as u
from Tasks.week8.cross_match import *
from HW.HW3 import f_to_m
from Tasks.week12.bad_data import circle

def qso_match(ra, dec, radius, seed=5341):
    """
    Given a circle's center and radius, cross matches
    sweep files and a QSO data file to find targets within
    the circle

    Parameters
    ----------
    ra : float
        The right ascension coordinate for the center of the circle

    dec : float
        The declination coordinate for the center of the circle

    radius : float
        The radius of the circle

    seed : int
        The seed to set the randomly pulled subset of objects to be classed
        as stars. Allows for reproduceability

    Returns
    -------
    psf_qsos : sweep file data array
        Sweep file data corresponding to confirmed QSOs in the
        given circle

    psf_star : sweep file data array
        A randomly selected subset of a sweep file data array
        that defines a 'star'.
    """

    #CTE Defining a circle
    circle_coos = {'RA':circle(ra, dec, radius)[0], 'DEC':circle(ra, dec, radius)[1]}
    coos = sc(ra*u.deg, dec*u.deg)

    #CTE Finding which sweep files that circle lies in and reading in their data
    sweep = in_sweeps('/d/scratch/ASTR5160/data/legacysurvey/dr9/south/sweep/9.0', circle_coos)
    sweep_stack = np.hstack([fits.open(i, memmap=True)[1].data for i in sweep])
    sweep_stack_coos = sc(sweep_stack['RA']*u.deg, sweep_stack['DEC']*u.deg)

    #CTE Limiting sweep data to objects that lie within the defined circle
    sweep_stack_coo_lim = sweep_stack[(coos.separation(sweep_stack_coos)) < radius*u.deg]

    #CTE Further limiting sweep data to objects with r-magnitude < 20
    sweep_r_lim = sweep_stack_coo_lim[np.where(f_to_m(sweep_stack_coo_lim['FLUX_R']) < 20)[0]]

    #CTE Taking subset of those objects that have 'TYPE' == 'PSF'
    psf_objs = sweep_r_lim[np.where(sweep_r_lim['TYPE'] == b'PSF')[0]]
    psf_coos = sc(psf_objs['RA']*u.deg, psf_objs['DEC']*u.deg)

    #CTE Reading in QSO data
    qso_file = '/d/scratch/ASTR5160/week10/qsos-ra180-dec30-rad3.fits'
    qso_dat = Table.read(qso_file)
    qso_coos = sc(qso_dat['RA']*u.deg, qso_dat['DEC']*u.deg)

    #CTE Coordinate matching QSO objects from PSF objects
    qso_match = qso_coos.search_around_sky(psf_coos, 0.5*u.arcsec)
    psf_qsos = psf_objs[qso_coos.search_around_sky(psf_coos, 0.5*u.arcsec)[0]]
    
    #CTE Pulling random subset of PSF objects to be a star archetype
    #CTE Setting the seed so that there is no risk of having nans (which kNN has
    #CTE issues dealing with)
    np.random.seed(seed)
    psf_star_ins = np.random.choice(np.arange(0, np.size(psf_objs)+1, 1), 275)
    psf_star = psf_objs[psf_star_ins]
    
    return psf_qsos, psf_star

def get_colors(tab):
    """
    Finds g-z and r-W1 colors from a given table

    Parameters
    ----------
    tab : sweep file data array or Table
        The table to get colors from. Must have atleast the following columns
        (all in units of nanomaggies):
        'FLUX_R', 'FLUX_G', 'FLUX_Z', 'FLUX_W1'

    Returns
    -------
    gz : list of floats
        The g-z colors for the given table objects

    rw1: list of floats
        The r-w1 colors for the given table objects
    """    

    r = f_to_m(tab['FLUX_R'])
    g = f_to_m(tab['FLUX_G'])
    z = f_to_m(tab['FLUX_Z'])
    w1 = f_to_m(tab['FLUX_W1'])
    
    gz = g - z
    rw1 = r - w1
    
    return gz, rw1

def build_archetypes(dat1, dat2):
    """
    Given two lists of data, reshapes them into a full data array
    readable by a kNN algorithm. Assigns an identity marker of "0" to dat1
    objects and "1" to dat2 objects.

    Parameters
    ----------
    dat1 : list of floats
        One list of data to go into the algorithm

    dat2 : list of floats
        The second list of data to go into the algorithm

    Returns
    -------
    all_dat : array of floats
        A full array of all the given data, reshaped into a 2x(len(dat1)+len(dat2))
        array

    dat_identites : array of ints
        An array of ones and zeros assigning identity markers to each object.
    """ 

    all_dat = np.reshape(np.concatenate((dat1, dat2)), 
                         (2, len(dat1[0])+len(dat2[0]))).T

    dat_identities = np.array(len(dat1[0])*[0] + len(dat2[0])*[1])
    
    return all_dat, dat_identities

def mk_mock(n, dat):
    """
    Creates n randomly generated mock data points 
    between the minimum and maximum of a data set

    Parameters
    ----------
    n : int
        The number of random points to generate

    dat : array of floats
        The model data to base the mock data off of.

    Returns
    -------
    mock_dat : array of floats
        The generated mock data.
    """

    mock_dat = []

    for i in range(2):
        col_min = np.min(dat[..., i])
        col_max = np.max(dat[..., i])
        mock_meas = np.random.random(n)*(col_max - col_min) + col_min
        mock_dat.append(mock_meas)
        
    mock_dat = np.reshape(mock_dat, (2, n)).T
        
    return mock_dat

def qso_ml_wrap():
    """
    Applies a kNN machine learning algorithm to
    a set of QSO and Star data, generating
    100000 random mock points to classify objects based
    on their g-z and r-W1 colors to test the algorithm.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """ 
    
    #CTE Getting QSO and Star archetypes established
    psf_qsos, psf_star = qso_match(180, 30, 3, 5341)
    
    psf_star_gz, psf_star_rw1 = get_colors(psf_star)
    psf_qsos_gz, psf_qsos_rw1 = get_colors(psf_qsos)

    star_dat = [psf_star_gz, psf_star_rw1]
    qso_dat = [psf_qsos_gz, psf_qsos_rw1]

    all_dat, dat_identities = build_archetypes(qso_dat, star_dat)

    #CTE Initiating kNN
    knn = neighbors.KNeighborsClassifier(n_neighbors=1)
    knn.fit(all_dat, dat_identities)

    #CTE Creating 100000 mock points
    mock_dat = mk_mock(100000, all_dat)

    #CTE Classifying the mock points
    mock_dat_class = knn.predict(mock_dat)

    #CTE Plotting mock data set g-z colors as a function of
    #CTE r-W1 colors and color coding based on classification
    colors = ['navy', 'red']
    names = ['Quasar', 'Star']

    fig, ax = plt.subplots()
    
    for i in range(2):
        datclass = mock_dat_class == i
        ax.scatter(mock_dat[datclass, 0], 
                   mock_dat[datclass, 1], label=names[i], color=colors[i], s = 1)

    ax.set_xlabel('(r - W1)', fontsize = 15)
    ax.set_ylabel('(g - z)', fontsize = 15)

    ax.legend(loc='best')
    
    save_show(fig, 'qso-ml.png')

if __name__ == '__main__':
    #CTE Until another comment says otherwise, this code
    #CTE was adapted from code provided by Adam Myers (ADM)

    iris = load_iris()
    colnames = "sepal_length sepal_width petal_length petal_width" 

    fig, ax = plt.subplots(1, 1, figsize=(8,6))
    for i in range(3):
        target_class = iris.target == i
        ax.scatter(iris.data[target_class, 0], iris.data[target_class, 1], s=90, label=iris.target_names[i])
        ax.tick_params(labelsize=14)
        ax.set_xlabel("Sepal length (cm)", size=14)
        ax.set_ylabel("Sepal width (cm)", size=14)
        ax.legend(prop={'size': 14})

    save_show(fig, 'sepal-lw.png')

    # ADM use the k-nearest neighbors algorithm (k-NN)
    # ADM with distances only to the nearest neighbor.
    knn = neighbors.KNeighborsClassifier(n_neighbors=1)
    knn.fit(iris.data, iris.target)
    mock_data = [5, 4, 1, 0]
    mock_data = [6, 3, 4, 1]

    # ADM let's map out the entire sepal_length/sepal_width space!
    # ADM (I'm restricting to just sepal_length and sepal_width as
    # ADM it's easier to picture a 2-D space than a 4-D space).
    n = 100000
    mock_data = []
    # ADM normally I don't condone appending to empty lists, but here
    # ADM I want to explicitly illustrate which columns I'm working
    # ADM on. This won't be a slow append, as it's only two columns.
    for i in range(2):
        col_min = np.min(iris.data[..., i])
        col_max = np.max(iris.data[..., i])
        # ADM generate random points in the space corresponding to the
        # ADM iris measurement of interest.
        mock_meas = np.random.random(n)*(col_max - col_min) + col_min
        mock_data.append(mock_meas)
    # ADM we now have a list of n*2 measurements, 
    # ADM but we want an array of 2 columns and n rows.
    mock_data = np.reshape(mock_data, (2, n)).T

    # ADM classify using the k-NN "black box"
    # ADM trained on the real-world iris data.
    # ADM but only use 2 columns as a simple illustration.
    knn = neighbors.KNeighborsClassifier(n_neighbors=1)
    knn.fit(iris.data[..., :2], iris.target)
    mock_target_class = knn.predict(mock_data)
    # ADM again, this for loop isn't strictly necessary, but
    # ADM it's a clear way to print the information to screen.

    # ADM let's plot the sepal_length, sepal_width space of the k-NN classifier.
    fig, ax = plt.subplots(1, 1, figsize=(8,6))
    for i in range(3):
        target_class = mock_target_class == i
        ax.scatter(mock_data[target_class, 0], mock_data[target_class, 1], s=10, label=iris.target_names[i])
        ax.tick_params(labelsize=14)
        ax.set_xlabel("Sepal length (cm)", size=14)
        ax.set_ylabel("Sepal width (cm)", size=14)
        ax.legend(prop={'size': 14})
    save_show(fig, 'sepal-knn.png')

    virginica_perc = np.size(np.where(mock_target_class == 2))/np.size(mock_target_class) * 100
    print(f'{virginica_perc:.3f}% of the mock irises were classified as virginica iris by the k-NN algorithm')

    #CTE This is the end of the adapted code.
    qso_ml_wrap()