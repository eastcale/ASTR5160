#!/usr/bin/env python3
from Tasks.week10.classif import *
from astropy.table import QTable, Table, Column, vstack
from astropy.coordinates import SkyCoord as sc
import astropy.units as u
import numpy as np
from Modules.plots import save_show
from HW.HW3 import f_to_m
from scipy import stats
import time

def get_mags(tab):
	"""
	Given a table that has (at least) the same columns as
	a Legacy Survey DR9 sweep file, collects the following
	(corrected for MW transmission) magnitudes in a dictionary:
	g, r, z, W1, W2, W3, and W4.

	Parameters
	----------
	tab : astropy.table Table or similar
		A table of some kind, that has (at least) the same columns
		as a Legacy Survey DR9 sweep file to collect magnitude information.

	Returns
	-------
	mag_dict : dictionary
		A dictionary containing transmission corrected magnitudes
		for every object in the given table.
	"""

	#CTE Collecting flux and transmission column names
	fluxes = ['FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FLUX_W2', 'FLUX_W3', 'FLUX_W4']
	trans = ['MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z', 'MW_TRANSMISSION_W1', 'MW_TRANSMISSION_W2', 'MW_TRANSMISSION_W3', 'MW_TRANSMISSION_W4']

	#CTE Creating strings with just the band names (G, R, Z etc.)
	bands = [i.split('_')[-1] for i in fluxes]

	#CTE Correcting fluxes (in nanomaggies) for MW transmission
	corr_flux = [tab[i]/tab[j] for i, j in zip(fluxes, trans)]

	#CTE Converting fluxes to magnitudes
	mags = f_to_m(np.stack(corr_flux))
	
	#CTE Compiling results in a dictionary
	mag_dict = {f"{j}_MAG":mags[i] for i, j in enumerate(bands)}

	return mag_dict

def get_cols(tab):
	"""
	Given a table that has (at least) the same columns
	as a Legacy Survey DR9 sweep file, collects predetermined
	colors in a dictionary

	Parameters
	----------
	tab : astropy.table Table or similar
		A table of some kind, that has (at least) the same columns
		as a Legacy Survey DR9 sweep file to collect colors from.

	Returns
	------
	col_dict : dictionary
		A dictionary containing predetermined colors for every
		object in the given table.
	"""

	#CTE Getting magnitudes from the table
	mags = get_mags(tab)

	#CTE Collecting colors
	col_dict = {'r-W1':mags['R_MAG']-mags['W1_MAG'], 'g-z':mags['G_MAG']-mags['Z_MAG'], 
	'W1-W2':mags['W1_MAG']-mags['W2_MAG'], 'W3-W4':mags['W3_MAG']-mags['W4_MAG']}

	return col_dict

def flag_cuts(tab):
	"""
	Given a table that has (at least) the same columns as
	a Legacy Survey DR9 sweep file, goes through a series of flag
	(or flag adjacent) cuts that have been pre-selected to be quasar markers.

	Parameters
	----------
	tab : astropy.table Table or similar
		A table of some kind, that has (at least) the same columns
		as a Legacy Survey DR9 sweep file to begin cutting.

	Returns
	------
	fc : boolean array
		A boolean array corresponding to objects that 
		satisfy the selected quasar markers.

	"""

	#CTE When reading in multiple sweep files, the way I stack makes some columns read as byte strings
	#CTE rather than Unicode strings. This line checks which instance we're dealing with and 
	#CTE sets them accordingly (specifically with the TYPE column so I can index for only == 'PSF')
	#CTE regardless of whether or not it is a byte string or Unicode string
	uc_byte_check = np.char.decode(tab['TYPE']) if isinstance(tab['TYPE'][0], bytes) else tab['TYPE']

	#CTE Limiting to r magnitudes brighter than 19
	#CTE Limiting to only PSF objects
	#CTE Limiting ALLMASK and WISEMASK to 0 (QSOs tend to have these values)
	fc = (f_to_m(tab['FLUX_R']) < 19) & (uc_byte_check == 'PSF') & (tab['ALLMASK_G'] == 0) & (
		tab['ALLMASK_R'] == 0) & (tab['ALLMASK_Z'] == 0) & (tab['WISEMASK_W1'] == 0) & (
		(tab['WISEMASK_W2'] == 0))

	return fc

def rW1_gz_cut(tab):
	"""
	Given a table that has (at least) the same columns as
	a Legacy Survey DR9 sweep file, makes two color-color cuts
	in r-W1 vs. g-z space based on pre-screened conditions
	for quasar candidates. 
	
	Parameters
	----------
	tab : astropy.table Table or similar
		A table of some kind, that has (at least) the same columns
		as a Legacy Survey DR9 sweep file to cut for quasars,

	Returns
	------
	cols : dictionary
		A dictionary of colors determined by the get_cols() function
		for every object in 'tab'

	fc : boolean array
		A boolean array corresponding to objects in 'tab'
		that satisfy the flag cuts determined by the
		flag_cuts() function.

	cc : boolean array
		A boolean array corresponding to the subset of
		'True' values in 'fc' that satisfy the cuts
		in r-W1 vs. g-z space.
	"""

	#CTE Getting colors for every object
	cols = get_cols(tab)

	#CTE Making flag cuts for every object
	fc = flag_cuts(tab)

	#CTE Separating data in to x and y axes
	x = cols['g-z'][fc]
	y = cols['r-W1'][fc]

	#CTE Comparing data to two predetermined line cuts
	#CTE and keeping all objects that are above BOTH lines
	cc = (y > line(x, 1, -1)) & (y > line(x, -2.5, -0.5))

	return cols, fc, cc

def W1W2_W3W4_cut(tab):
	"""
	Given a table that has (at least) the same columns as
	a Legacy Survey DR9 sweep file, makes a color-color
	cut in W3-W4 vs. W1-W2 space based on pre-screened
	conditions for quasar candidates

	Parameters
	----------
	tab : astropy.table Table or similar
		A table of some kind, that has (at least) the same columns as
		a Legacy Survey DR9 sweep file to make cuts on.

	Returns
	------
	cols : dictionary
		A dictionary of colors determined by the get_cols() function
		for every object in 'tab'

	fc : boolean array
		A boolean array corresponding to objects that satisfy
		the flag cuts determined by the flag_cuts() function

	cc1 : boolean array
		A boolean array corresponding to the subset of 'True'
		values in 'fc' that satisfy the cuts in r-W1 vs. g-z space
		determined by the rW1_gz_cut() function.

	cc2 : boolean array
		A boolean array corresponding to the subset of 'True'
		values in 'cc1' that satisfy the cuts in W3-W4 vs. W1-W2
		space. 
	"""

	#CTE Getting colors, flag cuts array, and the initial color-color cuts
	cols, fc, cc1 = rW1_gz_cut(tab)

	#CTE Separating data into x and y axes
	x = cols['W1-W2'][fc][cc1]
	y = cols['W3-W4'][fc][cc1]

	#CTE Comparing data to a predetermined line cut
	#CTE and keeping objects that are above the line.
	cc2 = y > line(x, -17, -4)

	return cols, fc, cc1, cc2

def splendid_function_v1(tab, area = None):
	"""
	Given a table that has (at least) the same columns
	as a Legacy Survey DR9 sweep file, identifies any quasar
	candidates in the table based on a series of flag and color-color cuts.
	Will also print to screen how many candidates were found in total.
	Optionally, an area (in square degrees) can be passed to determine
	how many quasar candidates have been found per square degree.
	If this is the case, the function will instead print to screen how many
	quasar candidaes were found as well as how many were found per square
	degree.
	
	Parameters
	----------
	tab : astropy.table Table or similar
		A table of some kind, that has (at least) the same columns
		as a Legacy Survey DR9 sweep file to identify quasar candidates in.

	area : None or int or float
		Default is None. If a value is given, it must be an integer or float
		corresponding to the area in square degrees spanned by the objects
		in 'tab'.

	Returns
	-------
	classif : boolean array
		A boolean array corresponding to object in 'tab' that satisfied
		all flag and color-color cuts and are likely quasar candidates.
	"""

	#CTE Starting the clock to see how long this function takes
	start_time = time.time()

	#CTE Going through all flag and color-color cuts
	cols, fc, cc1, cc2 = W1W2_W3W4_cut(tab)

	#CTE Creating a boolean array of False that is the same
	#CTE length as 'tab'
	class_mask = np.zeros_like(fc, dtype=bool)

	#CTE Setting indices where all cuts were satisfied to True
	fc_in = np.where(fc)[0]
	cc1_in = fc_in[cc1]
	cc2_in = cc1_in[cc2]
	class_mask[cc2_in] = True

	classif = class_mask

	#CTE Printing results
	if area != None:
		print(f'{sum(classif)} QSO candidates found ({sum(classif)/area} per square degree).')
	else:
		print(f'{sum(classif)} QSO candidates found.')

	#CTE Stopping the clock
	end_time = time.time()

	#CTE Printing function run time.
	print(f'Total Runtime: {end_time - start_time}.')

	return classif

if __name__ == '__main__':
	#CTE Setting file path for qso file.
	qsos = Table.read('/d/scratch/ASTR5160/week10/qsos-ra180-dec30-rad3.fits')
	qso_coos = sc(qsos['RA']*u.deg, qsos['DEC']*u.deg)

	#CTE Setting the directory for sweep files and making a sweep stack
	sweeps_dir = '/d/scratch/ASTR5160/data/legacysurvey/dr9/south/sweep/9.0'
	qso_insw = in_sweeps(sweeps_dir, qsos)
	print('Making sweep stack...')
	sw_stack1 = np.hstack([fits.open(i, memmap=True)[1].data for i in qso_insw])
	print('Sweep stack made...')

	#CTE Running the cuts function to test it
	test = splendid_function_v1(sw_stack1, 173)