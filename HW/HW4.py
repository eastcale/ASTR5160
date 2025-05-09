#!/usr/bin/env python3

def splendid_function(inp, area = None):
	"""
	Given access to data that has (at least) the same columns
	as a Legacy Survey DR9 sweep file, identifies any quasar
	candidates in the data based on a series of flag and color-color cuts.
	Will also print to screen how many candidates were found in total.
	Optionally, an area (in square degrees) can be passed to determine
	how many quasar candidates have been found per square degree.
	If this is the case, the function will instead print to screen how many quasar
	candidates were found as well as how many were found per square degree.
	
	Parameters
	----------
	inp : str or astropy.table Table
		If given a string, the string must be the full directory path to a .fits file
		containing (at least) the same columns as a Legacy Survey DR9 sweep file.

		If given a table, the table must contain (at least) the same columns
		as a Legacy Survey DR9 sweep file.

	area : None or float
		Default is None. If a value is given, it must be an integer or float
		corresponding to the area in square degrees spanned by the objects
		in the given data.

	Returns
	-------
	classif : boolean array
		A boolean array corresponding to objects in the given data
		that satisfied all flag and color-color cuts and are likely quasar candidates.
	"""

	#CTE Importing packages within the function to allow
	#CTE full of the function use within an independent python shell
	print('Importing packages...')
	import time
	#CTE Starting the timer to see how long this function takes
	start_time = time.time()

	import numpy as np
	from astropy.io import fits
	from HW.HW3 import f_to_m
	from HW.color_cut import get_cols, get_mags, flag_cuts, rW1_gz_cut, W1W2_W3W4_cut
	from Tasks.week10.classif import line

	#CTE Reading in data
	print('Reading in data...')
	if type(inp) == str:
		#CTE If a directory path is given, open the data and take only relevant columns
		with fits.open(inp, memmap=True) as dat:
			fluxes = ['FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FLUX_W2', 'FLUX_W3', 'FLUX_W4']
			trans = ['MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z', 'MW_TRANSMISSION_W1', 
			'MW_TRANSMISSION_W2', 'MW_TRANSMISSION_W3', 'MW_TRANSMISSION_W4']
			flags = ['TYPE', 'ALLMASK_G', 'ALLMASK_R', 'ALLMASK_Z', 'WISEMASK_W1', 'WISEMASK_W2']

			columns = fluxes + trans + flags
			tab = {i : dat[1].data[i] for i in columns}
	else:
		#CTE If a table is already given, proceed.
		tab = inp

	#CTE Making all flag and color-color cuts.
	print('Making cuts...')
	cols, fc, cc1, cc2 = W1W2_W3W4_cut(tab)

	print('Compiling classifications...')
	#CTE Creating a boolean array of False that is the same length
	#CTE as the given data
	class_mask = np.zeros_like(fc, dtype=bool)

	#CTE Setting indices that satisfy all cuts to True
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

	#CTE Printing the function run time
	print(f'Total Runtime: {end_time - start_time}.')

	return classif

if __name__ == '__main__':
	import argparse

	parser = argparse.ArgumentParser(description = 'This module takes as an input either\
	an astropy table with (at least) the same columns as a Legacy Survey DR9 sweep file\
	or the full directory path to such a file. The user can optionally give the total area in square\
	degrees that the survey covers as well. With this information, the module goes through a series of flag\
	and color-color cuts to classify any QSO candidates within the file. The module will print to screen\
	how many QSO candidates have been found. If an area is given, the module\
	will also print how many QSO candidates are found per square degree.')
	parser.add_argument('path_or_table', 
	help='What data would you like to use? Provide either an astropy table with (at least)\
	the same columns as a Legacy Survey DR9 sweep file, or the full directory path\
	to such a file.')
	parser.add_argument('--area', '-i', default = None, type=float,
	help='What is the area of the sky (in square degrees) that your data covers?\
	The module is functional without this argument, but will not provide a QSO area density if no area is given.')
	args = parser.parse_args()

	classification = splendid_function(args.path_or_table, args.area)
	print(classification)
