#!/usr/bin/env python3

def splendid_function(inp, area = None):
	"""
	
	Parameters
	----------

	Returns
	-------

	"""

	print('Importing packages...')
	import time
	start_time = time.time()

	import numpy as np
	from astropy.io import fits
	from HW.HW3 import f_to_m
	from HW.color_cut import get_cols, get_mags, flag_cuts, rW1_gz_cut, W1W2_W3W4_cut
	from Tasks.week10.classif import line

	print('Reading in data...')
	if type(inp) == str:
		with fits.open(inp, memmap=True) as dat:
			fluxes = ['FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FLUX_W2', 'FLUX_W3', 'FLUX_W4']
			trans = ['MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z', 'MW_TRANSMISSION_W1', 
			'MW_TRANSMISSION_W2', 'MW_TRANSMISSION_W3', 'MW_TRANSMISSION_W4']
			flags = ['TYPE', 'ALLMASK_G', 'ALLMASK_R', 'ALLMASK_Z', 'WISEMASK_W1', 'WISEMASK_W2']

			columns = fluxes + trans + flags
			tab = {i : dat[1].data[i] for i in columns}
	else:
		tab = inp

	print('Making cuts...')
	cols, fc, cc1, cc2 = W1W2_W3W4_cut(tab)

	print('Compiling classifications...')
	class_mask = np.zeros_like(fc, dtype=bool)
	fc_in = np.where(fc)[0]
	cc1_in = fc_in[cc1]
	cc2_in = cc1_in[cc2]
	class_mask[cc2_in] = True

	classif = class_mask

	if area != None:
		print(f'{sum(classif)} QSO candidates found ({sum(classif)/area} per square degree).')
	else:
		print(f'{sum(classif)} QSO candidates found.')

	end_time = time.time()

	print(f'Total Runtime: {end_time - start_time}.')

	return classif

if __name__ == '__main__':
	import argparse

	parser = argparse.ArgumentParser(description = 'This module takes as an input either\
	an astropy table with (at least) the same columns as a Legacy Survey DR9 sweep file\
	or the full directory path to such a file. The user can optionally give the total area in square\
	degrees that the survey covers as well. With this information, the module goes through a series of flag\
	and color-color cuts to classify any QSO candidates within the file. If an area is given, the module\
	will print how many QSO candidates are found per square degree.')
	parser.add_argument('path_or_table', 
	help='What data would you like to use? Provide either an astropy table with (at least)\
	the same columns as a Legacy Survey DR9 sweep file, or the full directory path\
	to such a file.')
	parser.add_argument('--area', '-i', default = None, type=float,
	help='What is the area of the sky (in square degrees) that your data covers?\
	The module is functional without this argument, but will not provide a QSO area density if no area is given.')
	args = parser.parse_args()

	classification = splendid_function(args.path_or_table, args.area)
