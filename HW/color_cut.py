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
	Parameters
	----------
	tab : 

	Returns
	-------

	"""

	fluxes = ['FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FLUX_W2', 'FLUX_W3', 'FLUX_W4']
	trans = ['MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z', 'MW_TRANSMISSION_W1', 'MW_TRANSMISSION_W2', 'MW_TRANSMISSION_W3', 'MW_TRANSMISSION_W4']

	bands = [i.split('_')[-1] for i in fluxes]
	ext_corr_flux = [tab[i]/tab[j] for i, j in zip(fluxes, trans)]
	mags = f_to_m(np.stack(ext_corr_flux))
	
	mag_dict = {f"{j}_MAG":mags[i] for i, j in enumerate(bands)}

	return mag_dict

def get_cols(tab):
	"""
	Parameters
	----------
	tab : 

	Returns
	------

	"""

	mags = get_mags(tab)

	col_dict = {'r-W1':mags['R_MAG']-mags['W1_MAG'], 'g-z':mags['G_MAG']-mags['Z_MAG'], 
	'W1-W2':mags['W1_MAG']-mags['W2_MAG'], 'W3-W4':mags['W3_MAG']-mags['W4_MAG']}

	return col_dict

def flag_cuts(tab):
	"""

	Parameters
	----------

	Returns
	------

	"""

	#CTE When reading in multiple sweep files, the way I stack makes some columns read as byte strings
	#CTE rather than Unicode strings. This line checks which instance we're dealing with and 
	#CTE sets them accordingly (specifically with the TYPE column so I can index for only == 'PSF')
	#CTE regardless of whether or not it is a byte string or Unicode string
	uc_byte_check = np.char.decode(tab['TYPE']) if isinstance(tab['TYPE'][0], bytes) else tab['TYPE']

	fc = (f_to_m(tab['FLUX_R']) < 19) & (uc_byte_check == 'PSF') & (tab['ALLMASK_G'] == 0) & (
		tab['ALLMASK_R'] == 0) & (tab['ALLMASK_Z'] == 0) & (tab['WISEMASK_W1'] == 0) & (
		(tab['WISEMASK_W2'] == 0))

	return fc

def rW1_gz_cut(tab):
	"""
	
	Parameters
	----------

	Returns
	------

	"""

	cols = get_cols(tab)
	fc = flag_cuts(tab)

	x = cols['g-z'][fc]
	y = cols['r-W1'][fc]

	cc = (y > line(x, 1, -1)) & (y > line(x, -2.5, -0.5))

	return cols, fc, cc

def W1W2_W3W4_cut(tab):
	"""

	Parameters
	----------

	Returns
	------

	"""

	cols, fc, cc1 = rW1_gz_cut(tab)

	x = cols['W1-W2'][fc][cc1]
	y = cols['W3-W4'][fc][cc1]

	cc2 = y > line(x, -17, -4)

	return cols, fc, cc1, cc2

def splendid_function_v1(tab, area = None):
	"""
	
	Parameters
	----------

	Returns
	-------

	"""
	start_time = time.time()

	cols, fc, cc1, cc2 = W1W2_W3W4_cut(tab)

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