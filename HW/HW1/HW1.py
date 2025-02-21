import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord as sc
from astropy.coordinates import EarthLocation as el
from astropy.coordinates import AltAz
from astropy.time import Time as ti
import calendar
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('month', help='What month are you observing in? Enter an int from 1-12.')
args = parser.parse_args()

path = open('/d/scratch/ASTR5160/week4/HW1quasarfile.txt')

dat = path.readlines()

def to_radec(dat):
	"""
	Given a list of strings in the form
	given by open(.../HW1quasarfile.txt).readlines(),
	rewrites them as proper right ascension and declination
	coordinates readable by astropy.SkyCoords.

	Parameters
	----------
	dat : list of str
		The list of strings returned from 
		open(.../HW1quasarfile.txt).readlines()

	Returns
	-------
	ra_array : array of str
		An array of right ascension coordinates
		in the form '00h00m00.00s'

	dec_array: array of str
		An array of declination coordinates in the form
		'00d00m00.0s'
	"""

	ra_list = [] #CTE Empty list to append to in the loop
	dec_list = []

	#CTE Indexing the proper slices of the string for each unit
	for i in dat:
		h = i[:2]
		m = i[2:4]
		s = i[4:9]

		sign = i[9]
		d = i[10:12]
		a_m = i[12:14]
		a_s = i[14:18]

		ra = f'{h}h{m}m{s}s' #CTE Formatting coordinates in the right form
		dec = f'{sign}{d}d{a_m}m{a_s}s'

		ra_list.append(ra)
		dec_list.append(dec)

	ra_array = np.array(ra_list)
	dec_array = np.array(dec_list)

	return ra_array, dec_array

def make_days(month, year):
	"""
	Given a month and a year, creates an array
	of dates for that month (accounts for leap years)

	Parameters
	----------
	month : int
		The month, in number form, to get the dates for
		(must be 1-12)

	year : int
		The specific year to do this function in.

	Returns
	-------
	day : array of int
		An array of dates from 1-28/29/30/31 depending
		on the month and year provided
	"""

	if month in (11, 4, 6, 9):
		day = np.arange(1, 31, 1)
	elif month in (1, 3, 5, 7, 8, 10, 12):
		day = np.arange(1, 32, 1)
	elif month == 2:
		if calendar.isleap(year) == True:
			day = np.arange(1, 30, 1)
		if calendar.isleap(year) == False:
			day = np.arange(1, 29, 1)

	return day

def max_alt(loc, coordinates, dates):
	"""
	Given a place, list of objects, and list of times,
	finds the object that has the highest altitude (and thus, lowest
	airmass) at the given location for each time.

	Parameters
	----------
	loc : astropy.EarthLocation object
		The location observations are taken at.

	coordinates : array of astropy.SkyCoord objects
		The array of coordinates for each object.

	dates : array astropy.time.Time objects
		Array of desired observation times.

	Returns
	-------
	airmass_arr : array of floats
		An array of airmasses at the altitude of the highest target in 'coordinates'
		on each time in 'dates'

	max_targ_arr : array of ints
		An array of indices corresponding to the object with the highest altitude
		for each time in 'dates'
	"""

	max_targ = []
	airmass = []

	for i in dates:
		altaz = AltAz(location = loc, obstime = i) #CTE Setting location and obseration time
		coos_altaz = coordinates.transform_to(altaz) #CTE Converting given coordinates to AltAz at the given time and place
		max_altitude_index = np.where((coos_altaz.alt.deg == max(coos_altaz.alt.deg)))[0][0] #CTE Getting the index of the object with
																												  #CTE the highest altitude
		air_mass = 1/(np.cos(np.radians(90-coos_altaz.alt.deg[max_altitude_index]))) #CTE Calculating the airmass of the object with highest alt
		max_targ.append(max_altitude_index)
		airmass.append(air_mass)

	max_targ_arr = np.array(max_targ)
	airmass_arr = np.array(airmass)

	return airmass_arr, max_targ_arr

ras, decs = to_radec(dat) #CTE Converting the coords from the .txt file into proper coordinates

coos = sc(ra = ras, dec = decs, unit=(u.hourangle, u.deg)) #CTE Generating a coord obj for each coordinate

kp = el.of_site('kpno') #CTE Setting observation location to Kitt Peak

year = ti.now().ymdhms[0] #CTE Getting the current year

month = int(args.month) #CTE Month is provided by user

day = make_days(month, year) #CTE Based on month and year, gets the proper dates

times = []

#CTE Since UTC is 7 hrs ahead of MST, setting obsveration times to 11 pm MST
#CTE is actually setting them to 6 am UTC on the next morning.
for i in day:
	if i != day[-1]:
		times.append(ti(f'{year}-{month}-{i+1}T06:00:00'))
	elif i == day[-1]:
		if month == 12:
			month2 = 1
		else:
			month2 = month + 1
		times.append(ti(f'{year}-{month2}-01T06:00:00'))

times = np.array(times) #CTE Turning 'times' list into an array to do math with
times_tab = np.array([i-7*u.hour for i in times]) #CTE Converting times to MST to print to table

airmass, max_targ = max_alt(kp, coos, times)

#CTE Evaluating each of the relevant parameters at the max altitude (lowest airmass) value.
q_coos = [dat[i] for i in max_targ]
ra = coos.ra.deg[max_targ]
dec = coos.dec.deg[max_targ]
airmasses = np.copy(airmass)

#CTE Making a data table
tab = {'Date' : times_tab, 'Quasar Coordinates (hms dms)' : q_coos, 'RA (deg)' : ra, 'DEC (deg)' : dec, 'Airmass' : airmasses}
df = pd.DataFrame(tab)
print(df)
df.to_csv(f'q-obs-{month}.csv', index=False)