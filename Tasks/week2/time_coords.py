from astropy import units as u
from astropy.coordinates import SkyCoord as sc
from astropy.coordinates import EarthLocation as el
from astropy.coordinates import Angle as ang
from astropy.coordinates import AltAz
from astropy.time import Time as ti
import numpy as np
import math

def to_decdeg(d, m, s):
	"""
	Given an angle in degrees, arcminutes,
	and arcseconds, converts to decimal degrees.

	Parameters
	----------
	d : float
		The degree value for the given angle.
	
	m : float
		The arcminutes value for the given angle.

	s : float
		The arcseconds value for the given angle.

	Returns
	-------
	ddeg : float
		The resulting angle in decimal degrees.
	"""
	
	ddeg = (d)+(m/60)+(s/3600)

	return ddeg

def time_to(current_time, time):
	"""
	Given the current time,
	finds the number of hours until a desired time

	Parameters
	----------
	current_time : astropy.time.Time
		The current time given by 
		astropy.time.Time.now()

	time : list of floats
		The desired time in the format 
		[hr, min, sec]

	Returns
	-------
	tot : float
		The number of hours until the desired time
	
	"""
	hr_to = (time[0] - current_time.ymdhms[3])
	min_to = (time[1] - current_time.ymdhms[4])
	sec_to = (time[2] - current_time.ymdhms[5])
	
	tot = ((hr_to) + (min_to/60) + (sec_to/3600))*u.hour
	
	return tot

coo = sc('05 55 10.3054 +07 24 25.4304', unit=(u.hourangle, u.deg)) #CTE Betelgeuse coordinates from Simbad

ra_dec_hand = 15*to_decdeg(5, 55, 10.3054) #CTE Converting to decimal degrees by hand to compare
decl_dec_hand = to_decdeg(7, 24, 25.4304)

today = ti.now() #CTE Today's date and time

days = np.arange(math.floor(today.mjd) - 5, math.floor(today.mjd) + 5, 1) #CTE Generating the MJDs for the 5 days before today
									  #CTE and the 5 days after today

wiro = el(lat = to_decdeg(41, 5, 49)*u.deg, lon = -1*to_decdeg(105, 58, 33)*u.deg, height = 2943) #CTE Geographical 													  #CTE coordinates of WIRO

obj = sc('12 00 00 +30 00 00', unit=(u.hourangle, u.deg)) #CTE Our given object's coordinates

utc_offset = -7 * u.hour #CTE UTC is 7 hours ahead of MST

current_time = ti.now() - utc_offset #CTE The current time in MST

time_to_11 = time_to(current_time, [23, 0, 0]) #CTE Finding the number hours from now to 11 pm

obs_time = np.array([current_time + time_to_11, current_time + time_to_11 + 30*u.day]) #CTE Setting observation times for 
										       #CTE 11 pm tonight as well as 11 pm
										       #CTE a month from now.

obj_altaz = obj.transform_to(AltAz(obstime=obs_time, location=wiro)) #CTE Finding Alt/Az coords. for our object

airmass = obj_altaz.secz #CTE .secz calculates the airmass where z is 90-alt

#CTE A bunch of print statements to summarize results
print("RA	RA Unit		DEC	DEC Unit \n-----------------------------------------------\n{}	Hours		{}	Decimal Degrees\n{}	Dec. Deg.	{}	Dec. Deg\n-----------------------------------------------\n                  By Hand\n-----------------------------------------------\n{}	Dec. Deg.	{}	Dec. Deg.\n-----------------------------------------------\nBy Hand Checks Out! :)".format(round(coo.ra.hour, 4), round(coo.dec.degree, 4), round(coo.ra.degree, 4), round(coo.dec.degree, 4), round(ra_dec_hand, 4), round(decl_dec_hand, 4)))

print("\nToday is {}\nThe Julian Date is {}\nThe Modified Julian Date is {}\nMJD+2400000.5={}; Checks out!".format(today, today.jd, today.mjd, today.mjd+2400000.5))

print("\nThe MJDs for 5 days before today to 5 days after today are:\n{}".format(days))

print("\nThe airmass of our source on {} is {}.\nThe airmass of our source on {} is {}.".format(obs_time[0], airmass[0], obs_time[1], airmass[1]))
