import numpy as np

def to_cart(ra, dec):
    """
    Converts given Equatorial Coordinates, converts to Cartesian Coordinates

    Parameters
    ----------
    ra : int or float
    	The right ascension (in hours) for the Equatorial coordinates to be converted

    dec: int or float
    	The right ascension (in degrees) for the Equatorial coordinates to be converted

    Returns
    -------
	xyz : list of floats
		The Cartesian Coordinates corresponding to the given Equatorial Coordinates
    """

    ra_rad = np.radians(15 * ra) #CTE Converting the RA to degrees, then to radians
    dec_rad = np.radians(dec) #CTE Converting the DEC to radians

    x = np.cos(ra_rad)*np.cos(dec_rad)
    y = np.sin(ra_rad)*np.cos(dec_rad)
    z = np.sin(dec_rad)

    xyz = [x, y, z]
    
    return xyz

def ra_cap(ra):
	"""
    Creates a vector 4-array (x, y, z, 1-h) for a spherical
    cap bounded by the given right ascension.

    Parameters
    ----------
    ra : int or float
    	The right ascension (in hours) bounding the spherical cap.

    Returns
    -------
	ra_bound : list of floats
		The resulting vector 4-array (x, y, z, 1-h) describing the spherical cap
		bounded by the given right ascension.
    """

	dec = 0

	xyz = to_cart(ra+90/15, dec) #CTE RA in hours + 90 deg / 15 hr/deg

	ra_bound = [xyz[0], xyz[1], xyz[2], 1]

	return ra_bound

def dec_cap(dec):
	"""
    Creates a vector 4-array (x, y, z, 1-h) for a spherical
    cap bounded by the given declination.

    Parameters
    ----------
    dec : int or float
    	The declination (in degrees) bounding the spherical cap.

    Returns
    -------
	dec_bound : list of floats
		The resulting vector 4-array (x, y, z, 1-h) describing the spherical cap
		bounded by the given declination.
    """

	ra = 0
	dec1 = 90

	xyz = to_cart(ra, dec1)

	h = 1-np.sin(np.radians(dec))

	dec_bound = [xyz[0], xyz[1], xyz[2], h]

	return dec_bound

def sph_cap(ra, dec, theta):
	"""
    Creates a vector 4-array (x, y, z, 1-h) for a spherical
    cap representing a circular field drawn on the surface
    of a sphere at the given Equitorial coordinates, with a
    radius of theta.

    Parameters
    ----------
    ra : int or float
    	The right ascension (in hours) where the cap is.

    dec : int or float
    	The declination (in hours) where the cap is.

    theta : int or float
    	The radius (in degrees) of the circular field.

    Returns
    -------
	cap : list of floats
		The resulting vector 4-array (x, y, z, 1-h) describing the spherical cap.
    """

	xyz = to_cart(ra, dec)

	h = 1-np.cos(np.radians(theta))

	cap = [xyz[0], xyz[1], xyz[2], h]

	return cap

def cap_to_str(cap, digits):
	"""
    Converts a given vector 4-array to a string with rounded terms.

    Parameters
    ----------
    cap : list of floats
    	The vector 4-array to convert into a string.

    digits : int
		The number of decimal places to round each term of the array to

    Returns
    -------
	capstr : str
		A string with each term of 'cap' rounded to the desired number of
		digits and separated by spaces.
    """

	capstr = f'{round(cap[0], digits)} {round(cap[1], digits)} {round(cap[2], digits)} {round(cap[3], digits)}'

	return capstr

racap = ra_cap(5)
deccap = dec_cap(36) #CTE Result is [6.123233995736766e-17, 0.0, 1.0, 0.41221474770752686] but am told it should be
							 #CTE [0, 0, 1, 0.4122...]. Rounding?

deccap[0] = 0.0 #CTE Manually setting deccap[0] to zero to match the notes

cap = sph_cap(5, 36, 1)

with open('caps.txt', 'w') as file:
	file.write(f'1 polygons\npolygon 1 ( 3 caps, 1 weight, 0 pixel, 0 str):\n')
	file.write(f' {cap_to_str(racap, 11)}\n')
	file.write(f' {cap_to_str(deccap, 17)}\n')
	file.write(f' {cap_to_str(cap, 9)}')