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

    #CTE Converting the RA to degrees, then to radians
    ra_rad = np.radians(15 * ra)

    #CTE Converting the DEC to radians
    dec_rad = np.radians(dec) 

    x = np.cos(ra_rad)*np.cos(dec_rad)
    y = np.sin(ra_rad)*np.cos(dec_rad)
    z = np.sin(dec_rad)

    xyz = [x, y, z]
    
    return xyz

def ra_cap(ra, pm=1):
	"""
    Creates a vector 4-array (x, y, z, 1-h) for a spherical
    cap bounded by the given right ascension.

    Parameters
    ----------
    ra : int or float
    	The right ascension (in hours) bounding the spherical cap.

    pm : int
        Options are 1 or -1. If 1, gives a standard vector for the spherical
        cap. If -1, gives a vector equal in magnitude but opposite in direction
        to the standard vector.

    Returns
    -------
	ra_bound : list of floats
		The resulting vector 4-array (x, y, z, 1-h) describing the spherical cap
		bounded by the given right ascension.
    """

	dec = 0

    #CTE RA in hours + 90 deg / 15 hr/deg
	xyz = to_cart(ra+90/15, dec) 

	ra_bound = [xyz[0], xyz[1], xyz[2], pm]

	return ra_bound

def dec_cap(dec, pm=1):
	"""
    Creates a vector 4-array (x, y, z, 1-h) for a spherical
    cap bounded by the given declination.

    Parameters
    ----------
    dec : int or float
    	The declination (in degrees) bounding the spherical cap.

    pm : int
        Options are 1 or -1. If 1, gives a standard vector for the spherical
        cap. If -1, gives a vector equal in magnitude but opposite in direction
        to the standard vector.

    Returns
    -------
	dec_bound : list of floats
		The resulting vector 4-array (x, y, z, 1-h) describing the spherical cap
		bounded by the given declination.
    """

	ra = 0
	dec1 = 90

	xyz = to_cart(ra, dec1)

	h = pm*(1-np.sin(np.radians(dec)))

	dec_bound = [xyz[0], xyz[1], xyz[2], h]

	return dec_bound

def sph_cap(ra, dec, theta, pm=1):
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

    pm : int
        Options are 1 or -1. If 1, gives a standard vector for the spherical
        cap. If -1, gives a vector equal in magnitude but opposite in direction
        to the standard vector.
        
    Returns
    -------
	cap : list of floats
		The resulting vector 4-array (x, y, z, 1-h) describing the spherical cap.
    """

	xyz = to_cart(ra, dec)

	h = pm*(1-np.cos(np.radians(theta)))

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

def to_poly(file_name, n_poly, n_cap, c_list, c1 = 'None', c2 = 'None', area = 'n'):
    """
    Writes a text file for n_poly polygons in a pymangle readable format. 
    Each polygon has a corresponding index in n_cap to specify the number of caps.
    c_list is a list that is n_poly lists long, each one with the corresponding
    number of caps.

    Parameters
    ----------
    file_name : str
        The name of the file to be written
    
    n_poly : int
        The number of polygons to be written to the file

    n_cap : list of ints
        The number of caps for each polygon.

    c_list : list of lists of ints
        For i in n_cap, each term in c_list has i number of caps,
        given by the cap_to_str function. ie., for n_cap = [1, 2], 
        c_list = [[cap1], [cap2, cap3]]

    c1 : list of floats
        The coordinates of the bottom left corner of the bounded shape
        (right now, only rectangular shapes are accepted) in the form:
        [ra (hours), dec (deg)]. Default is 'None' is no input is given.

    c2 : list of floats
        The coordinates of the upper right corner of the bounded shape
        (right now, only rectangular shapes are accepted) in the form:
        [ra (hours), dec (deg)]. Default is 'None' is no input is given.

    area : str
        Options are 'y' or 'n'. If 'y', then the area of the bounded shape
        is calculated in steradians and put into the polygon file. If 'n',
        then the area is given to be 0 steradians. Default is 'n'.

    Returns
    -------
    None
    """

    #CTE Determining the area of the polygon
    if area == 'n':
        ster = 0
    elif area == 'y' and c1 != 'None' and c2 != 'None':
        da =  np.radians(15 * c2[0]) - np.radians(15 * c1[0])
        dd =  np.sin(np.radians(c2[1])) - np.sin(np.radians(c1[1]))

        ster = round((da * dd), 8)

    #CTE Writing the polygon parameters in a pymangle readable format.
    with open(file_name, 'w') as file:
        file.write(f'{n_poly} polygons\n')
        for i in range(n_poly):
            file.write(f'polygon {i+1} ( {n_cap[i]} caps, 1 weight, 0 pixel, {ster} str):\n')
            for j in range(n_cap[i]):
                file.write(f' {c_list[i][j]}\n')

if __name__ == '__main__':
    racap = ra_cap(5)

    #CTE Result is [6.123233995736766e-17, 0.0, 1.0, 0.41221474770752686] but am told it should be                             
    #CTE [0, 0, 1, 0.4122...]. Rounding? VVV
    deccap = dec_cap(36)

    #CTE Manually setting deccap[0] to zero to match the notes 
    deccap[0] = 0.0 

    cap = sph_cap(5, 36, 1)

    ra_cap_str = cap_to_str(racap, 11)
    dec_cap_str = cap_to_str(deccap, 17)
    cap_str = cap_to_str(cap, 9)

    c_list = [[ra_cap_str, dec_cap_str, cap_str]]

    to_poly('caps.txt', 1, [3], c_list)