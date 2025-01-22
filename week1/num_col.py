# Importing necessary packages
import os
import numpy as np

# Generating array of points from 1 to 10
one_ten=np.arange(1, 11, 1)

def write(filename):
	"""
	Write two columns each ranging from 1-10 into a given file.
	
	Parameters
	----------
	filename : string
		The name of the file that will be used
	
	Returns
    	-------
    	File
    		A file containing two columns, each ranging from 1-10
	"""
	f=open(filename, 'x')
	for i, j in zip(one_ten, one_ten):
		f.write('{}	{}\n'.format(i, j))
	f.close()

if os.path.exists('num-col.txt'): # If the file already exists, delete it 
				  # before writing to a new one
	os.remove('num-col.txt')
	write('num-col.txt')
else:
	write('num-col.txt') # If the file does not exist, make one and write to it
