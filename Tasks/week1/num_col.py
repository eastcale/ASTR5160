#CTE Importing necessary packages
import os
import numpy as np
import subprocess

#CTE Generating array of points from 1 to 10
one_ten=np.arange(1, 11, 1)

def write(file_name):
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
	f=open(file_name, 'x')
	for i, j in zip(one_ten, one_ten):
		f.write('{}	{}\n'.format(i, j))
	f.close()

def write_show(file_name):
	"""
	Used to save the output of write(). Opens in gedit if user says yes.
	
	Parameters
	----------
	file_name : string
		The desired name of the file to be saved.
		
	tab : astropy.table.Table
		The table to be written to a fits file
	"""
	if os.path.exists(file_name): 
		yn = input("Overwrite {} (y/n)? ".format(file_name))
		
		if yn == 'y':	
			os.remove(file_name)
			write(file_name)
		elif yn == 'n':
			pass
	else:
		write(file_name)
		
	dis = input("Open {} (y/n)?".format(file_name))
	
	if dis == "y":
		subprocess.Popen(['gedit', file_name])
	else:
		pass
	
write_show('num-col.txt')
