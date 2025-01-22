#CTE Importing necessary packages
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image as im #CTE This package is used to open the image 
			    #CTE with any viewer available
import os
import num_col #CTE Will run num_col.py before running this script

dat = np.loadtxt('num-col.txt', unpack=True) #CTE Reading in the result of num_col.py

fig, ax=plt.subplots()

ax.scatter(dat[0], dat[1], color='gold', marker='x')
ax.plot(dat[0], dat[1], color='red')


if os.path.exists('test-plot.png'): #CTE If the resulting image exists, delete it and save a new one
	os.remove('test-plot.png')
	fig.savefig('test-plot.png')
else:
	fig.savefig('test-plot.png') #CTE If the image does not exist, save a new one

im.open('test-plot.png').show() #CTE Use whatever image viewer is available to view the plot

