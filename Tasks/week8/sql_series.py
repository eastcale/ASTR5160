#!/usr/bin/env python3

import sys
import matplotlib.pyplot as plt
sys.path.insert(0, "/d/users/caleb/Py_Modules")
from plots import min_max, tick_array, ticks, save_show, ticks2
from matplotlib_inline.backend_inline import set_matplotlib_formats
set_matplotlib_formats('svg')
import numpy as np
import pandas as pd

#CTE Example SQL Search to get ID, RA, DEC, and G mag for objects within 2' of (300, -1)
# select
#     p.ObjID, p.ra, p.dec, p.g
# from
#     photoObj p, dbo.fGetNearbyObjEq(300,-1,2) n
# where
#     p.objID = n.objID

dat = pd.read_csv('/dor1/caleb/Classes/AstrTechII/Notes/week8/week8-code/sql_gmag.csv')
df = pd.DataFrame(dat)
ID, ra, dec, g = np.array(df['ObjID']), np.array(df['ra']), np.array(df['dec']), np.array(df['g'])

def set_bin(x, xmax, xmin):
	ind = np.where((x >= xmin) & (x < xmax))
	return ind

bins = np.arange(np.floor(min(g)), np.floor(max(g)), 1)

index_dict = dict(zip([f'bin{i+1}' for i in range(len(bins)-1)], 
	[set_bin(g, bins[i+1], bins[i]) for i in range(len(bins)-1)])) 

coos_dict = dict(zip([f'{i} coos' for i in index_dict.keys()], 
	[[ra[index_dict[i]], dec[index_dict[i]]] for i in index_dict.keys()]))

x_ticks = [300-0.04, 300+0.04, 0.01, 0.002]
y_ticks = [-1-0.04, -1+0.04, 0.01, 0.002]

fig, ax = plt.subplots()

size = np.flip(np.linspace(0.001, 200, np.size(bins)))

for i, j in zip(coos_dict.keys(), size):
	ax.scatter(coos_dict[i][0], coos_dict[i][1], color='navy', s = j, alpha = 0.5, edgecolors='none')

ax.set_aspect('equal')
ax.invert_xaxis()

save_show(fig, 'sql-plot-g.svg')
