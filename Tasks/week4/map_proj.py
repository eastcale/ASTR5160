import sys
sys.path.insert(0, "/d/users/caleb/Py_Modules")
from plots import min_max, tick_array, ticks, save_show
import matplotlib.pyplot as plt
import numpy as np
from numpy.random import random
import dustmaps
from dustmaps.config import config
from dustmaps.sfd import SFDQuery
from astropy import units as u
from astropy.coordinates import SkyCoord as sc
from astropy import wcs
dustdir = "/d/scratch/ASTR5160/data/dust/v0_1/maps"
config["data_dir"] = dustdir
sfd = SFDQuery()

ra = 2*np.pi*(random(10000)-0.5)
dec = np.arcsin(1.-random(10000)*2.)

x_lab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']

#CTE Plotting ra/dec in a flat grid
fig_f, ax_f = plt.subplots()

ax_f.scatter(ra, dec, color='orange', s=2)
ticks(ax_f, ra, dec, 'none', 'none')
ax_f.grid(color='k', linestyle='solid', linewidth=0.5)

save_show(fig_f, 'flat-proj.png')

#CTE Plotting ra/dec in an Aitoff projection
fig_a = plt.figure()
ax_a = fig_a.add_subplot(111, projection = 'aitoff')

ax_a.scatter(ra, dec, color='orange', s=2)
ax_a.grid(color='blue', linestyle='dashed', linewidth=1.5)

ax_a.set_xticklabels(x_lab, weight=800)

save_show(fig_a, 'aitoff-proj.png')

#CTE Plotting ra/dec in a Lambert projection
fig_l = plt.figure()
ax_l = fig_l.add_subplot(111, projection = 'lambert')

ax_l.scatter(ra, dec, color='orange', s=2)
ax_l.grid(color='blue', linestyle='dashed', linewidth=1.5)

ax_l.set_xticklabels(x_lab, weight=800)

save_show(fig_l, 'lambert-proj.png')

#CTE Generating a meshgrid of ra/dec points
ra_grid, dec_grid = np.meshgrid(np.arange(0.5, 359.5, 1), np.arange(-89.5, 89.5, 1))

coos = sc(ra=ra_grid*u.deg, dec=dec_grid*u.deg)

ebvs = sfd(coos) #CTE Finding E(B-V) for each coordinate in the grid

w = wcs.WCS(naxis=2)
w.wcs.ctype = ["RA---AIT", "DEC--AIT"]
x, y = w.wcs_world2pix(ra_grid*u.deg, dec_grid*u.deg, 1) #CTE Converting ra_grid/dec_grid
														 #CTE to Cartesian coordinates

#CTE Plotting ra_grid/dec_grid in an Aitoff projection
fig_a1 = plt.figure()
ax_a1 = fig_a1.add_subplot(111, projection='aitoff')

ax_a1.scatter(ra_grid, dec_grid, color='orange', s=0.1)
ax_a1.grid(color='blue', linestyle='dashed', linewidth=1.5)

ax_a1.set_xticklabels(x_lab, weight=800)

save_show(fig_a1, 'aitoff-proj-radec.png')

#CTE Plotting Cartesian coordinates in an Aitoff Projection
fig_a2 = plt.figure()
ax_a2 = fig_a2.add_subplot(111, projection='aitoff')

ax_a2.scatter(x, y, color='orange', s=0.1)
ax_a2.grid(color='blue', linestyle='dashed', linewidth=1.5)

ax_a2.set_xticklabels(x_lab, weight=800)

save_show(fig_a2, 'aitoff-proj-xy.png')

#CTE Creating a contour plot of E(B-V) at each point; levels are
#CTE finnicky, but E(B-V) values make sense.
fig_a3 = plt.figure()
ax_a3 = fig_a3.add_subplot(111, projection='aitoff')

conts_a3 = ax_a3.contourf(x, y, ebvs, levels=10, vmax=42, vmin=21)
fig_a3.colorbar(conts_a3)

ax_a3.set_xticklabels(x_lab, weight=800)

ax_a3.grid(color='k', alpha=0.5)

save_show(fig_a3, 'aitoff-contour.png')
