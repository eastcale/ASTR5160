import sys
import os
import numpy as np
import pymangle
lap_path = '/Users/calebeastlund/Desktop/'
dep_path = '/d/dor1/caleb/Classes/AstrTechII/Notes/week5/week5-code/'
if os.path.exists(lap_path):
	sys.path.insert(0, lap_path)
elif os.path.exists(dep_path):
	sys.path.insert(0, dep_path)
	sys.path.insert(0, "/d/users/caleb/Py_Modules")
from plots import min_max, tick_array, ticks, save_show
from sph_caps import sph_cap, cap_to_str, to_poly
import matplotlib.pyplot as plt

#CTE FYI FOR ADAM: I made signifcant changes to my sph_caps.py
#CTE file from week5 in light of this assignment. Many of the functions I
#CTE call have been altered since they were last evaluated.

#CTE Also, since I am not defining any functions, do I still need
#CTE an 'if __name__ == '__main__'' statement? 

#CTE Converting ra to hourangle 
ra_cap1 = 76/15
ra_cap2 = 75/15 
dec_cap1 = 36
dec_cap2 = 35
r_cap1 = 5
r_cap2 = 5

#CTE Generating the vector for the two caps
cap1 = sph_cap(ra_cap1, dec_cap1, r_cap1)
cap2 = sph_cap(ra_cap2, dec_cap2, r_cap2)

#CTE Converting caps to pymangle readable format 
cap1_str = cap_to_str(cap1, 8)
cap2_str = cap_to_str(cap2, 8)

p1 = [[cap1_str, cap2_str]]
p2 = [[cap1_str], [cap2_str]]

to_poly('intersection.ply', 1, [2], p1)
to_poly('bothcaps.ply', 2, [1, 1], p2)

#CTE Creating a mask for each polygon
minter = pymangle.Mangle('intersection.ply')
mboth = pymangle.Mangle('bothcaps.ply')

#CTE Generating 10,000 random points for each of the 2 polygons
inter_rand_ra, inter_rand_dec = minter.genrand(10000) 
both_rand_ra, both_rand_dec = mboth.genrand(10000)

#CTE Plotting each of the 2 polygons
fig, ax = plt.subplots()

ax.scatter(both_rand_ra, both_rand_dec, color='navy', s=2, alpha=0.5)
ax.scatter(inter_rand_ra, inter_rand_dec, color='red', s=2, alpha=0.5)

save_show(fig, 'inter-both.png')
#CTE 'bothcaps.ply' covers the entire area of both caps
#CTE 'intersection.ply' covers only the section of each cap that intersects with the other

#CTE Flipping the constraint on cap1
cap1_flip = sph_cap(ra_cap1, dec_cap1, r_cap1, -1)

#CTE Converting new cap to pymangle readable format
cap1_str_flip = cap_to_str(cap1_flip, 8)

p1_flip1 =  [[cap1_str_flip, cap2_str]]

to_poly('intersection-flip1.ply', 1, [2], p1_flip1)

#CTE Creating new mask like minter but with the const. on cap1 flipped
mflip1 = pymangle.Mangle('intersection-flip1.ply')

#CTE Generating 10,000 random ponts for mflip1
flip1_rand_ra, flip1_rand_dec = mflip1.genrand(10000)

#CTE Plotting/comparing 'intersection.ply' and 'intersection-flip1.ply'
fig, ax = plt.subplots()

ax.scatter(inter_rand_ra, inter_rand_dec, color='navy', s=2, alpha=0.5)
ax.scatter(flip1_rand_ra, flip1_rand_dec, color='red', s=2, alpha=0.5)

save_show(fig, 'inter-flip1.png')
#CTE Inverting cap1 and making an intersection mask covers only the section
#CTE of cap1 that does NOT intersect with cap2

#CTE Flipping the constraint on cap2
cap2_flip = sph_cap(ra_cap2, dec_cap2, r_cap2, -1)

#CTE Converting the new cap2 to pymangle readable format
cap2_str_flip = cap_to_str(cap2_flip, 8)

p1_flip2 = [[cap1_str, cap2_str_flip]]

to_poly('intersection-flip2.ply', 1, [2], p1_flip2)

#CTE Creating new mask like minter but with the const. on cap2 flipped
mflip2 = pymangle.Mangle('intersection-flip2.ply')

flip2_rand_ra, flip2_rand_dec = mflip2.genrand(10000)

#CTE Plotting/comparing 'intersection.ply', 'intersection-flip1.ply' and 'intersection-flip2.ply'
fig, ax = plt.subplots()

ax.scatter(flip2_rand_ra, flip2_rand_dec, color='green', s=2, alpha=0.5)
ax.scatter(inter_rand_ra, inter_rand_dec, color='navy', s=2, alpha=0.5)
ax.scatter(flip1_rand_ra, flip1_rand_dec, color='red', s=2, alpha=0.5)

save_show(fig, 'inter-flip12.png')
#CTE Inverting cap2 and making an intersection mask covers only the section
#CTE of cap2 that does NOT intersect with cap1

#CTE Making a new polygon with both cap1 and cap2 to inverted
p1_flip12 = [[cap1_str_flip, cap2_str_flip]]

to_poly('flip12.ply', 1, [2], p1_flip12)

#CTE Creating a mask with both caps inverted
mflip12 = pymangle.Mangle('flip12.ply')

flip12_rand_ra, flip12_rand_dec = mflip12.genrand(1000000)

#CTE Plotting the mask with both caps inverted
fig, ax = plt.subplots()

ax.scatter(flip12_rand_ra, flip12_rand_dec, color='navy', s=2, alpha=0.5)

save_show(fig, 'flip12.png')
#CTE Inverting both cap1 and cap2 and making a polygon with 2 caps
#CTE covers all area that is NOT in cap1 OR cap2