import numpy as np
import astropy as ap
import matplotlib.pyplot as plt
from matplotlib import colors
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.wcs import utils
from astropy import units as u

plt.rcParams.update({'figure.max_open_warning': 0})

# initializing all directories
import os, subprocess
from os.path import isdir, isfile, join

f = 'new-image.fits'
try:
	file = fits.open(f)
except Exception as e:
	print(f)
hdr = file[0].header
img = file[0].data

# img_star_rotated = rotate(img, a)


# object id from fits file - inconsistent naming --> frustrating
# obj_id = hdr["OBJECT"][:-2].replace('_', ' ')
# object id from directory name --> string splicing
# obj_id = f.split('_')
# obj_id = obj_id[0][2:] + ' ' + obj_id[1]
# if '2016 GE1' not in obj_id: continue
# if '2015 VH65' not in obj_id: continue
# if not ('2016 GE1' in obj_id and '70o13' in f): continue

plt.figure()
plt.title(f)
plt.imshow(img, cmap='gray', vmin=1800, vmax=20000)

# WCS stuff
w = WCS(hdr)
c = SkyCoord(f'{hdr["CRVAL1"]} {hdr["CRVAL2"]}', unit=(u.deg, u.deg))
# c = SkyCoord(f'{obj[7]} {obj[8]}', unit=(u.deg, u.deg))
target_x, target_y = np.round(utils.skycoord_to_pixel(c, w))

# print(target_x, target_y, f'{obj[7]} {obj[8]}')


# trail_centroid = [trail_start[0], int((trail_start[1]+trail_end[1])/2 + .5)]

# # add offset from skycoord_to_pixel
# x_offset = (-target_x + trail_centroid[0])
# y_offset = (-target_y + trail_centroid[1]) 

# plt.plot(target_x , target_y , 'r+')


frame_center = np.array(img.shape)/2
f_center = utils.pixel_to_skycoord(frame_center[0], frame_center[1], w)

frame_top = [0,img.shape[1]/2]
f_top     = utils.pixel_to_skycoord(frame_top[0], frame_top[1], w)

frame_left = [img.shape[0]/2, 0]
f_left     = utils.pixel_to_skycoord(frame_left[0], frame_left[1], w)

d_dec = f_center.separation(f_top).deg * 1
d_ra  = f_center.separation(f_left).deg * 1

print(f_center.ra.deg, f_center.dec.deg)
print(d_ra, d_dec)

# d_ra = 1
# d_dec = 1

# args = ['./refcat', f'{f_center.ra.deg}', f'{f_center.dec.deg}', '-rect', f'{d_ra},{d_dec}', '-dir 00_m_16/']
args_str = f'./refcat {f_center.ra.deg} {f_center.dec.deg} -rect {d_ra},{d_dec} -dir 00_m_16/'
# args_str = f'./refcat {f_center.ra.deg} {f_center.dec.deg} -rad 1.0 -dir 00_m_16/'

print(args_str)

# RA, Dec, g, r, i, z, J, cyan, orange.
stars = np.array(os.popen(args_str).read().split('\n')[:-1])
print()

refcat = []
for i in stars:
	refcat.append(np.array(i.split(), dtype=float))
refcat = np.array(refcat)
# print(refcat.shape)

my_cat = np.loadtxt('catalog.cat')
my_coords = SkyCoord(ra=my_cat[:,3] *u.degree, dec=my_cat[:,4]*u.degree)
print(my_coords.shape)
plt.scatter(my_cat[:,1], my_cat[:,2], c='red', label='sextractor')


refcat_ra_dec = SkyCoord(ra=refcat[:,0]*u.degree, dec=refcat[:,1]*u.degree, frame='fk5')
print(refcat_ra_dec.shape)
refcat_x, refcat_y = np.round(utils.skycoord_to_pixel(refcat_ra_dec, w))

idx, d2d, d3d = my_coords.match_to_catalog_sky(refcat_ra_dec)

print(idx)

# plt.scatter(my_cat[:,1], my_cat[:,2], label='SExtractor')



# inv_filter = np.where()



# plt.scatter(refcat_x, refcat_y, color='red', label='all refcat stars')
plt.scatter(refcat_x[idx] , refcat_y[idx], color='blue', label='matched refcat')

print(np.max(d2d.deg))

print(refcat[idx])

plt.legend()


plt.show()