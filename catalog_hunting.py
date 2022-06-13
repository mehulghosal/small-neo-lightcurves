import numpy as np
import astropy as ap
import matplotlib.pyplot as plt
from matplotlib import colors
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.wcs import utils
from astropy import units as u
from scipy.ndimage import rotate
from magic_star import point_rotation

star_params = np.loadtxt('star_parameters.csv', delimiter=',', dtype=object, skiprows=1)

L = np.array(star_params[:,1], dtype=float)
a = np.array(star_params[:,2], dtype=float)

star_params = star_params[np.where(L>0)]
L = np.array(star_params[:,1], dtype=float)
a = np.array(star_params[:,2], dtype=float)

# to get absolute mag and orbital information
from astroquery.jplhorizons import Horizons

plt.rcParams.update({'figure.max_open_warning': 0})

# initializing all directories
import os, subprocess
from os.path import isdir, isfile, join
directory = './'	
dir_names = [directory+f+'/' for f in os.listdir(directory) if isdir(join(directory,f))] 

# output = open('output_rates.csv', 'w+')

input_file = np.loadtxt('input.csv', dtype=object, skiprows=1, usecols=(i for i in range(25)), delimiter=',')

mins = {'g':100, 'r': 150, 'i': 250}

for d in dir_names:
	file_names = [d+f for f in os.listdir(d) if isfile(join(d,f))]
	# stars      = file_names[]
	if 'GE1' not in d: continue
	star_index = -1
	for i in range(len(star_params)):
		# pass
		if d.split('_')[1] in star_params[i,0]:
			a = float(star_params[i,2])

	for f in file_names:

		# if '66' not in f: continue

		try:
			file = fits.open(f)
		except Exception as e:
			# print(f)
			continue
		print(f)
		hdr = file[0].header
		img = file[0].data

		img_star_rotated = rotate(img, a)


		# object id from fits file - inconsistent naming --> frustrating
		obj_id = hdr["OBJECT"][:-2].replace('_', ' ')
		# object id from directory name --> string splicing
		obj_id = f.split('_')
		obj_id = obj_id[0][2:] + ' ' + obj_id[1]
		# if '2016 GE1' not in obj_id: continue
		# if '2015 VH65' not in obj_id: continue
		# if not ('2016 GE1' in obj_id and '70o13' in f): continue

		plt.figure()
		plt.title(f)

		obj_rows = input_file[np.where(input_file[:,1]==obj_id),:][0]

		try:
			obj = obj_rows[np.where(obj_rows[:,0]==f.split('/')[-1])][0]
			trail_start = np.array(obj[-4:-2], dtype=int)
			trail_end	= np.array(obj[-2:], dtype=int)
			c = SkyCoord(f'{obj[7]} {obj[8]}', unit=(u.deg, u.deg))
		except Exception as e:
			# print(f,obj[-4:-2],obj[-2:])
			plt.close()
			continue

		angle = -1*np.arctan2(trail_end[0]-trail_start[0], trail_end[1]-trail_start[1]) * 180/np.pi
		img_rotated = rotate(img, angle)
		plt.imshow(img, cmap='gray', norm=colors.LogNorm(vmin=mins[hdr['FILTER'][0]]))

		# trail_start = point_rotation(trail_start[0], trail_start[1], a, img, img_star_rotated)
		# trail_end   = point_rotation(trail_end  [0], trail_end  [1], a, img, img_star_rotated)
		# plt.plot([trail_start[0], trail_end[0]], [trail_start[1], trail_end[1]], marker='*')

		# WCS stuff
		w = WCS(hdr)
		# c = SkyCoord(f'{hdr["CRVAL1"]} {hdr["CRVAL2"]}', unit=(u.deg, u.deg))
		
		# c = SkyCoord(f'{hdr["CRVAL1"]} {hdr["CRVAL2"]}', unit=(u.deg, u.deg))

		target_x, target_y = np.round(utils.skycoord_to_pixel(c, w))

		trail_centroid = [trail_start[0], int((trail_start[1]+trail_end[1])/2 + .5)]

		# # add offset from skycoord_to_pixel
		x_offset = (-target_x + trail_centroid[0])
		y_offset = (-target_y + trail_centroid[1]) 

		plt.plot(target_x + x_offset, target_y + y_offset, 'r+')
		plt.plot(target_x , target_y , 'g+')


		# args = ['./refcat', f'{f_center.ra.deg}', f'{f_center.dec.deg}', '-rect', f'{d_ra},{d_dec}', '-dir 00_m_16/']
		# args_str = f'./refcat {f_center.ra.deg} {f_center.dec.deg} -rect {d_ra},{d_dec} -dir 00_m_16/'
		args_str = f'./refcat {c.ra.deg} {c.dec.deg} -rad 0.5 -dir 00_m_16/'

		print(args_str)

		# RA, Dec, g, r, i, z, J, cyan, orange.
		ref_stars = np.array(os.popen(args_str).read().split('\n')[:-1])

		# print(stars)
		# print(stars.stderr)
		# print(stars)
		refcat = []
		for i in ref_stars:
			refcat.append(np.array(i.split(), dtype=float))
		refcat = np.array(refcat)

		refcat_ra_dec      = SkyCoord(ra=refcat[:,0]*u.degree, dec=refcat[:,1]*u.degree, frame='fk5')
		refcat_x, refcat_y = np.round(utils.skycoord_to_pixel(refcat_ra_dec, w))


		try:
			fitted_stars = np.loadtxt(f'{f[:-4]}_params.txt')
		except Exception as e:
			print(e, f)
			continue

		# ast_fit = fitted_stars[0]
		str_fit = fitted_stars[:]

		star_x = fitted_stars[:,4]
		star_y = fitted_stars[:,5]

		flux   = fitted_stars[:,-1]
		mag    = -2.5 * np.log10(flux)
		# print('mag', mag)

		# print(star_params[:,0,5:/])
		# print(np.where(d.split('_')[1] == star_params[:,0][5:]))

		# rotated_x = []
		# rotated_y = []
		# for i in range(len(refcat_x)):
		# 	rot_x, rot_y = point_rotation(refcat_x[i], refcat_y[i], a, img, img_star_rotated)

		# print(our_catalog)
		# print(refcat_ra_dec)
		
		# # rotating back to (master?) frame

		a+=angle

		if a<0: 
			a *= -np.pi/180
			m = img.shape[0] * np.abs(np.sin(a))
			star_x_rot =  (star_x -m) * np.cos(a) + star_y * np.sin(a) 
			star_y_rot = -(star_x -m) * np.sin(a) + star_y * np.cos(a)

		elif a>0:
			a *= -np.pi/180
			m = img.shape[1] * np.abs(np.sin(a))
			star_x_rot =  (star_x) * np.cos(a) + (star_y -m) * np.sin(a)
			star_y_rot = -(star_x) * np.sin(a) + (star_y -m) * np.cos(a)


		
		our_catalog = utils.pixel_to_skycoord(star_x_rot, star_y_rot, w)
		# our_catalog = utils.pixel_to_skycoord(star_x_rot, star_y_rot, w)
		# print(our_catalog)


		# plt.scatter(star_x_rot + x_offset, star_y_rot + y_offset, label='rotated fitted stars')
		# plt.scatter(star_x_rot[idX], star_y_rot[idX], label='rotated fitted stars')

		idx, d2d, d3d = our_catalog.match_to_catalog_sky(refcat_ra_dec, nthneighbor=1)
		idX, d2D, d3D = refcat_ra_dec.match_to_catalog_sky(our_catalog, nthneighbor=1)

		max_sep = 10.0 * u.arcsec
		sep_constraint = d2d < max_sep
		# refcat_matches = 

		# print(idX)
		# print(refcat_ra_dec[idx[sep_constraint]])
		# print((d2D.deg)*3600)
		# print(len(sep_constraint[0]))

		plt.scatter(refcat_x[idx] , refcat_y[idx], label='refcat stars')
		# plt.scatter(refcat_x , refcat_y, label='refcat stars')
		# plt.scatter(refcat_x[sep_constraint], refcat_y[sep_constraint], label='refcat stars')

		plt.scatter(star_x_rot, star_y_rot, label='rotated fitted stars')

		# plt.scatter(star_x, star_y, label='unrotated')
		plt.legend()
		plt.xlim((0, img.shape[1]))
		plt.ylim((img.shape[0], 0))

		# fig = plt.figure()
		# ax  = fig.add_subplot('111', projection='aitoff')


		# # ax.scatter(refcat_ra_dec.ra.radians, refcat_ra_dec.dec.radians)
		# # ax.scatter()


		# if True: break

	# if True: break

	plt.show()
# output.close()