import numpy as np
import astropy as ap
import matplotlib.pyplot as plt
from matplotlib import colors
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.wcs import utils
from astropy import units as u


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
	if 'VH1' not in d: continue

	for f in file_names:
		try:
			file = fits.open(f)
		except Exception as e:
			print(f)
			continue
		hdr = file[0].header
		img = file[0].data


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
		plt.imshow(img, cmap='gray', norm=colors.LogNorm(vmin=mins[hdr['FILTER'][0]]))

		obj_rows = input_file[np.where(input_file[:,1]==obj_id),:][0]

		try:
			obj = obj_rows[np.where(obj_rows[:,0]==f.split('/')[-1])][0]
			trail_start = np.array(obj[-4:-2], dtype=int)
			trail_end	= np.array(obj[-2:], dtype=int)
		except Exception as e:
			# print(f,obj[-4:-2],obj[-2:])
			plt.close()
			continue

		plt.plot([trail_start[0], trail_end[0]], [trail_start[1], trail_end[1]], marker='*')

		# WCS stuff
		w = WCS(hdr)
		c = SkyCoord(f'{hdr["CRVAL1"]} {hdr["CRVAL2"]}', unit=(u.deg, u.deg))
		# c = SkyCoord(f'{obj[7]} {obj[8]}', unit=(u.deg, u.deg))
		target_x, target_y = np.round(utils.skycoord_to_pixel(c, w))

		print(target_x, target_y, f'{obj[7]} {obj[8]}')
		

		trail_centroid = [trail_start[0], int((trail_start[1]+trail_end[1])/2 + .5)]

		# add offset from skycoord_to_pixel
		x_offset = (-target_x+trail_centroid[0])
		y_offset = (-target_y + trail_centroid[1]) 

		# plt.plot(target_x + x_offset, target_y + y_offset, 'r+')
		plt.plot(target_x + x_offset, target_y + y_offset, 'r+')

		# ast_sky_x, ast_sky_y = np.round(utils.pixel_to_skycoord(trail_middle[0], trail_middle[1], w))

		a = utils.pixel_to_skycoord(trail_centroid[0], trail_centroid[1], w)


		frame_center = np.array(img.shape)/2
		f_center = utils.pixel_to_skycoord(frame_center[0], frame_center[1], w)

		frame_top = [0,img.shape[1]/2]
		f_top     = utils.pixel_to_skycoord(frame_top[0], frame_top[1], w)

		frame_left = [img.shape[0]/2, 0]
		f_left     = utils.pixel_to_skycoord(frame_left[0], frame_left[1], w)

		d_dec = f_center.separation(f_top).deg * 4
		d_ra  = f_center.separation(f_left).deg * 4

		print(f_center.ra.deg, f_center.dec.deg)
		print(d_ra, d_dec)

		args = ['./refcat', f'{f_center.ra.deg}', f'{f_center.dec.deg}', '-rect', f'{d_ra},{d_dec}', '-dir 00_m_16/']
		# args_str = f'./refcat {f_center.ra.deg} {f_center.dec.deg} -rect {d_ra},{d_dec} -dir 00_m_16/'
		args_str = f'./refcat {f_center.ra.deg} {f_center.dec.deg} -rad 1.0 -dir 00_m_16/'

		print(args_str)

		# RA, Dec, g, r, i, z, J, cyan, orange.
		stars = np.array(os.popen(args_str).read().split('\n')[:-1])
		print()
		# print(stars)
		# print(stars.stderr)
		# print(stars)
		refcat = []
		for i in stars:
			refcat.append(np.array(i.split(), dtype=float))
		refcat = np.array(refcat)
		print(refcat.shape)


		try:
			fitted_stars = np.loadtxt(f'{f[:-4]}_params.txt')
		except Exception as e:
			print(e, f)
			continue

		ast_fit = fitted_stars[0]
		str_fit = fitted_stars[1:]

		our_catalog = utils.pixel_to_skycoord(fitted_stars[:,-2], fitted_stars[:,-1], w)
		print(our_catalog)
		# for i in our_catalog: print(i)


		refcat_ra_dec = SkyCoord(ra=refcat[:,0]*u.degree, dec=refcat[:,1]*u.degree, frame='fk5')
		print(refcat_ra_dec)
		refcat_x, refcat_y = np.round(utils.skycoord_to_pixel(refcat_ra_dec, w))

		idx, d2d, d3d = our_catalog.match_to_catalog_sky(refcat_ra_dec)

		print(idx)
		print(len(refcat_ra_dec))



		# inv_filter = np.where()



		plt.scatter(refcat_x, refcat_y, color='red', label='all refcat stars')
		plt.scatter(refcat_x[idx] , refcat_y[idx], color='blue', label='matched refcat')

		plt.scatter(fitted_stars[:,-2], fitted_stars[:,-1], color='green', label='fitted stars')
		print(np.max(d2d.deg))

		# obj = Horizons(id=obj_id, location='sun', epochs=str(hdr["MJD-OBS"]+2400000.5)) # adding 2.4 million to convert from MJD to JD
		# ob1 = Horizons(id=obj_id, location='568@399', epochs=str(hdr["MJD-OBS"]+2400000.5)) # for RA and DEC rates
		# ele = obj.elements()
		# eph = ob1.ephemerides()

		# print(f'{f} mjd: {hdr["MJD-OBS"]}. filter: {hdr["FILTER"]}')
		# print(f"{ele['H'][0]}  {ele['e'][0]}  {ele['a'][0]}.")



		# to_write = f"{f.split('/')[-1]}, {obj_id}, {ele['H'][0]}, {ele['a'][0]}, {ele['e'][0]}, {ele['incl'][0]}, {hdr['MJD-OBS']+2400000.5}, {hdr['RA_DEG']}, {hdr['DEC_DEG']}, {hdr['FILTER']}, {eph['RA_rate'][0]/60}, {eph['DEC_rate'][0]/60} \n"
		# to_write = f"{eph['DEC_rate'][0]/60} \n"

		# output.writelines(to_write)
		plt.legend()
		if True: break

	# if True: break

	plt.show()
# output.close()