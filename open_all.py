import numpy as np
import astropy as ap
import matplotlib.pyplot as plt
from matplotlib import colors
from astropy.io import fits
from scipy.ndimage import rotate
from astropy.wcs import WCS
from magic_star import point_rotation
# to get absolute mag and orbital information
from astroquery.jplhorizons import Horizons

plt.rcParams.update({'figure.max_open_warning': 0})

# initializing all directories
import os
from os.path import isdir, isfile, join
directory = './'	
dir_names = [directory+f+'/' for f in os.listdir(directory) if isdir(join(directory,f))] 

# output = open('output_rates.csv', 'w+')

input_file = np.loadtxt('input.csv', dtype=object, skiprows=1, usecols=(i for i in range(25)), delimiter=',')
i=0
mins = {'g':100, 'r': 150, 'i': 250}
start_times = []
for d in dir_names:
	file_names = [d+f for f in os.listdir(d) if isfile(join(d,f))]
	if not 'LT1_2016_06_07' in d: continue

	for f in file_names:
		try:
			file = fits.open(f)
		except Exception as e:
			print(f)
			continue
		hdr = file[0].header
		img = file[0].data

		# for experimenting, lol
		start_times.append(float(hdr['MJDATE']))

		# object id from fits file - inconsistent naming --> frustrating
		obj_id = hdr["OBJECT"][:-2].replace('_', ' ')
		# object id from directory name --> string splicing
		obj_id = f.split('_')
		obj_id = obj_id[0][2:] + ' ' + obj_id[1]
		# if '2016 LT1' not in obj_id: continue
		# if '2015 VH65' not in obj_id: continue
		# if not ('2016 GE1' in obj_id and '70o13' in f): continue

		plt.figure()
		plt.title(f)
		plt.imshow(img, cmap='gray', norm=colors.LogNorm(vmin=mins[hdr['FILTER'][0]]))

		obj_rows = input_file[np.where(input_file[:,1]==obj_id),:][0]

		# try:
		# 	obj = obj_rows[np.where(obj_rows[:,0]==f.split('/')[-1])][0]
		# 	trail_start = np.array(obj[-4:-2], dtype=int)
		# 	trail_end	= np.array(obj[-2:], dtype=int)
		# except Exception as e:
		# 	print(f,'skipped')
		# 	plt.close()
		# 	continue

		# plt.scatter([trail_start[0], trail_end[0]], [trail_start[1], trail_end[1]], marker='*', label='2016 LT1 trail')
		# # plt.legend()

		# angle = -1*np.arctan2(trail_end[0]-trail_start[0], trail_end[1]-trail_start[1]) * 180/np.pi
		# img_rotated = rotate(img, angle)
		# # ax[0].imshow(img_rotated, cmap='gray', norm=colors.LogNorm(vmin=mins[hdr['FILTER'][0]]))

		# trail_start = np.array(point_rotation(trail_start[0], trail_start[1], angle, img, img_rotated), dtype=int)
		# trail_end	= np.array(point_rotation(trail_end[0]  , trail_end[1]  , angle, img, img_rotated), dtype=int)

		# plt.figure()
		# plt.title(f)
		# plt.imshow(img_rotated, cmap='gray', norm=colors.LogNorm(vmin=mins[hdr['FILTER'][0]]))
		# plt.scatter([trail_start[0], trail_end[0]], [trail_start[1], trail_end[1]], marker='*', label='2016 LT1 trail rotated')

		# plt.legend()
		

		# obj = Horizons(id=obj_id, location='sun', epochs=str(hdr["MJD-OBS"]+2400000.5)) # adding 2.4 million to convert from MJD to JD
		# ob1 = Horizons(id=obj_id, location='568@399', epochs=str(hdr["MJD-OBS"]+2400000.5)) # for RA and DEC rates
		# ele = obj.elements()
		# eph = ob1.ephemerides()

		# print(f'{f} mjd: {hdr["MJD-OBS"]}. filter: {hdr["FILTER"]}')
		# print(f"{ele['H'][0]}  {ele['e'][0]}  {ele['a'][0]}.")



		# to_write = f"{f.split('/')[-1]}, {obj_id}, {ele['H'][0]}, {ele['a'][0]}, {ele['e'][0]}, {ele['incl'][0]}, {hdr['MJD-OBS']+2400000.5}, {hdr['RA_DEG']}, {hdr['DEC_DEG']}, {hdr['FILTER']}, {eph['RA_rate'][0]/60}, {eph['DEC_rate'][0]/60} \n"
		# to_write = f"{eph['DEC_rate'][0]/60} \n"

		# output.writelines(to_write)

		# if True: break

	# plt.figure()
	# plt.plot(np.sort(start_times))
	i+=1
	plt.show()
# output.close()
print(i)