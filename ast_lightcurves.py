import numpy as np
import astropy as ap
import matplotlib.pyplot as plt
from matplotlib import colors
from astropy.io import fits
from scipy.ndimage import rotate
from scipy.stats import mode
from scipy.optimize import curve_fit
from astropy.wcs import WCS, utils
from magic_star import point_rotation, reverse_rotation
from astropy.coordinates import SkyCoord
from astropy import units as u

# to get absolute mag and orbital information
from astroquery.jplhorizons import Horizons

plt.rcParams.update({'figure.max_open_warning': 0})

paperheight = 10
paperwidth = 13
margin = 1.0

fontsize_standard = 28

# initializing all directories
import os
from os.path import isdir, isfile, join
directory = './'	
dir_names = [directory+f+'/' for f in os.listdir(directory) if isdir(join(directory,f))] 

filt_ind = {'g':21, 'r': 25, 'i': 29}
mins = {'g':100, 'r': 150, 'i': 250}

for d in dir_names:
	lc_dirs = [d+f for f in os.listdir(d) if isdir(join(d,f))] 
	if not 'XD169' in d: continue

	times , mags , mags_err  = [] , [] , [] 
	zps , zps_err = [] , []
	fig_combined, ax_combined = plt.subplots(figsize=((paperwidth*1.15) - 2 * margin, (paperheight*1.15) - 2 * margin))


	img_name = ''	

	for ld in lc_dirs :

		if 'data/LCLIST' in ld or 'git' in ld: continue
		# if not ('51' in ld): continue

		# print(ld)

		lc_files = [join(ld,f) for f in os.listdir(ld) if isfile(join(ld,f))]

		# print( lc_files )

		
		# print(lc_files)

		for f in lc_files :

			if not 'lightcurve_asteroid' in f: continue
			if '59o' in f: continue
			
			# if not ( '66' in f or '67' in f or '68' in f or '69' in f or '70o' in f or '71o' in f or '71o' in f)  : continue


			fits_name = ('/'.join(f.split('/')[:-1]) + '.flt')
			# if 'on' in f : fits_name = ('/'.join(f.split('/')[:-1]) + '.fits')
			# else: continue

			try:
				fits_file = fits.open(fits_name)
				img_name  = f.split('/')[2].split('o')[0]
				print(fits_name)
			except Exception as e:
				print(f'NO FITS FILE FOUND: {fits_name}')
				continue

			# print(join(ld,img_name)+'_calibrated_lightcurve.txt')



			hdr = fits_file[0].header
			img = fits_file[0].data

			try:
				t , flux , flux_err, _ , _ = np.loadtxt(f , unpack=True)
			except Exception as e:
				print(e)
				t , flux , flux_err = np.loadtxt(f , unpack=True)

			# t -= t[0]
			# t *= 24*3600

			flux = np.abs(flux)

			print(hdr['FILTER'])
			# fig , ax = plt.subplots(figsize=((paperwidth*1.15) - 2 * margin, (paperheight*1.15) - 2 * margin))
			# ax.errorbar ( t , flux , flux_err , fmt='s' , markerfacecolor='blue' , markeredgecolor='black' , ecolor='black' , capthick=2 , markersize=7 , capsize=3  )
			# ax.set_xlim(0,60)

			mag = -2.5 * np.log10 ( flux )
			mag_err = np.abs(1.0875 * flux_err / flux )

			# fig_mag, ax_mag = plt.subplots(figsize=((paperwidth*1.15) - 2 * margin, (paperheight*1.15) - 2 * margin))
			# ax_mag.errorbar ( t , mag , mag_err , fmt='s' , markerfacecolor='blue' , markeredgecolor='black' , ecolor='black' , capthick=2 , markersize=7 , capsize=3   )
			# ax_mag.set_xlim(0,60)

			zp , zp_err = np.loadtxt(d+ img_name+'.zp' , unpack=True)

			fig_cal , ax_cal = plt.subplots(figsize=((paperwidth*1.15) - 2 * margin, (paperheight*1.15) - 2 * margin))
			ax_cal.errorbar ( t , mag + zp , (mag_err**2 + zp_err**2)**.5 , fmt='s' , markerfacecolor='blue' , markeredgecolor='black' , ecolor='black' , capthick=2 , markersize=7 , capsize=3   )
			# ax_cal.set_xlim(0,60)

			# ax_combined.errorbar(t , mag + zp , (mag_err**2 + zp_err**2)**.5 , fmt='s' , markerfacecolor='blue' , markeredgecolor='black' , ecolor='black' , capthick=2 , markersize=7 , capsize=3   )

			times.append(t)
			mags.append (mag + zp )
			mags_err.append((mag_err**2 + zp_err**2)**.5)
			zps.append(zp)
			zps_err.append(zp_err)

			np.savetxt ( join(ld,img_name)+'_calibrated_lightcurve.txt' , np.vstack([t , mag+zp , (mag_err**2 + zp_err**2)**.5 ]).T , header=f'zp={zp}')

			# plt.show()



	# t_0  = np.min([x[0] for x in times ])

	times = np.hstack(times) 
	mags  = np.hstack(mags)
	mags_err  = np.hstack(mags_err)
	# print(times)

	t_0 = np.min(times)
	# if True: break

	# for i in range(len(times)):
	ax_combined.errorbar((times - t_0)*24*3600 , mags, mags_err , fmt='s' , markerfacecolor='blue' , markeredgecolor='black' , ecolor='black' , capthick=2 , markersize=7 , capsize=3   )

	print( )
	np.savetxt ( d + 'corrected_timeseries.txt' , np.vstack ([ times , mags , mags_err]).T , header=' '.join(lc_dirs) )

	plt.show()
