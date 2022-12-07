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
from astroquery.jplhorizons import Horizons

def line ( x , m , b): return m * x + b

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

r_mags = []
r_errs = []
g_mags = []
g_errs = []
i_mags = []
i_errs = []


for d in dir_names:
	lc_dirs = [d+f for f in os.listdir(d) if isdir(join(d,f))] 
	if not 'GE1' in d: continue

	times , mags , mags_err  = [] , [] , [] 
	fig_combined, ax_combined = plt.subplots(figsize=((paperwidth*1.15) - 2 * margin, (paperheight*1.15) - 2 * margin))


	img_name = ''	

	for ld in lc_dirs :

		if 'data/LCLIST' in ld or 'git' in ld: continue

		lc_files = [join(ld,f) for f in os.listdir(ld) if isfile(join(ld,f))]

		for f in lc_files :

			if not 'calibrated_lightcurve' in f: continue
			
			fits_name = ('/'.join(f.split('/')[:-1]) + '.flt')

			try:
				fits_file = fits.open(fits_name)
				img_name  = f.split('/')[2].split('o')[0]
				print(fits_name)
			except Exception as e:
				print(f'NO FITS FILE FOUND: {fits_name}')
				continue

			hdr = fits_file[0].header
			img = fits_file[0].data

			img_filter = hdr['FILTER'][0]
			print(img_filter)


			time , mag , mag_err = np.loadtxt ( f , unpack=True)

			param , param_cov = curve_fit (line , time , mag , sigma=mag_err , absolute_sigma=True)
			print(param)

			fig, ax = plt.subplots(figsize=((paperwidth*1.15) - 2 * margin, (paperheight*1.15) - 2 * margin))
			ax.errorbar ( time - np.min(time) , mag , mag_err , fmt='s' , markerfacecolor='blue' , markeredgecolor='black' , ecolor='black' , capthick=2 , markersize=7 , capsize=3   )

			ax.plot ( time - np.min(time) , line(time, *param) , color='red' )

			corr_mag = mag - line(time, *param) 
			# corr_err = mag_err - line(time, *param)
			corr_err = mag_err

			ax.errorbar ( time - np.min(time) , corr_mag  , corr_err , fmt='s' , markerfacecolor='red' , markeredgecolor='black' , ecolor='black' , capthick=2 , markersize=7 , capsize=3   )

			if img_filter == 'r': 
				r_mags.append(mag)
				r_errs.append(mag_err)
			elif img_filter == 'g':
				g_mags.append(mag)
				g_errs.append(mag_err)
			elif img_filter == 'i':
				i_mags.append(mag)
				i_errs.append(mag_err)

			np.savetxt ( join(ld,img_name)+'_flat_lightcurve.txt' , np.vstack([time , corr_mag , corr_err ]).T )

			times.append(time)
			mags.append(corr_mag)
			mags_err.append(corr_err)

	times = np.hstack(times) 
	mags  = np.hstack(mags)
	mags_err  = np.hstack(mags_err)

	r_mags = np.hstack(r_mags)
	r_errs = np.hstack(r_errs)
	g_mags = np.hstack(g_mags)
	g_errs = np.hstack(g_errs)
	i_mags = np.hstack(i_mags)
	i_errs = np.hstack(i_errs)

	t_0 = np.min(times)

	r_mean = np.average ( r_mags , weights=1/r_errs**2 )
	g_mean = np.average ( g_mags , weights=1/g_errs**2 )
	i_mean = np.average ( i_mags , weights=1/i_errs**2 )
	print(f'g-r={g_mean-r_mean:.2f}')
	print(f'r-i={r_mean-i_mean:.2f}')

	fig1 , ax1 = plt.subplots( figsize=((paperwidth*1.15) - 2 * margin, (paperheight*1.15) - 2 * margin) )
	ax1.errorbar((times - t_0)*24*3600 , mags + r_mean, mags_err , fmt='s' , markerfacecolor='blue' , markeredgecolor='black' , ecolor='black' , capthick=2 , markersize=7 , capsize=3   )
	np.savetxt ( d + 'flat_norm_timeseries.txt' , np.vstack ([ times , mags + r_mean , mags_err]).T , header=' '.join(lc_dirs) )
	plt.show()