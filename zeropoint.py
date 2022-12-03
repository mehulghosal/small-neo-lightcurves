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

def line( x , m , b): return m * x + b

def line_one ( x , b): return x + b

def quadratic ( x , a , b , c) : return a * x ** 2 + b * x + c

def exponential ( x , A , b , c): return A * np.exp(x * b) + c

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
	if not 'EV84' in d: continue
	
	fig, ax = plt.subplots(figsize=((paperwidth*1.15) - 2 * margin, (paperheight*1.15) - 2 * margin))
	mags , mags_err , ref_mag , ref_mag_err = [] , [] , [] , [] 

	img_name = ''


	for ld in lc_dirs :

		if 'data/LCLIST' in ld or 'git' in ld: continue
		# if not ('66' in ld): continue

		# print(ld)

		lc_files = [join(ld,f) for f in os.listdir(ld) if isfile(join(ld,f))]

		# print( lc_files )

		
		# print(lc_files)

		for f in lc_files :

			if not 'ref' in f: continue
			
			if not '15o' in f: continue


			fits_name = ('/'.join(f.split('/')[:-1]) + '.flt')
			if 'on' in f : fits_name = ('/'.join(f.split('/')[:-1]) + '.fits')
			else: continue


			try:
				fits_file = fits.open(fits_name)
				img_name  = f.split('/')[2].split('on')[0]
				print(fits_name)
			except Exception as e:
				print(f'NO FITS FILE FOUND: {fits_name}')
				continue

			hdr = fits_file[0].header
			img = fits_file[0].data

			# exp_time   = float(hdr['EXPMEAS'])
			# gain       = float(hdr['GAIN'])
			# rd_noise   = float(hdr['RDNOISE']) 

			try:
				star_id , mag , mag_err , g , g_err , r , r_err , i , i_err = np.loadtxt ( f , skiprows=1 , unpack=True)
			except Exception as e: 
				print( f'NO REF STARS: {fits_name}')
				continue


			ref , ref_err = 0 , 0
			filt = hdr['FILTER']
			if 'g' in filt:
				ref , ref_err = g , g_err
			elif 'r' in filt:
				ref , ref_err = r , r_err
			elif 'i' in filt: 
				ref , ref_err = i , i_err

			try:
				for i in range(len(star_id)):
					mags.append(mag[i])
					mags_err.append(mag_err[i])
					ref_mag.append(ref[i])
					ref_mag_err.append(ref_err[i])
			except Exception as e:
				mags.append(mag)
				mags_err.append(mag_err)
				ref_mag.append(ref)
				ref_mag_err.append(ref_err)
			
	mags = np.array(mags).flatten()
	# print(mags)
	mags_err = np.array(mags_err).flatten()
	ref_mag = np.array(ref_mag).flatten()
	ref_mag_err = np.array(ref_mag_err).flatten()

	# print( mags , mags.shape)

	ax.errorbar ( mags , ref_mag , ref_mag_err , mags_err , fmt='s' , markerfacecolor='blue' , markeredgecolor='black' , ecolor='black' , capthick=2 , markersize=7 , capsize=3  )

	try:
		param , param_cov = curve_fit ( line_one , mags[mags < 1e10] , ref_mag[mags < 1e10] , sigma=ref_mag_err[mags < 1e10] , absolute_sigma=True )
		print( 'line slope one: ',  param[0] , np.diag(param_cov)[0]**.5)
	except Exception as e:
		print(e)
		print(mags)
		print(ref_mag)
		print(ref_mag_err)
		# continue

	ax.plot (mags , line_one(mags , *param) , label=f'y=x + {param[0]:.1f}' , color='blue')

	s = 100

	outliers_filter = np.where( (ref_mag > line_one(mags , *param) - s*np.diag(param_cov)[0]**.5) & (ref_mag < line_one(mags , *param) + s*np.diag(param_cov)[0]**.5) )
	# print(outliers_filter)

	ax.errorbar ( mags[outliers_filter] , ref_mag[outliers_filter] , ref_mag_err[outliers_filter] , mags_err[outliers_filter] , fmt='s' , markerfacecolor='red' , markeredgecolor='black' , ecolor='black' , capthick=2 , markersize=7 , capsize=3  )
	param_1 , param_cov_1 = curve_fit ( line_one , mags [outliers_filter], ref_mag[outliers_filter] , sigma=ref_mag_err[outliers_filter] , absolute_sigma=True )
	print( 'line, outliers removed: ',  param[0] , np.diag(param_cov)[0]**.5)

	ax.plot (mags , line_one(mags , *param_1) , label=f'y=x + {param_1[0]:.1f}' , color='red')

	plt.show()

	print(d+ img_name+'.zp')

	np.savetxt ( d+ img_name+'.zp' , [param[0] , np.diag(param_cov)[0]**.5] )


	# if True: break

