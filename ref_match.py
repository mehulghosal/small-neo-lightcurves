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

filt_ind = {'g':2, 'r': 3, 'i': 4}
mins = {'g':100, 'r': 150, 'i': 250}


for d in dir_names:
	lc_dirs = [d+f for f in os.listdir(d) if isdir(join(d,f))] 

	for ld in lc_dirs :

		if 'data/LCLIST' in ld or 'git' in ld: continue
		# print(ld)

		lc_files = [join(ld,f) for f in os.listdir(ld) if isfile(join(ld,f))]
		# print(os.listdir(ld))
		# print( lc_files )

		for f in lc_files :
			if not 'star_params' in f: continue
			if not 'GE1' in f: continue

			fits_name = ('/'.join(f.split('/')[:-1]) + '.flt')

			try:
				fits_file = fits.open(fits_name)
				print(fits_name)
			except Exception as e:
				print(f'NO FITS FILE FOUND: {fits_name}')
				continue

			hdr = fits_file[0].header
			img = fits_file[0].data

			exp_time   = float(hdr['EXPMEAS'])
			gain       = float(hdr['GAIN'])
			rd_noise   = float(hdr['RDNOISE']) 

			star_id , ra , dec , s , L , A , b , x , y , a , flux = np.loadtxt ( f , skiprows=1 , unpack=True )

			fig , ax = plt.subplots()
			ax.set_title(fits_name)
			ax.imshow(img, cmap='gray', norm=colors.LogNorm(vmin=mins[hdr['FILTER'][0]]))

			# transforming coordinates back to original frame
			x_0 , y_0 = reverse_rotation ( x , y , np.mean(a) , img)

			# transforming ra dec to x y
			w = WCS(hdr)
			ra_dec = SkyCoord ( ra=ra*u.degree , dec=dec*u.degree )
			# x_ , y_ = utils.skycoord_to_pixel ( SkyCoord(ra=ra*u.degree , dec=dec*u.degree) , w )


			# refcat magic
			# RA, Dec, g, r, i, z, J, cyan, orange.
			refcat = []

			args_str = f'./refcat {np.mean(ra_dec.ra.deg)} {np.mean(ra_dec.dec.deg)} -rad 1 -dir 00_m_16/'
			ref_stars = np.array(os.popen(args_str).read().split('\n')[:-1])
			for j in ref_stars:
				refcat.append(np.array(j.split(), dtype=float))
			refcat = np.array(refcat)
			# print(len(refcat))
			# end refcat magic

			ref_ra_dec = SkyCoord ( ra=refcat[:,0]*u.degree , dec=refcat[:,1]*u.degree)

			idx , d2d , d3d = ra_dec.match_to_catalog_sky (ref_ra_dec , nthneighbor=1)

			# fig_1 , ax1 = plt.subplots()
			hist , bins  = np.histogram(d2d.arcsec , bins=np.linspace(0 , 200 , 100) , range=[0,200])
			# ax1.set_xlabel('offset in arcsec')

			# basically converting bins --> integers so we find the mode. digitize gives me the index of which bin each d2d goes into
			# to be fair this is from stack overflow and it might be sketchy and untested
			binsd = bins[np.digitize ( d2d.arcsec , bins , right=False )-1] 
			# print(binsd)

			bins_mode = mode ( binsd  )[0]
			mode_err  = (bins_mode ** .5 )/2
			print(f'Mode offset: {bins_mode[0]} +/- {mode_err}')
			
			# now we constrain the offsets by +/- 1" around mode offset
			dist_filter = np.where ( (d2d.arcsec <= bins_mode + mode_err) & (d2d.arcsec >= bins_mode - mode_err)  )
			print(f'+/- 1 sigma from mode offsets in arcsec: {d2d.arcsec[dist_filter]}')

			matches = ref_ra_dec[idx[dist_filter]]

			ref_x , ref_y = utils.skycoord_to_pixel(matches , w)

			ax.scatter ( ref_x , ref_y , label='ref matches'  )
			ax.scatter ( x_0[dist_filter], y_0[dist_filter] , label='My stars matches')
			# ax.scatter ( x_0, y_0 , label='All my stars')

			# ax.scatter ( x_0[idx_[dist_filter_]] , y_0[idx_[dist_filter_]] , label='My stars')


			instrumental_mag = -2.5*np.log10(flux)
			ref_mag =  refcat[idx[dist_filter]] [:,filt_ind[hdr['FILTER'][0]]]

			fig_cal , ax_cal = plt.subplots( )
			ax_cal .scatter ( instrumental_mag[dist_filter] , ref_mag )

			param , param_cov = curve_fit ( line , instrumental_mag[dist_filter] , ref_mag )
			print('general line: ', param , np.diag(param_cov)**.5)
			ax_cal.plot (instrumental_mag[dist_filter] , line(instrumental_mag[dist_filter] , *param) , label=f'y={param[0]}x + {param[1]}')

			param_one , one_cov = curve_fit ( line_one , instrumental_mag[dist_filter] , ref_mag )
			print( 'line slope one: ',  param_one , np.diag(one_cov)**.5)
			ax_cal.plot (instrumental_mag[dist_filter] , line_one(instrumental_mag[dist_filter] , *param_one) , label=f'y=x + {param[0]}')


			ax.legend()
			ax_cal.legend()
			plt.tight_layout()
			print()



	plt.show()


		# if True: break

