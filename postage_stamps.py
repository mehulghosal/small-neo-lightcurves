import numpy as np
import astropy as ap
import matplotlib.pyplot as plt
from matplotlib import colors
from astropy.io import fits
from scipy.ndimage import rotate
from astropy.wcs import WCS, utils
from magic_star import point_rotation, reverse_rotation, trail_view
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

			print(img.shape)

			exp_time   = float(hdr['EXPMEAS'])
			gain       = float(hdr['GAIN'])
			rd_noise   = float(hdr['RDNOISE']) 

			star_id , ra , dec , s , L , A , b , x , y , a , flux = np.loadtxt ( f , skiprows=1 , unpack=True )


			img_rotated = rotate(img , a[0])

			fig , ax = plt.subplots()
			ax.set_title(fits_name)
			ax.imshow(img_rotated, cmap='gray', norm=colors.LogNorm(vmin=mins[hdr['FILTER'][0]]))
			ax.scatter ( x , y )

			for i in range( len(star_id) ):
				plt.figure()
				trail = trail_view ( img_rotated , s[i]*1.2 , L[i]*1.2 , A[i] , b[i] , x[i] , y[i] )
				plt.imshow( trail , )
				plt.title(f'star #{i}, s={s[i]:.1f}, L={L[i]:.1f} , x={x[i]:.1f} , y={y[i]:.1f}')
				plt.tight_layout()

			plt.show()
		if True: break