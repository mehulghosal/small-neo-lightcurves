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
	lc_dirs = [d+f for f in os.listdir(d) if isdir(join(d,f))]

	for ld in lc_dirs :

		if 'data/LCLIST' in ld or 'git' in ld: continue
		print(ld)

		lc_files = [join(ld,f) for f in os.listdir(ld) if isfile(join(ld,f))]
		# print(os.listdir(ld))
		# print( lc_files )

		for f in lc_files :
			if 'asteroid' in f: continue
			print(f)
			arr  = np.loadtxt ( f )
			L    = len(arr)//3
			mjd  = arr[:L]
			flux = arr[L:2*L]
			f_er = arr[2*L:3*L]

			to_write = np.array ( [mjd , flux , f_er ] ).T
			print( to_write.shape )
			np.savetxt ( f , to_write)
		# 	if True: break
		# if True: break

	# if not 'LT1_2016_06_07' in d: continue