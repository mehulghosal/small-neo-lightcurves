import warnings, subprocess, sys
import numpy as np
import astropy as ap

from astropy.time import Time
from astropy.table import Table
from astropy.timeseries import LombScargle
from astropy.timeseries import TimeSeries
from astropy.io import fits
from scipy.ndimage import rotate
from scipy.special import erf
from astropy.wcs import WCS
from astropy.wcs import utils
from astropy.coordinates import SkyCoord
from astropy import units as u

from astropy.utils.exceptions import AstropyWarning

import matplotlib.pyplot as plt


# initializing all directories
import os
from os.path import isdir, isfile, join
directory = './'	
dir_names = [directory+f+'/' for f in os.listdir(directory) if isdir(join(directory,f))] 



# next thirty lines are from magic_star.py : just iterating through the directories and finding fits files

if __name__ == '__main__':
	
	for d in dir_names:
		file_names = [d+f for f in os.listdir(d) if isfile(join(d,f))]
		yea = False

		# if '00_m_16' in d or 'small_asteroid_lightcurve' in d: continue
		# if not ('HD3' in d or 'VH65' in d or 'XD169' in d or 'XR169' in d or 'CE3' in d or 'CG18' in d or 'CW30' in d or 'EN156' in d or 'JB' in d or 'LT1' in d): continue
		if not ('EN156' in d or 'EL157' in d or 'FF14' in d or 'GE1' in d) : continue

		start_times = []
		lightcurves = []
		errors      = []

		for f in file_names:

			try:
				file = fits.open(f)
				print(f)
			except Exception as e:
				print(f , e)
				continue

			output_dir_name = f'{f[:-4]}/'
			if 'on' in f:
				output_dir_name = f'{f[:-5]}/'
			if not isdir(output_dir_name):
				os.mkdir(output_dir_name)

			sex_output_name = str(output_dir_name+'sex.cat')

			sex = subprocess.run(['sex', f, '-DETECT_MINAREA', str(200), '-CATALOG_NAME' , sex_output_name ], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

			img = file[0].data
			# fig, ax = plt.subplots()
			# ax.imshow(img , vmin=250 , vmax=1000 , cmap='gray')

			sex_output = np.loadtxt ( sex_output_name , skiprows=9 )

			x = sex_output [ : , 5]
			y = sex_output [ : , 6]

			# ax.scatter ( x , y) 
			# plt.show()

		# 	if True: break
		# if True: break