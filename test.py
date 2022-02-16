import numpy as np
import astropy as ap
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib import colors

from astropy.wcs import WCS
from astropy.wcs import utils
from astropy.coordinates import SkyCoord
from astropy import units as u

from astroquery.jplhorizons import Horizons

plt.rcParams.update({'figure.max_open_warning': 0})

import os
from os.path import isfile, join
directory = './2016_GE1_2016_04_04_UTC/'
# directory = './2016_EV84_2016_03_12_2016_UTC/'
file_names = [directory+f for f in os.listdir(directory) if isfile(join(directory,f))] 
mins = {'g':75, 'r': 175, 'i': 230}


for f in file_names:
	try:
		file = fits.open(f)
	except Exception as e:
		continue
	
	hdr = file[0].header
	img = file[0].data

	# if 'o22' not in f: continue
	plt.figure()
	# plt.imshow(img, cmap='magma', vmin=0, vmax=750)
	plt.imshow(img, cmap='magma', norm=colors.LogNorm(vmin=mins[hdr['FILTER'][0]]))
	# plt.imshow(np.log(img), cmap='gray', vmin=4.6)


	plt.title(f)
	

plt.show()