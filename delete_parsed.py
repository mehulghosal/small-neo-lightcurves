import numpy as np
import astropy as ap
from astropy.io import fits

import os
from os.path import isdir, isfile, join

directory = './'
dir_names = [directory+f+'/' for f in os.listdir(directory) if isdir(join(directory,f))] 

for d in dir_names:
	if '00_m_16' in d or 'pycache' in d or 'small_asteroid_lightcurve' in d : continue
	file_names = [d+f for f in os.listdir(d) if isfile(join(d,f))] 

	for f in file_names :

		try:
			file = fits.open(f)
			# print(f)
		except Exception as e: 
			# print(e , f)
			continue

		if 'on' in f: 
			print( f )
			os.remove(f)
