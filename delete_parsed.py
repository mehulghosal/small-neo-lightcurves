import numpy as np
import astropy as ap
from astropy.io import fits

import os, shutil
from os.path import isdir, isfile, join

directory = './'
dir_names = [directory+f+'/' for f in os.listdir(directory) if isdir(join(directory,f))] 

for d in dir_names:
	if '00_m_16' in d or 'pycache' in d or 'small_asteroid_lightcurve' in d : continue
	d_names = [d+f for f in os.listdir(d) if isdir(join(d,f))] 

	for d2 in d_names :

		# print(d2)

		if d2[-2] == 'o':

			print(d2)
			# shutil.rmtree(d2)

		# if 'on' in f: 
		# 	print( f )
		# 	os.remove(f)


