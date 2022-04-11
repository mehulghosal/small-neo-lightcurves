import numpy as np
import matplotlib.pyplot as plt
import astropy as ap
from astropy.timeseries import LombScargle
from astropy.io import fits
import os
from os.path import isdir, isfile, join

directory  = './'
dir_names = [directory+f+'/' for f in os.listdir(directory) if isdir(join(directory,f))] 

for d in dir_names:
	if 'GE1' not in d: continue
	file_names = [d+f for f in os.listdir(d) if isfile(join(d,f))]

	lc_files = []
	for f in file_names:
		# if ('lightcurve' not in f) and ('.txt' in f) and ('params' not in f) :
		if ('lightcurve' in f):
			lc_files.append(f)

	for f in lc_files:
		# lightcurve = np.loadtxt(f)
		# for lc in lightcurve:
		# 	plt.figure()
		# 	plt.scatter(np.arange(0,len(lc),1), lc)
		# 	plt.title(f)
		lightcurve = np.loadtxt(f)
		t = lightcurve[0]
		i = lightcurve[1]
		print(len(t))

		plt.figure()
		plt.scatter(t, i)
	plt.show()