import numpy as np
import matplotlib.pyplot as plt
import astropy as ap
from astropy.timeseries import LombScargle
from astropy.io import fits
import os
from os.path import isfile, join

directory  = './output/'
file_names = np.array([directory+f for f in os.listdir(directory) if isfile(join(directory,f))], dtype=str)
print(file_names)

lightcurves = []
indices     = []

for i in file_names:
	if 'indices' in i:
		indices.append(i)
	else: lightcurves.append(i)


print(len(lightcurves), len(indices))

for i in range(len(lightcurves)):
	lc = lightcurves[i]
	f = np.loadtxt(lc)
	time = f[:,0]
	intensity = f[:,1]

	# index_array = np.loadtxt(indices[i])

	rolling_ind = 0
	if False:
		for index in index_array:
			start = int(rolling_ind)
			end   = int(start + index - 1)
			rolling_ind += i
			plt.figure()
			plt.scatter(time[start:end], intensity[start:end])
			plt.title(lc)
	else:
		plt.figure()
		plt.scatter(time, intensity)
		plt.title(lc)
	plt.show()