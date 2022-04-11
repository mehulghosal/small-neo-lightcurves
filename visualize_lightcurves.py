import numpy as np
import matplotlib.pyplot as plt
import astropy as ap
from astropy.timeseries import LombScargle
from astropy.io import fits
import os
from os.path import isfile, join

directory  = './output/'
file_names = np.array([directory+f for f in os.listdir(directory) if isfile(join(directory,f))], dtype=str)

lightcurves = []
indices     = []

for i in file_names:
	if 'indices' in i:
		indices.append(i)
	else: lightcurves.append(i)


print(len(lightcurves), len(indices))

for i in range(len(indices)):
	lc = lightcurves[i]
	f = np.loadtxt(lc)
	time = f[:,0]
	intensity = f[:,1]

	index_array = np.loadtxt(indices[i])

	rolling_ind = 0
	p = []
	if True:
	# if 'GE1' in lc:
		for index in index_array:
			start = int(rolling_ind)
			end   = int(start + index - 1)
			rolling_ind += i
			frequency, power = LombScargle(time[start:end], intensity[start:end]).autopower()
			peak_frequency   = frequency[np.argmax(power)]
			peak_period      = 1/peak_frequency * 24 * 3600
			print('peak period: ', peak_period )
			p.append(peak_period)
			# plt.figure()
			# plt.scatter(time[start:end], intensity[start:end])
			# plt.title(lc)
		print(np.median(p), np.mean(p))
		frequency, period = LombScargle(time, intensity).autopower()
		peak_frequency   = frequency[np.argmax(power)]
		peak_period      = 1/peak_frequency * 24 * 3600
		print('peak period: ', peak_period )
		plt.figure()
		plt.scatter(time, intensity)
		plt.title(lc)
		print()
	plt.show()