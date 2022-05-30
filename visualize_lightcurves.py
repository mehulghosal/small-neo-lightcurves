import numpy as np
import matplotlib.pyplot as plt
import astropy as ap
from astropy.timeseries import LombScargle
from astropy.io import fits
import os
from os.path import isfile, join, isdir
from magic_star import bin_lightcurve, periodogram, fold_lightcurve

directory = './'	
dir_names = [directory+f+'/' for f in os.listdir(directory) if isdir(join(directory,f))] 

# directory  = './output/'
# file_names = np.array([directory+f for f in os.listdir(directory) if isfile(join(directory,f))], dtype=str)
count=0
if __name__ == '__main__':
	
	for d in dir_names:
		file_names = [d+f for f in os.listdir(d) if isfile(join(d,f))]
		yea = False

		if  not ('GE1' in d): continue
		# if not f_name in d: continue
		#if not ('2015_TG24' in d or '2016_NM15' in d or '2015_VH1' in d): continue

		for f in file_names:
			if '_lightcurve' not in f:
				continue
			# if count==5: continue
			# if not count==1: continue
			count+=1

			lc_file = np.loadtxt(f)

			# plt.figure()
			# plt.scatter(lc_file[:, 0], lc_file[:,1])

			plt.figure()
			time   = lc_file[:,0]
			errs   = lc_file[:,2]
			binned = bin_lightcurve(lc_file[:,1], int(len(lc_file)), np.nanmedian)

			binned = -2.5*np.log10(binned)
			# binned = lc_file[:,1]
			# plt.scatter(np.linspace(0, 60, len(binned)), binned)
			plt.scatter(time, binned)
			# plt.errorbar(np.linspace(0, 60, len(binned)), binned, yerr=errs)
			# plt.scatter(lc)
			plt.title(f)

			# period, power, peak_period = periodogram(np.arange(len(binned)), binned, num_maxes = 20)
			period, power, peak_period = periodogram(time, binned, num_maxes = 100)
			# plt.figure ()
			# plt.plot(period, power)
			# plt.title(f)
			# plt.xlim((0, 100))

			print('peak_period: ', peak_period[np.where(peak_period>1)][:5])
			actual_peak = peak_period[:] *2
			for T in peak_period:
				if T>14: 
					actual_peak=T
					break
			print(actual_peak)

			# folded_lightcurves, phase = fold_lightcurve(binned, time, actual_peak, exp_time=len(binned))
			# print('phase', phase)
			# plt.figure()
			# count_lc = 0
			# for lc in folded_lightcurves:
			# 	plt.scatter(np.arange(len(lc)), lc, label=count_lc)
			# 	count_lc+=1
			# plt.title(f)
			# plt.legend()

			# WORKS!!!
			folded_lightcurve = fold_lightcurve(binned, time, actual_peak, exp_time=time[-1]-time[0])




			# print( 'time: ', folded_lightcurve.time)
			# print(type(folded_lightcurve['data']))
			# print(folded_lightcurve.colnames )

			phase = np.array(folded_lightcurve['time'].value)
			data  = np.array(folded_lightcurve['data'])
			print(phase)
			print(data)

			plt.figure()
			# plt.scatter( np.array(folded_lightcurve['time']) , np.array(folded_lightcurve['data']) )
			plt.scatter( phase, data )

		plt.show()