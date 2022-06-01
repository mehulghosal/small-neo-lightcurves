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

		if  not ('CD31' in d): continue
		# if not f_name in d: continue
		#if not ('2015_TG24' in d or '2016_NM15' in d or '2015_VH1' in d): continue

		for f in file_names:
			if '_lightcurve' not in f:
				continue
			# if count==5: continue
			count+=1
			# if not count==1: continue


			lc_file = np.loadtxt(f)

			time   = lc_file[:,0] 
			errs   = lc_file[:,2]
			binned = bin_lightcurve(lc_file[:,1], int(len(lc_file)), np.nanmedian)

			binned = -2.5*np.log10(binned)
			
			fig = plt.figure()
			ax1 = fig.add_subplot(111)
			ax2 = ax1.twiny()
			ax1.scatter(time, binned)
			ax1.set_xlabel('mjd')
			# ax2.plot(np.linspace(0, 60, 94), binned)
			ax2.set_xticks(np.linspace(0,60,6))
			ax2.set_xbound(lower=0, upper=60)
			ax2.set_xticklabels(np.linspace(0,60,6))
			# ax2.cla()
			tit = ax1.set_title(f)
			tit.set_y(1.1)
			fig.subplots_adjust(top=0.85)
			ax2.set_xlabel('seconds')
			

			# period, power, peak_period = periodogram(np.arange(len(binned)), binned, num_maxes = 20)
			period, power, peak_period = periodogram(time, binned, num_maxes = 100)
			plt.figure ()
			plt.plot(period, power)
			plt.title(f)
			plt.xlim((0, 100))

			print('peak_period: ', peak_period[np.where(peak_period>1)][:5])
			actual_peak = peak_period[:5] * 1
			for T in peak_period:
				if T>28: 
					actual_peak=T
					break
			print(actual_peak)
			# actual_peak = [15]


			# WORKS!!!
			# folded_lcs = []
			# for period in actual_peak:
			# 	print('period: ', period)
			# 	phase, data = fold_lightcurve(binned, time, period, exp_time=time[-1]-time[0])
			# 	folded_lcs.append(data)
			# 	print(phase)
			# 	print(data)

			# 	plt.figure()
			# 	# plt.scatter( np.array(folded_lightcurve['time']) , np.array(folded_lightcurve['data']) )
			# 	plt.scatter( phase, data )
			# 	plt.title(period)


			phase, data = fold_lightcurve(binned, time, actual_peak, exp_time=time[-1]-time[0])


			plt.figure()
			plt.scatter(phase*24*3600, data)
			# plt.scatter(np.linspace(-actual_peak, actual_peak, len(data)), data)
			plt.xlabel('seconds')
			plt.ylabel('Mag')
			plt.title (f'Folded on {actual_peak:.2f} seconds')

		plt.show()