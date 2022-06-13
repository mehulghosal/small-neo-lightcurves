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
		lightcurves = []
		times       = []

		if  not ('GG1' in d): continue
		# if not f_name in d: continue
		#if not ('2015_TG24' in d or '2016_NM15' in d or '2015_VH1' in d): continue

		for f in file_names:
			if '_lightcurve' not in f:
				continue
			# if count==5: continue
			count+=1
			# if not count==5: continue

			lc_file = np.loadtxt(f)

			time_  = lc_file[:,0]
			times.append(time_)

			errs   = lc_file[:,2]

			b      = 1.5
			binned = bin_lightcurve(lc_file[:,1], int(len(lc_file)/b), np.median)
			errs_b = bin_lightcurve(errs**2, int(len(lc_file)/b), np.sum) ** .5

			mag_b  = -2.5*np.log10(binned)
			mag_er = errs_b * (1.08574 / binned) ** .5 

			# binned = -2.5*np.log10(binned)
			lightcurves.append(binned)
			time   = np.linspace(time_[0], time_[-1], len(binned))
			
			fig = plt.figure()
			ax1 = fig.add_subplot(111)
			ax2 = ax1.twiny()
			# ax1.scatter(time, binned)
			# ax1.errorbar(time, binned, yerr=errs)
			# time_b = np.linspace(time[0], time[-1], len(binned))
			# ax1.errorbar(time, binned, errs_b, fmt='b.', linewidth=0, elinewidth=2)
			# ax1.errorbar(time, mag_b, mag_er, fmt='b.', linewidth=0, elinewidth=2)
			# ax1.scatter   (time, binned)
			ax1.scatter   (time, mag_b)

			ax1.set_xlabel('mjd')
			# ax2.plot(np.linspace(0, 60, 94), binned)
			ax2.set_xticks(np.linspace(0,60,6))
			ax2.set_xbound(lower=0, upper=60)
			ax2.set_xticklabels(np.linspace(0,60,6))
			# ax2.cla()
			tit = ax1.set_title('2016 GE1 lightcurve')
			tit.set_y(1.1)
			fig.subplots_adjust(top=0.85)
			ax2.set_xlabel('seconds')
			

			# period, power, peak_period = periodogram(np.arange(len(binned)), binned, num_maxes = 20)
			# period, power, peak_period = periodogram(time, binned, num_maxes = 100)
			# period, power, peak_period = periodogram(time, binned, num_maxes = 100, err=errs)
			period, power, peak_period = periodogram(time, mag_b, num_maxes = 100, err=errs_b)
			plt.figure ()
			plt.plot(period, power)
			plt.title('Lomb-Scargle periodogram')
			plt.xlim((2, 50))

			print('peak_period: ', peak_period[np.where(peak_period>1)][:5])
			actual_peak = peak_period[:15] * 1
			for T in peak_period:
				if T>15 and T<17: 
					actual_peak=T*2
					break
			print(actual_peak)
			# plt.axvline(x=actual_peak/2, color='g', linestyle='-', label=f'peak period: {actual_peak/2:.2f} s ')
			plt.legend()
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


			# phase, data = fold_lightcurve(mag_b, time, actual_peak, exp_time=time[-1]-time[0])


			# plt.figure()
			# plt.scatter(phase*24*3600+actual_peak/2, data)
			# # plt.scatter(np.linspace(-actual_peak, actual_peak, len(data)), data)
			# plt.xlabel('seconds')
			# plt.ylabel('Mag')
			# plt.title (f'Lightcurve folded on {actual_peak:.2f} seconds')

		# print(np.array(times))
		# ind = np.argsort(np.array(times, dtype=object)[0])
		# print(ind)
		# lightcurves = np.array(lightcurves)[ind]
		# combined_lc = np.concatenate(lightcurves)
		# combined_T  = np.concatenate(np.array(times)[ind])
		# plt.figure()
		# plt.scatter(combined_T, combined_lc)

		plt.show()