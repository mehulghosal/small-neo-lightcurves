import numpy as np
import matplotlib.pyplot as plt
import astropy as ap
from astropy.timeseries import LombScargle
from astropy.io import fits
import os
from os.path import isfile, join, isdir
from testing import bin_lightcurve, periodogram, fold_lightcurve

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
		img_filters = []

		if  not ('GE1' in d): continue

		# if not f_name in d: continue
		#if not ('2015_TG24' in d or '2016_NM15' in d or '2015_VH1' in d): continue

		for f in file_names:

			if ('_lightcurve' not in f)  : continue
			print()

			fits_file  = fits.open(f[:-15] + '.flt')
			img_filter = fits_file[0].header['FILTER']
			img_filters.append(img_filter[0])
			print(img_filter)
			# if ('70o13' not in f) and  ('72o13' not in f) : continue
			# if 		: continue 
			# if count==5: continue
			count+=1

			lc_file = np.loadtxt(f)

			time_  = lc_file[:,0]
			times.append(time_)

			errs   = lc_file[:,2]

			# b      = 1.1
			# binned = bin_lightcurve(lc_file[:,1], int(len(lc_file)/b))
			# errs_b = bin_lightcurve(errs**2, int(len(lc_file)/b), np.sum) ** .5
			binned = lc_file[:,1]

			mag_b  = -2.5*np.log10(binned)
			mag_er = 1.08574 * errs / binned

			zp_filename = f[:-15]+'_zeropoint.txt'
			# print(zp_filename)
			zeropoint = np.loadtxt( zp_filename )
			zp_err    = zeropoint[1] 
			zeropoint = zeropoint[0]
			# zeropoint = 0
			# print(zeropoint)

			mag_cal = mag_b + zeropoint
			mag_cal_err = (mag_er**2 + zp_err**2)**.5

			# binned = -2.5*np.log10(binned)
			lightcurves.append(mag_cal)
			time   = np.linspace(time_[0], time_[-1], len(binned))
			
			fig = plt.figure()
			ax1 = fig.add_subplot(111)
			ax2 = ax1.twiny()
			# ax1.scatter(time, binned)
			# ax1.errorbar(time, binned, yerr=errs)
			# time_b = np.linspace(time[0], time[-1], len(binned))
			# ax1.errorbar(time, binned, errs, fmt='b.', linewidth=0, elinewidth=2)
			# ax1.errorbar(time, mag_b, mag_er, fmt='b.', linewidth=0, elinewidth=2)
			# ax1.scatter   (time, binned)
			# ax1.errorbar  (time, binned, errs, elinewidth=2)
			# ax1.scatter   (time_, mag_cal)
			ax1.errorbar  ( time_ , mag_cal , mag_cal_err , elinewidth=2 )

			ax1.set_xlabel('mjd')
			# ax2.plot(np.linspace(0, 60, 94), binned)
			ax2.set_xticks(np.linspace(0,60,6))
			ax2.set_xbound(lower=0, upper=60)
			ax2.set_xticklabels(np.linspace(0,60,6))
			# ax2.cla()
			# tit = ax1.set_title('2016 GE1 lightcurve')
			tit = ax1.set_title(f)
			tit.set_y(1.1)
			fig.subplots_adjust(top=0.85)
			ax2.set_xlabel('seconds')
			

			# period, power, peak_period = periodogram(np.arange(len(binned)), binned, num_maxes = 20)
			period, power, peak_period = periodogram(time_, mag_cal, num_maxes = 100)
			# period, power, peak_period = periodogram(time, binned, num_maxes = 100, err=errs)
			# period, power, peak_period = periodogram(time, mag_b, num_maxes = 100, err=errs_b)
			fig_per, ax_per = plt.subplots()
			ax_per.plot(period, power)
			ax_per.set_title('Lomb-Scargle periodogram')
			ax_per.set_xlim((2, 50))

			actual_peak = peak_period[:30] 
			for T in peak_period:
				if T>16 and T<19: 
					actual_peak=T*2
					break
			# actual_peak = 32.49
			print('actual peak [s]: ', actual_peak)

			# phase, data = fold_lightcurve( time, mag_cal, actual_peak)
			# # phase, data = fold_lightcurve( time_, mag_cal, actual_peak , exp_time=time[-1]-time[0])

			# fig_fold, ax_fold = plt.subplots()

			# ax_fold.scatter(phase*24*3600+actual_peak/2, data)
			# # ax_fold.scatter(phase, data)
			# # plt.scatter(np.linspace(-actual_peak, actual_peak, len(data)), data)
			# ax_fold.set_xlabel('seconds')
			# ax_fold.set_ylabel('Mag')
			# ax_fold.set_title (f'Lightcurve folded on {actual_peak:.2f} seconds')
			# if True: break

		
		times = np.array(times, dtype=object)
		img_filters = np.array(img_filters, dtype=object)

		start_times = []
		for t in times: start_times.append(t[0])
		ind = np.argsort(start_times)

		lightcurves = np.array(lightcurves, dtype=object)
		norms       = np.array([np.median(i) for i in lightcurves])
		
		img_filters = img_filters[ind]
		print(img_filters)

		combined_lc = np.concatenate(lightcurves[ind])
		# combined_lc = np.concatenate(lightcurves[ind]/norms[ind])
		combined_T  = np.concatenate(times[ind])

		fig_com, ax_com = plt.subplots()

		# print(combined_T.shape, combined_lc.shape)
		ax_com.scatter(combined_T, combined_lc)
		ax_com.set_title('combined lightcurves')


		# combined_period, combined_power, combined_peak = periodogram( combined_T , combined_lc , num_maxes = 50 )
		# fig_c_per, ax_c_per = plt.subplots()
		# ax_c_per.plot(combined_period , combined_power)
		# print(combined_peak[combined_peak>10])

		# P = combined_peak[combined_peak>10] [0] *2
		# print(P)
		# phase , data = fold_lightcurve(combined_T , combined_lc , P)
		# fig_fold , ax_fold = plt.subplots()
		# ax_fold.scatter(phase, data)

		plt.show()