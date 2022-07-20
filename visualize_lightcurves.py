import numpy as np
import matplotlib.pyplot as plt
import astropy as ap
from astropy.timeseries import LombScargle
from astropy.io import fits
import os
from os.path import isfile, join, isdir
from magic_star import bin_lightcurve, periodogram, fold_lightcurve, normalize_lightcurves, periodogram_xo
from gp_test import gaussian_process

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

		if not ('LT1' in d): continue
		# count += 1
		# if   ('NM15' in d): continue

		# if not f_name in d: continue
		#if not ('2015_TG24' in d or '2016_NM15' in d or '2015_VH1' in d): continue

		for f in file_names:

			if ('_lightcurve' not in f)  : continue
			print()

			fits_file  = fits.open(f[:-15] + '.flt')
			img_filter = fits_file[0].header['FILTER']
			img_filters.append(img_filter[0])
			print(f)
			print(img_filter)
			# if ('70o13' not in f) and  ('72o13' not in f) : continue
			# if 		: continue 
			# if count==5: continue
			# count+=1

			lc_file = np.loadtxt(f)

			time_  = lc_file[:,0]
			times.append(time_)

			errs   = lc_file[:,2]

			binned = lc_file[:,1]

			mag_b  = -2.5*np.log10(binned)
			mag_er = 1.08574 * errs / binned

			zp_filename = f[:-15]+'_zeropoint.txt'
			# print(zp_filename)
			# zeropoint = np.loadtxt( zp_filename )
			# zp_err    = zeropoint[1] 
			# zeropoint = zeropoint[0]
			zeropoint,zp_err = 0,0

			# print(zeropoint)

			mag_cal = mag_b + zeropoint
			mag_cal_err = (mag_er**2 + zp_err**2)**.5

			# binned = -2.5*np.log10(binned)
			lightcurves.append(mag_cal)
			time   = np.linspace(time_[0], time_[-1], len(binned))
			
			fig = plt.figure()
			ax1 = fig.add_subplot(111)
			ax2 = ax1.twiny()

			ax1.scatter ( time_ , mag_cal - np.mean(mag_cal) , color='blue' , s=5 )
			# ax1.errorbar  ( time_ , mag_cal - np.mean(mag_cal) , mag_cal_err , fmt='b.' , elinewidth=2 , linewidth=0 )

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
			
			# period, power, peak_period = periodogram(time_, mag_cal, num_maxes = 100 , err = mag_cal_err)
			period, power, peak_period = periodogram(time_, mag_cal, num_maxes = 150 )

			# print(period, power)

			# period, power, peak_period = periodogram_xo(time_ , mag_cal , num_maxes = 100 , err = mag_cal_err , min_period=1e-6 , max_period=.5)

			fig_per, ax_per = plt.subplots()
			ax_per.plot(period, power)
			# ax_per.set_title('Lomb-Scargle periodogram')
			ax_per.set_xlim((2, 50))
			# ax_per.set_ylim((0, 1))

			actual_peak = 0
			for T in peak_period:
				if (T>14 and T<16) or ( T>28 and T<40 ) :  
					actual_peak=T*2
					break
			# actual_peak = 32.49
			# actual_peak = peak_period * 3600
			# print
			print('actual peak [s]: ', actual_peak)

			# map_soln = gaussian_process( time_ , mag_cal , mag_cal_err , actual_peak )
			# print(map_soln)
			# ax1.plot(time_ , map_soln["pred"], label='GP model')

			# '''

			# phase, data = fold_lightcurve( time, mag_cal, actual_peak)
			# phase, data, errs_folded = fold_lightcurve( time_ , mag_cal , actual_peak , mag_cal_err )

			# fig_fold, ax_fold = plt.subplots()

			# # ax_fold.scatter(phase*24*3600+actual_peak/2, data)
			# ax_fold.errorbar(phase*24*3600+actual_peak/2, data - np.median(mag_cal), errs_folded, fmt='b.' , elinewidth=2 , linewidth=0 )

			# # ax_fold.scatter(phase, data)
			# # plt.scatter(np.linspace(-actual_peak, actual_peak, len(data)), data)
			# ax_fold.set_xlabel('seconds')
			# ax_fold.set_ylabel('Mag')
			# ax_fold.set_title (f'Lightcurve folded on {actual_peak:.2f} seconds')
			# if True: break
			'''

		# '''
		
		# plt.show()
		times       = np.array(times, dtype=object)
		img_filters = np.array(img_filters, dtype=object)
		lightcurves = np.array(lightcurves, dtype=object)

		normed_lcs , norms = normalize_lightcurves(lightcurves) 

		start_times = [t[0] for t in times]
		# for t in times: start_times.append(t[0])
		ind = np.argsort(start_times)

		# print(np.array(norms)[ind])
		normed_lcs = np.array(normed_lcs)[ind]
		norms      = np.array(norms     )[ind]

		img_filters = img_filters[ind]

		combined_lc = np.concatenate(lightcurves[ind])
		normed_lcs  = np.concatenate(normed_lcs[ind])
		# combined_lc = np.concatenate(lightcurves[ind]/norms[ind])
		combined_T  = np.concatenate(times[ind])


		fig_com, ax_com = plt.subplots()
		fig_com_norm , ax_com_norm = plt.subplots()

		# print(combined_T.shape, combined_lc.shape)
		# ax_com.scatter(combined_T, combined_lc)
		# ax_com.set_title('combined, calibrated lightcurves')
		ax_com.scatter(combined_T , combined_lc)
		ax_com_norm.scatter(combined_T , normed_lcs)
		# ax_com_norm.set_title('combined, normalized lightcurves')
		# '''


		combined_period, combined_power, combined_peak = periodogram( combined_T , combined_lc , num_maxes = 30 )
		fig_c_per, ax_c_per = plt.subplots()
		ax_c_per.plot(combined_period , combined_power)
		ax_c_per.set_xlim((2, 1000))

		# P = combined_peak[combined_peak>10] [0] *2
		print(combined_peak[combined_peak>10])
		# actual_peak = 0
		# '''

		# fig_fold , ax_fold = plt.subplots()
		for T in combined_peak:
			if T>100 and T<400: 
				actual_peak=T*2
				print('combined peak period [s] : ' , actual_peak)
				phase , data, errs = fold_lightcurve(combined_T , normed_lcs , actual_peak)
				plt.figure()
				plt.title(actual_peak)
				plt.errorbar(phase * 24*3600 + actual_peak/2, data, errs, fmt='b.' , elinewidth=2 , linewidth=0 )

		# '''
		# plt.show()