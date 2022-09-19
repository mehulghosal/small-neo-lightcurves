import numpy as np
import astropy as ap
import matplotlib.pyplot as plt
from matplotlib import colors
from astropy.io import fits
from scipy.ndimage import rotate
from astropy.wcs import WCS
from magic_star import point_rotation
# to get absolute mag and orbital information
from astroquery.jplhorizons import Horizons
from PyAstronomy.pyTiming import pyPeriod
from magic_star import periodogram

plt.rcParams.update({'figure.max_open_warning': 0})

paperheight = 10
paperwidth = 13
margin = 1.0

fontsize_standard = 28

def bin_lightcurve ( time , flux , flux_err , n=2 ):
	t = []
	f = []
	f_e = []
	for i in range( 0 , len(time)-n , n ):
		t.append(time[i])
		f.append ( np.sum (flux[i:i+n]) )
		f_e.append( (np.sum ( flux_err[i:i+n]**2 ))**.5 )
	return np.array(t) , np.array(f) , np.array(f_e)

# initializing all directories
import os
from os.path import isdir, isfile, join
directory = './'	
dir_names = [directory+f+'/' for f in os.listdir(directory) if isdir(join(directory,f))] 

# output = open('output_rates.csv', 'w+')

input_file = np.loadtxt('input.csv', dtype=object, skiprows=1, usecols=(i for i in range(25)), delimiter=',')

mins = {'g':100, 'r': 150, 'i': 250}
start_times = []
for d in dir_names:
	lc_dirs = [d+f for f in os.listdir(d) if isdir(join(d,f))] 

	for ld in lc_dirs :

		if 'data/LCLIST' in ld or 'git' in ld: continue
		# print(ld)

		lc_files = [join(ld,f) for f in os.listdir(ld) if isfile(join(ld,f))]
		# print(os.listdir(ld))
		# print( lc_files )

		for f in lc_files :
			if not 'lightcurve_asteroid' in f: continue
			if not 'FF14' in f: continue
			try :
				jd, flux, flux_err , mag, mag_err = np.loadtxt(f , unpack=True, skiprows=1)
				# except:
				# 	print(f'failed {f}')
				# 	continue
				print()
				print(f)

				if not isdir (ld+'/lightcurve_figs/'): os.mkdir (ld+'/lightcurve_figs/')
				if not isdir (ld+'/periodogram_figs/'): os.mkdir (ld+'/periodogram_figs/')

				fapLevels = np.array([ 0.35, 0.05, 0.005])
				seconds_from_start = (jd - jd[0]) * 24 * 3600
				periods = np.linspace (0 , np.max(seconds_from_start) , 1000)

				# seconds_from_start = np.linspace ( 0 , 60 , len(mag))
				fig_lc , ax_lc = plt.subplots(figsize=((paperwidth*1.15) - 2 * margin, (paperheight*1.15) - 2 * margin))
				unbinned = ax_lc.errorbar( seconds_from_start , mag - np.mean(mag) , mag_err , fmt='s' , markerfacecolor='blue' , markeredgecolor='black' , ecolor='black' , capthick=2 , markersize=7 , capsize=3 )
				ax_lc.set_ylim([-1,1])
				ax_lc.invert_yaxis()
				ax_lc.set_xlabel('Seconds from start')
				ax_lc.set_ylabel('Instrumental Magnitude')
				plt.tight_layout()
				plt.savefig(ld+'/lightcurve_figs/asteroid_lightcurve-bin1.png')

				# clp = pyPeriod.Gls(( seconds_from_start , flux , flux_err ) , freq=1/periods , verbose=True )
				clp = pyPeriod.Gls(( seconds_from_start , flux , flux_err )  , verbose=True )
				plevel = clp.powerLevel(fapLevels)

				print(1/clp.freq)



				fig_p , ax_p = plt.subplots(figsize=((paperwidth*1.15) - 2 * margin, (paperheight*1.15) - 2 * margin))
				ax_p.plot( 1/clp.freq  , clp.power, color='grey' , label='PyAstronomy GLs')
				# ax_p.plot( period , power )
				ax_p.plot([min(1/clp.freq), max(1/clp.freq)], [plevel[0]]*2, '-.', c='black',lw=3,label =r'$1\sigma$')
				ax_p.plot([min(1/clp.freq), max(1/clp.freq)], [plevel[1]]*2, ':', c='black',lw=3,label =r'$2\sigma$')
				ax_p.plot([min(1/clp.freq), max(1/clp.freq)], [plevel[2]]*2, '..', c='black',lw=3,label =r'$3\sigma$')

				ax_p.legend()
				ax_p.set_xlabel('Period [s]')
				# ax_p.set_xlim([0 , 60])
				ax_p.set_xscale('log')
				ax_p.set_ylim([0 , 1])
				plt.tight_layout()
				plt.savefig(ld+'/periodogram_figs/periodogram_bin1.png')

				print( f'PyAstronomy GLs best period: {1/clp.freq[np.argmax(clp.power)]}s'  )
				plt.show()

				if True: break


				t_2 , lc_2 , lc_err_2 = bin_lightcurve ( seconds_from_start , flux , flux_err  , 2)
				mag_2 = -2.5 * np.log10 ( lc_2 )
				mag_err_2 = 1.08574 * lc_err_2 / lc_2 
				fig_lc_2 , ax_lc_2 = plt.subplots(figsize=((paperwidth*1.15) - 2 * margin, (paperheight*1.15) - 2 * margin))
				bin_2 = ax_lc_2.errorbar( t_2 , mag_2 - np.mean(mag_2) , mag_err_2 , fmt='s' , markerfacecolor='blue' , markeredgecolor='black' , ecolor='black' , capthick=2 , markersize=7 , capsize=3 )
				
				ax_lc_2.set_ylim([-1,1])
				ax_lc_2.invert_yaxis()
				ax_lc_2.set_xlabel('Seconds from start')
				ax_lc_2.set_ylabel('Instrumental Magnitude')
				plt.tight_layout()
				plt.savefig(ld+'/lightcurve_figs/asteroid_lightcurve-bin2.png')

				clp2 = pyPeriod.Gls(( t_2 , lc_2 , lc_err_2 ) , freq=1/periods)
				plevel2 = clp2.powerLevel(fapLevels)

				fig_p2 , ax_p2 = plt.subplots(figsize=((paperwidth*1.15) - 2 * margin, (paperheight*1.15) - 2 * margin))
				ax_p2.plot( 1/clp2.freq  , clp2.power, color='grey' , label='PyAstronomy GLs')
				# ax_p.plot( period , power )
				ax_p2.plot([min(1/clp2.freq), max(1/clp2.freq)], [plevel2[0]]*2, '-.', c='black',lw=3,label =r'$2\sigma$')
				ax_p2.plot([min(1/clp2.freq), max(1/clp2.freq)], [plevel2[1]]*2, ':', c='black',lw=3,label =r'$3\sigma$')
				ax_p2.legend()
				ax_p2.set_xlabel('Period [s]')
				ax_p2.set_xlim([0 , 60])
				ax_p2.set_ylim([0 , 1])
				plt.tight_layout()
				plt.savefig(ld+'/periodogram_figs/periodogram_bin2.png')
				print( f'PyAstronomy GLs best period: {1/clp2.freq[np.argmax(clp2.power)]}s'  )


				t_3 , lc_3 , lc_err_3 = bin_lightcurve ( seconds_from_start , flux , flux_err  , 3)
				mag_3 = -2.5 * np.log10 ( lc_3 )
				mag_err_3 = 1.08574 * lc_err_3 / lc_3 
				fig_lc_3 , ax_lc_3 = plt.subplots(figsize=((paperwidth*1.15) - 2 * margin, (paperheight*1.15) - 2 * margin))
				bin_3 = ax_lc_3.errorbar( t_3 , mag_3 - np.mean(mag_3) , mag_err_3 , fmt='s' , markerfacecolor='blue' , markeredgecolor='black' , ecolor='black' , capthick=2 , markersize=7 , capsize=3 )
				
				ax_lc_3.set_ylim([-1,1])
				ax_lc_3.invert_yaxis()
				ax_lc_3.set_xlabel('Seconds from start')
				ax_lc_3.set_ylabel('Instrumental Magnitude')
				plt.tight_layout()
				plt.savefig(ld+'/lightcurve_figs/asteroid_lightcurve-bin3.png')


				clp3 = pyPeriod.Gls(( t_3 , lc_3 , lc_err_3 ) , freq=1/periods)
				plevel3 = clp3.powerLevel(fapLevels)

				fig_p4 , ax_p4 = plt.subplots(figsize=((paperwidth*1.15) - 2 * margin, (paperheight*1.15) - 2 * margin))
				ax_p4.plot( 1/clp3.freq  , clp3.power, color='grey' , label='PyAstronomy GLs')
				# ax_p.plot( period , power )
				ax_p4.plot([min(1/clp3.freq), max(1/clp3.freq)], [plevel3[0]]*2, '-.', c='black',lw=3,label =r'$2\sigma$')
				ax_p4.plot([min(1/clp3.freq), max(1/clp3.freq)], [plevel3[1]]*2, ':', c='black',lw=3,label =r'$3\sigma$')
				ax_p4.legend()
				ax_p4.set_xlabel('Period [s]')
				ax_p4.set_xlim([0 , 60])
				ax_p4.set_ylim([0 , 1])
				plt.tight_layout()
				plt.savefig(ld+'/periodogram_figs/periodogram_bin3.png')
				print( f'PyAstronomy GLs best period: {1/clp3.freq[np.argmax(clp3.power)]}s'  )


				t_4 , lc_4 , lc_err_4 = bin_lightcurve ( seconds_from_start , flux , flux_err  , 4)
				mag_4 = -2.5 * np.log10 ( lc_4 )
				mag_err_4 = 1.08574 * lc_err_4 / lc_4 
				fig_lc_4 , ax_lc_4 = plt.subplots(figsize=((paperwidth*1.15) - 2 * margin, (paperheight*1.15) - 2 * margin))
				bin_4 = ax_lc_4.errorbar( t_4 , mag_4-np.mean(mag_4) , mag_err_4 , fmt='s' , markerfacecolor='blue' , markeredgecolor='black' , ecolor='black' , capthick=2 , markersize=7 , capsize=3 )
				ax_lc_4.set_ylim([-1,1])
				ax_lc_4.invert_yaxis()
				ax_lc_4.set_xlabel('Seconds from start')
				ax_lc_4.set_ylabel('Instrumental Magnitude')
				plt.tight_layout()
				plt.savefig(ld+'/lightcurve_figs/asteroid_lightcurve-bin4.png')

				clp4 = pyPeriod.Gls(( t_4 , lc_4 , lc_err_4 ) , freq=1/periods)
				plevel4 = clp4.powerLevel(fapLevels)

				fig_p5 , ax_p5 = plt.subplots(figsize=((paperwidth*1.15) - 2 * margin, (paperheight*1.15) - 2 * margin))
				ax_p5.plot( 1/clp4.freq  , clp4.power, color='grey' , label='PyAstronomy GLs')
				# ax_p.plot( period , power )
				ax_p5.plot([min(1/clp4.freq), max(1/clp4.freq)], [plevel4[0]]*2, '-.', c='black',lw=3,label =r'$2\sigma$')
				ax_p5.plot([min(1/clp4.freq), max(1/clp4.freq)], [plevel4[1]]*2, ':', c='black',lw=3,label =r'$3\sigma$')
				ax_p5.legend()
				ax_p5.set_xlabel('Period [s]')
				ax_p5.set_xlim([0 , 60])
				ax_p5.set_ylim([0 , 1])
				plt.tight_layout()
				plt.savefig(ld+'/periodogram_figs/periodogram_bin4.png')
				print( f'PyAstronomy GLs best period: {1/clp4.freq[np.argmax(clp4.power)]}s'  )


				t_5 , lc_5 , lc_err_5 = bin_lightcurve ( seconds_from_start , flux , flux_err  , 5)
				mag_5 = -2.5 * np.log10 ( lc_5 )
				mag_err_5 = 1.08574 * lc_err_5 / lc_5 
				fig_lc_5 , ax_lc_5 = plt.subplots(figsize=((paperwidth*1.15) - 2 * margin, (paperheight*1.15) - 2 * margin))
				bin_5 = ax_lc_5.errorbar( t_5 , mag_5-np.mean(mag_5) , mag_err_5 , fmt='s' , markerfacecolor='blue' , markeredgecolor='black' , ecolor='black' , capthick=2 , markersize=7 , capsize=3 )
				ax_lc_5.set_ylim([-1,1])
				ax_lc_5.invert_yaxis()
				ax_lc_5.set_xlabel('Seconds from start')
				ax_lc_5.set_ylabel('Instrumental Magnitude')
				plt.tight_layout()
				plt.savefig(ld+'/lightcurve_figs/asteroid_lightcurve-bin5.png')

				clp5 = pyPeriod.Gls(( t_5 , lc_5 , lc_err_5 ) , freq=1/periods)
				plevel5 = clp5.powerLevel(fapLevels)

				fig_p6, ax_p6 = plt.subplots(figsize=((paperwidth*1.15) - 2 * margin, (paperheight*1.15) - 2 * margin))
				ax_p6.plot( 1/clp5.freq  , clp5.power, color='grey' , label='PyAstronomy GLs')
				# ax_p.plot( period , power )
				ax_p6.plot([min(1/clp5.freq), max(1/clp5.freq)], [plevel5[0]]*2, '-.', c='black',lw=3,label =r'$2\sigma$')
				ax_p6.plot([min(1/clp5.freq), max(1/clp5.freq)], [plevel5[1]]*2, ':', c='black',lw=3,label =r'$3\sigma$')
				ax_p6.legend()
				ax_p6.set_xlabel('Period [s]')
				ax_p6.set_xlim([0 , 60])
				ax_p6.set_ylim([0 , 1])
				plt.tight_layout()
				plt.savefig(ld+'/periodogram_figs/periodogram_bin5.png')
				print( f'PyAstronomy GLs best period: {1/clp5.freq[np.argmax(clp5.power)]}s'  )

				plt.show()
			except:
				continue
		if True: break
			