import numpy as np
import astropy as ap
import matplotlib.pyplot as plt
from matplotlib import colors
from astropy.io import fits
from scipy.ndimage import rotate
from scipy.stats import mode
from scipy.optimize import curve_fit
from astropy.wcs import WCS, utils
from magic_star import point_rotation, reverse_rotation
from astropy.coordinates import SkyCoord
from astropy import units as u
from astroquery.jplhorizons import Horizons
from PyAstronomy.pyTiming import pyPeriod

plt.rcParams.update({'font.size': 22})

def mag2flux ( mag , mag_err ): 
	f = 10 ** ( mag/ (-2.5) )
	return f , mag_err * f / 1.0875

def flux2mag ( flux , flux_err ): 
	return -2.5 * np.log10 ( flux ) , 1.0875 * flux_err / flux

def bin_lightcurve_ ( time , flux , flux_err , n=2 ):

	flux = mag2flux ( flux , flux_err )

	t = []
	f = []
	f_e = []
	for i in range( 0 , len(time)-n , n ):
		t.append(time[i])
		f.append ( np.sum (flux[i:i+n]) )
		f_e.append( (np.sum ( flux_err[i:i+n]**2 ))**.5 )
	return (np.array(t) , *flux2mag(np.array(f), np.array(f_e)) )

def bin_lightcurve ( time , mag , mag_err , n=2 ):

	t = []
	f = []
	f_e = []
	for i in range( 0 , len(time)-n , n ):
		t.append(time[i])
		f.append ( np.average ( mag[i:i+n] , weights=1/mag_err[i:i+n]**2 ) )
		f_e.append( (np.sum ( mag_err[i:i+n]**2 ))**.5 )
	return np.array(t) , np.array(f), np.array(f_e)  

def line ( x , m , b): return m * x + b

plt.rcParams.update({'figure.max_open_warning': 0})

paperheight = 10
paperwidth = 13
margin = 1.0

fontsize_standard = 28

# initializing all directories
import os
from os.path import isdir, isfile, join
directory = './'	
dir_names = [directory+f+'/' for f in os.listdir(directory) if isdir(join(directory,f))] 

filt_ind = {'g':21, 'r': 25, 'i': 29}
mins = {'g':100, 'r': 150, 'i': 250}

r_mags = []
r_errs = []
r_time = []
g_mags = []
g_errs = []
g_time = []
i_mags = []
i_errs = []
i_time = []

for d in dir_names:
	lc_dirs = [d+f for f in os.listdir(d) if isdir(join(d,f))] 
	if not 'FF14' in d: continue

	times , mags , mags_err, uncor  = [] , [] , [] , []
	fig_combined, ax_combined = plt.subplots(figsize=((paperwidth*1.15) - 2 * margin, (paperheight*1.15) - 2 * margin))


	img_name = ''	

	for ld in lc_dirs :

		if 'data/LCLIST' in ld or 'git' in ld: continue

		lc_files = [join(ld,f) for f in os.listdir(ld) if isfile(join(ld,f))]

		for f in lc_files :

			if not 'calibrated_lightcurve.txt' in f: continue
			if '52o' in f  : continue
			
			fits_name = ('/'.join(f.split('/')[:-1]) + '.flt')

			try:
				fits_file = fits.open(fits_name)
				img_name  = f.split('/')[2].split('o')[0]
				print(fits_name)
			except Exception as e:
				print(f'NO FITS FILE FOUND: {fits_name}')
				continue

			hdr = fits_file[0].header
			img = fits_file[0].data

			img_filter = hdr['FILTER'][0]
			print(img_filter)


			time , mag , mag_err = np.loadtxt ( f , unpack=True)

			time , mag , mag_err = time[np.isfinite(mag) * np.isfinite(mag_err)] , mag[np.isfinite(mag) * np.isfinite(mag_err)] , mag_err[np.isfinite(mag) * np.isfinite(mag_err)]

			a = np.average ( mag , weights=1/mag_err**2 )
			# e = 1 / np.sum ( 1 / mag_err**2 )
			# print(e)
			e = .3
			outliers = np.where ( (mag < a + e) & (mag > a - e)  )
			time , mag , mag_err = time[outliers] , mag[outliers] , mag_err[outliers] , 

			# time , mag , mag_err = bin_lightcurve ( time , mag , mag_err , n=3)

			# print(time)
			# print(mag)
			# print(mag_err)

			param , param_cov = curve_fit (line , time , mag , sigma=mag_err , absolute_sigma=True)
			# print(param , param_cov)

			seconds_from_start = (time - np.min(time)) * 24 * 3600

			fig, ax = plt.subplots(figsize=((paperwidth*1.15) - 2 * margin, (paperheight*1.15) - 2 * margin))
			ax.errorbar ( seconds_from_start , mag , mag_err , fmt='s' , markerfacecolor='blue' , markeredgecolor='black' , ecolor='black' , capthick=2 , markersize=7 , capsize=3   )

			# ax.plot ( time - np.min(time) , line(time, *param) , color='red' )
			ax.invert_yaxis()
			ax.set_xlabel('Seconds from start')
			ax.set_ylabel('Calibrated Magnitude')
			ax.set_xlim([-1 , 61])
			# ax.set_title(f)
			plt.tight_layout()

			corr_mag = mag - line(time, *param) 
			# corr_err = mag_err - line(time, *param)
			corr_err = mag_err

			

			# ax.errorbar ( time - np.min(time) , corr_mag  , corr_err , fmt='s' , markerfacecolor='red' , markeredgecolor='black' , ecolor='black' , capthick=2 , markersize=7 , capsize=3   )
			
			fapLevels = np.array([ 0.32, 0.05, 0.005, 5e-5])
			clp = pyPeriod.Gls(( seconds_from_start , corr_mag , corr_err )  , verbose=False )
			plevel = clp.powerLevel(fapLevels)

			fig_p , ax_p = plt.subplots(figsize=((paperwidth*1.15) - 2 * margin, (paperheight*1.15) - 2 * margin))
			ax_p.plot( 1/clp.freq  , clp.power, color='black' , label='PyAstronomy GLs')
			# ax_p.plot( period , power )
			ax_p.plot([min(1/clp.freq), max(1/clp.freq)], [plevel[0]]*2, '-.', c='black',lw=3,label =r'$1\sigma$')
			# ax_p.plot([min(1/clp.freq), max(1/clp.freq)], [plevel[1]]*2, ':' , c='black',lw=3,label =r'$2\sigma$')
			ax_p.plot([min(1/clp.freq), max(1/clp.freq)], [plevel[2]]*2, ':', c='black'  ,lw=3,label =r'$3\sigma$')
			ax_p.plot([min(1/clp.freq), max(1/clp.freq)], [plevel[3]]*2, '--', c='red'  ,lw=3,label =r'$5\sigma$')

			ax_p.legend()
			ax_p.set_xlabel('Period [s]')
			ax_p.set_xlim([min(1/clp.freq) , 60])
			# ax_p.set_xscale('log')
			ax_p.set_ylim([0 , 1])
			plt.tight_layout()
			# plt.savefig(join(ld,img_name)+'/periodogram_figs/periodogram.png')

			print( f'PyAstronomy GLs best period: {1/clp.freq[np.argmax(clp.power)]}s'  )

			if img_filter == 'r': 
				r_mags.append(mag)
				r_errs.append(mag_err)
				r_time.append(time)
			elif img_filter == 'g':
				g_mags.append(mag)
				g_errs.append(mag_err)
				g_time.append(time)
			elif img_filter == 'i':
				i_mags.append(mag)
				i_errs.append(mag_err)
				i_time.append(time)

			np.savetxt ( join(ld,img_name)+'_flat_lightcurve.txt' , np.vstack([time , corr_mag , corr_err ]).T )
			
			fig_dir = '/'.join(f.split('/')[:-2])+'/figs/'
			if not isdir(fig_dir):
				os.mkdir(fig_dir)

			print(fig_dir + img_name + '_elixir_calibrated_lightcurve.png')
			fig_p.savefig(fig_dir + img_name + '_elixir_calibrated_ls.png')
			fig.savefig(fig_dir + img_name + '_elixir_calibrated_lightcurve.png')

			# np.savetxt ( join(ld,img_name)+'_flat_lightcurve.txt' , np.vstack([time , corr_mag , corr_err ]).T )
			# if True: break

			times.append(time)
			mags.append(corr_mag)
			mags_err.append(corr_err)
			uncor.append(mag)

			print()

	times = np.hstack(times) 
	mags  = np.hstack(mags)
	mags_err  = np.hstack(mags_err)
	uncor = np.hstack(uncor)

	r_mags = np.hstack(r_mags)
	r_errs = np.hstack(r_errs)
	r_time = np.hstack(r_time)
	g_mags = np.hstack(g_mags)
	g_errs = np.hstack(g_errs)
	g_time = np.hstack(g_time)
	i_mags = np.hstack(i_mags)
	i_errs = np.hstack(i_errs)
	i_time = np.hstack(i_time)

	t_0 = np.min(times)

	r_mean = np.average ( r_mags , weights=1/r_errs**2 )
	g_mean = np.average ( g_mags , weights=1/g_errs**2 )
	i_mean = np.average ( i_mags , weights=1/i_errs**2 )
	print(f'g-r={g_mean-r_mean:.2f}')
	print(f'r-i={r_mean-i_mean:.2f}')

	seconds_from_start = (times - t_0)*24*3600

	fig_cal , ax_cal = plt.subplots( figsize=((paperwidth*1.15) - 2 * margin, (paperheight*1.15) - 2 * margin) )
	ax_cal.errorbar( seconds_from_start , uncor, mags_err , fmt='s' , markerfacecolor='blue' , markeredgecolor='black' , ecolor='black' , capthick=2 , markersize=7 , capsize=3   )
	ax_cal.hlines(g_mean, 0, max(seconds_from_start) , colors='red' , linestyle='--' , label=f'g mean={g_mean:.2f}')
	ax_cal.hlines(r_mean, 0, max(seconds_from_start) , colors='red' , linestyle=':' , label=f'r mean={r_mean:.2f}')
	ax_cal.hlines(i_mean, 0, max(seconds_from_start) , colors='red' , linestyle='-.' , label=f'i mean={i_mean:.2f}')
	ax_cal.legend()
	ax_cal.set_xlim(-1,max(seconds_from_start)+1)
	ax_cal.set_xticks(np.arange(0, max(seconds_from_start) , 60))
	ax_cal.set_xlabel('Seconds from exposure start')
	ax_cal.set_ylabel('Calibrated Magnitude')
	plt.gca().invert_yaxis()
	plt.tight_layout()

	fig1 , ax1 = plt.subplots( figsize=((paperwidth*1.15) - 2 * margin, (paperheight*1.15) - 2 * margin) )
	ax1.errorbar( seconds_from_start , mags + r_mean, mags_err , fmt='s' , markerfacecolor='blue' , markeredgecolor='black' , ecolor='black' , capthick=2 , markersize=7 , capsize=3   )
	np.savetxt ( d + 'flat_norm_timeseries.txt' , np.vstack ([ times , mags + r_mean , mags_err]).T , header=' '.join(lc_dirs) )


	fapLevels = np.array([ 0.35, 0.05, 0.005])
	clp = pyPeriod.Gls(( seconds_from_start , mags + r_mean , mags_err )  , verbose=False )
	plevel = clp.powerLevel(fapLevels)

	fig_a , ax_a = plt.subplots(figsize=((paperwidth*1.15) - 2 * margin, (paperheight*1.15) - 2 * margin))
	ax_a.plot( 1/clp.freq  , clp.power, color='black' , label='PyAstronomy GLs')
	# ax_p.plot( period , power )
	ax_a.plot([min(1/clp.freq), max(1/clp.freq)], [plevel[0]]*2, '-.', c='black',lw=3,label =r'$1\sigma$')
	ax_a.plot([min(1/clp.freq), max(1/clp.freq)], [plevel[1]]*2, ':' , c='black',lw=3,label =r'$2\sigma$')
	ax_a.plot([min(1/clp.freq), max(1/clp.freq)], [plevel[2]]*2, '--', c='red'  ,lw=3,label =r'$3\sigma$')

	ax_a.set_title('combined LS')
	ax_a.legend()
	ax_a.set_xlabel('Period [s]')
	ax_a.set_xlim([min(1/clp.freq) , max(seconds_from_start)/2])
	ax_a.set_xscale('log')
	ax_a.set_ylim([0 , 1])
	plt.tight_layout()
	print( f'PyAstronomy GLs best period: {1/clp.freq[np.argmax(clp.power)]}s'  )

	# plt.show()



	clp_r    = pyPeriod.Gls(( (r_time-min(r_time))*24*3600 , r_mags + r_mean , r_errs )  , verbose=False )
	plevel_r = clp.powerLevel(fapLevels)

	fig_r , ax_r = plt.subplots(figsize=((paperwidth*1.15) - 2 * margin, (paperheight*1.15) - 2 * margin))
	ax_r.plot( 1/clp_r.freq  , clp_r.power, color='black' , label='PyAstronomy GLs')
	# ax_p.plot( period , power )
	ax_r.plot([min(1/clp_r.freq), max(1/clp_r.freq)], [plevel_r[0]]*2, '-.', c='black',lw=3,label =r'$1\sigma$')
	ax_r.plot([min(1/clp_r.freq), max(1/clp_r.freq)], [plevel_r[1]]*2, ':' , c='black',lw=3,label =r'$2\sigma$')
	ax_r.plot([min(1/clp_r.freq), max(1/clp_r.freq)], [plevel_r[2]]*2, '--', c='red'  ,lw=3,label =r'$3\sigma$')

	ax_r.set_title('r LS')
	ax_r.legend()
	ax_r.set_xlabel('Period [s]')
	ax_r.set_xlim([min(1/clp_r.freq) , max(1/clp_r.freq)])
	ax_r.set_xscale('log')
	ax_r.set_ylim([0 , 1])
	plt.tight_layout()
	print( f'PyAstronomy GLs best period: {1/clp_r.freq[np.argmax(clp_r.power)]}s for r'  )

	clp_i    = pyPeriod.Gls(( (i_time-min(i_time))*24*3600 , i_mags + r_mean , i_errs )  , verbose=False )
	plevel_i = clp.powerLevel(fapLevels)

	fig_i , ax_i = plt.subplots(figsize=((paperwidth*1.15) - 2 * margin, (paperheight*1.15) - 2 * margin))
	ax_i.plot( 1/clp_i.freq  , clp_i.power, color='black' , label='PyAstronomy GLs')
	# ax_p.plot( period , power )
	ax_i.plot([min(1/clp_i.freq), max(1/clp_i.freq)], [plevel_i[0]]*2, '-.', c='black',lw=3,label =r'$1\sigma$')
	ax_i.plot([min(1/clp_i.freq), max(1/clp_i.freq)], [plevel_i[1]]*2, ':' , c='black',lw=3,label =r'$2\sigma$')
	ax_i.plot([min(1/clp_i.freq), max(1/clp_i.freq)], [plevel_i[2]]*2, '--', c='red'  ,lw=3,label =r'$3\sigma$')

	ax_i.set_title('i LS')
	ax_i.legend()
	ax_i.set_xlabel('Period [s]')
	ax_i.set_xlim([min(1/clp_i.freq) , max(1/clp_i.freq)])
	ax_i.set_xscale('log')
	ax_i.set_ylim([0 , 1])
	plt.tight_layout()
	print( f'PyAstronomy GLs best period: {1/clp_i.freq[np.argmax(clp_i.power)]}s for i'  )


	clp_g    = pyPeriod.Gls(( (g_time-min(g_time))*24*3600 , g_mags + r_mean , g_errs )  , verbose=False )
	plevel_g = clp.powerLevel(fapLevels)

	fig_g , ax_g = plt.subplots(figsize=((paperwidth*1.15) - 2 * margin, (paperheight*1.15) - 2 * margin))
	ax_g.plot( 1/clp_g.freq  , clp_g.power, color='grey' , label='PyAstronomy GLs')
	# ax_p.plot( period , power )
	ax_g.plot([min(1/clp_g.freq), max(1/clp_g.freq)], [plevel_g[0]]*2, '-.', c='black',lw=3,label =r'$1\sigma$')
	ax_g.plot([min(1/clp_g.freq), max(1/clp_g.freq)], [plevel_g[1]]*2, ':' , c='black',lw=3,label =r'$2\sigma$')
	ax_g.plot([min(1/clp_g.freq), max(1/clp_g.freq)], [plevel_g[2]]*2, '--', c='red'  ,lw=3,label =r'$3\sigma$')

	ax_g.set_title('g LS')
	ax_g.legend()
	ax_g.set_xlabel('Period [s]')
	ax_g.set_xlim([min(1/clp_g.freq) , max(1/clp_g.freq)])
	ax_g.set_xscale('log')
	ax_g.set_ylim([0 , 1])
	plt.tight_layout()
	print( f'PyAstronomy GLs best period: {1/clp_g.freq[np.argmax(clp_g.power)]}s for g'  )




	print()
	# plt.savefig(join(ld,img_name)+'/periodogram_figs/periodogram.png')




	plt.show()