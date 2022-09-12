import os, subprocess
from os.path import isdir, isfile, join
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import astropy as ap
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS, utils
from astropy.timeseries import LombScargle
from astropy.coordinates import SkyCoord
from scipy.ndimage import rotate
from scipy.optimize import curve_fit
from magic_star import take_lightcurve, point_rotation, reverse_rotation
from debugging import display_streak


def linear_function(x , m , b):
	return x * m + b

def line_slope_one(x , b):
	return x + b

directory  = './'
dir_names = [directory+f+'/' for f in os.listdir(directory) if isdir(join(directory,f))] 
mins = {'g':100, 'r': 150, 'i': 250}

for d in dir_names:
	if 'GE1' not in d: continue
	file_names = [d+f for f in os.listdir(d) if isfile(join(d,f))]

	for f in file_names:
	# if ('lightcurve' not in f) and ('.txt' in f) and ('params' not in f) :
		if ('_params.txt' not in f) : continue
		# if ('71o13' in f) : continue
		# if : continue

		star_params = np.loadtxt(f)
		print(len(star_params))
		# LOADING STAR PARAMS
		star_angle = star_params[0, -2]
		print(star_angle)

		star_length = star_params[:,1]
		star_s      = star_params[:,0]
		star_fwhm   = star_s * 2.355
		star_flux   = star_params[:,-1]

		centroid_x = star_params[:,-4]
		centroid_y = star_params[:,-3]

		trail_start_y = centroid_y - star_length/2
		trail_end_y   = centroid_y + star_length/2

		# GETTING IMAGE 
		print(f, f[:-11]+".flt")

		fits_img = fits.open(f[:-11]+".flt")
		hdr 	 = fits_img[0].header
		img 	 = fits_img[0].data

		cen_x_r , cen_y_r  = [] , []

		for i in range(len(centroid_y)):
			c_x = centroid_x[i]
			c_y = centroid_y[i]
			c_x_r , c_y_r = reverse_rotation(c_x, c_y, star_angle, img)
			cen_x_r.append(c_x_r)
			cen_y_r.append(c_y_r)
		cen_x_r = np.array(cen_x_r)
		cen_y_r = np.array(cen_y_r)


		inst_mag = -2.5 * np.log10(star_flux)
		mag_filter = np.where(inst_mag <= 0)
		inst_mag = inst_mag[mag_filter]



		cen_y_r = cen_y_r [ mag_filter ]
		cen_x_r = cen_x_r [ mag_filter ]


		w = WCS(hdr)
		c = SkyCoord(f'{hdr["CRVAL1"]} {hdr["CRVAL2"]}', unit=(u.deg, u.deg))
		fit_ra_dec = utils.pixel_to_skycoord(cen_x_r , cen_y_r , w )

		img_star_rotated = rotate(img, star_angle)

		fig_lc, ax_lc = plt.subplots(3,5)
		ax_lc[0,0].set_title(f)

		binning = 50

		sum_lc = np.zeros((binning))

		for i in range(9):
			lc = take_lightcurve(img_star_rotated, [centroid_x[i], trail_start_y[i]], [centroid_x[i], trail_end_y[i]], fwhm=star_fwhm[i], binning=binning)[0]
			# print(lc.shape)
			sum_lc += lc
			# print(len(lc))
			# ax_lc[i%3, i%5].scatter(np.arange( len(lc) ), lc)
			ax_lc[i%3, i%5].imshow(display_streak(img_star_rotated, star_s[i] , star_length[i] , star_angle , 0 , centroid_x[i] , centroid_y[i] ))


		fig, ax = plt.subplots()
		ax.scatter(np.arange(binning), sum_lc/np.median(sum_lc))

		
		# args_str = f'./refcat {c.ra.deg} {c.dec.deg} -rad 0.5 -dir 00_m_16/'
		args_str = f'./refcat {c.ra.deg} {c.dec.deg} -rect 0.25,0.25 -dir 00_m_16/'

		# RA, Dec, g, r, i, z, J, cyan, orange.
		ref_stars = np.array(os.popen(args_str).read().split('\n')[:-1])
		refcat = []
		for i in ref_stars:
			refcat.append(np.array(i.split(), dtype=float))
		refcat = np.array(refcat)

		ref_mag = []

		img_filter = hdr['FILTER'][0]
		
		if   img_filter == 'g': ref_mag = refcat[:,2]
		elif img_filter == 'r': ref_mag = refcat[:,3]
		elif img_filter == 'i': ref_mag = refcat[:,4]

		refcat_ra_dec      = SkyCoord(ra=refcat[:,0]*u.degree, dec=refcat[:,1]*u.degree, frame='fk5')
		refcat_x, refcat_y = utils.skycoord_to_pixel(refcat_ra_dec, w)

		# constraining to image dimensions
		image_dim     = np.where((refcat_x > 0) & (refcat_x < img.shape[1]) & (refcat_y > 0) & (refcat_y < img.shape[0]) )
		refcat_x      = refcat_x[image_dim] #0, img.shape[1]
		refcat_y      = refcat_y[image_dim]
		refcat_ra_dec = refcat_ra_dec[image_dim]
		ref_mag       = ref_mag[image_dim]


		fig_unr, ax_unr = plt.subplots()
		ax_unr.imshow(img, cmap='gray', norm=colors.LogNorm(vmin=mins[hdr['FILTER'][0]]))
		# ax_unr.scatter(refcat_x, refcat_y, label='refcat')
		ax_unr.set_xlim((0, img.shape[1]))
		ax_unr.set_ylim((img.shape[0], 0))
		
		# idx, d2d, d3d = fit_ra_dec.match_to_catalog_sky(refcat_ra_dec, nthneighbor=1)
		idx, d2d, d3d = refcat_ra_dec.match_to_catalog_sky(fit_ra_dec, nthneighbor=1)
		print(d2d.arcsec)

		dist_filter = np.where(d2d.arcsec < 75)
		# idx = idx[dist_filter]
		# fit_ra_dec = fit_ra_dec[dist_filter]
		
		ax_unr.scatter(refcat_x     , refcat_y     , label='refcat')
		ax_unr.scatter(cen_x_r[idx[dist_filter]] , cen_y_r[idx[dist_filter]] , label='fitted')
		ax_unr.legend()


		fig_mag, ax_mag = plt.subplots()
		# ax_mag.scatter(inst_mag[dist_filter], ref_mag[idx[dist_filter]])
		ax_mag.scatter(inst_mag[idx], ref_mag)

		# print(ref_mag[idx[dist_filter]] )

		cal_fit, cal_fit_cov = curve_fit ( line_slope_one , inst_mag[idx[dist_filter]] , ref_mag[dist_filter] )
		# line_label = f'M = {cal_fit[0]}*m + {cal_fit[1]}'
		# line_label = f'M = m + {cal_fit[0]}'
		print(cal_fit) 
		ax_mag.plot( inst_mag[idx[dist_filter]] , line_slope_one( inst_mag[idx[dist_filter]] , *cal_fit )  )
		# ax_mag.legend()
		# print(cal_fit)
		# print('fit (1 sigma) errors : ' , np.diag(cal_fit_cov) **.5)
		ax_mag.set_title(f'ZP = {cal_fit} +/- {np.diag(cal_fit_cov) **.5}')
		# print(f'{f[:-11]}_zeropoint.txt')

		# if True:
		# 	np.savetxt(f'{f[:-11]}_zeropoint.txt' , np.vstack([cal_fit , np.diag(cal_fit_cov) **.5 ]) )

		if True: break



	
	plt.show()