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
from magic_star import take_lightcurve, point_rotation
from debugging import display_streak

def reverse_rotation(star_x, star_y, a, img):
	a *= -np.pi/180
	if a>0: 
		m = img.shape[0] * np.abs(np.sin(a))
		star_x_rot =  (star_x -m) * np.cos(a) + star_y * np.sin(a) 
		star_y_rot = -(star_x -m) * np.sin(a) + star_y * np.cos(a)

	elif a<0:
		# a *= -np.pi/180
		m = img.shape[1] * np.abs(np.sin(a))
		star_x_rot =  (star_x) * np.cos(a) + (star_y -m) * np.sin(a)
		star_y_rot = -(star_x) * np.sin(a) + (star_y -m) * np.cos(a)
	return star_x_rot, star_y_rot

def linear_function(x , m , b):
	return x * m + b

directory  = './'
dir_names = [directory+f+'/' for f in os.listdir(directory) if isdir(join(directory,f))] 
mins = {'g':100, 'r': 150, 'i': 250}

for d in dir_names:
	if 'FF14' not in d: continue
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

		img_star_rotated = rotate(img, star_angle)

		fig_lc, ax_lc = plt.subplots(3,5)
		ax_lc[0,0].set_title(f)

		binning = 50

		sum_lc = np.zeros((binning))

		for i in range(10):
			lc = take_lightcurve(img_star_rotated, [centroid_x[i], trail_start_y[i]], [centroid_x[i], trail_end_y[i]], fwhm=star_fwhm[i], binning=binning)[0]
			# print(lc.shape)
			sum_lc += lc
			# print(len(lc))
			# ax_lc[i%3, i%5].scatter(np.arange( len(lc) ), lc)
			#  s, L, a, b, x_0, y_0, width=1
			ax_lc[i%3, i%5].imshow(display_streak(img_star_rotated, star_s[i] , star_length[i] , star_angle , 0 , centroid_x[i] , centroid_y[i] ))


		fig, ax = plt.subplots()
		ax.scatter(np.arange(binning), sum_lc/np.median(sum_lc))

		w = WCS(hdr)
		c = SkyCoord(f'{hdr["CRVAL1"]} {hdr["CRVAL2"]}', unit=(u.deg, u.deg))

		# args_str = f'./refcat {c.ra.deg} {c.dec.deg} -rad 0.5 -dir 00_m_16/'
		args_str = f'./refcat {c.ra.deg} {c.dec.deg} -rect 0.25,0.25 -dir 00_m_16/'
		# print(args_str)

		# RA, Dec, g, r, i, z, J, cyan, orange.
		ref_stars = np.array(os.popen(args_str).read().split('\n')[:-1])
		refcat = []
		for i in ref_stars:
			refcat.append(np.array(i.split(), dtype=float))
		refcat = np.array(refcat)

		refcat_ra_dec      = SkyCoord(ra=refcat[:,0]*u.degree, dec=refcat[:,1]*u.degree, frame='fk5')
		refcat_x, refcat_y = utils.skycoord_to_pixel(refcat_ra_dec, w)

		# constraining to image dimensions
		image_dim     = np.where((refcat_x > 0) & (refcat_x < img.shape[1]) & (refcat_y > 0) & (refcat_y < img.shape[0]) )
		refcat_x      = refcat_x[image_dim] #0, img.shape[1]
		refcat_y      = refcat_y[image_dim]
		refcat_ra_dec = refcat_ra_dec[image_dim]

		# refcat_x      = np.delete(refcat_x , 1)
		# refcat_y      = np.delete(refcat_y , 1)
		# refcat_ra_dec = np.delete(refcat_ra_dec , 1)


		# print(refcat_x, refcat_y)


		fit_ra_dec = utils.pixel_to_skycoord(cen_x_r , cen_y_r , w )

		fig_unr, ax_unr = plt.subplots()
		ax_unr.imshow(img, cmap='gray', norm=colors.LogNorm(vmin=mins[hdr['FILTER'][0]]))
		# ax_unr.scatter(refcat_x, refcat_y, label='refcat')
		ax_unr.set_xlim((0, img.shape[1]))
		ax_unr.set_ylim((img.shape[0], 0))
		

		# print(fit_ra_dec)

		idx, d2d, d3d = refcat_ra_dec.match_to_catalog_sky(fit_ra_dec, nthneighbor=1)

		
		ax_unr.scatter(refcat_x     , refcat_y     , label='refcat')
		ax_unr.scatter(cen_x_r[idx] , cen_y_r[idx] , label='fitted')
		ax_unr.legend()



		# print(refcat)
		# refcat_g = refcat[:,2]
		# refcat_r = refcat[:,3]
		# refcat_i = refcat[:,4]
		ref_mag = []

		img_filter = hdr['FILTER'][0]
		
		if   img_filter == 'g': ref_mag = refcat[:,2]
		elif img_filter == 'r': ref_mag = refcat[:,3]
		elif img_filter == 'i': ref_mag = refcat[:,4]

		ref_mag = ref_mag[image_dim]
		# print(ref_mag.shape)

		# ref_mag = np.delete(ref_mag , 1)

		inst_mag = -2.5 * np.log10(star_flux)


		# print(inst_mag[idx])

		fig_mag, ax_mag = plt.subplots()
		ax_mag.scatter(inst_mag[idx], ref_mag)


		cal_fit, cal_fit_cov = curve_fit ( linear_function , inst_mag[idx] , ref_mag )
		line_label = f'M = {cal_fit[0]}*m + {cal_fit[1]}'
		ax_mag.plot( inst_mag[idx] , linear_function( inst_mag[idx] , *cal_fit ) , label=line_label )
		ax_mag.legend()
		# print(cal_fit)
		print('fit (1 sigma) errors : ' , np.diag(cal_fit_cov) **.5)
		# print(f'{f[:-11]}_zeropoint.txt')

		np.savetxt(f'{f[:-11]}_zeropoint.txt'    , np.vstack([cal_fit , np.diag(cal_fit_cov) **.5 ]) )




		# if True: break



	
	plt.show()