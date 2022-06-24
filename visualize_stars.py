import os, subprocess
from os.path import isdir, isfile, join
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import astropy as ap
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs import utils
from astropy.timeseries import LombScargle
from astropy.coordinates import SkyCoord
from scipy.ndimage import rotate
from magic_star import take_lightcurve, point_rotation
from debugging import display_streak

def reverse_rotation():
	pass

directory  = './'
dir_names = [directory+f+'/' for f in os.listdir(directory) if isdir(join(directory,f))] 
mins = {'g':100, 'r': 150, 'i': 250}

for d in dir_names:
	if 'FF14' not in d: continue
	file_names = [d+f for f in os.listdir(d) if isfile(join(d,f))]

	for f in file_names:
	# if ('lightcurve' not in f) and ('.txt' in f) and ('params' not in f) :
		if ('params.txt' not in f): continue
		# if '66o13' not in f: continue

		star_params = np.loadtxt(f)
		print(len(star_params))
		# LOADING STAR PARAMS
		star_angle = star_params[0, -2]

		star_length = star_params[:,1]
		star_s      = star_params[:,0]
		star_fwhm   = star_s * 2.355
		star_flux   = star_params[:,-1]

		centroid_x = star_params[:,4]
		centroid_y = star_params[:,5]

		trail_start_y = centroid_y - star_length/2
		trail_end_y   = centroid_y + star_length/2

		# GETTING IMAGE 
		print(f, f[:-11]+".flt")

		fits_img = fits.open(f[:-11]+".flt")
		hdr 	 = fits_img[0].header
		img 	 = fits_img[0].data

		img_star_rotated = rotate(img, star_angle)

		fig_lc, ax_lc = plt.subplots(3,5)
		ax_lc[0,0].set_title(f)

		binning = 50

		sum_lc = np.zeros((binning))

		for i in range(15):
			lc = take_lightcurve(img_star_rotated, [centroid_x[i], trail_start_y[i]], [centroid_x[i], trail_end_y[i]], fwhm=star_fwhm[i], binning=binning)[0]
			# print(lc.shape)
			sum_lc += lc
			# print(len(lc))
			# ax_lc[i%3, i%5].scatter(np.arange( len(lc) ), lc)
			#  s, L, a, b, x_0, y_0, width=1
			ax_lc[i%3, i%5].imshow(display_streak(img_star_rotated, star_s[i] , star_length[i] , star_angle , 0 , centroid_x[i] , centroid_y[i] ))


		fig, ax = plt.subplots()
		ax.scatter(np.arange(binning), sum_lc/np.median(sum_lc))

		fig_im, ax_im = plt.subplots()
		ax_im.set_title(f)
		ax_im.imshow(img_star_rotated, cmap='gray', norm=colors.LogNorm(vmin=mins[hdr['FILTER'][0]]))
		ax_im.scatter(centroid_x, centroid_y)

		w = WCS(hdr)
		c = SkyCoord(f'{hdr["CRVAL1"]} {hdr["CRVAL2"]}', unit=(u.deg, u.deg))

		# args_str = f'./refcat {c.ra.deg} {c.dec.deg} -rad 0.5 -dir 00_m_16/'
		args_str = f'./refcat {c.ra.deg} {c.dec.deg} -rect 0.25,0.25 -dir 00_m_16/'
		print(args_str)

		# RA, Dec, g, r, i, z, J, cyan, orange.
		ref_stars = np.array(os.popen(args_str).read().split('\n')[:-1])
		refcat = []
		for i in ref_stars:
			refcat.append(np.array(i.split(), dtype=float))
		refcat = np.array(refcat)

		refcat_ra_dec      = SkyCoord(ra=refcat[:,0]*u.degree, dec=refcat[:,1]*u.degree, frame='fk5')
		refcat_x, refcat_y = utils.skycoord_to_pixel(refcat_ra_dec, w)
		# print(refcat_x, refcat_y)

		refcat_x_, refcat_y_ = [], []
		for i in range(len(refcat_x)):
			rot_x , rot_y = point_rotation(refcat_x[i], refcat_y[i], star_angle, img, img_star_rotated)
			refcat_x_.append(rot_x)
			refcat_y_.append(rot_y)

		ax_im.scatter(refcat_x_, refcat_y_)
		ax_im.set_xlim((0, img_star_rotated.shape[1]))
		ax_im.set_ylim((img_star_rotated.shape[0], 0))


		fig_unr, ax_unr = plt.subplots()
		ax_unr.imshow(img, cmap='gray', norm=colors.LogNorm(vmin=mins[hdr['FILTER'][0]]))
		ax_unr.scatter(refcat_x, refcat_y)
		ax_unr.set_xlim((0, img.shape[1]))
		ax_unr.set_ylim((img.shape[0], 0))

		# fitted_cat = SkyCoord()

		# idx, d2d, d3d = our_catalog.match_to_catalog_sky(refcat_ra_dec, nthneighbor=1)

		if True: break



	
	plt.show()