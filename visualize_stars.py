import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import astropy as ap
from astropy.timeseries import LombScargle
from astropy.io import fits
from os.path import isdir, isfile, join
from scipy.ndimage import rotate
from magic_star import take_lightcurve, point_rotation

directory  = './'
dir_names = [directory+f+'/' for f in os.listdir(directory) if isdir(join(directory,f))] 
mins = {'g':100, 'r': 150, 'i': 250}


for d in dir_names:
	if 'GE1' not in d: continue
	file_names = [d+f for f in os.listdir(d) if isfile(join(d,f))]

	for f in file_names:
	# if ('lightcurve' not in f) and ('.txt' in f) and ('params' not in f) :
		if ('params.txt' not in f): continue

		star_params = np.loadtxt(f)
		print(len(star_params))
		# LOADING STAR PARAMS
		star_angle = star_params[0, -2]

		star_length = star_params[:,1]
		star_fwhm   = star_params[:,0] * 2.355
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
		# for i in range(10):
		# 	ax_lc[i%2, i%5].

		sum_lc = np.zeros((100))

		for i in range(15):
			lc = take_lightcurve(img_star_rotated, [centroid_x[i], trail_start_y[i]], [centroid_x[i], trail_end_y[i]], fwhm=star_fwhm[i], binning=100)[0]
			# print(lc.shape)
			sum_lc += lc
			# print(len(lc))
			ax_lc[i%3, i%5].scatter(np.arange( len(lc) ), lc)


		fig_im, ax_im = plt.subplots()
		ax_im.set_title(f)
		ax_im.imshow(img_star_rotated, cmap='gray', norm=colors.LogNorm(vmin=mins[hdr['FILTER'][0]]))
		ax_im.scatter(centroid_x, centroid_y)

		fig, ax = plt.subplots()
		ax.scatter(np.arange(100), sum_lc/np.median(sum_lc))


		if True: break



	
	plt.show()