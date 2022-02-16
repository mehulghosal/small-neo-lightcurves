import warnings, subprocess, sep
import numpy as np
import astropy as ap
import matplotlib.pyplot as plt
from matplotlib import colors
from astropy.io import fits
from scipy.ndimage import rotate
from scipy.special import erf
from astropy.wcs import WCS
from astropy.wcs import utils
from astropy.coordinates import SkyCoord
from astropy import units as u

from astropy.utils.exceptions import AstropyWarning

plt.rcParams.update({'figure.max_open_warning': 0})
warnings.simplefilter('ignore', AstropyWarning)

# initializing all directories
import os
from os.path import isdir, isfile, join
directory = './'	
dir_names = [directory+f+'/' for f in os.listdir(directory) if isdir(join(directory,f))] 

input_file = np.loadtxt('input.csv', dtype=object, skiprows=1, usecols=(i for i in range(25)), delimiter=',')

mins = {'g':100, 'r': 150, 'i': 250}

from scipy.optimize import curve_fit

# rotate points by angle a [degrees]
def point_rotation(x,y,a,img,img_rot):
	a = -a * np.pi/180
	x_0, y_0 = 0, 0
	x_0_, y_0_ = img.shape[0]*np.abs(np.sin(a)), img.shape[1]*np.abs(np.sin(a))
	x_, y_ = int((x-x_0)*np.cos(a) - (y-y_0)*np.sin(a)), int((x-x_0)*np.sin(a) + (y-y_0)*np.cos(a))
	# to account for direction of rotation
	if a>0: x_+= int(x_0_)
	elif a<0: y_+= int(y_0_)

	return x_, y_

def model(x, s, m, a, c, b, d):
	return c*np.exp(-.5* ((x-m)/s)**2) + a*x + b + d*x**2

def quadratic(x, a, b, c, d, e):
	return a*x**2 + b*x + c + d*x**3 + e*x**4

img_rot, centroid = 0, 0
count = 0
# Veres 2012 eq 3
# r = [x,y], s = sigma, L is length, a is angle, b is background noise (holding constant for now)
# img_rot and centroid are not fitting variables - want to pass these in as constants; centroid = [x,y]
def trail_model(r, s, L, a, b):
	global img_rot, centroid, count
	img = rotate(img_rot, a)
	centroid_x, centroid_y = point_rotation(centroid[0], centroid[1], a, img_rot, img)
	print(r, count)
	x, y = point_rotation(r[0], r[1], a, img_rot, img)
	# x, y = r[0], r[1]
	x-=centroid_x
	y-=centroid_y

	trail = img[centroid_x-s:centroid_x+s, centroid_y-L, centroid_y+L]

	flux = np.sum(trail)


	return flux/(L/2 * s * (8 * np.pi)**.5) * np.exp(-((x+y)**2 )/(2*s**2)) * erf((x+y + L/4)/ (s*2**.5)) - erf((x+y - L/4)/ (s*2**.5)) + b






for d in dir_names:
	file_names = [d+f for f in os.listdir(d) if isfile(join(d,f))]

	for f in file_names:
		try:
			file = fits.open(f)
		except Exception as e:
			print(f)
			continue
		hdr = file[0].header
		img = file[0].data


		# object id from directory name --> string splicing
		obj_id = f.split('_')
		obj_id = obj_id[0][2:] + ' ' + obj_id[1]

		if not ('2016 GE1' in obj_id and '70o13' in f): continue
		# if '2016 CD3' not in obj_id: continue

		# plt.figure()
		fig, ax = plt.subplots(1,3)
		ax[0].set_title(f)
		# ax[0].imshow(img, cmap='gray', norm=colors.LogNorm(vmin=mins[hdr['FILTER'][0]]))

		obj_rows = input_file[np.where(input_file[:,1]==obj_id),:][0]
		
		
		try:
			obj = obj_rows[np.where(obj_rows[:,0]==f.split('/')[-1])][0]
			trail_start = np.array(obj[-4:-2], dtype=int)
			trail_end	= np.array(obj[-2:], dtype=int)
			# trail_start = np.array([1196, 3980])
			# trail_end = np.array([1303, 4175])
		except Exception as e:
			print(f,obj[-4:-2],obj[-2:])
			plt.close()
			continue

		
		angle = -1*np.arctan2(trail_end[0]-trail_start[0], trail_end[1]-trail_start[1]) * 180/np.pi
		img_rotated = rotate(img, angle)
		# ax[0].imshow(img_rotated, cmap='gray', norm=colors.LogNorm(vmin=mins[hdr['FILTER'][0]]))

		trail_start = np.array(point_rotation(trail_start[0], trail_start[1], angle, img, img_rotated), dtype=int)
		trail_end	= np.array(point_rotation(trail_end[0]  , trail_end[1]  , angle, img, img_rotated), dtype=int)
		trail_length = trail_end[1] - trail_start[1]

		# assuming vertical streaks for drawing rectangles and moving down 
		obj_width = 25
		
		obj_rect = img_rotated[trail_start[1]:trail_end[1], trail_start[0]-obj_width:trail_start[0]+obj_width]
		# ax[0].imshow(obj_rect, cmap='gray', norm=colors.LogNorm(vmin=300))

		col_sums = np.sum(obj_rect, axis=0)
		print(col_sums.shape)
		# col_sums /= np.max(col_sums)
		rect_width = np.arange(0, 2*obj_width, 1)
		param_vals, param_covs = curve_fit(model, rect_width, col_sums, p0=[3, obj_width, .03, 60000, 20000, -3])

		fwhm = int(param_vals[0] * 2.355)
		# fwhm = 6
		print(param_vals)
		print(np.diag(param_covs))

		ax[2].scatter(rect_width, col_sums, label='column sums')
		ax[2].plot(rect_width, model(rect_width, *param_vals), label='model fit')
		ax[2].legend()

		centroid_deviation = -obj_width + param_vals[1] # if negative, trail is to the left, if positive, trail to right
		height_correction = int((trail_end[1] - trail_start[1]) * .2 + .5) # 20% more rows above and below to get some sky 

		# correcting trail start/end
		trail_start[0] += int(centroid_deviation+.5)
		trail_end[0]   += int(centroid_deviation+.5)
		trail_start[1] -= height_correction
		trail_end[1]   += height_correction

		print(trail_end-trail_start)
		# asteroid trail length in 70o13 is 101 tall
		ax[0].plot([trail_start[0], trail_end[0]], [trail_start[1], trail_end[1]], marker='*')


		obj_width = 1*fwhm
		sky_width = 2*fwhm
		obj_rect = img_rotated[trail_start[1]:trail_end[1], trail_start[0]-obj_width:trail_start[0]+obj_width]


		sky_left  = img_rotated[trail_start[1]:trail_end[1], trail_start[0]-obj_width-sky_width:trail_start[0]-obj_width]
		sky_right = img_rotated[trail_start[1]:trail_end[1], trail_start[0]+obj_width:trail_start[0]+obj_width+sky_width]

		obj_row_sums = np.array([np.sum(i) for i in obj_rect])
		sky_left_row_sum  = np.array([np.sum(i) for i in sky_left ])
		sky_right_row_sum = np.array([np.sum(i) for i in sky_right])
		sky_row_avg = (sky_right_row_sum+sky_left_row_sum)/(sky_right.shape[1]+sky_left.shape[1])

		obj_minus_sky = obj_row_sums - sky_row_avg * obj_rect.shape[1]

		ax[0].imshow(img_rotated, cmap='gray', norm=colors.LogNorm(vmin=np.median(sky_row_avg)))
		# ax[0].imshow(obj_rect, cmap='gray', norm=colors.LogNorm(vmin=np.median(sky_row_avg)))

		sigma_row = obj_minus_sky + (len(obj_row_sums)) * (sky_row_avg + hdr['RDNOISE']**2) + (len(obj_row_sums))**2 * sky_row_avg**.5
		sigma_row = sigma_row ** .5

		# x = np.arange(0, 101, 101/len(obj_row_sums))
		x = np.arange(0, len(obj_row_sums), 1)
		# ax[1].plot(x, obj_minus_sky)

		param_vals, param_covs = curve_fit(quadratic, x, obj_minus_sky, sigma=sigma_row)
		ax[1].plot(x, quadratic(x, *param_vals))
		# ax[1].set_ylim([np.min(obj_minus_sky),25000])
		print(param_vals)
		print(np.diag(param_covs))
		ax[1].errorbar(x, obj_minus_sky, yerr = sigma_row, fmt='r', capsize=3, linewidth=2, elinewidth=1)

		# WCS stuff
		w = WCS(hdr)
		c = SkyCoord(f'{obj[7]} {obj[8]}', unit=(u.deg, u.deg))
		target_x, target_y = np.round(utils.skycoord_to_pixel(c, w))
		
		# ax[0].plot(target_x, target_y, 'b+')
		
		# https://sep.readthedocs.io/en/v0.4.x/api/sep.extract.html
		'''
		background = sep.Background(img)
		back = background.back()
		background.subfrom(img)
		print(np.median(background))
		# sep.set_sub_object_limit(100)
		sextractor = sep.extract(img, 3, err=back)
		print(len(sextractor))
		for s in sextractor:
			x, y = point_rotation(s[7], s[8], angle, img, img_rotated)
			# x, y = s[7], s[8]
			ax[0].plot(x, y, 'b+')
		'''

		# source extractor test
		sex = subprocess.run(['sex', f, '-DETECT_MINAREA', str(trail_length*fwhm)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
		sex_output = np.loadtxt('test.cat', skiprows=9)
		print(sex_output.shape)
		star_x = sex_output[:,5]
		star_y = sex_output[:,6]

		star_x_min = sex_output[:,1]
		star_y_min = sex_output[:,2]
		star_x_max = sex_output[:,3]
		star_y_max = sex_output[:,4]
		for i in range(len(star_x)):
			star_x[i], star_y[i] = point_rotation(star_x[i], star_y[i], angle, img, img_rotated)
			star_x_min[i], star_y_min[i] = point_rotation(star_x_min[i], star_y_min[i], angle, img, img_rotated)
			star_x_max[i], star_y_max[i] = point_rotation(star_x_max[i], star_y_max[i], angle, img, img_rotated)
			# print(star_x[i], star_y[i])

		bad_stars = np.where((star_x<trail_length) | (star_x>img_rotated.shape[1]-trail_length) | (star_y<trail_length) | (star_y>img_rotated.shape[0]-trail_length))
		print(bad_stars)
		star_x = np.delete(star_x, bad_stars, 0)
		star_y = np.delete(star_y, bad_stars, 0)
		star_x_min = np.delete(star_x_min, bad_stars, 0)
		star_y_min = np.delete(star_y_min, bad_stars, 0)
		star_x_max = np.delete(star_x_max, bad_stars, 0)
		star_y_max = np.delete(star_y_max, bad_stars, 0)
		# global centroid

		coords = np.indices(img_rotated.shape).T.reshape((img_rotated.shape[0], img_rotated.shape[1], 2)).flatten().reshape((img_rotated.shape[0] * img_rotated.shape[1], 2))
		coords = list(zip(coords[:,0], coords[:,1]))
		img_rot = img_rotated
		flattened_img = img_rotated.flatten()
		print(flattened_img.shape)
		print(coords[0], coords[-1])
		print()
		for i in range(len(star_x)):
			centroid = star_x[i], star_y[i]
		
			p, p_cov = curve_fit(trail_model, coords, flattened_img, p0=[3, trail_length, 0, np.mean(sky_row_avg)])
			print(p)
			print(cov)




		ax[0].scatter(star_x, star_y, c='orange', s=2, label='centroid')
		ax[0].scatter(star_x_min, star_y_min, c='green', s=2, label='mins')
		ax[0].scatter(star_x_max, star_y_max, c='purple', s=2, label='maxes')
		ax[0].legend()


		# ax[0].legend()
	
	
		
	plt.show()
	# if True: break
# output.close()