import warnings, subprocess
import numpy as np
import astropy as ap
import matplotlib.pyplot as plt
from matplotlib import colors
from astropy.io import fits
from scipy.ndimage import rotate
from scipy.optimize import curve_fit
from scipy.optimize import OptimizeWarning
from scipy.stats import norm
# from astropy.io.ascii import sextractor
from astropy.wcs import WCS
from astropy.wcs import utils
from astropy.coordinates import SkyCoord
from astropy import units as u

from astropy.utils.exceptions import AstropyWarning


plt.rcParams.update({'figure.max_open_warning': 0})
warnings.simplefilter('ignore', AstropyWarning)
warnings.simplefilter('ignore', OptimizeWarning)


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


# def star_box_model(x, t_1, t_2, a,b,c,d,e,f,g, a_1,b_1,c_1,d_1,e_1,f_1,g_1):
# def star_box_model(x, t_1, t_2, a,b,c,d,e,f,g,back):
def star_box_model(x, t_1, t_2, a, w, d, c, a_1, w_1, d_1, c_1):
	# r[np.where((x <= t_2) & (x >=t_1))] = a*x**6+b*x**5+c*x**4+d*x**3+e*x**2+f*x+g
	# r[np.where((x >= t_2) | (x <=t_1))] = a_1*x**6+b_1*x**5+c_1*x**4+d_1*x**3+e_1*x**2+f_1*x+g_1
	# return np.piecewise(x, [(x <= t_2) & (x >=t_1), (x >= t_2) | (x <=t_1)], [lambda x: a*x**6+b*x**5+c*x**4+d*x**3+e*x**2+f*x+g, lambda x: a_1*x**6+b_1*x**5+c_1*x**4+d_1*x**3+e_1*x**2+f_1*x+g_1])
	return np.piecewise(x, [(x <= t_2) & (x >=t_1), (x >= t_2) | (x <=t_1)], [lambda x: a*np.sin(w*x+d) + c, lambda x: a_1*np.sin(w_1*x+d_1)+c_1])

img_rot, centroid, = 0, 0
count = 0


# initializing all directories
import os
from os.path import isdir, isfile, join
directory = './'	
dir_names = [directory+f+'/' for f in os.listdir(directory) if isdir(join(directory,f))] 

input_file = np.loadtxt('input.csv', dtype=object, skiprows=1, usecols=(i for i in range(25)), delimiter=',')

mins = {'g':100, 'r': 150, 'i': 250}

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
		# print(np.median(img))
		

		# object id from directory name --> string splicing
		obj_id = f.split('_')
		obj_id = obj_id[0][2:] + ' ' + obj_id[1]

		if 'EN156' not in obj_id: continue
		# if '2016 CD31' not in obj_id: continue
		# if '2016 EN156' not in obj_id: continue

		# plt.figure()
		fig, ax = plt.subplots(1,3)
		ax[0].set_title(f)
		# plt.imshow(img, cmap='gray', norm=colors.LogNorm(vmin=mins[hdr['FILTER'][0]]))

		obj_rows = input_file[np.where(input_file[:,1]==obj_id),:][0]

		try:
			obj = obj_rows[np.where(obj_rows[:,0]==f.split('/')[-1])][0]
			trail_start = np.array(obj[-4:-2], dtype=int)
			trail_end	= np.array(obj[-2:], dtype=int)
		except Exception as e:
			print(f,obj[-4:-2],obj[-2:])
			plt.close()
			continue
		
		angle = -1*np.arctan2(trail_end[0]-trail_start[0], trail_end[1]-trail_start[1]) * 180/np.pi
		# absolutely love commenting out approximations :)
		# if np.abs(angle)<2: angle=0
		img_rotated = rotate(img, angle, reshape=True, axes=(0,1))
		
		# ax[0].imshow(img_rotated, cmap='gray', norm=colors.LogNorm(vmin=mins[hdr['FILTER'][0]]))
		
		print(f)
		print(angle)
		print(trail_start, trail_end)

		trail_start = np.array(point_rotation(trail_start[0], trail_start[1], angle, img, img_rotated), dtype=int)
		trail_end	= np.array(point_rotation(trail_end[0]  , trail_end[1]  , angle, img, img_rotated), dtype=int)
		trail_length = trail_end[1] - trail_start[1]

		# assuming vertical streaks for drawing rectangles and moving down 
		obj_width = 25
		right_width = obj_width
		left_width = obj_width
		
		obj_rect = img_rotated[trail_start[1]:trail_end[1], trail_start[0]-left_width:trail_start[0]+right_width]

		# to fit gaussians across sum of columns  -> mega signal, kinda what you'd expect of a point spread function, but column spread function
		col_sums = np.sum(obj_rect, axis=0)
		# col_sums /= np.max(col_sums)
		rect_width = np.arange(0, left_width+right_width, 1)
		param_vals, param_covs = curve_fit(model, rect_width, col_sums, p0=[3, obj_width, .03, 60000, 20000, -3])

		fwhm = int(param_vals[0] * 2.355)
		print(param_vals)
		print(np.diag(param_covs))

		centroid_deviation = -obj_width + param_vals[1] # if negative, trail is to the left, if positive, trail to right
		height_correction = int((trail_end[1] - trail_start[1]) * .2 + .5) # 20% more rows above and below to get some sky 

		# correcting trail start/end
		trail_start[0] += int(centroid_deviation+.5)
		trail_end[0]   += int(centroid_deviation+.5)
		trail_start[1] -= height_correction
		trail_end[1]   += height_correction

		trail_centroid = np.array([trail_start[0], np.mean([trail_start[1], trail_end[1]])])

		ax[0].plot([trail_start[0], trail_end[0]], [trail_start[1], trail_end[1]], marker='*', label='asteroid')
		print(centroid_deviation)

		# redefining obj rect to be centered based on corrected trail start/end
		# obj_rect = img_rotated[trail_start[1]:trail_end[1], trail_start[0]-3*fwhm:trail_start[0]+3*fwhm]
		# col_sums = np.sum(obj_rect, axis=0)
		# rect_width = np.arange(0, 6*fwhm, 1)
		# param_vals, param_covs = curve_fit(model, rect_width, col_sums, p0=[3, obj_width, .03, 60000, 20000, -3])

		ax[2].scatter(rect_width, col_sums, label='column sums')
		ax[2].plot(rect_width, model(rect_width, *param_vals), label='model fit')
		ax[2].legend()

		# redefining obj rectangle to be centered on corrected trail w/ width of fwhm
		right_width = int(fwhm)
		left_width  = int(fwhm)
		obj_rect = img_rotated[trail_start[1]-height_correction:trail_end[1]+height_correction, trail_start[0]-left_width:trail_start[0]+right_width]

		# total sky width of 4*fwhm
		sky_width = 5*right_width

		sky_left  = img_rotated[trail_start[1]-height_correction:trail_end[1]+height_correction, trail_start[0]-left_width-sky_width:trail_start[0]-left_width]
		sky_right = img_rotated[trail_start[1]-height_correction:trail_end[1]+height_correction, trail_start[0]+right_width:trail_start[0]+right_width+sky_width]

		obj_row_sums = np.array([np.sum(i) for i in obj_rect])
		sky_left_row_sum  = np.array([np.sum(i) for i in sky_left ])
		sky_right_row_sum = np.array([np.sum(i) for i in sky_right])
		sky_row_avg = (sky_right_row_sum+sky_left_row_sum)/(sky_right.shape[1]+sky_left.shape[1])

		obj_minus_sky = obj_row_sums - sky_row_avg * obj_rect.shape[1]

		sigma_row = obj_minus_sky + (len(obj_row_sums)) * (sky_row_avg + hdr['RDNOISE']**2) + (len(obj_row_sums))**2 * sky_row_avg**.5
		sigma_row = sigma_row ** .5

		x = np.arange(0, len(obj_row_sums), 1)
		ax[1].errorbar(x, obj_minus_sky, yerr = sigma_row, fmt='r', capsize=3, linewidth=2, elinewidth=1)

		# WCS stuff
		w = WCS(hdr)
		c = SkyCoord(f'{obj[7]} {obj[8]}', unit=(u.deg, u.deg))
		target_x, target_y = np.round(utils.skycoord_to_pixel(c, w))
		target_x, target_y = point_rotation(trail_start[0], trail_start[1], angle, img, img_rotated)


		# star stuff - SExtractor for detecting star trails
		sex = subprocess.run(['sex', f, '-DETECT_MINAREA', str(trail_length*fwhm)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
		# sex = subprocess.run(['sex', f, '-DETECT_MINAREA', str(trail_length*fwhm)])
		sex_output = np.loadtxt('test.cat', skiprows=9)
		print(sex_output.shape)
		star_x = sex_output[:,5]
		star_y = sex_output[:,6]

		star_x_min = sex_output[:,1]
		star_y_min = sex_output[:,2]
		star_x_max = sex_output[:,3]
		star_y_max = sex_output[:,4]

		dist_to_asteroid = []

		for i in range(len(star_x)):
			star_x[i], star_y[i] = point_rotation(star_x[i], star_y[i], angle, img, img_rotated)
			dist_to_asteroid.append((star_x[i] - trail_centroid[0])**2 + (star_y[i] - trail_centroid[1])**2)
			star_x_min[i], star_y_min[i] = point_rotation(star_x_min[i], star_y_min[i], angle, img, img_rotated)
			star_x_max[i], star_y_max[i] = point_rotation(star_x_max[i], star_y_max[i], angle, img, img_rotated)

		# num = len(dist_to_asteroid)
		num = 25
		dist_to_asteroid = np.array(dist_to_asteroid)
		dist_sorted = np.argsort(dist_to_asteroid)
		star_x = star_x[dist_sorted][:num]
		star_y = star_y[dist_sorted][:num]
		star_x_min = star_x_min[dist_sorted][:num]
		star_y_min = star_y_min[dist_sorted][:num]
		star_x_max = star_x_max[dist_sorted][:num]
		star_y_max = star_y_max[dist_sorted][:num]
		# filtering bad stars out: too close to edge to be useful; asteroid trail; getting only 15 nearest stars
		bad_stars = np.where((star_x<trail_length) | (star_x>img_rotated.shape[1]-trail_length) | (star_y<trail_length) | (star_y>img_rotated.shape[0]-trail_length))
		bad_stars = np.append(bad_stars, np.where((star_x<trail_start[0]+fwhm) & (star_x>trail_start[0]-fwhm) & (star_y<trail_end[1]) & (star_y>trail_start[1]))) # want to get rid of asteroid too
		bad_stars = np.append(bad_stars, np.where((star_x<trail_start[0]+fwhm) & (star_x>trail_start[0]-fwhm) & (star_y<trail_end[1]) & (star_y>trail_start[1]))) # want to get rid of asteroid too
		print(bad_stars)

		star_x = np.delete(star_x, bad_stars, 0)
		star_y = np.delete(star_y, bad_stars, 0)
		star_x_min = np.delete(star_x_min, bad_stars, 0)
		star_y_min = np.delete(star_y_min, bad_stars, 0)
		star_x_max = np.delete(star_x_max, bad_stars, 0)
		star_y_max = np.delete(star_y_max, bad_stars, 0)

		ax[0].scatter(star_x, star_y, c='orange', s=2, label='centroid')
		ax[0].imshow(img_rotated, cmap='gray', norm=colors.LogNorm(vmin=np.median(sky_row_avg))) #  setting min value to sky background median 


		
		# ax[0].legend()
	
	plt.show()
# output.close()