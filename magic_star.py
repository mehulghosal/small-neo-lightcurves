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
from scipy.optimize import least_squares

# rotate points by angle a [degrees]
def point_rotation(x,y,a,img,img_rot):
	a = -a * np.pi/180
	x_0, y_0 = 0, 0
	x_0_, y_0_ = img.shape[0]*np.abs(np.sin(a)), img.shape[1]*np.abs(np.sin(a))
	x_, y_ = np.array((x-x_0)*np.cos(a) - (y-y_0)*np.sin(a), dtype=int), np.array((x-x_0)*np.sin(a) + (y-y_0)*np.cos(a), dtype=int)
	# to account for direction of rotation
	if a>0: x_+= int(x_0_)
	elif a<0: y_+= int(y_0_)

	if x_<0: x_=0
	if y_<0: y_=0
	# if x_>img_rot.shape[0]: x_=img_rot.shape[0]
	# if y_>img_rot.shape[1]: y_=img_rot.shape[1]

	return x_, y_

def model(x, s, m, a, c, b, d):
	return c*np.exp(-.5* ((x-m)/s)**2) + a*x + b + d*x**2

def quadratic(x, a, b, c, d, e):
	return a*x**2 + b*x + c + d*x**3 + e*x**4


img_rot, centroid, = 0, 0
count = 0

# Veres 2012 eq 3
# r = [x,y], s = sigma, L is length, a is angle, b is background noise (holding constant for now)
# img_rot and centroid are not fitting variables - want to pass these in as constants; centroid = [x,y]
def trail_model(x, y, s, L, a, b, x_0, y_0):

	global img_rot, star_x_ext, star_y_ext, centroid
	
	# print(x, y)

	# trail = img_rot[int(c_y-L/2+0.5):int(c_y+L/2+.5) , int(c_x-s*2.355+.5): int(c_x+s*2.355+.5) ]
	trail = img_rot[int(star_y_ext[0]):int(star_y_ext[1]), int(star_x_ext[0]): int(star_x_ext[1])]

	flux   = np.sum(trail)
	a      = (a) * np.pi/180
	cosine = np.cos(a)
	sine   = np.sin(a)

	flux_term   = flux/(L * 2 * s * (2 * np.pi)**.5)
	exponential = np.exp( -(( (x-x_0)*sine + (y-y_0)*cosine )**2 ) / (2*s**2) )
	erf1 = erf(( (x-x_0) * cosine + (y-y_0) * sine + L/2) / (s*2**.5)) 
	erf2 = erf(( (x-x_0) * cosine + (y-y_0) * sine - L/2) / (s*2**.5))

	return flux_term * exponential * (erf1-erf2) + b

def draw_model(s, L, a, b, c_x, c_y):

	global img_rot, star_x_ext, star_y_ext, centroid

	# dont actually know if this meshgrid business works??? come back to this first if breaks
	# xx, yy = np.meshgrid(np.arange(0, img_rot.shape[1]), np.arange(0, img_rot.shape[0]))
	xx, yy = np.meshgrid( np.arange(0, img_rot.shape[1]), np.arange(0, img_rot.shape[0]) )

	model = trail_model(xx, yy, s, L, a, b, c_x, c_y)	#assuming this is 2FWHM wide and 2L tall

	# print(img.shape, rotate(img,-a).shape)
	return model
	

def residual(par):
	global img_rot, star_x_ext, star_y_ext, centroid
	s, L, a, b, c_x, c_y = par[0], par[1], par[2], par[3], par[4], par[5]
	
	model = draw_model(s, L, a, b, c_x, c_y)

	observed = img_rot[int(star_y_ext[0]):int(star_y_ext[1]),int( star_x_ext[0]): int(star_x_ext[1])]
	# observed_row_sums = np.array([np.sum(i) for i in observed])
	# observed_col_sums = np.sum(observed, axis=0)

	model_slice = model[int(star_y_ext[0]):int(star_y_ext[1]),int( star_x_ext[0]): int(star_x_ext[1])]
	# model_row_sums = np.array([np.sum(i) for i in model_slice])
	# model_col_sums = np.sum(model_slice, axis=0)

	# print(model.shape, observed.shape)
	# print()
	# print(model.shape)
	# print(observed.shape)

	r = np.sqrt(np.sum(model_slice-observed)**2)
	# r = np.sqrt( (model_row_sums-observed_row_sums)**2 ) 
	# r = np.sqrt((model_col_sums-observed_col_sums)**2)

	return r


for d in dir_names:
	file_names = [d+f for f in os.listdir(d) if isfile(join(d,f))]
	yea = False

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

		trail_centroid = np.array([trail_start[0], np.mean([trail_start[1], trail_end[1]])])

		print('trail length: ', trail_length)
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

		ax[0].imshow(img_rotated, cmap='gray', norm=colors.LogNorm(vmin=np.median(sky_row_avg))) #  setting min value to sky background median 

		sigma_row = obj_minus_sky + (len(obj_row_sums)) * (sky_row_avg + hdr['RDNOISE']**2) + (len(obj_row_sums))**2 * sky_row_avg**.5 # from magnier
		sigma_row = sigma_row ** .5

		# x = np.arange(0, 101, 101/len(obj_row_sums))
		x = np.arange(0, len(obj_row_sums), 1)
		# ax[1].plot(x, obj_minus_sky)

		param_vals, param_covs = curve_fit(quadratic, x, obj_minus_sky, sigma=sigma_row)
		# ax[1].plot(x, quadratic(x, *param_vals))
		# ax[1].set_ylim([np.min(obj_minus_sky),25000])
		print(param_vals)
		print(np.diag(param_covs))
		# UNCOMMENT LATER, maybe
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

		# source extractor !!
		sex = subprocess.run(['sex', f, '-DETECT_MINAREA', str(trail_length*fwhm)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
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
			# print(star_x[i], star_y[i])

		# filtering based on distance to asteroid
		num = len(dist_to_asteroid)
		dist_to_asteroid = np.array(dist_to_asteroid)
		dist_sorted = np.argsort(dist_to_asteroid)
		star_x = star_x[dist_sorted][:num]
		star_y = star_y[dist_sorted][:num]
		star_x_min = star_x_min[dist_sorted][:num]
		star_y_min = star_y_min[dist_sorted][:num]
		star_x_max = star_x_max[dist_sorted][:num]
		star_y_max = star_y_max[dist_sorted][:num]

		# filtering bad stars from sextractor
		bad_stars = np.where((star_x<trail_length) | (star_x>img_rotated.shape[1]-trail_length) | (star_y<trail_length) | (star_y>img_rotated.shape[0]-trail_length)) # too close to edge
		bad_stars = np.append(bad_stars, np.where((star_x<trail_start[0]+fwhm) & (star_x>trail_start[0]-fwhm) & (star_y<trail_end[1]) & (star_y>trail_start[1]))) # want to get rid of asteroid too
		bad_stars = np.append(bad_stars, np.where((star_x<trail_start[0]+fwhm) & (star_x>trail_start[0]-fwhm) & (star_y<trail_end[1]) & (star_y>trail_start[1]))) # want to get rid of asteroid too
		print(bad_stars)
		star_x = np.delete(star_x, bad_stars, 0)
		star_y = np.delete(star_y, bad_stars, 0)
		star_x_min = np.delete(star_x_min, bad_stars, 0)
		star_y_min = np.delete(star_y_min, bad_stars, 0)
		star_x_max = np.delete(star_x_max, bad_stars, 0)
		star_y_max = np.delete(star_y_max, bad_stars, 0)
		
		# f_stars, ax_stars = plt.subplots(2, 5)
		
		stars        = []
		trail_starts = []
		trail_ends   = []
		residuals    = []

		for i in range(len(star_x)):
		# for i in range(1):
			centroid = star_x[i], star_y[i]
			img_rot = img_rotated
			x_correction = (star_x_min[i] - star_x_max[i])*.20
			y_correction = (star_y_min[i] - star_y_max[i])*.20
			star_x_ext = int(star_x_min[i]-x_correction), int(star_x_max[i]+x_correction)
			star_y_ext = int(star_y_min[i]-y_correction), int(star_y_max[i]+y_correction)
			# print(centroid, star_x_ext, star_y_ext)

			L_0 = ((star_x_max-star_x_min)**2 + ((star_y_max-star_y_min)))**.5
			a_0 = np.arctan2( star_y_max-star_y_min , star_x_max-star_x_min) # in radians

			p0 = np.array([4, 227, 180-a_0[i]*180/np.pi, np.mean(sky_row_avg), centroid[0], centroid[1]])
			# p0 = np.array([2, L_0[i], 180-a_0[i]*180/np.pi, np.mean(sky_row_avg), centroid[0], centroid[1]])


			fit = least_squares(residual, p0, loss='linear', ftol=0.5, xtol=0.5, gtol=0.5, bounds=([1, trail_length/2, -180, 0, 0, 0],[10, trail_length*5, 180, 5e3, img_rotated.shape[1], img_rotated.shape[0]]))
			residuals.append([residual(p0), residual(fit.x)])

			print('p0:', p0)
			print('residual(p0): ', residuals[i][0])
			print('residual(fit params): ', residuals[i][1])

			# fit = least_squares(residual, p0, loss='linear', method='lm', ftol=0.5, xtol=0.5, gtol=0.5)
			# np.append(stars, fit.x, axis=0)
			# fit = least_squares(residual, p0, loss='linear')
			param = fit.x
			s, L, a, b, x_0, y_0 = fit.x[0], fit.x[1], fit.x[2], fit.x[3], fit.x[4], fit.x[5]
			
			# rotating frame to angle output by fit
			# img_star_rotated = rotate(img_rot, a)
			# x_0, y_0 = point_rotation(x_0, y_0, a, img_rot, img_star_rotated)
			
			# # have to reset global variables
			# centroid = x_0, y_0
			# img_rot  = img_star_rotated
			# star_max = point_rotation(star_x_ext[1], star_y_ext[1], a, img_rot, img_star_rotated)
			# star_min = point_rotation(star_x_ext[0], star_y_ext[0], a, img_rot, img_star_rotated)
			
			# star_x_ext = int(star_min[0]), int(star_max[0])
			# star_y_ext = int(star_min[1]), int(star_max[1])

			# p0 = [s, L, 20, b, x_0, y_0]
			# print(p0, end=' ')
			# fit = least_squares(residual, p0, loss='huber', ftol=.05, xtol=0.05, gtol=.05, bounds=([1, L-10, 0, 0, 0, 0],[10, L*3, 360, 5e3, img_star_rotated.shape[1], img_star_rotated.shape[0]]))
			# fit = least_squares(residual, p0, loss='linear', method='lm', ftol=0.05, xtol=0.05, gtol=0.05)

			print('fit', fit.x)
			print('fit success', fit.success)

			model = draw_model(*fit.x)
			# if i<10: 
			# 	ax_stars[int(i/5), i%5].set_title(str(centroid))
			# 	ax_stars[int(i/5), i%5].imshow(model)
			# 	a/x[3].imshow(model)
			
			# y coords flipped for wacky pyplot reasons
			# star_trail_start = np.array([x_0 - L/2 * np.cos(a*np.pi/180), y_0 + L/2 * np.sin(a*np.pi/180)])
			# star_trail_end   = np.array([x_0 + L/2 * np.cos(a*np.pi/180), y_0 - L/2 * np.sin(a*np.pi/180)])
			star_trail_start = np.array([centroid[0] - L/2 * np.cos(a*np.pi/180), centroid[1] + L/2 * np.sin(a*np.pi/180)])
			star_trail_end   = np.array([centroid[0] + L/2 * np.cos(a*np.pi/180), centroid[1] - L/2 * np.sin(a*np.pi/180)])

			trail_starts.append(star_trail_start)
			trail_ends  .append(star_trail_end  )

			stars.append(param)
			print(' ')


			# if fit.x[-2] == centroid[0] and fit.x[-1] == centroid[1]: print(True)
	
			# p, p_cov = curve_fit(trail_model, coords, flattened_img, p0=[3, trail_length, 0, np.mean(sky_row_avg)])
			
		stars 		 = np.array(stars)
		trail_starts = np.array(trail_starts)
		trail_ends   = np.array(trail_ends)
		print(stars.shape)

		s_std        = np.std(stars[:,0])
		length_std   = np.std(stars[:,1])
		angle_std    = np.std(stars[:,2])

		s_mean  	 = np.mean(stars[:,0])
		length_mean  = np.mean(stars[:,1])
		angle_mean   = np.mean(stars[:,2])

		# throwing away outliers
		threshold = 1 # sigmas

		star_filter  = np.where( (stars[:,0]<=s_mean+threshold*s_std) & (stars[:,0]>=s_mean-threshold*s_std) & (stars[:,1]<=length_mean+threshold*length_std) & (stars[:,1]>=length_mean-threshold*length_std) & (stars[:,2]<=angle_mean+threshold*angle_std) & (stars[:,2]>=angle_mean-threshold*angle_std) )
		stars        = stars       [star_filter]
		trail_starts = trail_starts[star_filter]
		trail_ends   = trail_ends  [star_filter]
		print(stars.shape)

		ax[0].plot([trail_starts[:,0], trail_ends[:,0]], [trail_starts[:,1], trail_ends[:,1]], 'y*', ms=3 )

		ax[0].scatter(star_x, star_y, c='orange', s=2, label='centroid')
		# ax[0].scatter(star_x_min, star_y_min, c='green', s=2, label='mins')
		# ax[0].scatter(star_x_max, star_y_max, c='purple', s=2, label='maxes')
		ax[0].legend()
		yea = True

		# star_residuals = open('star_Residuals.txt', 'w+')
		# star_residuals.writelines(str(residuals))

		# ax[0].legend()
	
	plt.show()
	if yea: break
# output.close()