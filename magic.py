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

# Veres 2012 eq 3
# r = [x,y], s = sigma, L is length, a is angle, b is background noise (holding constant for now)
# img_rot and centroid are not fitting variables - want to pass these in as constants; centroid = [x,y]
def trail_model(x, y, s, L, a, b_1, b_2, b_3, x_0, y_0):

	global img_rot, star_x_ext, star_y_ext, centroid
	
	L_but_longer = L*1.4
	s_but_wider  = s*1.8

	# trail = img_rot[int(c_y-L/2+0.5):int(c_y+L/2+.5) , int(c_x-s*2.355+.5): int(c_x+s*2.355+.5) ]
	trail = img_rot[int(y_0 - L_but_longer/2):int(y_0 + L_but_longer/2) , int(x_0 - s_but_wider*2.355):int(x_0 + s_but_wider*2.355)]
	# print( 'trail shape', trail.shape)

	flux   = np.sum(trail)
	a      = (a) * np.pi/180
	cosine = np.cos(a)
	sine   = np.sin(a)

	flux_term   = flux/(L * 2 * s * (2 * np.pi)**.5)
	exponential = np.exp( -(( (x-x_0)*sine + (y-y_0)*cosine )**2 ) / (2*s**2) )
	erf1 = erf(( (x-x_0) * cosine + (y-y_0) * sine + L/2) / (s*2**.5)) 
	erf2 = erf(( (x-x_0) * cosine + (y-y_0) * sine - L/2) / (s*2**.5))
	background = b_1 + b_2*x + b_3*y

	return flux_term * exponential * (erf1-erf2) + background

def draw_model(s, L, a, b_1, b_2, b_3, c_x, c_y):

	global img_rot, star_x_ext, star_y_ext, centroid

	# dont actually know if this meshgrid business works??? come back to this first if breaks
	# xx, yy = np.meshgrid(np.arange(0, img_rot.shape[1]), np.arange(0, img_rot.shape[0]))
	xx, yy = np.meshgrid( np.arange(0, img_rot.shape[1]), np.arange(0, img_rot.shape[0]) )

	model = trail_model(xx, yy, s, L, a, b_1, b_2, b_3, c_x, c_y)	#assuming this is 2FWHM wide and 2L tall

	# print(img.shape, rotate(img,-a).shape)
	return model
	

def residual(par):
	global img_rot, star_x_ext, star_y_ext, centroid
	s, L, a, b_1, b_2, b_3, x_0, y_0 = par[0], par[1], par[2], par[3], par[4], par[5], par[6], par[7]
	
	model = draw_model(s, L, a, b_1, b_2, b_3, x_0, y_0)

	# L_but_longer = L*1.2
	s_but_wider  = s

	box_y_width = np.abs(star_y_ext[1] - star_y_ext[0]) * 1.4

	# observed = img_rot[int(y_0 - L_but_longer/2):int(y_0 + L_but_longer/2) , int(x_0 - s_but_wider*2.355):int(x_0 + s_but_wider*2.355)]
	observed = img_rot[int(centroid[1] - box_y_width/2):int(centroid[1] + box_y_width/2) , int(centroid[0] - s_but_wider*2*2.355):int(centroid[0] + s_but_wider*2*2.355)]
	# observed_row_sums = np.array([np.sum(i) for i in observed])
	# observed_col_sums = np.sum(observed, axis=0)

	# model_slice = model[int(y_0 - L_but_longer/2):int(y_0 + L_but_longer/2) , int(x_0 - s_but_wider*2.355):int(x_0 + s_but_wider*2.355)]
	model_slice = model[int(centroid[1] - box_y_width/2):int(centroid[1] + box_y_width/2) , int(centroid[0] - s_but_wider*2*2.355):int(centroid[0] + s_but_wider*2*2.355)]
	# model_row_sums = np.array([np.sum(i) for i in model_slice])
	# model_col_sums = np.sum(model_slice, axis=0)

	# print(model_slice.shape, observed.shape)
	# print()
	# print(model.shape)
	# print(observed.shape)

	r = np.sqrt(np.sum(model_slice-observed)**2)
	# r = np.sqrt( (model_row_sums-observed_row_sums)**2 ) 
	# r = np.sqrt((model_col_sums-observed_col_sums)**2)

	return r


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

		# if '2016 GE1' not in obj_id: continue
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

		dist_to_asteroid = np.array(dist_to_asteroid)
		dist_sorted = np.argsort(dist_to_asteroid)
		star_x = star_x[dist_sorted][:10]
		star_y = star_y[dist_sorted][:10]
		star_x_min = star_x_min[dist_sorted][:10]
		star_y_min = star_y_min[dist_sorted][:10]
		star_x_max = star_x_max[dist_sorted][:10]
		star_y_max = star_y_max[dist_sorted][:10]
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



		L_0 = ((star_x_max-star_x_min)**2 + ((star_y_max-star_y_min)))**.5
		# to rotate image -- negative angle from vertical
		a_0 = -1*np.arctan2( star_x_max-star_x_min,  star_y_max-star_y_min) * 180/np.pi
		
		f_stars, ax_stars = plt.subplots(2, 5)
		
		stars        = []
		trail_starts = []
		img_stars    = []
		trail_ends   = []
		residuals    = []


		for i in range(len(star_x)):
		# for i in range(1):
			centroid = star_x[i], star_y[i]
			# img_rot = img_rotated
			# x_correction = (star_x_min[i] - star_x_max[i])*.10
			# y_correction = (star_y_min[i] - star_y_max[i])*.10
			x_correction, y_correction = 0,0
			star_x_ext = int(star_x_min[i]-x_correction), int(star_x_max[i]+x_correction)
			star_y_ext = int(star_y_min[i]-y_correction), int(star_y_max[i]+y_correction)
			# print(centroid, star_x_ext, star_y_ext)

			img_star_rotated = rotate(img_rotated, a_0[i])
			img_rot = img_star_rotated

			centroid = point_rotation(centroid[0], centroid[1], a_0[i], img_rotated, img_star_rotated)
			upper_left = point_rotation(star_x_min[i] , star_y_min[i] , a_0[i], img_rotated, img_star_rotated)
			lower_rite = point_rotation(star_x_max[i] , star_y_max[i] , a_0[i], img_rotated, img_star_rotated)

			star_x_ext = upper_left[0], lower_rite[0]
			star_y_ext = upper_left[1], lower_rite[1]

			p0 = np.array([3, 225, 90, np.mean(sky_row_avg), 0, 0, centroid[0], centroid[1]])
			# p0 = np.array([2, L_0[i], 180-a_0[i]*180/np.pi, np.mean(sky_row_avg), centroid[0], centroid[1]])

			param_bounds = ([1, L_0[i]/2, 90+2*a_0[i], 0, 0, 0, 0, 0], [10, L_0[i]*5, 90-2*a_0[i], 2e3, 500, 500, img_star_rotated.shape[1], img_star_rotated.shape[0] ])

			fit = least_squares(residual, p0, loss='huber', ftol=0.5, xtol=0.5, gtol=0.5, bounds=param_bounds)
			residuals.append([residual(p0), residual(fit.x)])

			print('p0:', p0)
			print('residual(p0): ', residuals[i][0])
			print('residual(fit params): ', residuals[i][1])

			param = p0
			s, L, a, b, x_0, y_0 = fit.x[0], fit.x[1], fit.x[2], fit.x[3], fit.x[4], fit.x[5]
			# s, L, a, b, x_0, y_0 = p0[0], p0[1], p0[2], p0[3], p0[4], p0[5]

			print('fit', fit.x)
			print('fit success', fit.success)
			star_model = draw_model(*fit.x)
			# L*=1.2
			# s*=1.2
			
				# ax[3].imshow(model)
			
			# y coords flipped for wacky pyplot reasons
			# star_trail_start = np.array([x_0 - L/2 * np.cos(a*np.pi/180), y_0 + L/2 * np.sin(a*np.pi/180)])
			# star_trail_end   = np.array([x_0 + L/2 * np.cos(a*np.pi/180), y_0 - L/2 * np.sin(a*np.pi/180)])
			
			# now rotating back to the asteroids reference
			centroid = star_x[i], star_y[i]
			a -= a_0[i]

			star_trail_start = np.array([centroid[0] - L/2 * np.cos(a*np.pi/180), centroid[1] + L/2 * np.sin(a*np.pi/180)])
			star_trail_end   = np.array([centroid[0] + L/2 * np.cos(a*np.pi/180), centroid[1] - L/2 * np.sin(a*np.pi/180)])

			trail_starts.append(star_trail_start)
			trail_ends  .append(star_trail_end  )
			img_stars   .append(img_star_rotated)
			stars       .append(param)
			
			print(' ')

			


			# if fit.x[-2] == centroid[0] and fit.x[-1] == centroid[1]: print(True)
	
			# p, p_cov = curve_fit(trail_model, coords, flattened_img, p0=[3, trail_length, 0, np.mean(sky_row_avg)])
			
		stars 		 = np.array(stars)
		trail_starts = np.array(trail_starts)
		trail_ends   = np.array(trail_ends)
		print('initially, ', stars.shape[0])

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
		print('filtering: ', stars.shape[0])

		ax[0].plot([trail_starts[:,0], trail_ends[:,0]], [trail_starts[:,1], trail_ends[:,1]], 'y*', ms=3 )

		ax[0].scatter(star_x, star_y, c='orange', s=2, label='centroid')
		# ax[0].scatter(star_x_min, star_y_min, c='green', s=2, label='mins')
		# ax[0].scatter(star_x_max, star_y_max, c='purple', s=2, label='maxes')
		ax[0].legend()

		# lightcurves of stars
		for i in range(len(stars)):
		# for i in range(1):
			trail_end     = np.array(point_rotation(trail_starts[i,0], trail_starts[i,1], a_0[i], img_rotated, img_star_rotated))
			trail_start   = np.array(point_rotation(trail_ends  [i,0], trail_ends  [i,1], a_0[i], img_rotated, img_star_rotated))
			print(trail_start, trail_end)

			fwhm = stars[i,0] * 2.355
			L = int(stars[i,1]*.2+.5)

			trail_start[1] -= L
			trail_end[1]   += L

			img_star_rotated = img_stars[i]

			str_width = int(1*fwhm)
			sky_width = int(2*fwhm)
			str_rect = img_star_rotated[trail_start[1]:trail_end[1], trail_start[0]-str_width:trail_start[0]+str_width]

			str_row_sums = np.array([np.sum(j) for j in str_rect])

			sky_left  = img_star_rotated[trail_start[1]:trail_end[1], trail_start[0]-str_width-sky_width:trail_start[0]-str_width]
			sky_right = img_star_rotated[trail_start[1]:trail_end[1], trail_start[0]+str_width:trail_start[0]+str_width+sky_width]

			sky_left_row_sum  = np.array([np.sum(i) for i in sky_left ])
			sky_right_row_sum = np.array([np.sum(i) for i in sky_right])
			sky_row_avg = (sky_right_row_sum+sky_left_row_sum)/(sky_right.shape[1]+sky_left.shape[1])

			str_minus_sky = str_row_sums - sky_row_avg * str_rect.shape[1]

			sigma_row = str_minus_sky + (len(str_row_sums)) * (sky_row_avg + hdr['RDNOISE']**2) + (len(str_row_sums))**2 * sky_row_avg**.5 # from magnier
			sigma_row = sigma_row ** .5
			
			x = np.arange(0, len(str_row_sums))

			normalize = np.max(str_minus_sky)
			# str_minus_sky /= normalize

			n_bins   = 5  # smooths over this many bins
			smoothed = np.array([np.sum(str_minus_sky[j:j+n_bins])/n_bins for j in range(0, len(str_minus_sky)-n_bins, n_bins)])
			x_smooth = np.arange(0, len(smoothed))

			print(str_minus_sky.shape)
			print(smoothed.shape)


			param_star, param_covs_star = curve_fit(star_box_model, x_smooth, smoothed, p0=[70, 290, 2000, 2*np.pi/10, 0, 17300, .01, 2*np.pi/20, 0, 0])
			print(i,param_star)


			if i<10 and i<stars.shape[0] : 
				ax_stars[int(i/5), i%5].set_title(str([star_x[i], star_y[i]]))
				# ax_stars[int(i/5), i%5].imshow(img_star_rotated[int(centroid[1] - L/2): int(centroid[1] + L/2), int(centroid[0] - s*2.355): int(centroid[0] + s*2.355)])		
				ax_stars[int(i/5), i%5].plot(x, str_minus_sky)
				# ax_stars[int(i/5), i%5].plot(x_smooth, smoothed)
				# ax_stars[int(i/5), i%5].plot(x, star_box_model(x,*param_star))
				# ax_stars[int(i/5), i%5].errorbar(x, str_minus_sky, yerr = sigma_row, fmt='o', capsize=3, linewidth=2, elinewidth=1)



		# ax[0].legend()
	
	plt.show()
# output.close()