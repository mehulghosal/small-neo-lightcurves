import warnings, subprocess, sys
import numpy as np
import astropy as ap

# import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.table import Table
from astropy.timeseries import LombScargle
from astropy.timeseries import TimeSeries
# from matplotlib import colors
from astropy.io import fits
from scipy.ndimage import rotate
from scipy.special import erf
from astropy.wcs import WCS
from astropy.wcs import utils
from astropy.coordinates import SkyCoord
from astropy import units as u

from astropy.utils.exceptions import AstropyWarning

try:
	f_name 		 = sys.argv[1]
	l_from_input = sys.argv[2]
	a_from_input = sys.argv[3]
	write_output = sys.argv[4]
except Exception as e:
	print(e)

# plt.rcParams.update({'figure.max_open_warning': 0})
warnings.simplefilter('ignore', AstropyWarning)

# initializing all directories
import os
from os.path import isdir, isfile, join
directory = './'	
dir_names = [directory+f+'/' for f in os.listdir(directory) if isdir(join(directory,f))] 

input_file = np.loadtxt('input.csv', dtype=object, skiprows=1, usecols=(i for i in range(25)), delimiter=',')

se_dir   = './SEoutput/'
se_files = [se_dir+f for f in os.listdir(se_dir) if isfile(join(se_dir,f))]

mins = {'g':100, 'r': 150, 'i': 250}

from scipy.optimize import curve_fit
from scipy.optimize import least_squares


"""
rotate points by angle a [degrees]
	origin (0,0) is to pleft of image. +x to the right, +y down

PARAMETERS
-----------

x 		: float
	ccd column pixel coordinate
y 		: float
	ccd row pixel coordinate
img 	: array
	original image rotated from
img_rot : array
	rotated image 

RETURNS
--------
x_ : float
	rotated CCD column pixel coordinate
y_ : float
	rotated CCD row pixel coordinate
"""
def point_rotation(x,y,a,img,img_rot):
	a = -a * np.pi/180
	x_0, y_0 = 0, 0
	x_0_, y_0_ = img.shape[0]*np.abs(np.sin(a)), img.shape[1]*np.abs(np.sin(a))
	x_, y_ = np.array((x-x_0)*np.cos(a) - (y-y_0)*np.sin(a), dtype=int), np.array((x-x_0)*np.sin(a) + (y-y_0)*np.cos(a), dtype=int)
	# to account for direction of rotation
	if   a>0: x_+= int(x_0_)
	elif a<0: y_+= int(y_0_)

	if x_<0: x_=0
	if y_<0: y_=0

	return x_, y_


'''
taking a lightcurve of streaked artifact in CCD image

PARAMETERS
----------
img 				: array, dtype=float 
	numpy nxm array representing CCD image
trail_start 		: array, dtype=float
	x, y CCD pixel coordinates of trail start (~top)
trail_end 			: array, dtype=float 
	x, y CCD pixel coordinates of trail end (~bottom)
fwhm 				: float
	(optional) FWHM of trailed Gaussian in CCD pixels; default = 4
b 					: float 
	(optional) average sky flux contribution per pixel; defaults to calculating sky from region around trail
height_correction 	: float
	(optional) pixels to extend lightcurve above/below trail start/stop; default = 0
display 			: bool
	(optional) will plt.show the lightcurve and mess the mojo up ; default = False
err 				: bool
	(optional) whether to return uncertainties on flux measurements ; default = False
binnning			: int
	(optional) will bin_lightcurve(lightcurve, binning), returned lightcurve will have len = binning
gain 				: float 
	(optional) e-/ADU for CCD image ; default = 1.6
rd_noise 			: float
	(optional) read noise in ADU of image ; default = 3
obj_width			: float
	(optional) width (in FWHM) either side of centroid to sum for object flux ; default = 1
sky_width 			: float 
	(optional) width (in FWHM) either side of object box to sum for sky flux ; default = 4
autotrim 			: bool 
	(optional) if True, will run obj_row_sums through curve_fit with another_box() !! doesnt do anything yet !!

RETURNS
--------

r : list
	to return uncertainties on measurement, err = True
		returns [fluxes[array(dtype=float)], uncertainties[array(dtype=float)], average sky measurement[float]]
	without uncertainties, returns: [ fluxes : array(dtype=float) ]
'''
def take_lightcurve(img, trail_start, trail_end, fwhm=4, b=None, height_correction=0, display=False, err=False, binning=None, gain=1.6, rd_noise=3, obj_width=1, sky_width=4, autotrim=False):
	obj_width = obj_width*fwhm
	sky_width = sky_width*fwhm

	trail_start_y = trail_start[1]
	trail_end_y   = trail_end  [1]

	if not height_correction == 0:
		trail_start_y -= height_correction
		trail_end_y   += height_correction

	obj_rect  = img[int(trail_start_y + .5 ):int(trail_end_y + .5 ), int(trail_start[0] - obj_width + .5):int(trail_end[0]+obj_width + .5)]

	sky_left  = img[int(trail_start_y + .5 ):int(trail_end_y + .5 ), int(trail_start[0] - obj_width - sky_width + .5) : int(trail_start[0] - obj_width + .5)]
	sky_right = img[int(trail_start_y + .5 ):int(trail_end_y + .5 ), int(trail_start[0] + obj_width + .5) : int(trail_start[0] + obj_width + sky_width + .5)]

	obj_row_sums      = np.array([np.sum(i) for i in obj_rect ])
	sky_left_row_sum  = np.array([np.sum(i) for i in sky_left ])
	sky_right_row_sum = np.array([np.sum(i) for i in sky_right])
	sky_row_sum       = sky_right_row_sum+sky_left_row_sum  # total sky counts

	if binning is not None:
		obj_row_sums = bin_lightcurve(obj_row_sums, binning)
		sky_row_sum  = bin_lightcurve(sky_row_sum , binning)

	sky_n_pixels      = sky_left_row_sum.size+sky_right_row_sum.size		# num sky pixels
	sky_row_avg       = sky_row_sum/sky_n_pixels							# sky counts/n_px

	if b is not None:
		sky_row_avg = b

	obj_minus_sky = obj_row_sums - sky_row_avg * obj_rect.shape[1]

	sigma_row = obj_minus_sky/gain + (obj_rect.shape[1]) * (sky_row_avg/gain + rd_noise**2) + (obj_rect.shape[1])**2 * (sky_row_sum**.5 / sky_n_pixels)**2 # from magnier
	sigma_row = sigma_row ** .5

	if display:
		plt.figure()
		t = np.arange(len(obj_minus_sky))
		plt.scatter(t, obj_minus_sky)

	r = []

	r.append(obj_minus_sky)
	
	if err: 
		r.append(sigma_row)
		r.append(sky_row_avg)
	return r

'''
another attempt at binning - this time extracting fractional pixel fluxes.

PARAMETERS
----------
lc 			 [array(dtype=float)]: original lightcurve to be rebinneddisplay
trail_length [int]				 : length we want our output lightcurve


RETURNS 
---------
array: 
	(shape 1xtrail_length) of same dtype as lc. 
	reorganizes bins and takes fractional pixel fluxes across adjacent pixels according to ratio between len(lc) and trail_length
'''
def bin_lightcurve(lc, trail_length):
	L = len(lc)
	length_ratio = L/trail_length

	binned = []
	for i in range(0, int(trail_length)):
		first = lc[ int( i * length_ratio ) ] * ( 1 - (i * length_ratio)%1  )
		j = (i+1) * length_ratio 
		secon = 0
		if not int(j) == len(lc): 
			# print(j)
			secon = lc[ int(j) ] * ( j % 1)
		third = np.sum( lc[ int( i * length_ratio + 1 ) : int(j) ] )
		s = first + secon + third
		# print(int(c+1), int(j), s, lc[int(c)] * (1-c%1), lc[int(j+1)] * (j%1))
		# c = j 
		binned.append(s)

	return np.array(binned)


'''
folding lightcurves on dominant period w/ astropy timeseries

PARAMETERS
----------
time 	   : array, dtype=float
	expects time values in mjd, same length as lightcurve
lightcurve : array, dtype=float
	lightcurve data values
period 	   : float 
	dominant period, 2x peak Lomb-Scargle period
exp_time   : float 
	doesn't change anything rn
phase 	   : float 
	doesn't change anythign rn

RETURNS
--------
phase : array  
	phase is the range of -period/2 to period/2
data  : array-like
	data is the folded lightcurve data
'''
def fold_lightcurve( time , lightcurve , period , exp_time=60 , phase=0 ):
	lightcurve_table = Table([Time(time, format='mjd'), lightcurve], names=('time', 'data'))
	ts = TimeSeries(data=lightcurve_table)
	folded_lc = ts.fold( period=period*u.second, normalize_phase=False)
	phase = np.array(folded_lc['time'].value)
	data  = np.array(folded_lc['data'])
	return phase, data

'''
lomb scargle periodogram of lightcurve

PARAMETERS
----------

time 	   : array, dtype=float 
	expects time values in mjd
lightcurve : array, dtype=float
	lightcurve data values, same length as time
num_maxes  : int		  
	(optional) number of maxima returned in (peak_period)
err 	   : array, dtype=float
	(optional) uncertainty on each value in lightcurve

RETURNS
----------
period      : array 
	range of periods
power       : array
	power associated with each period
peak_period : tuple
	num_maxes long tuple with dominant periods
'''
def periodogram(time, lightcurve, num_maxes=1, err=None):
	if err is None:
		frequency, power = LombScargle(time, lightcurve).autopower()
	else:
		frequency, power = LombScargle(time, lightcurve, err).autopower()
	period = 1/frequency * 24*3600
	# period = 1/frequency

	# peak_frequency = frequency[np.argmax(power)]
	peak_period    = period[(-power).argsort()[:num_maxes]]
	# print('peak period: ', peak_period )

	return period, power, (peak_period)


'''
perpendicular summing to lightcurves -- testing trailed point spread function by summing along columns
runs this summed/replicated PSF through curve_fit, "model" is  a Gaussian

PARAMETERS
----------
img 		: array, dtype=float 
	numpy nxm array representing CCD image
trail_start : array, dtype=float
	x, y CCD pixel coordinates of trail start (~top)
trail_end 	: array, dtype=float 
	x, y CCD pixel coordinates of trail end (~bottom)
obj_width 	: float
	(optional) CCD pixel columns to sample over ; default=25
display 	: bool
	(optional) !! doesn't do anything yet !!

RETURNS
----------
param_vals : array 
	curve_fit best fit parameters: [ s , m , a , c , b , d ]
param_covs : array
	covariance matrix of best fit parameteres. to get 1 sigma errors on params, np.sqrt(np.diag(param_covs))
obj_width  : tuple
	returns the same param?  forgot why i needed this lol
'''
def trail_spread_function(img, trail_start, trail_end, obj_width=25, display = False):
		
	obj_rect = img[int(trail_start[1] + .5):int(trail_end[1] + .5), int(trail_start[0]-obj_width + .5):int(trail_start[0]+obj_width + .5)]

	col_sums = np.sum(obj_rect, axis=0)
		
	# col_sums /= np.max(col_sums)
	rect_width = np.arange(0, 2*obj_width, 1)
	param_vals, param_covs = curve_fit( gaussian_1D , rect_width , col_sums , p0=[ 3 , obj_width , .03 , 60000 , 20000 , -3 ] )

	# ax[2].scatter(rect_width, col_sums, label='column sums')
	# ax[2].plot(rect_width, model(rect_width, *param_vals), label='model fit')
	# ax[2].legend()
	return param_vals, param_covs, obj_width

'''
one dimensional Gaussian function - used for curve_fit in trail_spread_function

PARAMETERS
----------
x : array
	independent axis
s : float
	Gaussian spread
m : float
	center of spread along x
a : float
	coefficient of linear term in background estimate
c : float
	vertical scaling of Gaussian 
b : float
	coefficient of constant term in background
d : float
	coefficient of quadratic term in background


RETURNS
----------
returns Gaussian model used for curve_fit 
'''
def gaussian_1D(x, s, m, a, c, b, d):
	return c*np.exp(-.5* ((x-m)/s)**2) + a*x + b + d*x**2


'''
one dimensional box function - used for nothing (hopefully) 
	~ potentially for fitting lightcurves to get start/endpoints 

PARAMETERS
----------
x   : array
	independent axis
t_1 : float
	coordinate along x of first boundry
t_2 : float
	coordinate along x of second boundry
s_1 : float
	value of lc where x < t_1
s_2 : float
	value of lc where t1 <= x <= t_2 
s_3 : float
	value of lc where x > t_2


RETURNS
----------
returns box model used for curve_fit 
'''
def another_box( x , t_1 , t_2 , s_1 , s_2 , s_3  ):
	r = np.zeros(x.shape)
	r[               : int(t_1 + .5) ] += s_1
	r[ int(t_1 + .5) : int(t_2 + .5) ] += s_2
	r[ int(t_2 + .5) :               ] += s_3
	return r

'''
one dimensional Fourier series function - used for nothing (hopefully) 
	~ potentially for fitting lightcurves to get pretties

PARAMETERS
----------
x      : array
	independent axis
params : array
	len(rows) is number of Fourier terms 
		a ~ vertical scaling of sine term
		b ~ frequency of sine term
		c ~ vertical offset of term

RETURNS
----------
returns Fourier model used for curve_fit 
'''
def fourier(x, *params):
    params = np.array(params).reshape(-1,3)
    a = params[:, 0]
    b = params[:, 1]
    c = params[:, 2]
    ret = a[0] * np.sin(np.pi / b[0] * x) + c[0]
    for deg in range(1, len(a)):
        ret += a[deg] * np.sin((deg+1) * np.pi / b[deg] * x) + c[deg]
    return ret

img_rot, centroid, = 0, 0
count = 0

'''

scipy.optimize.curve_fit attempt: 6/13/2022
2d trail function from Veres 2012, Gaussian convolved with a straight line


! obj type for now -- some typa np ndarray
coord [obj  ]: is this just meshgrid ? 
	what if i wrap the 
		xx, yy = np.meshgrid( np.arange(0, img_rot.shape[1]), np.arange(0, img_rot.shape[0]) ) 
	into [xx, yy]

	ORRRRR, i just give it img_rot.shape that draw_model uses,, and it creates the xx, yy
	maybe img_rot?

PARAMETERS
-----------
s   : float 
	Gaussian spread 
L   : float
	trail length
a   : float
	angle from positive horizontal
b   : float
	constant estimate of background flux
x_0 : float
	CCD pixel column number of trail centroid 
y_0 : float
	CCD pixel row number of trail centroid

RETURNS
-------- 
	trail spread function for curve_fit's pleasure, this is flattened

	parameters will be an array of 6 floats

	covariance matrix will be 5x5 matrix of cross correlation of parameters square root of diagonals gives uncertainties on fit parameters

'''
def trail_model_2d( coord , s , L , a , b , x_0 , y_0 ):

	global flux, img_rot

	model  = draw_model ( s , L , a, b , x_0 , y_0 )
	return model.flatten()


'''
actually doing the Veres 2012 eq 3 calculations for every x, y given
	in usage, have to explicitely define img_rot and flux, so these are the variables you expect!!

PARAMETERS
-----------
x 	: array, dtype=int
	array of CCD column coordinates in px
y 	: array, dtype=int
	array of CCD row coordinates in px
s   : float 
	Gaussian spread 
L   : float
	trail length
a   : float
	angle from positive horizontal
b_1 : float
	constant estimate of background flux
x_0 : float
	CCD pixel column number of trail centroid 
y_0 : float
	CCD pixel row number of trail centroid

RETURNS
-------- 
	2d numpy array -- same shape as img_rot. scalar background estimate is b_1, with trail drawn vertically at x_0, y_0
'''
def trail_model(x, y, s, L, a, b_1, x_0, y_0):

	global img_rot, flux
	
	# ok i think this needs to be > 1
	L_but_longer = L*1.3
	s_but_wider  = s*1.3

	# trail = img_rot[int(c_y-L/2+0.5):int(c_y+L/2+.5) , int(c_x-s*2.355+.5): int(c_x+s*2.355+.5) ]
	trail = img_rot[int(y_0 - L_but_longer/2):int(y_0 + L_but_longer/2+1) , int(x_0 - s_but_wider*2.355 + .5):int(x_0 + s_but_wider*2.355 + .5)]
	# print( 'trail shape', trail.shape)

	flux   = np.sum(trail)
	a      = (a) * np.pi/180
	cosine = np.cos(a)
	sine   = np.sin(a)

	flux_term   = flux/(L * 2 * s * (2 * np.pi)**.5)
	exponential = np.exp( -(( (x-x_0)*sine + (y-y_0)*cosine )**2 ) / (2*s**2) )
	erf1 = erf(( (x-x_0) * cosine + (y-y_0) * sine + L/2) / (s*2**.5)) 
	erf2 = erf(( (x-x_0) * cosine + (y-y_0) * sine - L/2) / (s*2**.5))
	background = b_1 

	return flux_term * exponential * (erf1-erf2) + background


'''
driver function for trail_model, but used in trail_model_2d to get the entire image.
	in usage, have to explicitely define img_rot

PARAMETERS
-----------
s   : float 
	Gaussian spread 
L   : float
	trail length
a   : float
	angle from positive horizontal
b_1 : float
	constant estimate of background flux
x_0 : float
	CCD pixel column number of trail centroid 
y_0 : float
	CCD pixel row number of trail centroid

RETURNS
-------- 
	2d numpy array -- same shape as img_rot
'''
def draw_model(s, L, a, b_1, c_x, c_y):

	global img_rot

	# dont actually know if this meshgrid business works??? come back to this first if breaks
	# xx, yy = np.meshgrid(np.arange(0, img_rot.shape[1]), np.arange(0, img_rot.shape[0]))
	xx, yy = np.meshgrid( np.arange(0, img_rot.shape[1]), np.arange(0, img_rot.shape[0]) )
	model = trail_model(xx, yy, s, L, a, b_1, c_x, c_y)	#assuming this is 2FWHM wide and 2L tall

	# print(img.shape, rotate(img,-a).shape)
	return model
	
'''

I guess this is where the shitshow begins i guess
	this will be less well commented for a while just because i haven't had the care to.
	so a lot of it will be ugly and inexplicable unless u ask me.

'''

if __name__ == '__main__':
	
	for d in dir_names:
		file_names = [d+f for f in os.listdir(d) if isfile(join(d,f))]
		yea = False

		# if not ('XD169' in d): continue
		if not f_name in d: continue
		#if not ('2015_TG24' in d or '2016_NM15' in d or '2015_VH1' in d): continue

		start_times = []
		lightcurves = []
		errors      = []

		for f in file_names:
			# if '06o13' not in f: continue
			try:
				file = fits.open(f)
				print(f)
			except Exception as e:
				print(f)
				continue
			hdr = file[0].header
			img = file[0].data

			exp_time   = float(hdr['EXPMEAS'])
			gain       = float(hdr['GAIN'])
			rd_noise   = float(hdr['RDNOISE'])
			# obs_filter = float(hdr['FILTE'])

			# object id from directory name --> string splicing
			obj_id = f.split('_')
			obj_id = obj_id[0][2:] + ' ' + obj_id[1]

			obj_rows = input_file[np.where(input_file[:,1]==obj_id),:][0]
			
			try:
				obj = obj_rows[np.where(obj_rows[:,0]==f.split('/')[-1])][0]
				trail_start = np.array(obj[-4:-2], dtype=int)
				trail_end	= np.array(obj[-2:], dtype=int)
				start_time  = float(obj[6])
			except Exception as e:
				# print(f,obj[-4:-2],obj[-2:])
				# plt.close()
				continue

			# global variable flux to capture the total flux of the trail
			flux = 0

			# NEGATIVE ANGLE OF ASTEROID TRAIL WRT HOME FRAME			
			angle       = -1*np.arctan2(trail_end[0]-trail_start[0], trail_end[1]-trail_start[1]) * 180/np.pi
			# IMG ROTATED TO ASTEROID TRAIL IS VERTICAL
			img_rotated = rotate(img, angle)

			ast_trail_start  = np.array(point_rotation(trail_start[0], trail_start[1], angle, img, img_rotated), dtype=int)
			ast_trail_end	 = np.array(point_rotation(trail_end  [0], trail_end  [1], angle, img, img_rotated), dtype=int)
			ast_trail_length = ast_trail_end[1] - ast_trail_start[1]

			# DOING TRAIL SPREAD TO GET FIRST APPROX FOR FWHM
			trail_spread, trail_spread_covs, trail_width = trail_spread_function(img_rotated, ast_trail_start, ast_trail_end, display=False)
			fwhm = int(trail_spread[0] * 2.355 + .5)
			# correcting trail start/end
			centroid_deviation  = trail_spread[1] - trail_width # if negative, trail is to the left, if positive, trail to right
			ast_trail_start[0] += int(centroid_deviation+.5)
			ast_trail_end  [0] += int(centroid_deviation+.5)

			trail_centroid = np.array([ast_trail_start[0], np.mean([ast_trail_start[1], ast_trail_end[1]])])

			# ASTEROID TRAIL FITTING
			img_rot     = img_rotated
			centroid    = trail_centroid
			# box_x_width = 30
			# box_y_width = ast_trail_length * 2

			p0           = np.array([trail_spread[0], ast_trail_length, 90, 200, trail_centroid[0], trail_centroid[1]])
			param_bounds = ([1, ast_trail_length/4, -180, 0, 0, 0], [15, 650, 180, 2e3, img_rotated.shape[1], img_rotated.shape[0] ])

			# TRAIL FITTING ATTEMPT WITH scipy.optimize.least_squares()
			# fit          = least_squares(residual, p0, loss='linear', ftol=0.05, xtol=0.05, gtol=0.05, bounds=param_bounds)

			# img_slice = img_rotated[int(centroid[1] - box_y_width/2 + .5):int(centroid[1] + box_y_width/2 + .5) , int(centroid[0] - box_x_width/2 + .5):int(centroid[0] + box_x_width/2 + .5)]

			# TRAIL FITTING WITH scipy.optimize.curve_fit()
			# ast_param , ast_param_cov = curve_fit(trail_model_2d, img_rotated, img_slice.flatten(), p0=p0)	
			#  display_streak(img_rot, 10, 300, 90, b, x_0, y_0, width=2)
			ast_param , ast_param_cov = curve_fit(trail_model_2d, img_rotated, img_rot.flatten(), p0=p0)	


			ast_flux = flux
			
			print('asteroid p0[  s , L , a , b , x_0 , y_0 ]: '            , p0)
			print('asteroid fit parameters [ s , L , a , b , x_0 , y_0 ]: ', ast_param)
			print('parameter uncertainties: '                              , np.sqrt(np.diag(ast_param_cov)))

			ast_fwhm			  = ast_param[0] * 2.355
			ast_trail_length	  = ast_param[1]
			# ast_height_correction = ast_trail_length * 0
			ast_height_correction = - int( ast_fwhm ) - 1
			# ast_height_correction = - ast_fwhm
			
			trail_centroid 		  = np.array([ast_param[4], ast_param[5]])

			ast_trail_start = np.array([trail_centroid[0] , trail_centroid[1] - ast_trail_length/2 ])
			ast_trail_end   = np.array([trail_centroid[0] , trail_centroid[1] + ast_trail_length/2 ])


			print('asteroid trail length: ', ast_trail_length)

			# asteroid trail length in 70o13 is 101 tall

			obj_minus_sky, sigma_row, sky_row_avg = take_lightcurve(img_rotated, ast_trail_start, ast_trail_end, fwhm=ast_fwhm, b=None, height_correction=ast_height_correction, display=False, err=True, gain=gain, rd_noise=rd_noise)
			
			# WCS stuff
			w = WCS(hdr)
			c = SkyCoord(f'{obj[7]} {obj[8]}', unit=(u.deg, u.deg))
			target_x, target_y = np.round(utils.skycoord_to_pixel(c, w))
			
			# source extractor !!
			# sex = subprocess.run(['sex', f, '-DETECT_MINAREA', str(trail_length*fwhm), '-CATALOG_NAME', '_'.join(f.split("/")[1:])[:-4] + '.cat'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
			try:
				se_index = [x for x in se_files if (f.split('/')[1] in x and f.split("/")[2].split(".")[0] in x)][0]
			except Exception as e:
				print(e)
				continue

			sex_output = np.loadtxt(se_index, skiprows=9)
			print('SExtractor found stars: ',sex_output.shape[0])
			star_x = sex_output[:,5]
			star_y = sex_output[:,6]

			dist_to_asteroid = []

			# to rotate to asteroid's reference -- SExtractor works with raw fits file data
			for i in range(len(star_x)):
				# star_x[i], star_y[i] = point_rotation(star_x[i], star_y[i], angle, img, img_rotated)
				# star_x_min[i], star_y_min[i] = point_rotation(star_x_min[i] , star_y_min[i] , angle, img, img_rotated)
				# star_x_max[i], star_y_max[i] = point_rotation(star_x_max[i] , star_y_max[i] , angle, img, img_rotated)

				dist_to_asteroid.append( ((star_x[i] - trail_centroid[0]) **2 + (star_y[i] - trail_centroid[1]) **2 ) **.5 )
				
			# filtering based on distance to asteroid
			dist_to_asteroid = np.array  (dist_to_asteroid)
			dist_sorted      = np.argsort(dist_to_asteroid)

			star_x     = star_x[dist_sorted]
			star_y     = star_y[dist_sorted]

			# filtering bad stars from sextractor
			bad_stars = np.where((star_x<ast_trail_length) | (star_x>img_rotated.shape[1]-ast_trail_length) | (star_y<ast_trail_length) | (star_y>img_rotated.shape[0]-ast_trail_length)) # too close to edge
			bad_stars = np.append(bad_stars, 0) # want to get rid of asteroid too
			# bad_stars = np.append(bad_stars, np.where((star_x<trail_start[0]+fwhm) & (star_x>trail_start[0]-fwhm) & (star_y<trail_end[1]) & (star_y>trail_start[1]))) # want to get rid of asteroid too
			print('filter on sextractor', len(bad_stars))
			star_x     = np.delete(star_x, bad_stars, 0)
			star_y     = np.delete(star_y, bad_stars, 0)
			
			l = float(l_from_input)
			a = float(a_from_input)
			
			cheats_on = False

			try:
				# star parametsrs
				cheat_codes = np.loadtxt(f'{f[:-4]}_params.txt')
				cheats_on   = True
			except Exception as e:
				print('Invalid cheat code: ', e, cheats_on)
			
			stars        = []
			trail_starts = []
			trail_ends   = []
			residuals    = []
			row_errs     = []
			row_flux     = []

			failed_log   = []

			i = 0
			while True:

				if i >= len(star_x) or i == 50: break

				# img_star_rotated = rotate(img, a)
				img_star_rotated = img

				# setting global variables for trail fitting
				# img_rot  = img_star_rotated 
				img_rot = img
				centroid = star_x[i], star_y[i]
				# centroid = point_rotation( centroid[0] , centroid[1] , a , img , img_star_rotated )

				str_p0       = np.array([3, l, a, np.mean(sky_row_avg), centroid[0], centroid[1]])
				param_bounds = ([1, l/2, -180, 0, 0, 0], [10, l*5, 180, 2e3, img_star_rotated.shape[1], img_star_rotated.shape[0] ])
				
				try:
					str_param, star_param_cov = curve_fit(trail_model_2d, img_star_rotated, img_star_rotated.flatten(), p0=str_p0)
				except Exception as e:
					print(e , f' LOL star fit failed , skipping trail number {i} for filname : {f}  ')
					failed_log.append(str_p0)
					continue
					
				s, L, A, b, x_0, y_0      = str_param[0], str_param[1], str_param[2], str_param[3], str_param[4], str_param[5]

				residual = np.sum(( trail_model_2d(0, *str_param) - img_star_rotated.flatten() ) ** 2 ) ** .5

				print('star parameters: '     , str_param)
				print('param uncertainties:, ', np.sqrt(np.diag(star_param_cov)))

				# capturing that global variable after the trail fit has converged
				str_flux = flux

				img_star_rotated = rotate(img, A)
				
				# keeping it rotated to star's reference, so don't actually need to go back to asteroid 
				x_0, y_0 = point_rotation( x_0 , y_0 , A , img , img_star_rotated )

				star_trail_start = np.array([x_0, y_0 - L/2 ])
				star_trail_end   = np.array([x_0, y_0 + L/2 ])

				fwhm = s * 2.355
				st_height_correction = int(ast_height_correction * L/ast_trail_length ) - 1
				# st_height_correction = - int(fwhm/2) - 1

				# str_minus_sky, sigma_row_star, str_sky_avg = take_lightcurve(img_star_rotated, star_trail_start, star_trail_end, binning=len(obj_minus_sky), fwhm=fwhm, height_correction=st_height_correction, display=False, err=True, gain=gain, rd_noise=rd_noise)
				str_minus_sky, sigma_row_star, str_sky_avg = take_lightcurve(img_star_rotated, star_trail_start, star_trail_end, fwhm=fwhm, display=False, err=True, gain=gain, rd_noise=rd_noise, height_correction=st_height_correction, binning=len(obj_minus_sky))

				norm = np.median(str_minus_sky)

				row_flux    .append(str_minus_sky /norm)
				row_errs    .append(sigma_row_star/norm)
				trail_starts.append(star_trail_start)
				trail_ends  .append(star_trail_end  )
				residuals   .append(residual)
				stars       .append(np.hstack((str_param, a, flux)))

				print(' ')
				i+=1
				
			row_flux = np.array(row_flux)
			row_errs = np.array(row_errs)

			stars 		 = np.array(stars)
			residuals    = np.array(residuals)
			trail_starts = np.array(trail_starts)
			trail_ends   = np.array(trail_ends)

			print('initially, ', stars.shape[0])

			s_std        = np.std(stars[:,0])
			length_std   = np.std(stars[:,1])
			angle_std    = np.std(stars[:,2])

			s_mean  	 = np.mean(stars[:,0])
			length_mean  = np.mean(stars[:,1])
			angle_mean   = np.mean(stars[:,2])

			# throwing away outliers, ig. 
			# TODO: fit more stars and increase threshold? 
			threshold = 2 # sigmas

			star_filter  = np.where( (stars[:,0]<=s_mean+threshold*s_std) & (stars[:,0]>=s_mean-threshold*s_std) & (stars[:,1]<=length_mean+threshold*length_std) & (stars[:,1]>=length_mean-threshold*length_std) & (stars[:,2]<=angle_mean+threshold*angle_std) & (stars[:,2]>=angle_mean-threshold*angle_std) )
			stars        = stars       [star_filter]
			trail_starts = trail_starts[star_filter]
			trail_ends   = trail_ends  [star_filter]
			residuals    = residuals   [star_filter]
			# total_flux   = total_flux  [star_filter]

			row_flux = row_flux[star_filter]
			print('filtering: ', stars.shape[0])

			# sorting by residuals from biiiig fit
			res_filter   = np.argsort(residuals)
			residuals    = residuals   [res_filter]
			stars        = stars       [res_filter]
			trail_starts = trail_starts[res_filter]
			trail_ends   = trail_ends  [res_filter]
			# total_flux   = total_flux  [res_filter]

			row_flux = row_flux[res_filter][:10]

			print('row_sums_smooth shape: ', row_flux.shape)

			row_avgs = np.nanmean(row_flux, axis=0)
		
			# norm = np.nanmean(row_avgs)
			# row_avgs/=norm

			#sky_corrected_lightcurve = obj_minus_sky[ast_start:ast_end] / row_avgs_smooth # this is the actual sky correction 
			sky_corrected_lightcurve = obj_minus_sky / row_avgs

			x = np.linspace(start_time, start_time + exp_time/(60*60*24), len(sky_corrected_lightcurve))

			directory_name = d.split('/')[1]

			if write_output == 'True':
				np.savetxt(f'{f[:-4]}_params.txt'    , stars )
				np.savetxt(f'{f[:-4]}_lightcurve.txt', np.array([ x , sky_corrected_lightcurve , sigma_row ]).T )

			print()
			file.close()
			if True: break

			# ax[0].legend()


		# plt.show()
	# output.close()

