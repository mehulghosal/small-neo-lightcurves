import matplotlib.pyplot as plt
import numpy as np
from magic_star import take_lightcurve
'''
DEBUGGING UTILITIES -- NOT FOR USE IN __MAIN__

	these are to be used in ipython test scipts with only three stars in row_sums_smooth

row_sums_smooth  [list or array-like: shape=(3, length)]: every row is a lightcurve, columns are row number/time whatever

returns fig, ax

'''

def plot_st_lcs(row_sums_smooth, display=False):
	fig, ax = plt.subplots(3, 1)
	for i in range(3):
		ax[i].plot(row_sums_smooth[i])
	if display: fig.show()
	return fig, ax


def plot_unbinned(img, trail_starts, trail_ends, display=False, height_correction=0):
	fig, ax  = plt.subplots(3, 1)
	for i in range(3):
		st_lc = take_lightcurve(img, trail_starts[i], trail_ends[i], height_correction=height_correction)[0]
		ax[i].plot(st_lc)
	if display: fig.show()
	return fig, ax

"""
returns cropped in rectangle displaying the trail from parameters
	usage:  > plt.imshow(display_streak(img_star_rotated, *stars[0]))

PARAMETERS
-----------
img   : array
	2d numpy array representing image we extract trail from 
s 	  : float
	Gaussian spread of trail
L 	  : float
	length of trail
a 	  : float
	angle of trail counterclockwise from positive x axis
b 	  : float
	scalar representation of background contribution 
x_0   : float
	CCD pixel column coordinate of centroid
y_0   : float
	CCD pixel row coordinate of centroid
width : int
	width in FWHM of returned view
height: int
	height in L of returned view

RETURNS
---------
obj_rect : array
	2 * width * s * 2.355 columns wide
	L * height rows tall
"""
def display_streak(img, s, L, a, b, x_0, y_0, width=1, height=1):
	obj_width  = width * s * 2.355
	obj_height = L * height
	obj_rect   = img[int(y_0 - obj_height/2 + .5) : int(y_0 + obj_height/2 + .5), int(x_0 - obj_width + .5) : int(x_0 + obj_width + .5) ]
	return obj_rect