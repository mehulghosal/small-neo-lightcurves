import pylab as plt
import sys
sys.path.insert(0, '/Users/bolin/NEO/Follow_up/APO_observing')
from apo_observing_functions import *
import matplotlib.pyplot as plt
#matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!afsfsdfasd
import string
import os
import sys
import numpy as np
import random
from decimal import *
import warnings
import re
#import pylab
import math
from matplotlib.widgets import Slider
from bisect import bisect_left
from bisect import bisect_right
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib.ticker import FuncFormatter
from scipy import stats
from sklearn.neighbors import KernelDensity
from statsmodels.nonparametric.kde import KDEUnivariate
from statsmodels.nonparametric.kernel_density import KDEMultivariate
import pyslalib.slalib as sla
import subprocess
import string
from math import log10, exp, sin, tan, sqrt, radians
import pyslalib.slalib as sla
import itertools
from astropy.time import Time
import matplotlib.pyplot as plt
plt.ion()
from scipy.interpolate import Akima1DInterpolator
from scipy.interpolate import interp1d
from apo_observing_functions import *
import zscale as z
from scipy.signal import savgol_filter
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

sed_paretheses =  'awk -F "[()]" \'{ for (i=2; i<NF; i+=2) print $i }\''

hours_to_deg = 15.0
years_to_months  = 12.
years_to_days = 365.
days_to_hours = 24.
hours_to_minutes = 60.
minutes_to_seconds = 60.
deg_to_arcsec = 3600.
djcal_precision = 5
astrores_seconds_rounding_format_RA = 2
astrores_seconds_rounding_format_Dec = 1
number_array_entries_for_hms = 3
number_coordinates = 2 #1 for RA and 1 for Dec
pipe = '|'
colon = ':'
epoch_cfht = '2000.0'
epoch_2000 = 2000.0
j_2000_mjd = '51544.5'
one = '1'
plus = '+'
minus = '-'
zero = '0'
blank = ''
space = ' '
dot = '.'
underscore = '_'
newline = '\n'
object_name_max_length_cfht = 20
RA_s_length_decimal = 3#string length for astrores formating
Jy_to_mJy = 1000.0
numbered_asteroid_string_max_length = 5
UT_to_HST = (10.0 / days_to_hours)
G = 6.67408e-11 #m^3 kg e^-1 s^-2
mass_sun_kg = 1.989e30
au_to_meters = 149597870700.0
seconds_to_days = 1.0 / (24.*3600.)

visir_constant = 36.0

#plotting variables
fontsize_standard = 14
cbar_label_pad = 17
plt.rcParams['axes.linewidth'] = 1
plt.rcParams['lines.linewidth'] = 1
plt.rcParams['patch.linewidth'] = 1
plt.rcParams['contour.negative_linestyle'] = 'solid'

font = {'weight' : 'bold',
        'size'   : fontsize_standard}
plt.rc('font', **font)
plt.rc('text', usetex=True)



#parula map
from matplotlib.colors import LinearSegmentedColormap

cm_data = [[0.2081, 0.1663, 0.5292], [0.2116238095, 0.1897809524, 0.5776761905],
 [0.212252381, 0.2137714286, 0.6269714286], [0.2081, 0.2386, 0.6770857143],
 [0.1959047619, 0.2644571429, 0.7279], [0.1707285714, 0.2919380952,
  0.779247619], [0.1252714286, 0.3242428571, 0.8302714286],
 [0.0591333333, 0.3598333333, 0.8683333333], [0.0116952381, 0.3875095238,
  0.8819571429], [0.0059571429, 0.4086142857, 0.8828428571],
 [0.0165142857, 0.4266, 0.8786333333], [0.032852381, 0.4430428571,
  0.8719571429], [0.0498142857, 0.4585714286, 0.8640571429],
 [0.0629333333, 0.4736904762, 0.8554380952], [0.0722666667, 0.4886666667,
  0.8467], [0.0779428571, 0.5039857143, 0.8383714286],
 [0.079347619, 0.5200238095, 0.8311809524], [0.0749428571, 0.5375428571,
  0.8262714286], [0.0640571429, 0.5569857143, 0.8239571429],
 [0.0487714286, 0.5772238095, 0.8228285714], [0.0343428571, 0.5965809524,
  0.819852381], [0.0265, 0.6137, 0.8135], [0.0238904762, 0.6286619048,
  0.8037619048], [0.0230904762, 0.6417857143, 0.7912666667],
 [0.0227714286, 0.6534857143, 0.7767571429], [0.0266619048, 0.6641952381,
  0.7607190476], [0.0383714286, 0.6742714286, 0.743552381],
 [0.0589714286, 0.6837571429, 0.7253857143],
 [0.0843, 0.6928333333, 0.7061666667], [0.1132952381, 0.7015, 0.6858571429],
 [0.1452714286, 0.7097571429, 0.6646285714], [0.1801333333, 0.7176571429,
  0.6424333333], [0.2178285714, 0.7250428571, 0.6192619048],
 [0.2586428571, 0.7317142857, 0.5954285714], [0.3021714286, 0.7376047619,
  0.5711857143], [0.3481666667, 0.7424333333, 0.5472666667],
 [0.3952571429, 0.7459, 0.5244428571], [0.4420095238, 0.7480809524,
  0.5033142857], [0.4871238095, 0.7490619048, 0.4839761905],
 [0.5300285714, 0.7491142857, 0.4661142857], [0.5708571429, 0.7485190476,
  0.4493904762], [0.609852381, 0.7473142857, 0.4336857143],
 [0.6473, 0.7456, 0.4188], [0.6834190476, 0.7434761905, 0.4044333333],
 [0.7184095238, 0.7411333333, 0.3904761905],
 [0.7524857143, 0.7384, 0.3768142857], [0.7858428571, 0.7355666667,
  0.3632714286], [0.8185047619, 0.7327333333, 0.3497904762],
 [0.8506571429, 0.7299, 0.3360285714], [0.8824333333, 0.7274333333, 0.3217],
 [0.9139333333, 0.7257857143, 0.3062761905], [0.9449571429, 0.7261142857,
  0.2886428571], [0.9738952381, 0.7313952381, 0.266647619],
 [0.9937714286, 0.7454571429, 0.240347619], [0.9990428571, 0.7653142857,
  0.2164142857], [0.9955333333, 0.7860571429, 0.196652381],
 [0.988, 0.8066, 0.1793666667], [0.9788571429, 0.8271428571, 0.1633142857],
 [0.9697, 0.8481380952, 0.147452381], [0.9625857143, 0.8705142857, 0.1309],
 [0.9588714286, 0.8949, 0.1132428571], [0.9598238095, 0.9218333333,
  0.0948380952], [0.9661, 0.9514428571, 0.0755333333],
 [0.9763, 0.9831, 0.0538]]

parula_map = LinearSegmentedColormap.from_list('parula', cm_data)

paperheight = 10
paperwidth = 13
margin = 1.0

fontsize_standard = 28

#lightcurve and PDM

#phase dispersion minimalization with detrended data

#test for mehul
from PyAstronomy.pyTiming import pyPDM
#lightcurve from GE1lc_errs.png plot in Mehul's July 26 email
DCTAPO_date_MJD, DCTAPO_mag, DCTAPO_mag_unc = np.loadtxt('/Users/bolin/NEO/Follow_up/APO_observing/mehul_test/GE_1_lightcurves/1917070o13_lightcurve.txt').T

#test plot
fig = plt.figure(figsize=((paperwidth*1.15) - 2 * margin, (paperheight*1.15) - 2 * margin))
line1 = plt.errorbar(x=DCTAPO_date_MJD, y=DCTAPO_mag, yerr=DCTAPO_mag_unc, fmt='s', markerfacecolor = 'blue', markeredgecolor = 'black', ecolor='black', capthick=2,markersize=7,capsize=3)
plt.show()

# Create artificial data with frequency = 3,
# period = 1/3
x = DCTAPO_date_MJD
y = DCTAPO_mag

# Get a ``scanner'', which defines the frequency interval to be checked.
# Alternatively, also periods could be used instead of frequency.
S = pyPDM.Scanner(minVal=0.0001, maxVal=.2, dVal=0.00001, mode="period")
P = pyPDM.PyPDM(x, y)
f1, t1 = P.pdmEquiBinCover(5, 3, S)

num_peak = 1.0
best_period = f1[np.where(t1== np.min(t1))]
print best_period
best_frequency = 1/best_period

fig = plt.figure(figsize=((paperwidth*1.15) - 2 * margin, (paperheight*1.15) - 2 * margin))
ax1 = fig.add_subplot(111)
ax1.semilogx(f1*24.*3600.,(t1),color='grey')
best_frequency1 = f1[np.where(t1== np.min(t1))]
plt.ylabel(r'$\mathrm{\Theta}$',fontsize=20)
plt.xlabel(r'$\mathrm{Lightcurve \; period \; (s)}$',fontsize=20)
ax1.axvline(best_period*24*3600., color='blue', linestyle='-',label =r'$\mathrm{Period:\; '+ str(np.round(best_period[0]*24*3600.,2))+'\;  s}$',linewidth=2.2)
ax1.legend(loc='lower right',prop={'size':16})
ax1.set_xlim(f1.min()*24.*3600.,f1.max()*24.*3600.)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.show()
plt.savefig('GE1lc_test_for_mehul_pdm.png')

#plot the lightcurve with the fitted lightcurve
time_unit_factor=24.*3600.
fig = plt.figure(figsize=((paperwidth*1.15) - 2 * margin, (paperheight*1.15) - 2 * margin))
ax1 = fig.add_subplot(111)
line1 = plt.errorbar(x=(DCTAPO_date_MJD-DCTAPO_date_MJD[0])*time_unit_factor, y=DCTAPO_mag, yerr=DCTAPO_mag_unc, fmt='s', markerfacecolor = 'blue', markeredgecolor = 'black', ecolor='black', capthick=2,markersize=7,capsize=3)
t = np.linspace(np.min(DCTAPO_date_MJD), np.max(DCTAPO_date_MJD),1000.)
Amplitude = (np.max(DCTAPO_mag) - np.min(DCTAPO_mag))
set_phase = np.pi*3.6
y = (Amplitude * 0.5* np.sin(2 * np.pi * t*num_peak*(24*3600)/((1.0/(best_frequency)) *24.*3600.) + set_phase))
y_offset = np.median(DCTAPO_mag)
ax1.plot((t-DCTAPO_date_MJD[0])*time_unit_factor, y+y_offset, color='black')
#plt.xlim(np.min(DCTAPO_date_MJD)-DCTAPO_date_MJD[0], np.max(DCTAPO_date_MJD)-DCTAPO_date_MJD[0])
MyFormatter2 = FuncFormatter(latex_ticks1)
ax1.axes.yaxis.set_major_formatter(MyFormatter2)
MyFormatter = FuncFormatter(latex_ticks1)
ax1.axes.xaxis.set_major_formatter(MyFormatter)
plt.xlabel(r'$\mathrm{Time \, (seconds \, after \, MJD \,' + str(np.round(t[0]-2400000.5,5)) + ')}$',fontsize=20)
plt.ylabel(r'$\mathrm{Flux}$',fontsize=20)
plt.gca().invert_yaxis()
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.show()
plt.savefig('GE1lc_lightcurve_for_mehul.png')

#lightcurve from EV84lc_errs.png plot in Mehul's July 26 email
DCTAPO_date_MJD, DCTAPO_mag, DCTAPO_mag_unc = np.loadtxt('/Users/bolin/NEO/Follow_up/APO_observing/mehul_test/EV84_lightcurves/1909610o22_lightcurve.txt').T
#test plot
fig = plt.figure(figsize=((paperwidth*1.15) - 2 * margin, (paperheight*1.15) - 2 * margin))
line1 = plt.errorbar(x=DCTAPO_date_MJD, y=DCTAPO_mag, yerr=DCTAPO_mag_unc, fmt='s', markerfacecolor = 'blue', markeredgecolor = 'black', ecolor='black', capthick=2,markersize=7,capsize=3)

import numpy
import matplotlib.pylab as plt
from PyAstronomy.pyTiming import pyPDM

# Create artificial data with frequency = 3,
# period = 1/3
x = DCTAPO_date_MJD
y = DCTAPO_mag

# Get a ``scanner'', which defines the frequency interval to be checked.
# Alternatively, also periods could be used instead of frequency.
S = pyPDM.Scanner(minVal=0.0001, maxVal=.2, dVal=0.00001, mode="period")
P = pyPDM.PyPDM(x, y)
f1, t1 = P.pdmEquiBinCover(5, 3, S)

num_peak = 1.0
best_period = f1[np.where(t1== np.min(t1))]
print best_period
best_frequency = 1/best_period

fig = plt.figure(figsize=((paperwidth*1.15) - 2 * margin, (paperheight*1.15) - 2 * margin))
ax1 = fig.add_subplot(111)
ax1.semilogx(f1*24.*3600.,(t1),color='grey')
best_frequency1 = f1[np.where(t1== np.min(t1))]
plt.ylabel(r'$\mathrm{\Theta}$',fontsize=20)
plt.xlabel(r'$\mathrm{Lightcurve \; period \; (s)}$',fontsize=20)
ax1.axvline(best_period*24*3600., color='blue', linestyle='-',label =r'$\mathrm{Period:\; '+ str(np.round(best_period[0]*24*3600.,2))+'\;  s}$',linewidth=2.2)
ax1.legend(loc='lower right',prop={'size':16})
ax1.set_xlim(f1.min()*24.*3600.,f1.max()*24.*3600.)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.savefig('EV84lc_test_for_mehul_pdm.png')

#plot the lightcurve with the fitted lightcurve
time_unit_factor=24.*3600.
fig = plt.figure(figsize=((paperwidth*1.15) - 2 * margin, (paperheight*1.15) - 2 * margin))
ax1 = fig.add_subplot(111)
line1 = plt.errorbar(x=(DCTAPO_date_MJD-DCTAPO_date_MJD[0])*time_unit_factor, y=DCTAPO_mag, yerr=DCTAPO_mag_unc, fmt='s', markerfacecolor = 'blue', markeredgecolor = 'black', ecolor='black', capthick=2,markersize=7,capsize=3)
t = np.linspace(np.min(DCTAPO_date_MJD), np.max(DCTAPO_date_MJD),1000.)
Amplitude = (np.max(DCTAPO_mag) - np.min(DCTAPO_mag))
set_phase = np.pi*2.2
y = (Amplitude * 0.5* np.sin(2 * np.pi * t*num_peak*(24*3600)/((1.0/(best_frequency)) *24.*3600.) + set_phase))
y_offset = np.median(DCTAPO_mag)
ax1.plot((t-DCTAPO_date_MJD[0])*time_unit_factor, y+y_offset, color='black')
#plt.xlim(np.min(DCTAPO_date_MJD)-DCTAPO_date_MJD[0], np.max(DCTAPO_date_MJD)-DCTAPO_date_MJD[0])
MyFormatter2 = FuncFormatter(latex_ticks1)
ax1.axes.yaxis.set_major_formatter(MyFormatter2)
MyFormatter = FuncFormatter(latex_ticks1)
ax1.axes.xaxis.set_major_formatter(MyFormatter)
plt.xlabel(r'$\mathrm{Time \, (seconds \, after \, MJD \,' + str(np.round(t[0]-2400000.5,5)) + ')}$',fontsize=20)
plt.ylabel(r'$\mathrm{Flux}$',fontsize=20)
plt.gca().invert_yaxis()
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.savefig('EV84lc_lightcurve_for_mehul.png')

#generate fake lightcurve data

N = 100
period_s = 25.
f = 1/period_s
A = 0.5
sig = 0.25
mag_err_array = np.ones(N) * sig

time = np.arange(float(N))
flux = A * np.sin(2. * np.pi*time*f)
# Adding the noise
flux += np.random.normal(0, sig, time.size)

#plot fake data
fig = plt.figure(figsize=((paperwidth*1.15) - 2 * margin, (paperheight*1.15) - 2 * margin))
line1 = plt.errorbar(x=time, y=flux, yerr=mag_err_array, fmt='s', markerfacecolor = 'blue', markeredgecolor = 'black', ecolor='black', capthick=2,markersize=7,capsize=3)
plt.xlabel('Time (s)')

#periodogram test fake data
from PyAstronomy.pyTiming import pyPeriod

'''
https://pyastronomy.readthedocs.io/en/latest/pyTimingDoc/pyPeriodDoc/gls.html
'''

clp = pyPeriod.Gls((time, flux, mag_err_array))

#false positive levels
# Define FAP levels of 5% and 0.5%
fapLevels = np.array([0.05, 0.005])
# Obtain the associated power thresholds
plevels = clp.powerLevel(fapLevels)

fig = plt.figure(figsize=((paperwidth*1.15) - 2 * margin, (paperheight*1.15) - 2 * margin))
ax1 = fig.add_subplot(111)
num_peak=1.0
periods = 1/clp.freq
best_period = periods[np.argmax(clp.power)]
line_width = 2.5
mult = 1.2
paperheight = 6.5*1.15
paperwidth = 9.5*1.15
margin = 0.5
ax1.semilogx(1/clp.freq,clp.power,color='grey')
ax1.axvline(best_period, color='blue', linestyle='-',label =r'$\mathrm{Period:\, '+ str(np.round((best_period),2))+'\;  s}$',linewidth=2.2)
plt.ylabel(r'$\mathrm{Power}$',fontsize=20)
plt.xlabel(r'$\mathrm{Lightcurve \; period \; (s)}$',fontsize=20)
ax1.set_xlim(periods.min(),periods.max())
plt.plot([min(periods), max(periods)], [plevels[0]]*2, '-.', c='black',lw=3,label =r'$2\sigma$')
plt.plot([min(periods), max(periods)], [plevels[1]]*2, ':', c='black',lw=3,label =r'$3\sigma$')
ax1.legend(loc='upper left',prop={'size':16})
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

#plot fake data and best sine wave
bestSine = clp.sinmod(time)

fig = plt.figure(figsize=((paperwidth*1.15) - 2 * margin, (paperheight*1.15) - 2 * margin))
line1 = plt.errorbar(x=time, y=flux, yerr=mag_err_array, fmt='s', markerfacecolor = 'blue', markeredgecolor = 'black', ecolor='black', capthick=2,markersize=7,capsize=3)
plt.plot(time, bestSine, "k--",lw=5)
plt.xlabel('Time (s)')

#folded data
fig = plt.figure(figsize=((paperwidth*1.15) - 2 * margin, (paperheight*1.15) - 2 * margin))
line1 = plt.errorbar(x=time/best_period-np.floor(time/best_period), y=flux, yerr=mag_err_array, fmt='s', markerfacecolor = 'blue', markeredgecolor = 'black', ecolor='black', capthick=2,markersize=7,capsize=3)
plt.xlabel('Phase')

#generate synthetic data based on Mehul's lightcurves

#test for mehul
#lightcurve from GE1lc_errs.png plot in Mehul's July 26 email
DCTAPO_date_MJD, DCTAPO_mag, DCTAPO_mag_unc = np.loadtxt('/Users/bolin/NEO/Follow_up/APO_observing/mehul_test/GE_1_lightcurves/1917070o13_lightcurve.txt').T

DCTAPO_date_seconds_from_start_s = (DCTAPO_date_MJD - DCTAPO_date_MJD[0]) * 3600. * 24.
DCTAPO_mag_norm = DCTAPO_mag-np.median(DCTAPO_mag)

#test plot
fig = plt.figure(figsize=((paperwidth*1.15) - 2 * margin, (paperheight*1.15) - 2 * margin))
line1 = plt.errorbar(x=DCTAPO_date_seconds_from_start_s, y=DCTAPO_mag_norm, yerr=DCTAPO_mag_unc, fmt='s', markerfacecolor = 'blue', markeredgecolor = 'black', ecolor='black', capthick=2,markersize=7,capsize=3)
plt.show()

period_s = 17.174283947263444
f = 1/period_s

A = np.mean([np.abs(np.max(DCTAPO_mag_norm)),np.abs(np.min(DCTAPO_mag_norm))])
sig = np.median(DCTAPO_mag_unc)

time = DCTAPO_date_seconds_from_start_s
mag_err_array = np.ones(len(time)) * sig

flux = A * np.sin(2. * np.pi*time*f)
# Adding the noise
flux += np.random.normal(0, sig, time.size)

#plot fake data
fig = plt.figure(figsize=((paperwidth*1.15) - 2 * margin, (paperheight*1.15) - 2 * margin))
line1 = plt.errorbar(x=time, y=flux, yerr=mag_err_array, fmt='s', markerfacecolor = 'blue', markeredgecolor = 'black', ecolor='black', capthick=2,markersize=7,capsize=3)
plt.xlabel('Time (s)')

#compare real data and fake data

#periodogram real data
from PyAstronomy.pyTiming import pyPeriod

clp = pyPeriod.Gls((DCTAPO_date_seconds_from_start_s, DCTAPO_mag_norm,DCTAPO_mag_unc))

#false positive levels
# Define FAP levels of 5 and 0.5%
fapLevels = np.array([0.05, 0.005])
# Obtain the associated power thresholds
plevels = clp.powerLevel(fapLevels)

fig = plt.figure(figsize=((paperwidth*1.15) - 2 * margin, (paperheight*1.15) - 2 * margin))
ax1 = fig.add_subplot(111)
num_peak=1.0
periods = 1/clp.freq
best_period = periods[np.argmax(clp.power)]
line_width = 2.5
mult = 1.2
paperheight = 6.5*1.15
paperwidth = 9.5*1.15
margin = 0.5
ax1.semilogx(1/clp.freq,clp.power,color='grey')
ax1.axvline(best_period, color='blue', linestyle='-',label =r'$\mathrm{Period:\, '+ str(np.round((best_period),2))+'\;  s}$',linewidth=2.2)
plt.ylabel(r'$\mathrm{Power}$',fontsize=20)
plt.xlabel(r'$\mathrm{Lightcurve \; period \; (s)}$',fontsize=20)
ax1.set_xlim(periods.min(),periods.max())
plt.plot([min(periods), max(periods)], [plevels[0]]*2, '-.', c='black',lw=3,label =r'$2\sigma$')
plt.plot([min(periods), max(periods)], [plevels[1]]*2, ':', c='black',lw=3,label =r'$3\sigma$')
ax1.legend(loc='upper left',prop={'size':16})
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

#plot fake data and best sine wave
bestSine = clp.sinmod(DCTAPO_date_seconds_from_start_s)

fig = plt.figure(figsize=((paperwidth*1.15) - 2 * margin, (paperheight*1.15) - 2 * margin))
line1 = plt.errorbar(x=DCTAPO_date_seconds_from_start_s, y=DCTAPO_mag_norm, yerr=DCTAPO_mag_unc, fmt='s', markerfacecolor = 'blue', markeredgecolor = 'black', ecolor='black', capthick=2,markersize=7,capsize=3)
plt.plot(DCTAPO_date_seconds_from_start_s, bestSine, "k--",lw=5)
plt.xlabel('Time (s)')

#folded data
fig = plt.figure(figsize=((paperwidth*1.15) - 2 * margin, (paperheight*1.15) - 2 * margin))
line1 = plt.errorbar(x=DCTAPO_date_seconds_from_start_s/best_period-np.floor(DCTAPO_date_seconds_from_start_s/best_period), y=DCTAPO_mag_norm, yerr=DCTAPO_mag_unc, fmt='s', markerfacecolor = 'blue', markeredgecolor = 'black', ecolor='black', capthick=2,markersize=7,capsize=3)
line1 = plt.errorbar(x=time/best_period-np.floor(time/best_period), y=flux, yerr=mag_err_array, fmt='s', markerfacecolor = 'red', markeredgecolor = 'black', ecolor='black', capthick=2,markersize=7,capsize=3)
plt.xlabel('Phase')
#amplitudes look ok



#lightcurve from GE1lc_errs.png plot in Mehul's July 26 email
DCTAPO_date_MJD, DCTAPO_mag, DCTAPO_mag_unc = np.loadtxt('/Users/bolin/NEO/Follow_up/APO_observing/test_lightcurve.txt').T

DCTAPO_date_seconds_from_start_s = DCTAPO_date_MJD
DCTAPO_mag_norm = DCTAPO_mag-np.median(DCTAPO_mag)

#test plot
fig = plt.figure(figsize=((paperwidth*1.15) - 2 * margin, (paperheight*1.15) - 2 * margin))
line1 = plt.errorbar(x=DCTAPO_date_seconds_from_start_s, y=DCTAPO_mag_norm, yerr=DCTAPO_mag_unc, fmt='s', markerfacecolor = 'blue', markeredgecolor = 'black', ecolor='black', capthick=2,markersize=7,capsize=3)
plt.show()

#heat map showing efficiency as function of period vs ampltidue
