#!/usr/bin/env python

#for small asteroid lightcurve project, estimates the efficiency of period determination of synthetic data with uncertainties compatible with real asteroids streak data determined over a range of lightcurve periods and amplitudes

#bryce bolin

import numpy as np
import argparse
from PyAstronomy.pyTiming import pyPeriod

'''
period_finding_efficiency_test_v1.py

Estimates period finding efficiency for streak asteroid lightcurve data as a function of period and amplitude 

Python modules required:

argparse
numpy
from PyAstronomy.pyTiming import pyPeriod
https://pyastronomy.readthedocs.io/en/latest/pyTimingDoc/pyPeriodDoc/gls.html


directions for use and sample execution:

python -i period_finding_efficiency_test_v1.py @if /Users/bolin/NEO/Follow_up/APO_observing/mehul_test/GE_1_lightcurves/1917070o13_lightcurve.txt @rn GE1_1917070o13 @per "9 130 16" @amp "0.09 1.21 0.16" @nt 200 @fa 0.0025 @ps .05 @as 0.1

output are two matricies:

One matrix has the efficiency of those that pass the x sigma/y probability test

The second matrix has the efficiency of those that pass the x sigma AND are within z % of the period w% similarity within the amplitude used to generate the synthetic data.

Currently only the second matrix is written and printed out.

Output is in .pyc files and as standard out

change log:

2022-Sep-14: first version producing output.

'''

def string_seperated_to_array_spaces(input_array,data_type_str):
    input_array = input_array.replace(" ",",")
    cols_index_storage = np.array([])
    starting_point = 0
    comma_count = 0
    for i in range(0, len(input_array)):
        test = input_array[starting_point:].find(',')
        #print i, test, starting_point
        if test != -1 and test!=0:
            cols_index_storage = np.append(cols_index_storage,float(input_array[starting_point:input_array[starting_point:].find(',')+starting_point]))
            starting_point = input_array[starting_point:].find(',')+starting_point
        if test == -1:
            cols_index_storage = np.append(cols_index_storage,float(input_array[starting_point:]))
            break
        if test == 0:
            starting_point += 1
            comma_count +=1
    return cols_index_storage.astype(data_type_str)


parser = argparse.ArgumentParser(prefix_chars='@')
parser.add_argument('@if', '@@in_file', help='read in mehuls data in format of mjd, flux, flux unc')
parser.add_argument('@rn', '@@run_name', help='name of test run, e.g., GE1_1917070o13')
parser.add_argument('@per', '@@period', help='period range in format lower bound, upper bound and bin width, e.g., 11 261 20, the bin width must span the a single binsize and be divisible by two such that the lower bound on the first bin is a non-zero integer or float')
parser.add_argument('@amp', '@@amplitude', help='amplitude range in format lower bound, upper bound and bin width, e.g., 0.11 2.21 0.2, the bin width must span the a single binsize and be divisible by two such that the lower bound on the first bin is a non-zero integer or float')
parser.add_argument('@nt', '@@numbertrials', help='number of trials per period, amplitude bin')
parser.add_argument('@fa', '@@falsealarm', help='false alarm probability expressed as a fraction, e.g., 0.005')
parser.add_argument('@ps', '@@periodsimilarity', help='period similarity to count as finding approximately the right period, can be 0.05')
parser.add_argument('@as', '@@amplitudesimilarity', help='amplitude similarity to count as finding approximately the right amplitude, can be 0.1')

args = parser.parse_args()

in_file = str(args.in_file)
run_name = str(args.run_name)
period_params_low_high_bin_size = string_seperated_to_array_spaces(args.period,'float')
amplitude_params_low_high_bin_size = string_seperated_to_array_spaces(args.amplitude,'float')
number_of_tests_per_bin = int(args.numbertrials)
false_alarm_probability_threshold = float(args.falsealarm)
period_similarity_threshold = float(args.periodsimilarity)
amplitude_similarity_threshold = float(args.amplitudesimilarity)

DCTAPO_date_MJD, DCTAPO_flux, DCTAPO_flux_unc = np.loadtxt(in_file).T

DCTAPO_date_seconds_from_start_s = (DCTAPO_date_MJD - DCTAPO_date_MJD[0]) * 3600. * 24.

time_s = DCTAPO_date_seconds_from_start_s

median_flux = np.median(DCTAPO_flux)

DCTAPO_flux_norm = DCTAPO_flux-median_flux

sig = np.median(DCTAPO_flux_unc)

flux_err_array = np.ones(len(time_s)) * sig

#number of tests per period and amplitude bin

period_low_s, period_high_s, period_bin_width_s = period_params_low_high_bin_size
amplitude_low_mag, amplitude_high_mag, amplitude_bin_width_mag = amplitude_params_low_high_bin_size

out_put_name = run_name + '_per_'+str(period_low_s)+'_'+str(period_high_s)+'_'+str(period_bin_width_s) +'_amp_' + str(amplitude_low_mag)+'_'+str(amplitude_high_mag)+'_'+str(amplitude_bin_width_mag)+'_num_trials_'+str(number_of_tests_per_bin)+'_false_alarm_prob_'+str(false_alarm_probability_threshold)+'_period_threshold_'+str(period_similarity_threshold)+'_amp_threshold_'+str(amplitude_similarity_threshold)

# ranges for period and ampltiude
#need to generate fake lightcurves with period and amplitude between upper and lower bounds for period and amplitude for each bin, starting at each bin center

#period bins
period_range_s = np.arange(period_low_s, period_high_s, period_bin_width_s)
period_range_low_s = period_range_s - (period_bin_width_s/2)
period_range_high_s = period_range_s + (period_bin_width_s/2)

#amplitude bins
amplitude_range_mag = np.arange(amplitude_low_mag, amplitude_high_mag, amplitude_bin_width_mag)
amplitude_range_low_mag = amplitude_range_mag - (amplitude_bin_width_mag/2)
amplitude_range_high_mag = amplitude_range_mag + (amplitude_bin_width_mag/2)

#matrices for storing the fraction passing 3 sigma test divided by total number of runs and fraction of runs passing 3 sigma test and period within 20% divided by the number of runs passing the 3 sigma test

#3 sigma divided by total
efficiency_period_vs_amplitude_3_sigma_pass = np.zeros(len(period_range_s) * len(amplitude_range_mag)).reshape(len(period_range_s),len(amplitude_range_mag))

#3 sigma + similarity divided by 3 sigma
efficiency_period_vs_amplitude_3_sigma_plus_pass_similar_period_pass = np.zeros(len(period_range_s) * len(amplitude_range_mag)).reshape(len(period_range_s),len(amplitude_range_mag))


#loop over range periods
for i in range(0,len(period_range_s)):
 #loop over range of amplitudes
 for j in range(0,len(amplitude_range_mag)):
  #execute number of tests
  for k in range(0,number_of_tests_per_bin):
   if k % 20 == 0:
    print '# period', i+1, 'of', len(period_range_s), 'amplitude' , j+1, 'of', len(amplitude_range_mag), 'run' , k+1, 'of', number_of_tests_per_bin
   #randomize period
   period_s = np.random.uniform(period_range_low_s[i],period_range_high_s[i])
   freq = 1/period_s
   #randomize amplitude
   amplitude_mag = np.random.uniform(amplitude_range_low_mag[j],amplitude_range_high_mag[j])
   flux_ratio = 10**(amplitude_mag/-2.5) #for lower peak
   A = ((1./flux_ratio) * median_flux) - median_flux
   #randomize phase
   phase = np.random.uniform(-2*np.pi,2*np.pi)
   #generate flux
   flux = (A * np.sin(2. * np.pi*(time_s-phase)*freq))
   #generate the noise
   noise_flux = np.random.normal(0, sig, time_s.size)
   #add the noise
   flux_noise_added = flux + noise_flux
   #do the generalized LS test
   clp = pyPeriod.Gls((time_s, flux_noise_added,flux_err_array))
   #false alarm probability (FAP) levels
   # Define FAP levels at x-sigma/ y%
   false_alarm_prob_level = np.array([false_alarm_probability_threshold])
   # Obtain the associated power thresholds
   false_alarm_power_level = clp.powerLevel(false_alarm_prob_level)
   #determine the best period
   periods_LS = 1/clp.freq
   max_arg = np.argmax(clp.power)
   best_power = clp.power[max_arg]
   best_period = periods_LS[max_arg]
   #false alarm test
   if best_power > false_alarm_power_level: #test if above 3 sigma
    efficiency_period_vs_amplitude_3_sigma_pass[i,j] += 1.0
    period_fraction = np.abs(best_period-period_s)/period_s
    if period_fraction<period_similarity_threshold: #test if above 3 sigma and within 20% of period
     #Add an amplitude check to the criteria as in Masiero et al. 2009
     bestSine = clp.sinmod(time_s)
     fit_sin_amplitude = np.max(bestSine)
     amplitude_fraction = np.abs(fit_sin_amplitude-A)/A
     if amplitude_fraction<amplitude_similarity_threshold:
        efficiency_period_vs_amplitude_3_sigma_plus_pass_similar_period_pass[i,j] += 1.0

efficiency_period_vs_amplitude_3_sigma_pass.dump(out_put_name+'_3_sigma.pyc')
efficiency_period_vs_amplitude_3_sigma_plus_pass_similar_period_pass.dump(out_put_name+'_3_sigma_and_threshold.pyc')
np.savetxt(out_put_name+'_3_sigma.txt', efficiency_period_vs_amplitude_3_sigma_pass)
np.savetxt(out_put_name+'_3_sigma_and_threshold.txt', efficiency_period_vs_amplitude_3_sigma_plus_pass_similar_period_pass)
