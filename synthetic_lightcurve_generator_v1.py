#!/usr/bin/env python

#for small asteroid lightcurve project

#bryce bolin


import numpy as np
import argparse


'''
synthetic_lightcurve_generator_v1.py

Generates synthetic lighcurves based off of the scatter in Mehul's data

Python modules required:

argparse
numpy

directions for use and sample execution:
python synthetic_lightcurve_generator_v1.py @if /Users/bolin/NEO/Follow_up/APO_observing/mehul_test/GE_1_lightcurves/1917070o13_lightcurve.txt @per 17.2 @amp 1000

'''
parser = argparse.ArgumentParser(prefix_chars='@')
parser.add_argument("@if", "@@in_file", help="read in mehul's data in format of mjd, flux, flux unc")
parser.add_argument("@per", "@@period", help="period in seconds")
parser.add_argument("@amp", "@@amplitude", help="amplitude of lightcurve in flux units")
args = parser.parse_args()


in_file = str(args.in_file)
period_s = float(args.period)
amplitude_flux = float(args.amplitude)


DCTAPO_date_MJD, DCTAPO_mag, DCTAPO_mag_unc = np.loadtxt(in_file).T

DCTAPO_date_seconds_from_start_s = (DCTAPO_date_MJD - DCTAPO_date_MJD[0]) * 3600. * 24.
DCTAPO_mag_norm = DCTAPO_mag-np.median(DCTAPO_mag)


period_s = 17.174283947263444
f = 1/period_s

time_difference_round_up_int_s = int(np.ceil(DCTAPO_date_seconds_from_start_s[-1]))
A = np.mean([np.abs(np.max(DCTAPO_mag_norm)),np.abs(np.min(DCTAPO_mag_norm))])
sig = np.median(DCTAPO_mag_unc)

time = DCTAPO_date_seconds_from_start_s
mag_err_array = np.ones(len(time)) * sig

flux = A * np.sin(2. * np.pi*time*f)
# Adding the noise
flux += np.random.normal(0, sig, time.size)

for i in range(0,len(flux)):
    print '%10.2f %10.2f %10.1f'%(time[i], flux[i], sig)
