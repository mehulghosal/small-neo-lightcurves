#!/usr/bin/env python

#for small asteroid lightcurve project, generates synthetic lightcurve data based off of real streak asteroid lightcurve data

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
python synthetic_lightcurve_generator_v1.py @if /Users/bolin/NEO/Follow_up/APO_observing/mehul_test/GE_1_lightcurves/1917070o13_lightcurve.txt @per 17.2 @amp 0.2

change log:

2022-Sep-13: changed amplitude unites from flux to mag. The out put is still in time_s, flux, flux unc

'''

parser = argparse.ArgumentParser(prefix_chars='@')
parser.add_argument("@if", "@@in_file", help="read in mehul's data in format of mjd, flux, flux unc")
parser.add_argument("@per", "@@period", help="period in seconds")
parser.add_argument("@amp", "@@amplitude", help="amplitude of lightcurve in mag units")
args = parser.parse_args()


in_file = str(args.in_file)
period_s = float(args.period)
amplitude_mag = float(args.amplitude)


DCTAPO_date_MJD, DCTAPO_flux, DCTAPO_flux_unc = np.loadtxt(in_file).T

DCTAPO_date_seconds_from_start_s = (DCTAPO_date_MJD - DCTAPO_date_MJD[0]) * 3600. * 24.
DCTAPO_flux_norm = DCTAPO_flux-np.median(DCTAPO_flux)


period_s = period_s
f = 1/period_s

median_flux = np.median(DCTAPO_flux)

flux_ratio = 10**(amplitude_mag/-2.5) #for lower peak

A = ((1./flux_ratio) * median_flux) - median_flux
sig = np.median(DCTAPO_flux_unc)

time_s = DCTAPO_date_seconds_from_start_s
mag_err_array = np.ones(len(time_s)) * sig
#generate random phase
phase = np.random.uniform(-2*np.pi,2*np.pi)
#generate flux
flux = A * np.sin(2. * np.pi*(time_s-phase)*f)
# Adding the noise
flux += np.random.normal(0, sig, time_s.size)

for i in range(0,len(flux)):
    print '%10.2f %10.2f %10.1f'%(time_s[i], flux[i], sig)
