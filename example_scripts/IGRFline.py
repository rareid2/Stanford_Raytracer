"""

function to generate magnetic field line from IGRF13 coefficients
using pyIGRF package and function adapted from MATLAB igrfline
author: Riley R
Apr. 2020

INPUT: coordinates as lat_start (deg N), lon_start (deg E), alt_start (m from center of Earth)
of starting point of magnetic field line

OUTPUT: lat, lon, alt list in deg N, deg E, m from center of Earth

"""
import numpy as np
import matplotlib.pyplot as plt
import pyIGRF
from raytracer_settings import *


# --------------------------- START FUNCTION ---------------------------

def IGRFline(lat_start, lon_start, alt_start):

    # seems to be good conditions to finish the field line
    nsteps = 1000
    distance = (alt_start/R_E)*R_E*4/1e3
    steplen = distance/np.abs(nsteps)

    # conevrt ray_datenum to year and fraction of the year - needed for pyIGRF
    sec2yr = 31536000
    day2sec = 86400

    sec_in_yr = day2sec*int(days_in_the_year) + milliseconds_day/(1e3)
    year_frac = ray_datenum.year + (sec_in_yr/sec2yr)

    # convert to km from surface of Earth for pyIGRF
    alt_start = (alt_start - R_E)/1e3

    # initialize
    lat = np.zeros(nsteps)
    lon = np.zeros(nsteps)
    alt = np.zeros(nsteps)
    lat[0] = lat_start
    lon[0] = lon_start
    alt[0] = alt_start

    for index in range(nsteps-1):

        # get magnetic field values
        B_out = pyIGRF.igrf_value(lat[index], lon[index], alt[index], year_frac)

        Bt, Bp, Br = B_out[3:6]

        Br = -Br

        B = np.hypot(Br, np.hypot(Bp, Bt))

        dr = Br / B
        dp = Bp / B
        dt = Bt / B

        # stop if field line hits Earth
        if alt[index] + steplen * dr < 0:
            print('hit earth')
            break
        # else update to next step
        else:
            alt[index + 1] = alt[index] + steplen * dr
            lon[index + 1] = lon[index] + steplen * dt / alt[index]
            lat[index + 1] = lat[index] + steplen * dp / (alt[index] * np.cos(np.deg2rad(lon[index])))

    # remove extra zeroes
    lat = [i for i in lat if i != 0]
    lon = [i for i in lon if i != 0]
    alt = [i for i in alt if i != 0]

    alt_out = [i * 1e3 + R_E for i in alt] # convert back to ECEF in m
    alt = alt_out

    return lat,lon,alt

# --------------------------- END FUNCTION ---------------------------

