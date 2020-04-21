"""

function to generate magnetic field line from IGRF13 coefficients
using pyIGRF package and function adapted from MATLAB igrfline
author: Riley R
Apr. 2020

INPUT: coordinates as lat (deg N), lon (deg E), alt (m from center of Earth)
of starting point of magnetic field line

OUTPUT: unit vector of direction of local magnetic field line in lla

"""

# only getting top half of field line - how to get other half? thought
# it would be done by changing phi direction ...

import numpy as np
import matplotlib.pyplot as plt
import pyIGRF
from raytracer_settings import *

# --------------------------- START FUNCTION ---------------------------

def IGRFdirection(lat, lon, alt):

    # convert to km from surface of Earth for pyIGRF
    alt = (alt - R_E) / 1e3

    # conevrt ray_datenum to year and fraction of the year - needed for pyIGRF
    sec2yr = 31536000
    day2sec = 86400

    sec_in_yr = day2sec*int(days_in_the_year) + milliseconds_day/(1e3)
    year_frac = ray_datenum.year + (sec_in_yr/sec2yr)

    # get magnetic field values
    B_out = pyIGRF.igrf_value(lat, lon, alt, year_frac)

    Bt, Bp, Br = B_out[3:6]
    Br = -Br

    return Br, Bp, Bt

# --------------------------- END FUNCTION ---------------------------