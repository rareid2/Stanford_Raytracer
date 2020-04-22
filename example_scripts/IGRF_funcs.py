"""

functions to generate magnetic field intensity and line from IGRF13 coefficients
using pyIGRF package and function adapted from MATLAB igrf package and
Dr. Austin Sousa's thesis work

author: Riley R
Apr. 2020

"""

# import packages
import numpy as np
import matplotlib.pyplot as plt
# for IGRF ODE
import pyIGRF
from scipy.integrate import ode
# for current time
from raytracer_settings import *
# Spacepy (for coordinate transforms)
from spacepy import coordinates as coord
from spacepy.time import Ticktock

import spacepy.time as spt
"""
function: IGRFdirection

INPUT: coordinates as lat (deg N), lon (deg E), alt (m from center of Earth)
of starting point of magnetic field line

OUTPUT: unit vector of direction of local magnetic field line in lla

"""

# --------------------------- START FUNCTION ---------------------------

def IGRFdirection(t, x, direction):

    # working in XZ plane
    # take in XZ coordinates in the cartesian form
    pos = coord.Coords([x[0], 0, x[2]], 'GEO', 'car')
    print(pos)

    # conevrt ray_datenum to year and fraction of the year - needed for pyIGRF
    sec2yr = 31536000
    day2sec = 86400

    sec_in_yr = day2sec*int(days_in_the_year) + milliseconds_day/(1e3)
    year_frac = ray_datenum.year + (sec_in_yr/sec2yr)

    # convert coordinates to lla
    pos.ticks = Ticktock(ray_datenum)  # add ticks
    pos_lla = pos.convert('GEO', 'sph')

    print(pos_lla)
    # get magnetic field values
    B_out = pyIGRF.igrf_value(pos_lla.lati, pos_lla.long, pos_lla.radi, year_frac)

    Bt, Bp, Br = B_out[3:6]
    Br = -Br

    Bdir = polar2cart(Br, Bt, Bp)
    print('Bdir', Bdir)

    return [Bdir[0][0], Bdir[2][0]]

# --------------------------- END FUNCTION ---------------------------

"""
function: IGRFline

INPUT: coordinates as lat_start (deg N), lon_start (deg E), alt_start (m from center of Earth)
of starting point of magnetic field line

OUTPUT: lat, lon, alt list in deg N, deg E, m from center of Earth

"""

# --------------------------- START FUNCTION ---------------------------

def IGRFline(p0, direction):

    x = []
    z = []
    dt = 0.1
    r = ode(IGRFdirection)
    r.set_integrator('vode')

    r.set_initial_value(p0, 0)
    r.set_f_params(direction)
    counts = 0
    while r.successful():
        r.integrate(r.t + dt)
        x.append(r.y[0])
        z.append(r.y[1])

        #         print r.y, counts
        counts += 1
        if np.linalg.norm(r.y) < 1:
            print
            "hit the earth!"
            break

        if counts > 500:
            print
            "max count!"
            break
    return x, z
    # seems to be good conditions to finish the field line
    #nsteps = 1000
    #distance = 100e3
    #steplen = distance / np.abs(nsteps)
    #distance = (alt_start/R_E)*R_E*4/1e3
    #steplen = distance/np.abs(nsteps)

    # convert to km from surface of Earth for pyIGRF
    #alt_start = (alt_start - R_E) / 1e3

    # conevrt ray_datenum to year and fraction of the year - needed for pyIGRF
    #sec2yr = 31536000
    #day2sec = 86400

    #sec_in_yr = day2sec*int(days_in_the_year) + milliseconds_day/(1e3)
    #year_frac = ray_datenum.year + (sec_in_yr/sec2yr)

    # initialize
    #lat = np.zeros(nsteps)
    """
    
    lon = np.zeros(nsteps)
    alt = np.zeros(nsteps)
    lat[0] = lat_start*D2R
    lon[0] = lon_start*D2R
    alt[0] = alt_start

    for index in range(nsteps-1):

        # get magnetic field values
        B_out = pyIGRF.igrf_value(lat[index]*R2D, lon[index]*R2D, alt[index], year_frac)

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
            lat[index + 1] = lat[index] + steplen * dt / alt[index]
            lon[index + 1] = lon[index] + steplen * dp / (alt[index] * np.cos(lat[index]))

    # remove extra zeroes
    lat = [i for i in lat if i != 0]
    lon = [i for i in lon if i != 0]
    alt = [i for i in alt if i != 0]

    alt_out = [i * 1e3 + R_E for i in alt] # convert back to ECEF in m
    lat_out = [i * R2D for i in lat]
    lon_out = [i * R2D for i in lon]

    return lat_out,lon_out,alt_out
"""
# --------------------------- END FUNCTION ---------------------------

"""
function: IGRFdirection

INPUT: coordinates as lat (deg N), lon (deg E), alt (m from center of Earth)
of starting point of magnetic field line

OUTPUT: unit vector of direction of local magnetic field line in lla



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

"""

import math

def polar2cart(r, theta, phi):
    return [
         r * math.sin(theta) * math.cos(phi),
         r * math.sin(theta) * math.sin(phi),
         r * math.cos(theta)
    ]