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
    distance = 100e3
    steplen = distance / np.abs(nsteps)
    #distance = (alt_start/R_E)*R_E*4/1e3
    #steplen = distance/np.abs(nsteps)

    # convert to km from surface of Earth for pyIGRF
    #alt_start = (alt_start - R_E) / 1e3
    #lat_start = lat_start * D2R # theta in sph
    #lon_start = lon_start * D2R # phi in sph

    # conevrt ray_datenum to year and fraction of the year - needed for pyIGRF
    sec2yr = 31536000
    day2sec = 86400

    sec_in_yr = day2sec*int(days_in_the_year) + milliseconds_day/(1e3)
    year_frac = ray_datenum.year + (sec_in_yr/sec2yr)

    # initialize
    lat = np.zeros(nsteps)
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
    print(alt)

    alt_out = [i * 1e3 + R_E for i in alt] # convert back to ECEF in m
    lat_out = [i * R2D for i in lat]
    lon_out = [i * R2D for i in lon]
    print(alt_out)

    return lat_out,lon_out,alt_out

# --------------------------- END FUNCTION ---------------------------

from spacepy import coordinates as coord
from spacepy.time import Ticktock

lat, lon, alt = IGRFline(0,0,3*R_E/1e3)
sph_coords = list(zip(alt, lat, lon))
cvals = coord.Coords(sph_coords, 'GEO', 'sph')
tvec_datetime = [ray_datenum + dt.timedelta(seconds=s) for s in range(len(sph_coords))]
cvals.ticks = Ticktock(tvec_datetime)  # add ticks
newcoord = cvals.convert('GEO', 'car')
#print(newcoord)

#lat, lon, alt = IGRFline(0,0,2*R_E)

#x = [r*np.sin(deg2rad) for r in lat if i != 0]


fig, ax = plt.subplots(1,1, sharex=True, sharey=True)
earth = plt.Circle((0, 0), 1, color='b', alpha=1)
ax.add_artist(earth)

plt.plot(newcoord.x/R_E, newcoord.z/R_E)
#plt.plot(newcoord.x/R_E, newcoord.y/R_E)
#plt.plot(newcoord.x/R_E, newcoord.z/R_E)
plt.show()