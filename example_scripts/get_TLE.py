"""

here is a function that will split up TLE's into an orbit
to plot nicely with ray tracer

get TLE's from http://celestrak.com/satcat/search.php

DSX NORAD ID: 44344
VPM NORAD ID: 45120

returns position of satellite for one orbit, with orbit starting at given time
in geocentric cartesian coordinates in meters

"""

# import packages
import numpy as np
import matplotlib.pyplot as plt

from skyfield.api import EarthSatellite, Topos, load, utc
import datetime as dt
from spacepy.coordinates import Coords
from spacepy.time import Ticktock

from raytracer_settings import *


# --------------------------------- START FUNCTION -------------------------------------
def get_TLE(line1,line2,sat_name):

    # load timescale UTC
    ts = load.timescale()

    # TLE form:
    # line1 = '1 44344U 19036F   20099.44261897 -.00000008 +00000-0 +00000-0 0  9998'
    # line2 = '2 44344 042.2458 098.1824 1975230 124.0282 256.3811 04.54371606013099'

    # find the satellite
    satellite = EarthSatellite(line1, line2, sat_name, ts)

    # find when TLE was generated - keep updated every 1-2 weeks
    print(sat_name, 'TLE is current as of:', satellite.epoch.utc_jpl())

    # grab time from ray_tracer settings
    datenum = ray_datenum.replace(tzinfo=utc)  # specifiy UTC time zone

    # calculate orbital period
    orbital_period = 3600*24/float(line2.rsplit(None, 1)[-1])

    # generate time vector
    t_pos = [datenum + dt.timedelta(seconds=s) for s in np.linspace(0,orbital_period)]
    t = ts.utc(t_pos)

    # find geocentric cartesian coordinates over orbit for satellite
    geocentric = satellite.at(t)
    # print(geocentric.position.km)

    return geocentric.position.m


"""
# visualize orbit in XZ plane
plt.plot(geocentric.position.km[0],geocentric.position.km[1])
plt.show()
"""


# --------------------------------- END FUNCTION -------------------------------------