"""

here is a function that will split up TLE's into an orbit
to plot nicely with ray tracer

get TLE's from http://celestrak.com/satcat/search.php

DSX NORAD ID: 44344
VPM NORAD ID: 45120

returns position of satellite for one orbit, with orbit starting at given time
<<<<<<< HEAD
=======
in geocentric cartesian coordinates in meters
>>>>>>> a40385af70d8379ef79a0df25c5c9dc386b09bc8

"""

# import packages
import numpy as np
import matplotlib.pyplot as plt
<<<<<<< HEAD
=======

>>>>>>> a40385af70d8379ef79a0df25c5c9dc386b09bc8
from skyfield.api import EarthSatellite, Topos, load, utc
import datetime as dt
from spacepy.coordinates import Coords
from spacepy.time import Ticktock
<<<<<<< HEAD
=======

>>>>>>> a40385af70d8379ef79a0df25c5c9dc386b09bc8
from raytracer_settings import *

# --------------------------------- START FUNCTION -------------------------------------
def get_TLE(line1, line2, sat_name):

    # load timescale UTC
    ts = load.timescale()

    # TLE form:
    # line1 = '1 44344U 19036F   20099.44261897 -.00000008 +00000-0 +00000-0 0  9998'
    # line2 = '2 44344 042.2458 098.1824 1975230 124.0282 256.3811 04.54371606013099'

    # find the satellite
    satellite = EarthSatellite(line1, line2)

    # find when TLE was generated - keep updated every 1-2 weeks
    print(sat_name, 'TLE is current as of:', satellite.epoch.utc_jpl())

    # grab time from ray_tracer settings
    datenum = ray_datenum.replace(tzinfo=utc)  # specifiy UTC time zone

    # calculate orbital period
    orbital_period = 3600*24/float(line2.rsplit(None, 1)[-1])

    # generate time vector
    t_pos = [datenum + dt.timedelta(seconds=s) for s in np.linspace(0,2*orbital_period)]
    t = ts.utc(t_pos)
    # find geocentric cartesian coordinates over orbit for satellite
    geocentric = satellite.at(t)

    return geocentric.position.m, t_pos
<<<<<<< HEAD

"""

#for testing
line1 = '1 44344U 19036F   20117.92375283 -.00000009 +00000-0 +00000-0 0  9998'
line2 = '2 44344 042.2535 091.4116 1974932 131.9492 246.6867 04.54371554013931'
x, t_pos = get_TLE(line1, line2, 'DSX')
# visualize orbit in XZ plane
plt.plot(x[0]/R_E,x[1]/R_E)
plt.show()
"""

=======
"""
# for testing
line1 = '1 44344U 19036F   20113.30349832 -.00000013 +00000-0 +00000-0 0  9992'
line2 = '2 44344 042.2517 093.1041 1975016 129.9693 249.1586 04.54371641013721'
x, t_pos = get_TLE(line1, line2, 'DSX')
# visualize orbit in XZ plane
plt.plot(x[0]/R_E,x[2]/R_E)
plt.show()
"""




>>>>>>> a40385af70d8379ef79a0df25c5c9dc386b09bc8
# --------------------------------- END FUNCTION -------------------------------------

# --------------------------------- START FUNCTION -------------------------------------
def get_pos(line1, line2, sat_name, datenum):

    # load timescale UTC
    ts = load.timescale()

<<<<<<< HEAD
    # find the satellite
    satellite = EarthSatellite(line1, line2)

    # grab time from ray_tracer settings
    datenum = datenum.replace(tzinfo=utc)  # specifiy UTC time zone
    t = ts.utc(datenum)

    # find satellite
    geocentric = satellite.at(t)

    return geocentric.position.m, t

    # this is GEI in meters cartesian
=======
    # TLE form:
    # line1 = '1 44344U 19036F   20099.44261897 -.00000008 +00000-0 +00000-0 0  9998'
    # line2 = '2 44344 042.2458 098.1824 1975230 124.0282 256.3811 04.54371606013099'

    # find the satellite
    satellite = EarthSatellite(line1, line2)

    # find when TLE was generated - keep updated every 1-2 weeks
    # print(sat_name, 'TLE is current as of:', satellite.epoch.utc_jpl())

    # grab time from ray_tracer settings
    datenum = datenum.replace(tzinfo=utc)  # specifiy UTC time zone
    t = ts.utc(datenum)
    # find geocentric cartesian coordinates over orbit for satellite
    geocentric = satellite.at(t)

    return geocentric.position.m, t
>>>>>>> a40385af70d8379ef79a0df25c5c9dc386b09bc8
