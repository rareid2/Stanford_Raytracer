"""

here is a script that will split up TLE's into orbits
to plot nicely with ray tracer

"""

# import packages
import numpy as np
import matplotlib.pyplot as plt
from skyfield.api import load, EarthSatellite
from skyfield.timelib import Time
import datetime as dt

from raytracer_settings import *

# set up TLE list
TLE_list = []

# updated 4/8/2020
DSX_TLE = '''1 44344U 19036F   20099.44261897 -.00000008 +00000-0 +00000-0 0  9998
2 44344 042.2458 098.1824 1975230 124.0282 256.3811 04.54371606013099'''
VPM_TLE = '''1 45120U 19071K   20099.39853967 +.00003361 +00000-0 +11326-3 0  9999
2 45120 051.6441 341.6999 0012197 157.2831 202.8693 15.33573896010262'''

TLE_list.append(DSX_TLE)
TLE_list.append(VPM_TLE)

# empty lists for each array
x_sat = []
y_sat = []
z_sat = []

# extract orbits for each sat
for TLE in TLE_list:
    L1, L2 = TLE.splitlines()

    # load = Loader('~/Documents/fishing/SkyData')
    ts   = load.timescale()

    Roadster = EarthSatellite(L1, L2)
    hours = np.arange(0, 6, 0.01)

    # get date and time from raytracer settings and convert for this package
    date_convert = dt.datetime(int(yearday[0:4]), 1, 1) + dt.timedelta(int(yearday[4:]) - 1)
    time = ts.utc(date_convert.year,date_convert.month, date_convert.day, hours)

    Rpos    = Roadster.at(time).position.m
    x_sat.append(Rpos[0])
    y_sat.append(Rpos[1])
    z_sat.append(Rpos[2])

# to visualize

plt.plot(x_sat[1]/R_E,z_sat[1]/R_E)
plt.show()

# end of script