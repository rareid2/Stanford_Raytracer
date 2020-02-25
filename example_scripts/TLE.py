import numpy as np
import matplotlib.pyplot as plt
from skyfield.api import load, EarthSatellite
from skyfield.timelib import Time
from spacepy import coordinates as coord
from spacepy.time import Ticktock
from datetime import datetime

R_E = 6371e3  # m

DSX_TLE = '''1 44344U 19036F   20051.25825015 -.00000031  00000-0  00000+0 0  9998
2 44344  42.2223 115.8421 1974926 103.3809 279.3131  4.54370801 10908'''
L1, L2 = DSX_TLE.splitlines()

# load = Loader('~/Documents/fishing/SkyData')
ts   = load.timescale()

Roadster = EarthSatellite(L1, L2)
hours = np.arange(0, 6, 0.01)
time = ts.utc(2020, 1, 1, hours)

Rpos    = Roadster.at(time).position.m
x_DSX, y_DSX, z_DSX = Rpos

plt.plot(x_DSX/(R_E),z_DSX/(R_E), c='g')
plt.show()
