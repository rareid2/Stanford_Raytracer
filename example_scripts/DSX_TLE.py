import numpy as np
import matplotlib.pyplot as plt
from skyfield.api import load, EarthSatellite
from skyfield.timelib import Time

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


n = 500
slope = []
for i in range(0, np.size(x_DSX)-1):
    slope.append((z_DSX[i+1]-z_DSX[i])/(x_DSX[i+1]-x_DSX[i]))
slope = np.array(slope)

x = np.linspace(0,3,10)
mytan = slope[n]*(x - (x_DSX[n])) + z_DSX[n]

plt.plot(x_DSX,z_DSX, 'g')
plt.plot(x_DSX[n],z_DSX[n], 'bo')
plt.plot(x, mytan, 'r')
plt.xlim([-3,3])
plt.ylim([-3,3])

plt.show()
theta = np.rad2deg(np.arctan(slope[n]))