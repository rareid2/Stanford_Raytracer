"""

functions to generate magnetic field intensity and line from IGRF13 coefficients
from Dr. Austin Sousa's thesis work
using a branched version of spacepy
author: Riley R
Apr. 2020

"""
import numpy as np
import spacepy.irbempy as irbem
import spacepy.coordinates as coord
from spacepy.time import Ticktock
import matplotlib.pyplot as plt
import datetime as dt
from scipy.integrate import ode
from raytracer_settings import *

"""
generate direction of local magnetic field line
inputs are the 3d position in cartesian geocentric coordinates
bmodel is defined as:
    Bmodel index number (from irbempy)
                - 0 = IGRF13
                - 1 = Eccentric tilted dipole
                - 2 = Jensen&Cain 1960
                - 3 = GSFC 12/66 updated to 1970
                - 4 = User-defined model (Default: Centred dipole + uniform [Dungey open model] )
                - 5 = Centred dipole
extfield = 0 (not sure what this does rn)
direction is either 1 or -1
B_dir can be called indivudally without using trace_fieldline_ODE
inputs to B_dir if called on its own: 
t = set at 0 if called indivuallay, else let ODE integrated call it
x = starting position 3 COMPONENT IN EARTH RADII GEOCENTRIC CART
"""
# need to inclide ray_datenum as an input here
def B_dir(t, x, bmodel, extfield, direction):
    pos = coord.Coords([x[0], x[1], x[2]], 'GEO', 'car')
    tv = Ticktock(ray_datenum)
    B = irbem.get_Bfield(tv, pos, extMag=extfield, options=[1, 0, 0, 0, bmodel], omnivals=None)
    Bmags = direction * B['Bvec'] / B['Blocal']
    return [Bmags[0][0], Bmags[0][1], Bmags[0][2]]

"""
trace field line uses ODE to trace along field line
call with initial position in Earth radii in
p0 = 3 COMPONENT XYZ position in eatrh radii in GEO cartesian
extfield = 0
bmodel = use IGRF = 0
direction = 1 or -1 
"""

def trace_fieldline_ODE(p0, bmodel, extfield, direction):

    x = []
    y = []
    z = []
    dt = 5e-2
    r = ode(B_dir)
    r.set_integrator('vode')

    r.set_initial_value(p0, 0)
    r.set_f_params(bmodel, extfield, direction)
    counts = 0
    while r.successful():
        r.integrate(r.t + dt)
        x.append(r.y[0])
        y.append(r.y[1])
        z.append(r.y[2])

        counts += 1
        if np.linalg.norm(r.y) < 1:
            # hit the earth
            break

        if counts > 500:
            print('max count')
            break
    return x, y, z

"""
#example call
startpoint = [1.5291777608361319, 1.5, -1.310595670385567]
direction = 1
x,y,z = trace_fieldline_ODE(startpoint,0,'0',direction)
plt.plot(x,z)
plt.show()
"""

def B_direasy(t, x):

    bmodel = 0       # IGRF13
    extfield = '0'
    direction = 1    # north

    # I have no idea why it seems like GEO is much more correct

    pos = coord.Coords([x[0], x[1], x[2]], 'GEO', 'car')
    tv = Ticktock(t, 'UTC')
    B = irbem.get_Bfield(tv, pos, extMag=extfield, options=[1, 0, 0, 0, bmodel], omnivals=None)
    Bmags = direction * B['Bvec'] / B['Blocal']
    return [Bmags[0][0], Bmags[0][1], Bmags[0][2]]

def findFootprints(t, x, hemis):

    # hemis valid cases are 'same', 'other', 'north' or 'south'

    bmodel = 0      # IGRF13
    extfield = '0'
    stopalt = 400   # km

    pos = coord.Coords([x[0], x[1], x[2]], 'GEO', 'car')
    tv = Ticktock(t, 'UTC')

    footpoint = irbem.find_footpoint(tv, pos, extMag=extfield, options=[1, 0, 0, 0, bmodel], hemi=hemis, alt=stopalt, omnivals=None)
    gdz = footpoint['loci']

    return (gdz)

"""
# passing in GEI coordinates cartesian?
startpoint = [3*R_E, 0, 0]
ray_datenum = dt.datetime(2020, 4, 29, 0, 0)
Bx, By, Bz = B_direasy(ray_datenum, startpoint)
dirB = np.reshape(np.array([Bx, By, Bz]), (1, 3))

footprint = findFootprints(ray_datenum, startpoint, 'north')
footprint.ticks = Ticktock(ray_datenum, 'UTC')
newcoord = footprint.convert('GEO', 'car')

thetalist = [0]
directions = []
# get bfield footprint from time and place pointing north at 400km
# rotate around direction of field line around x axis
for theta in thetalist:
    R = [[1, 0, 0], [0, np.cos(D2R * theta), - np.sin(D2R * theta)],
         [0, np.sin(D2R * theta), np.cos(D2R * theta)]]
    direction = np.matmul(dirB, np.reshape(np.array(R), (3, 3)))
    direction = direction / np.linalg.norm(direction)
    # add that normalized direction
    directions.append(np.squeeze(direction))

fig, ax = plt.subplots(1,1, sharex=True, sharey=True)
plt.arrow(startpoint[0]/R_E, startpoint[2]/R_E, directions[0][0], directions[0][2])
earth = plt.Circle((0, 0), 1, color='b', alpha=1, zorder=100)
ax.add_artist(earth)
plt.ylim([-5, 5])
plt.xlim([-5, 5])
plt.show()
"""