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

"""
generate direction of local magnetic field line
inputs are the 3d position in cartesian geocentric coordinates in earth radii
bmodel is defined as:
    Bmodel index number (from irbempy)
                - 0 = IGRF13
                - 1 = Eccentric tilted dipole
                - 2 = Jensen&Cain 1960
                - 3 = GSFC 12/66 updated to 1970
                - 4 = User-defined model (Default: Centred dipole + uniform [Dungey open model] )
                - 5 = Centred dipole
extfield = 0
direction is either 1 or -1
"""

def B_dir(t, x, bmodel, extfield, direction, ray_datenum):
    pos = coord.Coords([x[0], x[1], x[2]], 'GEO', 'car')
    tv = Ticktock(ray_datenum, 'UTC')
    B = irbem.get_Bfield(tv, pos, extMag=extfield, options=[1, 0, 0, 0, bmodel], omnivals=None)
    Bmags = direction * B['Bvec'] / B['Blocal']
    return [Bmags[0][0], Bmags[0][1], Bmags[0][2]]

def trace_fieldline_ODE(p0, bmodel, extfield, direction, ray_datenum):

    x = []
    y = []
    z = []
    dt = 5e-2
    r = ode(B_dir)
    r.set_integrator('vode')

    r.set_initial_value(p0, 0)
    r.set_f_params(bmodel, extfield, direction, ray_datenum)
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


# takes in datetime object in UTC and GEO car coordinates in earth raddi
# returns unit vec of B field direction in cartestian
def B_direasy(t, x, thatdir):

    bmodel = 0       # IGRF13
    extfield = '0'
    direction = thatdir 

    pos = coord.Coords([x[0], x[1], x[2]], 'GEO', 'car')
    tv = Ticktock(t, 'UTC')
    B = irbem.get_Bfield(tv, pos, extMag=extfield, options=[1, 0, 0, 0, bmodel], omnivals=None)
    Bmags = direction * B['Bvec'] / B['Blocal']
    
    return [Bmags[0][0], Bmags[0][1], Bmags[0][2]]


# takes in datetime object in UTC and GEO car coordinates in earth radii
# returns spherical coord (alt, lat, lon) in GDZ in km deg deg
def findFootprints(t, x, hemis):

    # hemis valid cases are 'same', 'other', 'north' or 'south'

    bmodel = 0      # IGRF13
    extfield = '0'
    stopalt = 450   # km

    pos = coord.Coords([x[0], x[1], x[2]], 'GEO', 'car')
    tv = Ticktock(t, 'UTC')

    footpoint = irbem.find_footpoint(tv, pos, extMag=extfield, options=[1, 0, 0, 0, bmodel], hemi=hemis, alt=stopalt, omnivals=None)
    gdz = footpoint['loci']

    return (gdz)
