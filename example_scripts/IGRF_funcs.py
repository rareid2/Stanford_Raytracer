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
from example_scripts.raytracer_settings import *


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

startpoint = [1.5291777608361319, 1.5, -1.310595670385567] # Start from the equator
direction = 1
x,y,z = trace_fieldline_ODE(startpoint,0,'0',direction)
plt.plot(x,z)
plt.show()
"""