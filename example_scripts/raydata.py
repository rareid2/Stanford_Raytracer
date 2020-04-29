
"""
here is a script that will call run_rays and save
"""

# import needed packages
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import datetime as dt
from example_scripts.raytracer_utils import readdump, read_rayfile, read_rayfiles, read_damp
from example_scripts.run_rays import run_rays
from example_scripts.raytracer_settings import *
from example_scripts.IGRF_funcs import B_dir, trace_fieldline_ODE
from spacepy import coordinates as coord
from spacepy.time import Ticktock
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm, ListedColormap, BoundaryNorm
from get_TLE import get_TLE

# --------------------------- Ray Tracing --------------------------
# define lists here - must be lists even if only one arg
# started at midnight on the 29th

# boomerang - 3.2? check  out ramp!

freq = [8.2e3]

orbitdata = np.genfromtxt('orbit_pos.txt')
timedata = np.genfromtxt('orbit_time.txt')

dsxdata = orbitdata[:,0:3]

# get when in apo/peri
apa = []
for i in range(len(dsxdata)):
    dat = dsxdata[i]
    x = dat[0] / R_E
    if x<-1 or x>1:
        apa.append(i)

dsxpositions = []
vpmpositions = []
raytime = []

for ap in apa:
    dsxpositions.append(orbitdata[ap][0:3])
    vpmpositions.append(orbitdata[ap][3:6])
    raytime.append(timedata[ap])
print('got all positions')

looplen = len(dsxpositions)

# iterate
for numcount in range(looplen):
    positions = [dsxpositions[numcount]]
    ray_datenum = raytime[numcount]

    thetalist = [0] # in deg

    # initialize - leave empty
    directions = []
    Blines = []

    for position in positions:
        # grab position
        startpoint = [position[0]/R_E, position[1]/R_E, position[2]/R_E]
        # get bfield direction
        Bx, By, Bz = B_dir(0, startpoint, 0, '0', 1)
        dirB = np.reshape(np.array([Bx, By, Bz]), (1,3))

        # grab full magnetic field line for plotting later
        Blines.append(trace_fieldline_ODE(startpoint, 0, '0', 1))
        Blines.append(trace_fieldline_ODE(startpoint, 0, '0', -1))

        # rotate around direction of field line around x axis
        for theta in thetalist:
            R = [ [1, 0, 0], [0, np.cos(D2R * theta), - np.sin(D2R * theta)],
                  [0, np.sin(D2R * theta), np.cos(D2R * theta)] ]
            direction = np.matmul(dirB, np.reshape(np.array(R), (3, 3)))
            direction = direction/np.linalg.norm(direction)
            # add that normalized direction
            directions.append(np.squeeze(direction))

    # run!
    run_rays(freq, positions, directions)

    # ---------------------- Load output directory -------------------------

    # Load all the rayfiles in the output directory
    file_titles = os.listdir(ray_out_dir)

    # create empty lists to fill with ray files and damp files
    raylist = []
    damplist = []

    for filename in file_titles:
        if '.ray' in filename:
            raylist += read_rayfile(os.path.join(ray_out_dir, filename))

    for filename in file_titles:
        if '.damp' in filename:
            damplist += read_damp(os.path.join(ray_out_dir, filename))

    # quick check: did the rays propagate?
    raylist = [checkray for checkray in raylist if not len(checkray["time"]) < 2]

    # abandon if not
    if raylist == []:
        sys.exit(0)

    # ------------------------ Coordinate Conversion --------------------------

    # convert to desired coordinate system into vector list rays
    rays = []
    for r in raylist:
        tmp_coords = coord.Coords(list(zip(r['pos'].x, r['pos'].y, r['pos'].z)), 'SM', 'car', units=['m', 'm', 'm'])
        tvec_datetime = [ray_datenum + dt.timedelta(seconds=s) for s in r['time']]
        tmp_coords.ticks = Ticktock(tvec_datetime)  # add ticks
        tmp_coords.sim_time = r['time']
        rays.append(tmp_coords)

    #initialize
    r_length = []
    rx = []
    ry = []
    rz = []

    for r in rays:
        rx.append(r.x / R_E)
        ry.append(r.y / R_E)
        rz.append(r.z / R_E)
        r_length.append(len(r))

    dlist = []
    for d in damplist:
        damp = d["damping"]
        damp = np.squeeze(np.array(damp))
        #if len(damp) < max(r_length):
        #    leftover = max(r_length) - len(damp)
        #    damp = np.concatenate((damp, np.zeros(int(leftover))), axis=0)
        dlist.append(damp)

    # ------------------------------- Saving ---------------------------------------

    # need to save the ray trajectories (or really just the end points)
    # damping
    # and the bfield line
    # also need the time tvec_datetime that goes with it...
    # amp of ray too?
    # what will I put into this plot? one axis will be time throughout the 2 week period
    # the next will be distance from ray end point to VPM throughout that time
    # so I'll need data points showing VPM's location for all time
    # and the ray location for the times I have them