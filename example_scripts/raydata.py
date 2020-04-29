
"""
here is a script that will call run_rays and save
"""

# import needed packages
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import datetime as dt
from dateutil import parser
from example_scripts.raytracer_utils import readdump, read_rayfile, read_rayfiles, read_damp
from example_scripts.run_rays import run_rays
from example_scripts.raytracer_settings import *
from example_scripts.IGRF_funcs import B_dir, trace_fieldline_ODE, findFootprints, B_direasy
from spacepy import coordinates as coord
from spacepy.time import Ticktock
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm, ListedColormap, BoundaryNorm
from get_TLE import get_TLE
import xlsxwriter

# --------------------------- Ray Tracing --------------------------
# define lists here - must be lists even if only one arg
# started at midnight on the 29th

# boomerang - 3.2? check  out ramp!

freq = [8.2e3]
thetalist = [0]  # in deg

# read in orbital data
# comes in as GCRS car coordinate system in m == ECI == GEI car??????
orbitdata = np.genfromtxt('orbit_pos.txt')
dsxdata = orbitdata[:,0:3]

# convert back to time data - time is in UTC
f = open('orbit_time.txt', 'r')
times = f.readlines()
timestr = [t.strip() for t in times]
timedata = [parser.parse(t) for t in timestr]

print('got data')

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
n = 150000

bfoots = []
workbook = xlsxwriter.Workbook('more_conjunctions.xlsx')
worksheet = workbook.add_worksheet()

s = 0
# iterate
for numcount in np.linspace(0,looplen-1,n):

    fig, ax = plt.subplots(1,1, sharex=True, sharey=True)
    lw = 2  # linewidth

    numcount = int(numcount)
    print(numcount)
    positions = [dsxpositions[numcount]]
    ray_datenum = raytime[numcount]

    # initialize - leave empty
    directions = []

    for position in positions:
        # grab position and find direction of local bfield
        startpoint = [position[0]/R_E, position[1]/R_E, position[2]/R_E]
        Bx, By, Bz = B_direasy(ray_datenum, startpoint)
        dirB = np.reshape(np.array([Bx, By, Bz]), (1, 3))

        # get bfield footprint from time and place pointing north at 400km
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
        tmp_coords.ticks = Ticktock(tvec_datetime, 'UTC')  # add ticks
        tmp_coords.sim_time = r['time']
        new_coords = tmp_coords.convert('GEO', 'car')
        rays.append(new_coords)

    #initialize
    rx = []
    ry = []
    rz = []

    for r in rays:
        rx.append(r.x / R_E)
        ry.append(r.y / R_E)
        rz.append(r.z / R_E)

    dlist = []
    for d in damplist:
        damp = d["damping"]
        damp = np.squeeze(np.array(damp))
        dlist.append(damp)

    def myplot(ax, xs, ys, zs, cmap):
        for x, y, z in zip(xs, ys, zs):
            plot = LineCollection([np.column_stack((x, y))], cmap=cmap, zorder=102)
            plot.set_array(z)
            ax.add_collection(plot)
        return plot

    line = myplot(ax, rx, rz, dlist, 'Reds')
    fig.colorbar(line, ax=ax, label = 'Normalized wave power')
    earth = plt.Circle((0, 0), 1, color='b', alpha=1, zorder=100)
    ax.add_artist(earth)
    ax.set_aspect('equal')
    max_lim = 4

    plt.xticks(np.arange(-max_lim, max_lim, step=1))
    plt.yticks(np.arange(-max_lim, max_lim, step=1))
    plt.xlabel('L (R$_E$)')
    plt.ylabel('L (R$_E$)')
    plt.xlim([-max_lim, max_lim])
    plt.ylim([-2.5, 2.5])

    savename = 'plots/ray_' + str(ray_datenum) + '.png'
    fig.savefig(savename, format='png')
    plt.close()
    #plt.show()

    # ------------------------------- Saving ---------------------------------------

    footprint = findFootprints(ray_datenum, startpoint, 'north')
    footprint.ticks = Ticktock(ray_datenum, 'UTC')
    footprint = footprint.convert('GEO', 'car')
    bfoots.append(footprint)
    for ir in range(len(rx[0])):
        worksheet.write(ir + 1, s * 5, rx[0][ir])
        worksheet.write(ir + 1, (s * 5) + 1, ry[0][ir])
        worksheet.write(ir + 1, (s * 5) + 2, rz[0][ir])
        worksheet.write(ir + 1, (s * 5) + 3, dlist[0][ir])

    s += 1

workbook.close()