
"""
here is a script that will call run_rays and plot the trajectory with
normalized wave power as a color scale
this is currently set for XZ coordinates in SM
ONLY WORKS FOR A SINGLE POSITION - CHANGE THE TIME FOR DIFF POSITIONS
"""

# import needed packages
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import datetime as dt
from raytracer_utils import readdump, read_rayfile, read_rayfiles, read_damp
from run_rays import run_rays
from raytracer_settings import *
from IGRF_funcs import B_dir, trace_fieldline_ODE, findFootprints, B_direasy
from spacepy import coordinates as coord
from spacepy.time import Ticktock
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm, ListedColormap, BoundaryNorm
from TLE_funcs import TLE2pos

# coordinate mania!
# TLES give us GEI, raytracer needs SM, and IGRF funcs needs GEO

# -------------------------------- SET TIME --------------------------------
# change time information here - use UTC -
year = 2020
month = 5
day = 17
hours = 12
minutes = 0
seconds = 0

ray_datenum = dt.datetime(year, month, day, hours, minutes, seconds)

# -------------------------------- GET POSITIONS --------------------------------
# these will be in ECI coordinates (GEI)
# DSX TLE
l11 = '1 44344U 19036F   20130.24435661 -.00000027 +00000-0 +00000+0 0  9994'
l21 = '2 44344 042.2583 086.8979 1974641 137.2296 239.9665 04.54371389014496'
# VPM TLE
l12 = '1 45120U 19071K   20132.49935632  .00001453  00000-0  55129-4 0  9998'
l22 = '2 45120  51.6416 181.7127 0011592 280.5137  79.4539 15.33820525015342'

lines1 = [l11, l12]
lines2 = [l21, l22]
satnames = ['DSX', 'VPM']

# get DSX and VPM positions for... 
plen = 1  # second
r, tvec = TLE2pos(lines1, lines2, satnames, plen, ray_datenum)

# convert to meters
dsx = [rpos*1e3 for rpos in r[0]]
vpm = [rpos*1e3 for rpos in r[1]]

# convert startpoint to SM for raytracer
dsxpos = coord.Coords(dsx, 'GEI', 'car', units=['m', 'm', 'm'])
dsxpos.ticks = Ticktock(ray_datenum, 'UTC') # add ticks
SM_dsx = dsxpos.convert('SM', 'car')

# convert vpm to SM for plotting
vpmpos = coord.Coords(vpm, 'GEI', 'car', units=['m', 'm', 'm'])
vpmpos.ticks = Ticktock(ray_datenum, 'UTC') # add ticks
SM_vpm = vpmpos.convert('SM', 'car')

# -------------------------------- DEFINE RAY DIRECTIONS --------------------------------
position = [float(SM_dsx.x), float(SM_dsx.y), float(SM_dsx.z)]
positions = []
freq = [18e3] # Hz
directions = []
thetalist = [0, 5, 10, 15, 20, 25, 30, 35, 45, -5, -10, -15, -20, -25, -30, -35, -45]  # in deg -- what angles to launch at? 

# grab position and find direction of local bfield
# convert to RE for bfield lib - SM is okay here
startpoint = [position[0]/R_E, position[1]/R_E, position[2]/R_E]
Bx, By, Bz = B_direasy(tvec, startpoint)
dirB = np.reshape(np.array([Bx, By, Bz]), (1, 3))

# rotate around direction of field line around x axis
for theta in thetalist:
    R = [ [1, 0, 0], [0, np.cos(D2R * theta), - np.sin(D2R * theta)],
        [0, np.sin(D2R * theta), np.cos(D2R * theta)] ]
    direction = np.matmul(dirB, np.reshape(np.array(R), (3, 3)))
    direction = direction/np.linalg.norm(direction)
    
    # add that normalized direction
    directions.append(np.squeeze(direction))

    # make sure position list matches direction list
    positions.append(position)
    
# -------------------------------- RUN RAYS --------------------------------
# convert for raytracer settings
days_in_the_year = ray_datenum.timetuple().tm_yday
days_in_the_year = format(days_in_the_year, '03d')

# yearday and miliseconds day are used by raytracer
yearday = str(year)+ str(days_in_the_year)   # YYYYDDD
milliseconds_day = hours*3.6e6 + minutes*6e4 + seconds*1e3

# run it!
run_rays(freq, positions, directions, yearday, milliseconds_day)

# -------------------------------- LOAD OUTPUT --------------------------------
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

# -------------------------------- CONVERT COORDINATES --------------------------------
# convert to desired coordinate system into vector list rays
rays = []
for r in raylist:
    tmp_coords = coord.Coords(list(zip(r['pos'].x, r['pos'].y, r['pos'].z)), 'SM', 'car', units=['m', 'm', 'm'])
    tvec_datetime = [ray_datenum + dt.timedelta(seconds=s) for s in r['time']]
    tmp_coords.ticks = Ticktock(tvec_datetime, 'UTC')  # add ticks
    tmp_coords.sim_time = r['time']
    #new_coords = tmp_coords.convert('GEO', 'car')
    #rays.append(new_coords)
    rays.append(tmp_coords)

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

# -------------------------------- PLOTTING --------------------------------
fig, ax = plt.subplots(1,1, sharex=True, sharey=True)
lw = 2  # linewidth

# for z axis rotation of ray, fieldlines, and sat positions
LLA_dsx = SM_dsx.convert('SM', 'sph')
th = LLA_dsx.long

# rotate and plot sat positions in SM RE
dsxx = startpoint[0] * np.cos(np.deg2rad(-th)) - startpoint[1] * np.sin(np.deg2rad(-th))
dsxy = startpoint[0] * np.sin(np.deg2rad(-th)) + startpoint[1] * np.cos(np.deg2rad(-th))
dsxz = startpoint[2]

vpmx = SM_vpm.x/R_E * np.cos(np.deg2rad(-th)) - SM_vpm.y/R_E * np.sin(np.deg2rad(-th))
vpmy = SM_vpm.x/R_E * np.sin(np.deg2rad(-th)) + SM_vpm.y/R_E * np.cos(np.deg2rad(-th))
vpmz = SM_vpm.z/R_E

plt.plot(dsxx, dsxz, '-go', zorder=105, label='DSX')
plt.plot(vpmx, vpmz, '-yo', zorder=106, label='VPM')

# rotate rays
rxr = rx * np.cos(np.deg2rad(-th)) - ry * np.sin(np.deg2rad(-th))
ryr = rx * np.sin(np.deg2rad(-th)) + ry * np.cos(np.deg2rad(-th))
rzr = rz

# create line segments for plotting
# one line for each ray
for p in range(len(rxr)):
    points = np.array([rxr[p], rzr[p]]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    line_segments = LineCollection(segments, cmap = 'Reds')
    ax.add_collection(line_segments)
    line_segments.set_array(dlist[p])
# add in color bar
axcb = fig.colorbar(line_segments, ax=ax, label = 'Normalized wave power')

# -------------------------------- EARTH AND IONO --------------------------------
earth = plt.Circle((0, 0), 1, color='b', alpha=0.5, zorder=100)
iono = plt.Circle((0, 0), (R_E + H_IONO) / R_E, color='g', alpha=0.5, zorder=99)
ax.add_artist(earth)
ax.add_artist(iono)

# -------------------------------- PLASMASPHERE --------------------------------
plasma_model_dump = os.path.join(ray_out_dir, 'model_dump_mode_1_XZ.dat')
d_xz = readdump(plasma_model_dump)
Ne_xz = d_xz['Ns'][0, :, :, :].squeeze().T * 1e-6
Ne_xz[np.isnan(Ne_xz)] = 0

# Axis spacing depends on how the modeldump was ran
psize = 10
px = np.linspace(-10, 10, 200)
py = np.linspace(-10, 10, 200)

# Colorbar limits (log space)
clims = [-2, 5]

# Plot background plasma (equatorial slice)
g = plt.pcolormesh(px, py, np.log(Ne_xz), cmap = 'twilight')
#fig.colorbar(g, ax=ax, orientation="horizontal", pad = 0.2, label= 'Plasmasphere density')

# ---------------------------------- BFIELD -----------------------------------

 # need to convert to GEO car
CAR_dsx = SM_dsx.convert('GEO', 'car')
bstart = [float(CAR_dsx.x)/R_E, float(CAR_dsx.y)/R_E, float(CAR_dsx.z)/R_E]

# (from IGRF13 model)
L_shells = [2, 3, 4]  # Field lines to draw in plane with DSX
for L in L_shells:
    Lx, Ly, Lz = trace_fieldline_ODE([L,0,0], 0, '0', 1, ray_datenum)
    plt.plot(Lx, Lz, color='b', linewidth=1, linestyle='dashed')
    Lx, Ly, Lz = trace_fieldline_ODE([L,0,0], 0, '0', -1, ray_datenum)
    plt.plot(Lx, Lz, color='b', linewidth=1, linestyle='dashed')
    Lx, Ly, Lz = trace_fieldline_ODE([-L,0,0], 0, '0', 1, ray_datenum)
    plt.plot(Lx, Lz, color='b', linewidth=1, linestyle='dashed')
    Lx, Ly, Lz = trace_fieldline_ODE([-L,0,0], 0, '0', -1, ray_datenum)
    plt.plot(Lx, Lz, color='b', linewidth=1, linestyle='dashed')
print('finished plotting field lines')

# plot field line from orbital position
Blines = []

Blines.append(trace_fieldline_ODE(bstart, 0, '0', 1, ray_datenum))
Blines.append(trace_fieldline_ODE(bstart, 0, '0', -1, ray_datenum))

for blinex, bliney, blinez in Blines:
    raytime = []

    # create list of the same times
    for m in range(int(len(blinex))):
        raytime.append(ray_datenum)
    
    # convert to SM coords
    bpos = np.column_stack((blinex, bliney, blinez))
    bpos = coord.Coords(bpos, 'GEO', 'car', units=['Re', 'Re', 'Re'])
    bpos.ticks = Ticktock(raytime, 'UTC') # add ticks
    SM_b = bpos.convert('SM', 'car')

    # rotate around z axis
    #LLA_dsx = SM_dsx.convert('SM', 'sph')
    #th = LLA_dsx.long
    brot_x = SM_b.x * np.cos(np.deg2rad(-th)) - SM_b.y * np.sin(np.deg2rad(-th))
    brot_y = SM_b.x * np.sin(np.deg2rad(-th)) + SM_b.y * np.cos(np.deg2rad(-th))
    brot_z = SM_b.z

    # plot
    plt.plot(brot_x, brot_z, color='r', linewidth=1, linestyle='dashed')

# -------------------------------- GET FOOTPRINT --------------------------------
# also in GEO car, so need to use bstart 
footprint = findFootprints(ray_datenum, bstart, 'north')
footprint.ticks = Ticktock(ray_datenum, 'UTC')
footprint = footprint.convert('SM', 'car')
# rotate around z axis as well
foot_x = footprint.x * np.cos(np.deg2rad(-th)) - footprint.y * np.sin(np.deg2rad(-th))
foot_y = footprint.x * np.sin(np.deg2rad(-th)) + footprint.y * np.cos(np.deg2rad(-th))
foot_z = footprint.z
plt.plot(foot_x, foot_z, '-ro', label='Bfield footpoint')

# -------------------------------- FORMATTING --------------------------------
ax.set_aspect('equal')
max_lim = 4

plt.xticks(np.arange(-max_lim, max_lim, step=1))
plt.yticks(np.arange(-max_lim, max_lim, step=1))
plt.xlabel('L (R$_E$)')
plt.ylabel('L (R$_E$)')
plt.xlim([-max_lim, max_lim])
plt.ylim([-2.5, 2.5])

mytitle = str(freq[0]/1e3) + 'kHz rays at ' + str(ray_datenum)
plt.title(mytitle)
ax.legend(loc = 'lower center', fontsize =  'x-small')

#savename = 'raytest.png'
#fig.savefig(savename, format='png')
#plt.close()
plt.show()
