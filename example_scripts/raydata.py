
"""
here is a script that will call run_rays and save
"""

# import needed packages
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import datetime as dt
from dateutil import parser
from raytracer_utils import readdump, read_rayfile, read_rayfiles, read_damp
from run_rays import run_rays
from raytracer_settings import *
from IGRF_funcs import B_dir, trace_fieldline_ODE, findFootprints, B_direasy
from spacepy import coordinates as coord
from spacepy.time import Ticktock
from TLE_funcs import TLE2posfast

# -------------------------------- SET TIME --------------------------------
# change time information here - use UTC -
year = 2020
month = 5
day = 17
hours = 0
minutes = 0
seconds = 0

ray_datenum = dt.datetime(year, month, day, hours, minutes, seconds)

# -------------------------------- GET POSITIONS --------------------------------
# DSX TLE
l11 = '1 44344U 19036F   20130.24435661 -.00000027 +00000-0 +00000+0 0  9994'
l21 = '2 44344 042.2583 086.8979 1974641 137.2296 239.9665 04.54371389014496'
# VPM TLE
l12 = '1 45120U 19071K   20132.49935632  .00001453  00000-0  55129-4 0  9998'
l22 = '2 45120  51.6416 181.7127 0011592 280.5137  79.4539 15.33820525015342'

lines1 = [l11, l12]
lines2 = [l21, l22]
satnames = ['DSX', 'VPM']

# get DSX and VPM positions for a day
r, tvec = TLE2posfast(lines1, lines2, satnames, 3, ray_datenum)

# convert to meters and SM coord
dsx = [rpos*1e3 for rpos in r[0]]
vpm = [rpos*1e3 for rpos in r[1]]

dray = []
dfoot = []

#dsxpos = coord.Coords(dsx, 'GDZ', 'car', units=['m', 'm', 'm'])
#dsxpos.ticks = Ticktock(ray_datenum, 'UTC') # add ticks
#SM_dsx = dsxpos.convert('SM', 'car')
#print(SM_dsx)

# -------------------------------- DEFINE RAY DIRECTIONS --------------------------------
positions = dsx
freq = [8.2e3] # Hz
directions = []
thetalist = [0]  # in deg -- what angles to launch at? 

i = 0
for position, rayt in zip(positions, tvec):
    print(position)
    # grab position and find direction of local bfield
    startpoint = [position[0]/R_E, position[1]/R_E, position[2]/R_E]
    Bx, By, Bz = B_direasy(rayt, startpoint)
    dirB = np.reshape(np.array([Bx, By, Bz]), (1, 3))

    # rotate around direction of field line around x axis
    for theta in thetalist:
        R = [ [1, 0, 0], [0, np.cos(D2R * theta), - np.sin(D2R * theta)],
            [0, np.sin(D2R * theta), np.cos(D2R * theta)] ]
        direction = np.matmul(dirB, np.reshape(np.array(R), (3, 3)))
        direction = direction/np.linalg.norm(direction)
        # add that normalized direction
        directions.append(np.squeeze(direction))
 
    # -------------------------------- RUN RAYS --------------------------------
    # convert for raytracer settings
    days_in_the_year = rayt.timetuple().tm_yday
    days_in_the_year = format(days_in_the_year, '03d')

    # yearday and miliseconds day are used by raytracer
    yearday = str(year)+ str(days_in_the_year)   # YYYYDDD
    milliseconds_day = hours*3.6e6 + minutes*6e4 + seconds*1e3

    # position is in GEO meters - is that correct? 
    run_rays(freq, [position], directions, yearday, milliseconds_day)

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

    # quick check: did the rays propagate?
    raylist = [checkray for checkray in raylist if not len(checkray["time"]) < 2]

    # abandon if not
    if raylist == []:
        sys.exit(0)

    # -------------------------------- CONVERT COORDINATES --------------------------------
    # convert to desired coordinate system into vector list rays
    rays = []
    for r in raylist:
        tmp_coords = coord.Coords(list(zip(r['pos'].x, r['pos'].y, r['pos'].z)), 'SM', 'car', units=['m', 'm', 'm'])
        tvec_datetime = [rayt + dt.timedelta(seconds=s) for s in r['time']]
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

    # -------------------------------- GET FOOTPRINT --------------------------------
    footprint = findFootprints(rayt, startpoint, 'north')
    footprint.ticks = Ticktock(rayt, 'UTC')
    footprint = footprint.convert('GEO', 'car')
    
    # -------------------------------- FIND DISTANCE --------------------------------
    dr = np.sqrt((rx[0][-1] - vpm[i][0]/R_E)**2 + (ry[0][-1] - vpm[i][1]/R_E)**2 + (rz[0][-1] - vpm[i][2]/R_E)**2)
    df = np.sqrt((footprint.x - vpm[i][0]/R_E)**2 + (footprint.y - vpm[i][1]/R_E)**2 + (footprint.z - vpm[i][2]/R_E)**2)

    dray.append(dr*R_E/1e3) # in km
    dfoot.append(df*R_E/1e3)  # in km
    print(i, 'dist is', dr*R_E/1e3, ' km')
    i+=1

# easy
fig, ax = plt.subplots(1, 1)
ax.plot(tvec, dray, label='field-aligned ray')
ax.plot(tvec, dfoot, label='fieldline footprint')

# formatting
for label in ax.get_xticklabels():
    label.set_rotation(40)
    label.set_horizontalalignment('right')

#ax.set_ylim([0,6000])
ax.set_xlabel('UTC Time')
ax.set_ylabel('Distance in km')
plt.legend()
ax.set_title('Distance from Launched Rays on DSX to VPM  \n for 8.2kHz field-aligned ray')
plt.savefig('distvstime.png')
plt.show()
