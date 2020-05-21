
"""
here is a script that will call run_rays and save
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
from TLE_funcs import TLE2pos
from haversine import haversine, Unit
from multiprocessing import Pool, cpu_count
import tempfile, shutil, time

# coordinate mania!
# TLES give us GEI car in km, raytracer needs SM car in m
# and IGRF funcs  B_dir, tracefieldline, and findfootprint need GEO car in RE
# and IGRF func Bdireasy need SM car in RE
# we want all in MAG SPH at end

# -------------------------------- SET TIME --------------------------------
# change time information here - use UTC -
year = 2020
month = 5
day = 28
hours = 0
minutes = 0
seconds = 0

ray_datenum = dt.datetime(year, month, day, hours, minutes, seconds)

# -------------------------------- GET POSITIONS --------------------------------
# these will be in ECI coordinates (GEI) in km
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
plen = 10  # second
r, tvec = TLE2pos(lines1, lines2, satnames, plen, ray_datenum)

# convert to meters
dsx = [rpos*1e3 for rpos in r[0]]
vpm = [rpos*1e3 for rpos in r[1]]

# convert startpoint to SM for raytracer
dsxpos = coord.Coords(dsx, 'GEI', 'car', units=['m', 'm', 'm'])
dsxpos.ticks = Ticktock(tvec, 'UTC') # add ticks
SM_dsx = dsxpos.convert('SM', 'car')

# convert vpm to MAG LLA for distance comp.
vpmpos = coord.Coords(vpm, 'GEI', 'car', units=['m', 'm', 'm'])
vpmpos.ticks = Ticktock(tvec, 'UTC') # add ticks
MAG_vpm = vpmpos.convert('MAG', 'sph')

# -------------------------------- DEFINE RAY DIRECTIONS --------------------------------
positions = np.column_stack((SM_dsx.x, SM_dsx.y, SM_dsx.z))

# need to check what this creates ^ 
freq = [8.2e3] # Hz
directions = []
thetalist = [0, 15, 20, 30, 35, 45, -15, -20, -30, -45]  # in deg -- what angles to launch at? 

rayperpos = int(len(thetalist))

def launchmanyrays(position, rayt):
#for position, rayt in zip(positions, tvec):

    # declare lists to return 
    dfoot = []
    dray = []

    finalpos = []
    # print('position is', position)
    # grab position and find direction of local bfield
    # convert to RE for bfield lib - SM is okay here
    startpoint = [position[0]/R_E, position[1]/R_E, position[2]/R_E]
    print('startpoint is', startpoint)
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
        # make sure same position for EACH direction
        finalpos.append(position)

    # -------------------------------- RUN RAYS --------------------------------
    # convert for raytracer settings
    days_in_the_year = rayt.timetuple().tm_yday
    days_in_the_year = format(days_in_the_year, '03d')

    # yearday and miliseconds day are used by raytracer
    yearday = str(year)+ str(days_in_the_year)   # YYYYDDD
    milliseconds_day = rayt.hour*3.6e6 + rayt.minute*6e4 + rayt.second*1e3

    # position is in GEO meters - is that correct?
    tmpdir = tempfile.mkdtemp() 
    run_rays(freq, finalpos, directions, yearday, milliseconds_day, tmpdir)

    # -------------------------------- LOAD OUTPUT --------------------------------
    # Load all the rayfiles in the output directory
    ray_out_dir = tmpdir
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
        dray.append(np.nan)
        dfoot.append(np.nan)
        return

    # -------------------------------- CONVERT COORDINATES --------------------------------
    # convert to desired coordinate system into vector list rays
    rays = []
    for r in raylist:
        tmp_coords = coord.Coords(list(zip(r['pos'].x, r['pos'].y, r['pos'].z)), 'SM', 'car', units=['m', 'm', 'm'])
        tvec_datetime = [rayt + dt.timedelta(seconds=s) for s in r['time']]
        tmp_coords.ticks = Ticktock(tvec_datetime, 'UTC')  # add ticks
        tmp_coords.sim_time = r['time']
        new_coords = tmp_coords.convert('MAG', 'sph')
        rays.append(new_coords)
        #rays.append(tmp_coords)

    #initialize
    rlat = []
    rlong = []
    rrad = []

    for r in rays:
        rlat.append(r.lati)
        rlong.append(r.long)
        rrad.append(r.radi)
    
    dlist = []
    for d in damplist:
        damp = d["damping"]
        damp = np.squeeze(np.array(damp))
        dlist.append(damp)

    # -------------------------------- GET FOOTPRINT --------------------------------
    # footprint needs GEO car
    dsxstart = coord.Coords([[position[0], position[1], position[2]]], 'SM', 'car', units=['m', 'm', 'm'])
    dsxstart.ticks = Ticktock(rayt, 'UTC') # add ticks
    GEO_dsx = dsxstart.convert('GEO', 'car')
    bstart = [float(GEO_dsx.x/R_E), float(GEO_dsx.y/R_E), float(GEO_dsx.z/R_E)]

    footprint = findFootprints(rayt, bstart, 'same')
    footprint.ticks = Ticktock(rayt, 'UTC')
    footprint = footprint.convert('MAG', 'sph')
    
    # -------------------------------- FIND DISTANCE --------------------------------
    # save a list for EACH?
    # or get average? let's try avg fist 
    foot = (footprint.lati, footprint.long)
    vpm = (MAG_vpm.lati[-1], MAG_vpm.long[-1])

    df = haversine(foot, vpm) # in km
    df = np.abs(df)
    #dfoot.append(df)

    # find ray distance
    dravg = []
    for raylat, raylong in zip(rlat, rlong):
        rayend = (raylat[-1], raylong[-1])
        dr = haversine(rayend, vpm) # in km
        dr = np.abs(dr)
        dravg.append(dr)

    # get average and save
    averageraydist = sum(dravg) / len(dravg)
    #dray.append(averageraydist)

    #print('dr is', dray)
    #print('df is', dfoot)
    
    # clear temp directories
    shutil.rmtree(tmpdir)

    return averageraydist, df, rayt

    #print(' \n \n \n  \n \n \n  \n \n \n  TIME IS:       \n \n \n', rayt, ' \n \n \n \n \n \n')

# parallelization
nmbrcores = cpu_count()
#nmbrcores = nmbrcores//2
lstarg = zip(positions, tvec)

start = time.time()
with Pool(nmbrcores) as p:
    results = p.starmap(launchmanyrays, lstarg)
    #print(results)

end = time.time()
print(f'time is with 2 cores {end-start}')

"""
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
"""