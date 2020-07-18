
"""
here is a script that will call run_rays and save
"""

# import needed packages
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import datetime as dt
from random import random
from raytracer_utils import readdump, read_rayfile, read_rayfiles, read_damp
from run_rays import run_rays
from raytracer_settings import *
from IGRF_funcs import B_dir, trace_fieldline_ODE, findFootprints, B_direasy
from spacepy import coordinates as coord
from spacepy.time import Ticktock
from TLE_funcs import TLE2pos
from haversine import haversine, Unit
from multiprocessing import Pool, cpu_count
import tempfile, shutil, time, pickle

# -------------------------------- EASY SETTINGS--------------------------------
# how long are we launching rays for?? IN SECONDS
plen = 30*60  # seconds
# how often to send the rays? for every second, set to 0
time_int = 100 
# how many rays per launch?
rayn = 1000
# what is the freq of the rays?
freq = [28e3] # Hz
# are we doing a MC sim?
MCsim = 1 # if not, the rays will just be random
crs_out = 'GEO' #what coord sys to save in? 
# where to save data (directory)? 
datadir = '/home/rileyannereid/workspace/SR-output/'
# do we want theta and phi output? 
anglefiles = 1 # 1 for yes
# finally... set start time!

# -------------------------------- SET TIME --------------------------------
# change time information here - use UTC -
year = 2020
month = 6
day = 12
hours = 8
minutes = 40
seconds = 0

ray_datenum = dt.datetime(year, month, day, hours, minutes, seconds)

# file setup
datadir = datadir + str(freq[0]/1e3) + 'kHz' + str(ray_datenum.month) + str(ray_datenum.day) + str(ray_datenum.year) + '/'
try:
    os.mkdir(datadir)
except OSError:
    print ("Creation of the directory %s failed" % datadir)
else:
    print ("Successfully created the directory %s" % datadir)

if MCsim == 1:
    datadir = datadir + str(freq[0]/1e3) + 'kHz' + str(ray_datenum.month) + str(ray_datenum.day) + str(ray_datenum.year) + str(ray_datenum.hour) + str(ray_datenum.minute) + '/'
    try:
        os.mkdir(datadir)
    except OSError:
        print ("Creation of the directory %s failed" % datadir)
    else:
        print ("Successfully created the directory %s" % datadir)
 
# -------------------------------- GET POSITIONS --------------------------------
r, tvec = TLE2pos(plen, ray_datenum)

# convert to meters
dsx = [rpos*1e3 for rpos in r[0]]
vpm = [rpos*1e3 for rpos in r[1]]

# convert startpoint to SM car for raytracer
GEIcar_dsx = coord.Coords(dsx, 'GEI', 'car', units=['m', 'm', 'm'])
GEIcar_dsx.ticks = Ticktock(tvec, 'UTC') # add ticks
SMcar_dsx = GEIcar_dsx.convert('SM', 'car')

# convert vpm to out coords sph
GEIcar_vpm = coord.Coords(vpm, 'GEI', 'car', units=['m', 'm', 'm'])
GEIcar_vpm.ticks = Ticktock(tvec, 'UTC') # add ticks
outsph_vpm = GEIcar_vpm.convert(crs_out, 'sph')

# -------------------------------- DEFINE RAY POS --------------------------------
dsxpositions = np.column_stack((SMcar_dsx.x, SMcar_dsx.y, SMcar_dsx.z))
vpmpositions = np.column_stack((outsph_vpm.radi, outsph_vpm.lati, outsph_vpm.long))

if time_int > 0:
    dsxpositions = dsxpositions[0::time_int]
    vpmpositions = vpmpositions[0::time_int]
    tvec = tvec[0::time_int]

# start running some loops!

def launchmanyrays(position, vpmpos, rayt):

    # grab position and find direction of local bfield
    SMcar_dsxpos = coord.Coords(position, 'SM', 'car', units=['m', 'm', 'm'])
    SMcar_dsxpos.ticks = Ticktock(rayt)
    GEOcar_dsx = SMcar_dsxpos.convert('GEO', 'car')

    # check with hemi we are in
    # lets just go to whichever hemisphere VPM is

    if vpmpos[1] > 0:
        dir = 1   # north
        dirstr = 'north'
    else:
        dir = -1  # south
        dirstr = 'south'

    Bstart = [float(GEOcar_dsx.x) / R_E, float(GEOcar_dsx.y) / R_E, float(GEOcar_dsx.z) / R_E]
    Bx, By, Bz = B_direasy(rayt, Bstart, dir)

    # convert to SM coordinates for raytracer
    dirB = np.reshape(np.array([Bx, By, Bz]), (1, 3))
    dirB = coord.Coords(dirB[0], 'GEO', 'car', units=['Re', 'Re', 'Re'])
    dirB.ticks = Ticktock(rayt, 'UTC') # add ticks
    SMsph_dirB = dirB.convert('SM', 'sph')

    # fill for raytracer call
    positions = []
    directions = []

    # -------------------------------- DEFINE RAY DIRECS --------------------------------
    thetalist = []
    philist = []
    if MCsim == 1:
        # generate random polar angles from a sin theta distribution
        for i in range(int(rayn)):
            xi = random.random()
            th = np.arccos(1-2*xi)
            pxi = random.random()
            thetalist.append(R2D*th)
        # generate random azimuth angles uniformly btwn 0 and pi (forward hemi)
        for i in range(int(rayn)):
            ph = np.pi * random.random()
            philist.append(R2D*ph)
    else: 
        # genereate random between 0 and pi
        for i in range(int(rayn)):
            th = np.pi * random.random()
            ph = np.pi * random.random()
            thetalist.append(R2D*th)
            philist.append(R2D*ph)

    # rotate to be defined w respect to B0
    thetalist = [th - 90 for th in thetalist]
    philist = [ph - 90 for ph in philist]

    # save those angles to parse later
    
    if anglefiles == 1:
        fname = datadir + str(freq[0]/1e3) + 'kHz' + str(rayt.month) + str(rayt.day) + str(rayt.year) + str(rayt.hour) + str(rayt.minute) + str(rayt.second) + 'thetalist.txt'
        with open(fname, "w") as outfile:
            outfile.write("\n".join(str(item) for item in thetalist))
        outfile.close()

        fname = datadir + str(freq[0]/1e3) + 'kHz' + str(rayt.month) + str(rayt.day) + str(rayt.year) + str(rayt.hour) + str(rayt.minute) + str(rayt.second) + 'philist.txt'
        with open(fname, "w") as outfile:
            outfile.write("\n".join(str(item) for item in philist))
        outfile.close()

    # rotate directions
    for theta, phi in zip(thetalist, philist):

        # increase (or decrease) polar angle
        newth = float(SMsph_dirB.lati) + theta
        newphi = float(SMsph_dirB.long) + phi
        Rot_dirB = [float(SMsph_dirB.radi), newth, newphi] 
        Rot_dirB = coord.Coords(Rot_dirB, 'SM', 'sph')
        Rot_dirB.ticks = Ticktock(rayt, 'UTC') # add ticks
        SMcar_dirB = Rot_dirB.convert('SM', 'car')

        direction = [float(SMcar_dirB.x), float(SMcar_dirB.y), float(SMcar_dirB.z)]
        
        # add the normalized direction (or zeros)
        directions.append(np.squeeze(direction))

        # make sure position list matches direction list
        positions.append(position)
   
    # -------------------------------- RUN RAYS --------------------------------
    # convert for raytracer settings
    days_in_the_year = rayt.timetuple().tm_yday
    days_in_the_year = format(days_in_the_year, '03d')

    # yearday and miliseconds day are used by raytracer
    yearday = str(year)+ str(days_in_the_year)   # YYYYDDD
    milliseconds_day = rayt.hour*3.6e6 + rayt.minute*6e4 + rayt.second*1e3

    # run it!
    tmpdir = tempfile.mkdtemp() 
    run_rays(freq, positions, directions, yearday, milliseconds_day, tmpdir)

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

    # -------------------------------- CONVERT COORDINATES --------------------------------
    # convert to desired coordinate system into vector list rays
    rays = []

    for r in raylist:
        tmp_coords = coord.Coords(list(zip(r['pos'].x, r['pos'].y, r['pos'].z)), 'SM', 'car', units=['m', 'm', 'm'])
        tvec_datetime = [rayt + dt.timedelta(seconds=s) for s in r['time']]
        tmp_coords.ticks = Ticktock(tvec_datetime, 'UTC')  # add ticks
        tmp_coords.sim_time = r['time']
        new_coords = tmp_coords.convert(crs_out, 'sph')
        rays.append(new_coords)

    #initialize
    rrad = []
    rlat = []
    rlong = []

    for r in rays:
        rrad.append(r.radi)
        rlat.append(r.lati)
        rlong.append(r.long)
        
    dlist = []
    for d in damplist:
        damp = d["damping"]
        damp = np.squeeze(np.array(damp))
        dlist.append(damp)

    # -------------------------------- GET FOOTPRINT --------------------------------
    # also in GEO car, so need to use bstart 
    GDZsph_foot = findFootprints(rayt, Bstart, dirstr)
    GDZsph_foot.units = ['km', 'deg', 'deg']
    GDZsph_foot.ticks = Ticktock(rayt, 'UTC')
    outsph_foot = GDZsph_foot.convert(crs_out, 'sph')

    # -------------------------------- CLEANUP --------------------------------
    foot = (float(outsph_foot.lati), float(outsph_foot.long))
    vpmf = (vpmpos[1], vpmpos[2])

    # find rayend points
    rayendls = []
    for raylat, raylong in zip(rlat, rlong):
        if len(raylat) > 2:
            rayend = (raylat[-1], raylong[-1])
            rayendls.append(rayend)
        else:
            zeros = (0,0)
            rayendls.append(zeros)

    # clear temp directories
    shutil.rmtree(tmpdir)

    rayt = rayt.strftime("%Y-%m-%d %H:%M:%S")
    print(rayt) # fix this

    return rayendls, vpmf, foot, rayt

# -------------------------------- PARALLELIZEEEE --------------------------------

# parallel
nmbrcores = cpu_count()
lstarg = zip(dsxpositions, vpmpositions, tvec)

start = time.time()
with Pool(nmbrcores) as p:
    results = p.starmap(launchmanyrays, lstarg)

end = time.time()
# print(f'time is with 2 cores {end-start}')

if MCsim == 1:
    addon = 'MCsim'
else:
    addon = ''

fname = datadir + str(freq[0]/1e3) + 'kHz' + str(ray_datenum.month) + str(ray_datenum.day) + str(ray_datenum.year) + str(ray_datenum.hour) + str(ray_datenum.minute) + addon + '.txt'
with open(fname, "w") as outfile:
    outfile.write("\n".join(str(item) for item in results))
outfile.close()

# -------------------------------- END --------------------------------