
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

# -------------------------------- SET TIME --------------------------------
# change time information here - use UTC -
year = 2020
month = 6
day = 11
hours = 21
minutes = 58
seconds = 0

ray_datenum = dt.datetime(year, month, day, hours, minutes, seconds)

# -------------------------------- GET POSITIONS --------------------------------
# these will be in ECI coordinates (GEI) in km
# last updated 6/8

# DSX TLE
l11 = '1 44344U 19036F   20159.06535339 -.00000040  00000-0  00000-0 0  9994'
l21 = '2 44344 042.2684 076.3406 1973825 149.5843 223.5230 04.54370887015809'
# VPM TLE
l12 = '1 45120U 19071K   20160.05856933  .00002378  00000-0  82577-4 0  9991'
l22 = '2 45120  51.6439  48.4744 0012296  22.9438 337.2091 15.33961651 19578'

lines1 = [l11, l12]
lines2 = [l21, l22]
satnames = ['DSX', 'VPM']

# get DSX and VPM positions for... 
plen = 1  # seconds
r, tvec = TLE2pos(lines1, lines2, satnames, plen, ray_datenum)

# convert to meters
dsx = [rpos*1e3 for rpos in r[0]]
vpm = [rpos*1e3 for rpos in r[1]]

# convert startpoint to SM car for raytracer
GEIcar_dsx = coord.Coords(dsx, 'GEI', 'car', units=['m', 'm', 'm'])
GEIcar_dsx.ticks = Ticktock(tvec, 'UTC') # add ticks
SMcar_dsx = GEIcar_dsx.convert('SM', 'car')

# convert vpm to MAG sph
GEIcar_vpm = coord.Coords(vpm, 'GEI', 'car', units=['m', 'm', 'm'])
GEIcar_vpm.ticks = Ticktock(tvec, 'UTC') # add ticks
MAGsph_vpm = GEIcar_vpm.convert('MAG', 'sph')

# -------------------------------- DEFINE RAY DIRECTIONS --------------------------------
dsxpositions = np.column_stack((SMcar_dsx.x, SMcar_dsx.y, SMcar_dsx.z))
vpmpositions = np.column_stack((MAGsph_vpm.radi, MAGsph_vpm.lati, MAGsph_vpm.long))

freq = [8.2e3] # Hz

# how many rays? 
rayn = 1e5
thetalist = []

# generate random angles from a sin theta distribution
for i in range(int(rayn)):
    xi = random.random()
    th = np.arccos(1-2*xi)
    pxi = random.random()
    thetalist.append(R2D*th)

# rotate to be defined w respect to B0
thetalist = [th - 90 for th in thetalist]

# save those angles to parse later
fname = str(freq[0]/1e3) + 'kray' + str(ray_datenum) + 'thetalist.txt'
with open(fname, "w") as outfile:
    outfile.write("\n".join(str(item) for item in thetalist))

outfile.close()

# start running rays!

def launchmanyrays(position, vpmpos, rayt):

    # grab position and find direction of local bfield
    SMcar_dsxpos = coord.Coords(position, 'SM', 'car', units=['m', 'm', 'm'])
    SMcar_dsxpos.ticks = Ticktock(rayt)
    GEOcar_dsx = SMcar_dsxpos.convert('GEO', 'car')
    GEOsph_dsx = SMcar_dsxpos.convert('GEO', 'sph')

    # check with hemi we are in
    if GEOsph_dsx.lati > 0:
        dir = 1   # north
    else:
        dir = -1  # south

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

    # rotate directions
    for theta in thetalist:

        # increase (or decrease) polar angle
        newth = float(SMsph_dirB.lati) + theta
        Rot_dirB = [float(SMsph_dirB.radi), newth, float(SMsph_dirB.long)] 
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

    # quick check: did the rays propagate?
    # raylist = [checkray for checkray in raylist if not len(checkray["time"]) < 2]
    #for r in raylist:
    #    if len(r['time']) < 2:
    #        dravg = 'nan'
    #        df = 'nan'
    #        #return dravg, df, rayt

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
    GDZsph_foot = findFootprints(rayt, Bstart, 'same')
    GDZsph_foot.units = ['km', 'deg', 'deg']
    GDZsph_foot.ticks = Ticktock(rayt, 'UTC')
    MAGsph_foot = GDZsph_foot.convert('MAG', 'sph')

    # -------------------------------- FIND DISTANCE --------------------------------
    foot = (float(MAGsph_foot.lati), float(MAGsph_foot.long))
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

    return rayendls, vpmf, foot, rayt

# -------------------------------- RUN --------------------------------

# parallel
nmbrcores = cpu_count()
lstarg = zip(dsxpositions, vpmpositions, tvec)

start = time.time()
with Pool(nmbrcores) as p:
    results = p.starmap(launchmanyrays, lstarg)

end = time.time()
# print(f'time is with 2 cores {end-start}')

fname = str(freq[0]/1e3) + 'kray' + str(ray_datenum.year) + str(ray_datenum.month) + str(ray_datenum.day) + str(ray_datenum.hour) + str(ray_datenum.minute) + 'MCsim.txt'
with open(fname, "w") as outfile:
    outfile.write("\n".join(str(item) for item in results))

outfile.close()

# -------------------------------- END --------------------------------