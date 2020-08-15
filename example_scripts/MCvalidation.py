
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

from run_model_dump import modeldump

import seaborn as sns
sns.set(style="whitegrid")


# start running some loops!

def launchmanyrays(position, vpmpos, rayt, opts):

    # unpack opts
    freq = opts[0] # has to be list format
    rayn = opts[1]
    datadir = opts[4]
    crs_out = opts[5]
    MCsim = opts[6]
    anglefiles = opts[7]

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
    year = rayt.year
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
        if r.radi[-1] < (minalt + 1e3): # get rid of rays that didnt make it
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
        rayend = (raylat[-1], raylong[-1])
        if len(raylat) > 2 and np.sign(vpmf[0]) == np.sign(rayend[0]):
            rayendls.append(rayend)
        else:
            rayendls.append(np.zeros(2))
    #avgdist = np.sum(dr) / len(dr)
    # clear temp directories
    shutil.rmtree(tmpdir)

    #rayt = rayt.strftime("%Y-%m-%d %H:%M:%S")
    #print(rayt) # fix this
    #print(rn)

    return rayendls


# -------------------------------- GET POSITIONS --------------------------------
def getpos(ray_datenum, opts):
    plen = opts[2]
    time_int = opts[3]
    crs_out = opts[5]

    r, tvec = TLE2pos(plen, ray_datenum, 1)

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
        
    return dsxpositions, vpmpositions, tvec

###############################################################################################

# just a test run
dates = [dt.datetime(2020,7,25,0,6,30)]

fs = [8.2e3]
bs = ['fullday' for i in range(len(dates))]

###############################################################################################

rayns = np.ones(int(1e5))

for cdate, cf, bsstr in zip(dates, fs, bs):
    year = cdate.year
    month = cdate.month
    day = cdate.day
    hours = cdate.hour
    minutes = cdate.minute
    seconds = cdate.second

    rayn = 1
    plen = 1  # seconds
    time_int = 0

    datadir = '/home/rileyannereid/workspace/SR-output/' + bsstr + '/'
    crs_out = 'GEO' #what coord sys to save in?
    MCsim = 1 # if not, the rays will just be random
    anglefiles = 1 # 1 for yes

    freq = [cf]
    ray_datenum = dt.datetime(year, month, day, hours, minutes, seconds) # get the raydatenum
    #modeldump(year, month, day, hours, minutes, seconds) # run model dump to update plasmasphere


    # file setup
    datadir = datadir + str(ray_datenum.month) + str(ray_datenum.day) + str(ray_datenum.year) + '/'

    #try:
    #    os.mkdir(datadir)
    #except OSError:
        #print ("Creation of the directory %s failed" % datadir)
    #else:
        #print ("Successfully created the directory %s" % datadir)

    if MCsim == 1:
        datadir = datadir + str(freq[0]/1e3) + 'kHz' + str(ray_datenum.month) + str(ray_datenum.day) + str(ray_datenum.year) + str(ray_datenum.hour) + str(ray_datenum.minute) + '/'
    #    try:
    #        os.mkdir(datadir)
        #except OSError:
            #print ("Creation of the directory %s failed" % datadir)
        #else:
            #print ("Successfully created the directory %s" % datadir)

    opts = [freq, rayn, plen, time_int, datadir, crs_out, MCsim, anglefiles]
    dsxpositions, vpmpositions, tvec = getpos(ray_datenum, opts)
    dsxpositions = dsxpositions[0]
    vpmpositions = vpmpositions[0]
    tvec = tvec[0]
    # run funcs 
"""   
dataout = []
for rni, rn in enumerate(rayns):
    print(rni)
    
    lolout = launchmanyrays(dsxpositions, vpmpositions, tvec, opts)
    dataout.append(lolout)

fname = datadir + str(8.2e3/1e3) + 'kHz' + str(ray_datenum.month) + str(ray_datenum.day) + str(ray_datenum.year) + str(ray_datenum.hour) + '6MCim.txt'
with open(fname, "w") as outfile:
    outfile.write("\n".join(str(item) for item in dataout))
outfile.close()
print('done')
"""

fname = datadir + str(8.2e3/1e3) + 'kHz' + str(ray_datenum.month) + str(ray_datenum.day) + str(ray_datenum.year) + str(ray_datenum.hour) + '6MCim.txt'

# open it
infile = open(fname,'r')
mydist = []

myvpm = (vpmpositions[1], vpmpositions[2])
for line in infile:
    #print(line)
    out = line
    out = out.replace('(', "")
    out = out.replace(')', "")
    out = out.replace('[', "")
    out = out.replace(']', "")
    #out = out.replace('nan', "0")
    out = out.replace(' ', "")

    out = out.rstrip()
    if not out: 
        continue
    #print(out)
    if out[0] == 'a':
        continue
    else:
        out = out.split(",")
        x = out
        raypos = (float(x[0]), float(x[1]))
        mydist.append(np.abs(haversine(raypos, myvpm)))


myx = [i for i in range(len(mydist))]
myy = []
csum = 0
for mi, md in enumerate(mydist):
    n = mi + 1
    csum += md
    myy.append(csum/n)


infile.close()

bb = 1
plt.plot(myx, myy)
plt.plot(myx, np.ones(len(myx)) * myy[-1]*(1 + bb/100), 'r--')
plt.plot(myx, np.ones(len(myx)) * myy[-1]*(1 - bb/100), 'r--')
plt.xlabel('number of rays')
plt.ylabel('avg dist to vpm [km]')
plt.show()

print('done')
