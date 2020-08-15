
# RUN IT

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

from launchmultiplerays import parallelrun, getpos

def readDSXlog(fnameDSX):
    dates = []
    fs = []
    bs = []

    infile = open(fnameDSX,'r')

    # goes thru line by line
    for line in infile:
        out = line
        year = 2020
        if out[0] == 'b':
            bs.append('burst')
        elif out[0] == 's':
            bs.append('survey')
        elif out[0] =='F':
            bs.append('survey')
            
        month = int(out[8])
        if out[10] == '0':
            day = int(out[11])
        else:
            day = int(out[10:12])

        if out[18] == '0':
            hour = int(out[19])
        else:
            hour = int(out[18:20])
        
        minute = int(out[21:23])
        dates.append(dt.datetime(year, month, day, hour, minute))

        if out[30] == '8':
            fs.append(8.2e3)
        elif out[30] == '2':
            fs.append(28e3)
        elif out[30] == 'H':
            fs.append(25e3)

    infile.close()
    return dates, fs, bs

###############################################################################################

def readconjlog(fnameconj):
    dates = []

    infile = open(fnameconj,'r')
    chdatel = 0

    # goes thru line by line
    for line in infile:
        out = line
        # skip if empty
        if not out.strip():
            continue

        year = 2020
            
        month = int(out[5:7])
        day = int(out[8:10])
        hour = int(out[11:13])
        minute = int(out[14:16])
        second = int(out[17:19])
        
        # convert to fractions of a day
        minphr = 60
        
        chmin = minute/minphr
        chdate = hour + chmin

        if np.abs(chdate - chdatel) > 30/minphr: # if more than 30 minutes apart, this is a new conjunction!
            dates.append(dt.datetime(year, month, day, hour, minute))
        chdatel = chdate

    infile.close()
    return dates


###############################################################################################

def fullweeksetup(weekstart, mo, clist):

    stdates = [dt.datetime(2020, mo, d, 0, 0) for d in range(weekstart, weekstart+8)]

    condtime = []
    fs = []
    bs = []

    if clist == 1: 
        for cdate in stdates:
            
            year = cdate.year
            month = cdate.month
            day = cdate.day
            hours = cdate.hour
            minutes = cdate.minute
            seconds = cdate.second
            ray_datenum = cdate

            datadir = '/home/rileyannereid/workspace/SR-output/' + 'fullday' + '/'
            datadir = datadir + str(ray_datenum.month) + str(ray_datenum.day) + str(ray_datenum.year) + '/'
            fname = datadir + str(ray_datenum.month) + str(ray_datenum.day) + str(ray_datenum.year) + 'conjlist.txt'
            newdates = readconjlog(fname)

            condtime.extend(newdates)
            for ii in newdates:
                fs.append(28e3)
                bs.append('fullday')
            for ii in newdates:
                fs.append(8.2e3)
                bs.append('fullday')
            condtime.extend(newdates)
        dates = condtime
    else:
        dates = stdates
        fs = [28e3 for d in range(weekstart, weekstart+8)] # arbitrary freq
        bs = ['fullday' for d in range(weekstart, weekstart+8)]

###############################################################################################

# just a test run
dates = [dt.datetime(2020,7,25,0,0,0)]

fs = [8.2e3]
bs = ['fullday' for i in range(len(dates))]

###############################################################################################


for cdate, cf, bsstr in zip(dates, fs, bs):
    year = cdate.year
    month = cdate.month
    day = cdate.day
    hours = cdate.hour
    minutes = cdate.minute
    seconds = cdate.second

    rayn = 100
    plen = 10*60  # seconds
    time_int = 100

    datadir = '/home/rileyannereid/workspace/SR-output/' + bsstr + '/'
    crs_out = 'GEO' #what coord sys to save in?
    MCsim = 1 # if not, the rays will just be random
    anglefiles = 1 # 1 for yes

    freq = [cf]
    ray_datenum = dt.datetime(year, month, day, hours, minutes, seconds) # get the raydatenum
    #modeldump(year, month, day, hours, minutes, seconds) # run model dump to update plasmasphere


    # file setup
    datadir = datadir + str(ray_datenum.month) + str(ray_datenum.day) + str(ray_datenum.year) + '/'

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

    opts = [freq, rayn, plen, time_int, datadir, crs_out, MCsim, anglefiles]

    # run funcs 
    dsxpositions, vpmpositions, tvec = getpos(ray_datenum, opts)
    parallelrun(dsxpositions, vpmpositions, tvec, opts)