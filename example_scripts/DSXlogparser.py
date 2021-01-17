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

###############################################################################################
def readDSXlog(fnameDSX):
    dates = []
    
    infile = open(fnameDSX,'r')

    # goes thru line by line
    for line in infile:
        out = line
        year = 2020
        day = int(out[0:2])

        if out[3] == 'S':
            month = 9
        else:
            month = 10

        hour = int(out[12:14])
        minute = int(out[15:17])

        dates.append(dt.datetime(year, month, day, hour, minute))
    infile.close()
    return dates

###############################################################################################

def readconjlog(fnameconj):
    condates = []
    chdatel = 0
    # goes thru line by line
    #for line in lines:
    with open(fnameconj) as f:
        for line in f: 
            out = line
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
                condates.append(dt.datetime(year, month, day, hour, minute))
            chdatel = chdate
    f.close()
    return condates

#dd = readconjlog('/home/rileyannereid/workspace/SR-output/fullday/8312020/8312020conjlist.txt')
#print(dd)
###############################################################################################

def readcloseconjlog(fnameconj):
    dates = []
    freq = []
    chdatel = 0

    # goes thru line by line
    for line in infile:
        out = line

        year = 2020
            
        month = int(out[5:7])
        day = int(out[8:10])
        hour = int(out[11:13])
        minute = int(out[14:16])
        second = int(out[17:19])
        cf = float(out[20:23])
        
        # convert to fractions of a day
        minphr = 60
        
        chmin = minute/minphr
        chdate = hour + chmin

        if np.abs(chdate - chdatel) > 30/minphr: # if more than 30 minutes apart, this is a new conjunction!
            dates.append(dt.datetime(year, month, day, hour, minute))
            freq.append(cf * 10**3)
        chdatel = chdate

    infile.close()
    return dates, freq

###############################################################################################
def readTNTlog(fnameTNT):
    #Timestamp[SUT] Duration StartFreq StopFreq +Vmax -vMax +FRes -FRes
    infile = open(fnameTNT,'r')

    tx_dict = {}
    tx_n = 0

    # goes thru line by line
    for line in infile:
        out = line
        if out[0]=='#':
            pass
        else:
            tx_n += 1
            tx_n_str = 'tx ' + str(tx_n)
            tx_dict[tx_n_str] = {}

            tx_dict[tx_n_str]['timestamp'] = dt.datetime.strptime(out[:23], "%Y-%m-%d %H:%M:%S.%f")
            tx_dict[tx_n_str]['duration'] = int(out[28:33])
            tx_dict[tx_n_str]['startfreq'] = int(out[37:45])
            tx_dict[tx_n_str]['stopfreq'] = int(out[46:51])
            
    infile.close()
    return tx_dict

###############################################################################################

def fullweeksetup(weekstart, mo, clist, weeklen):

    stdates = [dt.datetime(2020, mo, d, 0, 0) for d in range(weekstart, weekstart+weeklen)]
    condtime = []
    fs = []
    bs = []

    if clist == 1: # for specific conjs
        for cdate in stdates:
            ray_datenum = cdate
            datadir = '/home/rileyannereid/workspace/SR-output/' + 'fullday' + '/'
            datadir = datadir + str(ray_datenum.month) + str(ray_datenum.day) + str(ray_datenum.year) + '/'
            fname = datadir + str(ray_datenum.month) + str(ray_datenum.day) + str(ray_datenum.year) + 'conjlist.txt'
            newdates = readconjlog(fname) # for specific conjs

            condtime.extend(newdates)
            for ii in newdates:
                fs.append(28e3)
                bs.append('fullday')
            for ii in newdates:
                fs.append(8.2e3)
                bs.append('fullday')
            condtime.extend(newdates)

        dates = condtime
    if clist == 0:
        dates = stdates
        fs = [28e3 for d in range(weekstart, weekstart+weeklen)] # arbitrary freq
        bs = ['fullday' for d in range(weekstart, weekstart+weeklen)]
    
    return dates, fs, bs

###############################################################################################

def runsomerays(dates, fs, bs, runopts):
    for cdate, cf, bsstr in zip(dates, fs, bs):
        year = cdate.year
        month = cdate.month
        day = cdate.day
        hours = cdate.hour
        minutes = cdate.minute
        seconds = cdate.second

        rayn = runopts[0]
        plen = runopts[1] #seconds
        time_int = runopts[2]
        MCsim = runopts[3] # if not, rays are random
        anglefiles = runopts[4] # 1 for yes

        datadir = '/home/rileyannereid/workspace/SR-output/' + bsstr + '/'
        crs_out = 'GEO' #what coord sys to save in?

        freq = [cf]
        ray_datenum = dt.datetime(year, month, day, hours, minutes, seconds) # get the raydatenum
        modeldump(year, month, day, hours, minutes, seconds) # run model dump to update plasmasphere

        # file setup
        datadir = datadir + str(ray_datenum.month) + str(ray_datenum.day) + str(ray_datenum.year) + '/'

        try:
            os.mkdir(datadir)
        except OSError:
            pass
        else:
            pass

        if MCsim == 1:
            datadir = datadir + str(freq[0]/1e3) + 'kHz' + str(ray_datenum.month) + str(ray_datenum.day) + str(ray_datenum.year) + str(ray_datenum.hour) + str(ray_datenum.minute) + '/'
            try:
                os.mkdir(datadir)
            except OSError:
                pass
            else:
                pass

        opts = [freq, rayn, plen, time_int, datadir, crs_out, MCsim, anglefiles]

        # run funcs 
        dsxpositions, vpmpositions, tvec = getpos(ray_datenum, opts)
        parallelrun(dsxpositions, vpmpositions, tvec, opts)

###############################################################################################

# just a test run
#dates = [dt.datetime(2020,7,25,0,0,0)]

#fs = [8.2e3]
#bs = ['fullday' for i in range(len(dates))]
#runsomerays(dates, fs, bs, [1000, 10*60, 30, 1, 1])


#dates, fs, bs = readDSXlog('DSXlogs.txt')
# set opts
#rayn = 100
#plen = 30*60
#timeint = 100
#MCsim = 1
#angfiles = 1

#runopts = [rayn, plen, timeint, MCsim, angfiles]
#runsomerays(dates, fs, bs, runopts)

###############################################################################################