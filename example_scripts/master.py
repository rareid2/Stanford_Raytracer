# master prediction run for a week -- will get updated TLEs from time of running
import numpy as np
import datetime as dt
from haversine import haversine, Unit
from DSXlogparser import runsomerays, fullweeksetup, readDSXlog

import sys
sys.path.insert(1, '/home/rileyannereid/workspace/scratches/')
from rayoutput_parser import readrayout

###################################### STEP 1 ######################################
# check the Kp index!

# full week of predictions
startday = 21
mo = 9
weeklen = 2 # add 1 more than you think (8 days for a week)
clist = 0

dates, fs, bs = fullweeksetup(startday, mo, clist, weeklen)

# set opts
rayn = 1
plen = 86400
timeint = 100
MCsim = 0
angfiles = 0

runopts = [rayn, plen, timeint, MCsim, angfiles]
#runsomerays(dates, fs, bs, runopts)
print('finished running full day')

###################################### STEP 2 ######################################

# check the possible actual conjunctions 
pday = 0
for cdate, cf, bsstr in zip(dates, fs, bs):
    #  ---------- some file set up ----------
    datadir = '/home/rileyannereid/workspace/SR-output/' + bsstr + '/'
    ray_datenum = cdate # get the raydatenum
    rename = str(ray_datenum.month) + str(ray_datenum.day) + str(ray_datenum.year)
    fname = datadir + rename + '/' + rename + bsstr + '.txt'
    #  ---------- some file set up ----------

    # save the possible conjunctions
    day = ray_datenum.day
    if np.abs(day - pday) > 0:
        if pday!= 0: 
            outfile.close() 
        # open txt file
        fname2 = datadir + rename + '/' + rename + 'conjlist.txt'
        outfile = open(fname2,'w')

    # get the output data (all of it)
    dataout = readrayout(fname, 0)
    raylatsdict, raylonsdict, vpm, foot, tvec = dataout
    lastt = 0

    df = []
    for latiray, longray, vpmpos, footpos in zip(raylatsdict.values(), raylonsdict.values(), vpm, foot):
        df.append(np.abs(haversine(footpos, vpmpos)))

    for dffi, dff in enumerate(df):
        checkd = tvec[dffi].hour*3600 + tvec[dffi].minute*60 + tvec[dffi].second
        if dff < 800: # try to keep it loose enough to find opposite hemisphere cases -- look into this
            if np.abs(lastt - checkd) > 30*60: # new conj!
                outfile.write(str(tvec[dffi]))
                outfile.write("\n")
                lastt = checkd

    pday = day
print('finished conjuring the conjunctions')

###################################### STEP 3 ######################################

# FINALLY narrow down the output and run some longer ones!
clist = 1
dates, fs, bs = fullweeksetup(startday, mo, clist, weeklen)

# set opts
rayn = 5000
plen = 30
timeint = 30
MCsim = 1
angfiles = 1

dates = readDSXlog('DSXlogs.txt')

fsl = []
fsh = []
for dd in dates:
    fsl.append(28e3)
    fsh.append(8.2e3)

fs = fsl.extend(fsh)
bs = []

dates = dates.extend(dates)
for dd in dates:
    bs.append('additional')

print('est. run time: ', len(dates) * rayn * plen/timeint * (1/3000), ' min')

runopts = [rayn, plen, timeint, MCsim, angfiles] # FULL runs -- takes a WHILE ! 3-4hrs (start the night before or early in the morn!)
runsomerays(dates, fs, bs, runopts)

################################# MOVE TO PLOTTING #################################