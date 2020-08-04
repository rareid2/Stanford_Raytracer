
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

from haversine import haversine, Unit
from multiprocessing import Pool, cpu_count
import tempfile, shutil, time, pickle
import scipy
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as mcolors
from matplotlib import cm # Colormaps
import cartopy.crs as ccrs
import pandas as pd
import os

from run_model_dump import modeldump

# where do we ACTUALLY need to launch from? 
rayt = dt.datetime(2020, 6, 14, 6, 59)
GEOsph_wsmrpt = coord.Coords([R_E, 32.943896, -106.41965], 'GEO', 'sph', units=['m', 'deg', 'deg'])
GEOsph_wsmrpt.ticks = Ticktock(rayt, 'UTC')
GEOcar_wsmrpt = GEOsph_wsmrpt.convert('GEO', 'car')
GEOcar_wsmrpt = [float(GEOcar_wsmrpt.x/R_E), float(GEOcar_wsmrpt.y/R_E), float(GEOcar_wsmrpt.z/R_E)]

Blines = []
Blines.append(trace_fieldline_ODE(GEOcar_wsmrpt, 0, '0', -1, rayt))

crs_out = 'GEO'
for blinex, bliney, blinez in Blines:
    raytime = []
    # create list of the same times
    for m in range(int(len(blinex))):
        raytime.append(rayt)
    
    # convert coords
    bpos = np.column_stack((blinex, bliney, blinez))
    bpos = coord.Coords(bpos, 'GEO', 'car', units=['Re', 'Re', 'Re'])

    bpos.ticks = Ticktock(raytime, 'UTC') # add ticks
    outsph_bline = bpos.convert(crs_out, 'sph')

for outbr_ind, outbr in enumerate(outsph_bline.radi):
    if np.abs((outbr * R_E) - (1000e3 + R_E)) < 100e3:
        wsmrstart = outsph_bline[outbr_ind]
        break

#print(wsmrstart) # got it! 

# 2 times to look at 6-14-2020 6:59 and 6-14-2020 11:54

# settings
datadir = '/home/rileyannereid/workspace/SR-output/WSMR/'

# define a WSMR at 1000 km up traced along fieldline from WSMR
GEOsph_wsmrpt = coord.Coords([float(wsmrstart.radi*R_E), float(wsmrstart.lati), float(wsmrstart.long)], 'GEO', 'sph', units=['m', 'deg', 'deg'])
GEOsph_wsmrpt.ticks = Ticktock(rayt, 'UTC')
SMcar_wsmrpt = GEOsph_wsmrpt.convert('SM', 'car')
GEOcar_wsmrpt = GEOsph_wsmrpt.convert('GEO', 'car')
Bstart = [float(GEOcar_wsmrpt.x) / R_E, float(GEOcar_wsmrpt.y) / R_E, float(GEOcar_wsmrpt.z) / R_E]

K = np.array([[1,0], [0,1]])
m = np.array([0, 0]).reshape(2, 1)
d = 2
n = int(1e4)
z = np.random.multivariate_normal(mean=m.reshape(d,), cov=K, size=n)
y = np.transpose(z)
dx = y[0] * 1000e3
dy = y[1] * 1000e3

# actually correct way to do this ?
#R = np.array([[1, 0, 0], [0, np.cos(-D2R * GEOsph_wsmrpt.lati), -np.sin(-D2R * GEOsph_wsmrpt.lati)], [0, np.sin(-D2R * GEOsph_wsmrpt.lati), np.cos(-D2R * GEOsph_wsmrpt.lati)]])
#scolist = []
#for x, y in zip(dx, dy): 
#    scc = np.array([[x], [y], [0]])
#    sc = np.matmul(R, scc)
#    # now translate
#    sco = (float(sc[0]), float(sc[1]), float(sc[2] - GEOcar_wsmrpt.z))
#    scolist.append(sco)
##sc = coord.Coords(scolist, 'GEO', 'car', units = ['m', 'm', 'm'])
#sc.ticks = Ticktock([rayt for i in np.ones(len(dx))], 'UTC')
#ndone = sc.convert('GEO', 'sph')


newc = [(float(GEOcar_wsmrpt.x + sx), float(GEOcar_wsmrpt.y + sy), float(GEOcar_wsmrpt.z)) for sx, sy in zip(dx, dy)]
sc = coord.Coords(newc, 'GEO', 'car', units = ['m', 'm', 'm'])
sc.ticks = Ticktock([rayt for i in np.ones(len(dx))], 'UTC')
ndone = sc.convert('GEO', 'sph')
SMpos = sc.convert('SM', 'car')

vpm = (1,1)

freq = [26.1e3]


# PLOTTING THE LAUNCCCH
fig,ax=plt.subplots(1,1)
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
ax.hist2d(ndone.long, ndone.lati, bins=30, cmap = 'BuPu')
ax.set_extent([-130, -80, 5, 55])
ax.gridlines(color="black", linestyle="dotted", draw_labels=True)
plt.title('Launch from WSMR' + str(rayt))
#plt.show()
plt.savefig(datadir + str(freq[0]/1e3) + str(rayt) + '.png')
plt.close()

# launch rays adapted for this
modeldump(rayt.year, rayt.month, rayt.day, rayt.hour, rayt.minute, rayt.second)


dir = -1  # south
dirstr = 'south'

# fill for raytracer call
cleanpos = [(float(pos.x), float(pos.y), float(pos.z)) for pos in SMpos]
positions = cleanpos # fill from pdf
directions = [np.zeros(3) for i in range(n-1)]

# -------------------------------- RUN RAYS --------------------------------
# convert for raytracer settings
days_in_the_year = rayt.timetuple().tm_yday
days_in_the_year = format(days_in_the_year, '03d')

# yearday and miliseconds day are used by raytracer
yearday = str(rayt.year)+ str(days_in_the_year)   # YYYYDDD
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

for filename in file_titles:
    if '.ray' in filename:
        raylist += read_rayfile(os.path.join(ray_out_dir, filename))

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

# -------------------------------- GET FOOTPRINT --------------------------------
# also in GEO car, so need to use bstart 
GDZsph_foot = findFootprints(rayt, Bstart, dirstr)
GDZsph_foot.units = ['km', 'deg', 'deg']
GDZsph_foot.ticks = Ticktock(rayt, 'UTC')
outsph_foot = GDZsph_foot.convert(crs_out, 'sph')

# -------------------------------- CLEANUP --------------------------------
foot = (float(outsph_foot.lati), float(outsph_foot.long))

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

results = rayendls, vpm, foot, rayt

fname = datadir + str(freq[0]/1e3) + 'kHz' + str(rayt.month) + str(rayt.day) + str(rayt.year) + str(rayt.hour) + str(rayt.minute) + '.txt'
with open(fname, "w") as outfile:
    outfile.write("\n".join(str(item) for item in results))
outfile.close()
