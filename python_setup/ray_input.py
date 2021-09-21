import numpy as np
import datetime as dt
import os
import json

# import constants and functions from this folder
from constants import *
from convert_coords import convert2
from helper_funcs import get_yearmiliday,getBdir,parallel_run_rays,read_rayfile

# Opening JSON file
cwd = os.getcwd()
print(cwd)
f = open(cwd+'/config.json')
# returns JSON object as a dictionary
settings = json.load(f)

#  -------------------------- import settings from config file -----------------------------
raytracer_dir = settings['raytracer_dir']
rayfile_directory = settings['rayfile_directory']

minalt = settings["minalt"]+R_E
Kp = settings['Kp']
AE = settings['AE']
Pdyn = settings['Pdyn']
Dst = settings['Dst']
ByIMF = settings['ByIMF']
BzIMF = settings['BzIMF']
W = settings['W']
use_tsyg = settings['use_tsyg']
use_IGRF = settings['use_IGRF']

t_max = settings['t_max']
dt0 = settings['dt0']
dtmax = settings['dtmax']
root = settings['root']
fixedstep = settings['fixedstep']
maxerr = settings['maxerr']
maxsteps = settings['maxsteps']
mode = settings['mode']

if mode == 'GCPM':
    md = 6
if mode == 'DE':
    md = 7

# -------------------- process datetime info -----------------------------------
ray_datenum = dt.datetime.strptime(settings['date_time_iso_utc'],"%Y-%m-%dT%H:%M:%S.%fZ")
# can't be timezone naive
ray_datenum.replace(tzinfo=dt.timezone.utc)

# --------------- get each ray position, freq, direcion --------------------------
positions = []
freqs = []
directions = []

for r in settings['rays']:
    alt = r['alt']
    lat = r['lat']
    lon = r['lon']
    ray_start = convert2([[alt+R_E,lat,lon]], [ray_datenum],'GEO','sph',['m','deg','deg'], 'SM', 'car', ['m','m','m'])
    positions.append(ray_start[0])
    
    freq = r['freq']
    freqs.append(freq)

    alpha = r['alpha'] 
    phi = r['phi']

    direc = getBdir(ray_start, ray_datenum, rayfile_directory, [alpha], [phi], md, freq)
    directions.append(direc[0])

# run the rays
parallel_run_rays(ray_datenum, positions, directions, freqs, rayfile_directory, md)