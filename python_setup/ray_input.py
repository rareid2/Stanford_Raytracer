import numpy as np
import datetime as dt
import os
import json

# import constants and functions from this folder
from constants import *
from convert_coords import convert2
from helper_funcs import get_yearmiliday,getBdir,parallel_run_rays,read_rayfile

#  -------------------------- import settings from config file -----------------------------
rayfile_directory = settings['rayfile_directory']
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