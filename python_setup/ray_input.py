import numpy as np
import datetime as dt
import os
import json

# import constants and functions from this folder
from constants import *
from convert_coords import convert2
from helper_funcs import get_yearmiliday,getBdir_src,parallel_run_rays

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

# only need to run this once! not in parallel!
# Create directory for inputs/outputs if doesn't already exist
ray_out_dir = os.path.join(rayfile_directory,dt.datetime.strftime(ray_datenum, '%Y-%m-%d_%H_%M_%S'))
if not os.path.exists(ray_out_dir):
    os.makedirs(ray_out_dir)

# --------------- get each ray position, freq, direcion --------------------------
positions = []
freqs = []
directions = []
thetas = []
phis = []

for ri,r in enumerate(settings['rays']):
    alt = r['alt']
    lat = r['lat']
    lon = r['lon']
    ray_start = convert2([[alt+R_E,lat,lon]], [ray_datenum],'GEO','sph',['m','deg','deg'], 'SM', 'car', ['m','m','m'])
    positions.append(ray_start[0])
    
    freq = r['freq']
    freqs.append(freq)

    alpha = r['alpha'] 
    phi = r['phi']

    thetas.append(alpha)
    phis.append(phi)

# get directions as a batch
directions = getBdir_src(ray_datenum, positions, thetas, phis, rayfile_directory)

print('initialized all rays')

# run the rays
parallel_run_rays(ray_datenum, positions, directions, freqs, ray_out_dir, md)