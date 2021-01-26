import numpy as np
import datetime as dt
import os

from constants_settings import *
from coordinates import create_spc, convert_spc
from satellites import sat
from bfield import Bfieldinfo, trace_fieldline_ODE
from run_rays import single_run_rays, parallel_run_rays
from raytracer_utils import read_rayfile
from ray_plots import plotray2D, plotrefractivesurface

# example call to ray tracer!
# FIRST, navigate to constants_settings and make sure the settings are correct for the run

# let's look at a conjunction between DSX and VPM:
# use the datetime package to define the start time -- make sure to use UTC timezone
ray_datenum = dt.datetime(2020, 9, 14, 22, 55, tzinfo=dt.timezone.utc)

# first, we need the positions of the satellites -- use the sat class
dsx = sat()             # define a satellite object
dsx.catnmbr = 44344     # prove NORAD ID
dsx.time = ray_datenum  # set time
dsx.getTLE_ephem()      # get TLEs nearest to this time -- sometimes this will lag

# propagate the orbit! setting sec=0 will give you just the position at that time
dsx.propagatefromTLE(sec=0, orbit_dir='future', crs='SM', carsph='car', units=['m','m','m'])

# now we have ray start point in the correct coordinates (SM cartesian in m)
ray_start = dsx.pos

# next, define the direction of the ray!

# set up a bfield class -- use this to find footpoints, trace fieldlines, and get unit vec of the fieldline
bfield_ray_start = Bfieldinfo()
bfield_ray_start.time = ray_datenum
bfield_ray_start.pos = ray_start

# get the unit vector direction going north in SM coords -- use spherical for ease of rotation
bfield_ray_start.Bfield_direction(hemis='north', crs='SM', carsph='sph')
bfield_dir = bfield_ray_start.unit_vec

# we're going to change the polar angle - let's rotate to be 45 deg
# NOTE THAT THIS IS 1-2 degree OFF of what the ray tracer will read -- for exactly field aligned, set ray_start_dir = np.zeros(3)
alpha = 30
ray_start_dir_c = create_spc(cor_array=[float(bfield_dir.radi[0]), float(bfield_dir.lati[0]) + alpha, float(bfield_dir.long[0])], dt_array=ray_datenum, crs='SM', carsph='sph', units=['Re','deg','deg'])
# convert back to SM car -- keep in Re units
ray_start_dir = convert_spc(cvals=ray_start_dir_c, dt_array=ray_datenum, crs='SM', carsph='car', units=['Re','Re','Re'])

# two ways to run

# run at a single time -- use run_rays and input a list of positions, directions, and freqs (ALL SAME LENGTH)
# generates one input file and one output file with all rays in it

nrays = 1 # how many rays
freq = 25e3  # Hz
rayfile_directory = '/home/rileyannereid/workspace/SR-output/rayfiles'

# simplest call
positions = []
directions = []
freqs = []
for n in range(nrays):
    positions.append([float(ray_start.x[0]), float(ray_start.y[0]), float(ray_start.z[0])])
    directions.append([float(ray_start_dir[0].x), float(ray_start_dir[0].y), float(ray_start_dir[0].z)])
    freqs.append(freq)

# time to run is about 1 sec every 10 rays 
single_run_rays(ray_datenum, positions, directions, freqs, rayfile_directory)

# OR run parallel at different times -- use parallel_run_rays and input a list of times, and LIST OF LISTS with positions, 
# directions, and frequencies

# let's say we want to re-run this every 10 seconds in time for a minute
#tvec = [ray_datenum + dt.timedelta(seconds=i) for i in range(0,60,10)]
#positions_list = [positions for i in range(len(tvec))]
#directions_list = [directions for i in range(len(tvec))]
#freqs_list = [freqs for i in range(len(tvec))]
#directory_list = [rayfile_directory for i in range(len(tvec))]

#parallel_run_rays(tvec, positions_list, directions_list, freqs_list, directory_list)


# -------------------------------------------------------------------------
# that's it! let's look at output

# Load all the rayfiles in the output directory
ray_out_dir = rayfile_directory + '/'+ dt.datetime.strftime(ray_datenum, '%Y-%m-%d %H:%M:%S')
file_titles = os.listdir(ray_out_dir)

# create empty lists to fill with ray files and damp files
raylist = []

for filename in file_titles:
    if '.ray' in filename:
        raylist += read_rayfile(os.path.join(ray_out_dir, filename))

#plotray2D(ray_datenum, raylist, ray_out_dir, 'GEO', 'car', units=['Re','Re','Re'])
#plotrefractivesurface(ray_datenum, raylist[0])