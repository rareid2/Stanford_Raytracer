
"""
here is a script that will call run_rays and plot the trajectory with
normalized wave power as a color scale
ONLY WORKS FOR A SINGLE POSITION - CHANGE THE TIME FOR DIFF POSITIONS
plotting in MAG coordinates
"""

# import needed packages
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
import sys
import datetime as dt
from raytracer_utils import readdump, read_rayfile, read_rayfiles, read_damp
from run_rays import run_rays
from raytracer_settings import *
from IGRF_funcs import B_dir, trace_fieldline_ODE, findFootprints, B_direasy
from spacepy import coordinates as coord
from spacepy import irbempy
from spacepy.time import Ticktock
from TLE_funcs import TLE2pos
import tempfile

# coordinate mania!
# TLES give us GEI, raytracer needs SM, and IGRF funcs needs GEO

# -------------------------------- SET TIME --------------------------------
# change time information here - use UTC -
year = 2020
month = 5
day = 20
hours = 0
minutes = 0
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
r, tvec = TLE2pos(lines1, lines2, satnames, 1, ray_datenum)

# redefine time here -- more accurate
# ray_datenum = tvec[0]

# convert to meters
dsx = [rpos*1e3 for rpos in r[0]]
vpm = [rpos*1e3 for rpos in r[1]]

# only grab first one - weird bug fix with JD dates
dsx = [dsx[0]]
vpm = [vpm[0]]

# convert startpoint to SM car for raytracer
GEIcar_dsx = coord.Coords(dsx, 'GEI', 'car', units=['m', 'm', 'm'])
GEIcar_dsx.ticks = Ticktock(ray_datenum, 'UTC') # add ticks
SMcar_dsx = GEIcar_dsx.convert('SM', 'car')

# convert vpm to MAG sph for plotting
GEIcar_vpm = coord.Coords(vpm, 'GEI', 'car', units=['m', 'm', 'm'])
GEIcar_vpm.ticks = Ticktock(ray_datenum, 'UTC') # add ticks
MAGsph_vpm = GEIcar_vpm.convert('MAG', 'sph')

# -------------------------------- DEFINE RAY DIRECTIONS --------------------------------
# start position of raytracer
position = [float(SMcar_dsx.x), float(SMcar_dsx.y), float(SMcar_dsx.z)]

freq = [26e3] # Hz
thetalist = [45, 50, 55, 60, 65, 70, 75, 80, 85, 90, -45, -50, -55, -60, -65, -70, -75, -80, -85, -90] # in deg -- what angles to launch at? 
#thetalist = [0, 5, 10, 15, 25, 30, 35, 40, 45, -45, -40, -35, -30, -25, -20, -15, -10, -5]
#[45, 50, 55, 60, 65, 70, 75, 80, 85, 90, -45, -50, -55, -60, -65, -70, -75, -80, -85]

# grab position and find direction of local bfield
GEOcar_dsx = GEIcar_dsx.convert('GEO', 'car')

# check with hemi we are in
GEOsph_dsx = GEIcar_dsx.convert('GEO', 'sph')
if GEOsph_dsx.lati > 0:
    dir = 1   # north
else:
    dir = -1  # south

Bstart = [float(GEOcar_dsx.x)/R_E, float(GEOcar_dsx.y)/R_E, float(GEOcar_dsx.z)/R_E]
Bx, By, Bz = B_direasy(ray_datenum, Bstart, dir)

# convert to SM coordinates for raytracer
dirB = np.reshape(np.array([Bx, By, Bz]), (1, 3))
dirB = coord.Coords(dirB[0], 'GEO', 'car', units=['Re', 'Re', 'Re'])
dirB.ticks = Ticktock(ray_datenum, 'UTC') # add ticks
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
    Rot_dirB.ticks = Ticktock(ray_datenum, 'UTC') # add ticks
    SMcar_dirB = Rot_dirB.convert('SM', 'car')

    direction = [float(SMcar_dirB.x), float(SMcar_dirB.y), float(SMcar_dirB.z)]
    
    # add the normalized direction (or zeros)
    directions.append(np.squeeze(direction))

    # make sure position list matches direction list
    positions.append(position)

# -------------------------------- RUN RAYS --------------------------------
# convert for raytracer settings
days_in_the_year = ray_datenum.timetuple().tm_yday
days_in_the_year = format(days_in_the_year, '03d')

# yearday and miliseconds day are used by raytracer
yearday = str(year)+ str(days_in_the_year)   # YYYYDDD
milliseconds_day = hours*3.6e6 + minutes*6e4 + seconds*1e3

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
    tvec_datetime = [ray_datenum + dt.timedelta(seconds=s) for s in r['time']]
    tmp_coords.ticks = Ticktock(tvec_datetime, 'UTC')  # add ticks
    tmp_coords.sim_time = r['time']
    new_coords = tmp_coords.convert('MAG', 'sph')
    rays.append(new_coords)

dlist = []
for d in damplist:
    damp = d["damping"]
    damp = np.squeeze(np.array(damp))
    dlist.append(damp)
    print(damp)
    
# -------------------------------- PLOTTING --------------------------------
fig, ax = plt.subplots(1,1, sharex=True, sharey=True)
lw = 2  # linewidth

# rotate plot to be in plane of view
MAGsph_dsx = GEIcar_dsx.convert('MAG', 'sph')
th = -MAGsph_dsx.long

# rotate long to be at prime merid
Rot_dsx = coord.Coords([float(MAGsph_dsx.radi), float(MAGsph_dsx.lati), float(MAGsph_dsx.long + th)], 'MAG', 'sph', units=['m', 'deg', 'deg'])
Rot_vpm = coord.Coords([float(MAGsph_vpm.radi), float(MAGsph_vpm.lati), float(MAGsph_vpm.long + th)], 'MAG', 'sph', units=['m', 'deg', 'deg'])
MAGcar_dsx = Rot_dsx.convert('MAG', 'car')
MAGcar_vpm = Rot_vpm.convert('MAG', 'car')

# plot sat locations
plt.plot(MAGcar_dsx.x / R_E, MAGcar_dsx.z / R_E, '-go', zorder=105, label='DSX')
plt.plot(MAGcar_vpm.x / R_E, MAGcar_vpm.z / R_E, '-yo', zorder=104, label='VPM')

# rotate rays to be at prime merid also and plot
for r, d in zip(rays, dlist):
    rrad = []
    rlat = []
    rlon = []
    rrad.append(r.radi)
    rlat.append(r.lati)
    rlon.append(r.long)

    rrlon = [rl + th for rl in rlon]
    rcoords = [np.column_stack([rr, rl, rrl]) for rr, rl, rrl in zip(rrad, rlat, rrlon)]

    Rot_ray = coord.Coords(rcoords[0], 'MAG', 'sph', units=['m', 'deg', 'deg'])
    MAGcar_ray = Rot_ray.convert('MAG', 'car')

    if len(MAGcar_ray.x) > 1:
        plotp = ax.scatter(MAGcar_ray.x / R_E, MAGcar_ray.z / R_E, c=d, s = 1, cmap = 'Reds', vmin = 0, vmax = 1.5, zorder = 103)
        # plotp = ax.scatter(MAGcar_ray.x / R_E, MAGcar_ray.z / R_E, c = 'Red', s = 1, zorder = 103)
# add in color bar - will be just for the last ray, but bounds are set
plt.colorbar(plotp, label = 'Normalized wave power')

# -------------------------------- EARTH AND IONO --------------------------------
earth = plt.Circle((0, 0), 1, color='b', alpha=0.5, zorder=100)
iono = plt.Circle((0, 0), (R_E + H_IONO) / R_E, color='g', alpha=0.5, zorder=99)
ax.add_artist(earth)
ax.add_artist(iono)

# -------------------------------- PLASMASPHERE --------------------------------
# bad pls change this
path2plasma = '/Users/rileyannereid/workspace/Stanford_Raytracer/example_scripts/modeldumps/'
plasma_model_dump = os.path.join(path2plasma, 'model_dump_mode_1_XZ.dat')
d_xz = readdump(plasma_model_dump)
Ne_xz = d_xz['Ns'][0, :, :, :].squeeze().T * 1e-6
Ne_xz[np.isnan(Ne_xz)] = 0

# Axis spacing depends on how the modeldump was ran
psize = 10
px = np.linspace(-10, 10, 200)
py = np.linspace(-10, 10, 200)

# Colorbar limits (log space)
clims = [-2, 5]

# Plot background plasma (equatorial slice)
g = plt.pcolormesh(px, py, np.log(Ne_xz), cmap = 'twilight')
#fig.colorbar(g, ax=ax, orientation="horizontal", pad = 0.2, label= 'Plasmasphere density')

# ---------------------------------- BFIELD -----------------------------------

# (from IGRF13 model)
L_shells = [2, 3, 4]  # Field lines to draw
#why dont we convert here? 
for L in L_shells:
    Lx, Ly, Lz = trace_fieldline_ODE([L,0,0], 0, '0', 1, ray_datenum)
    plt.plot(Lx, Lz, color='b', linewidth=1, linestyle='dashed')
    Lx, Ly, Lz = trace_fieldline_ODE([L,0,0], 0, '0', -1, ray_datenum)
    plt.plot(Lx, Lz, color='b', linewidth=1, linestyle='dashed')
    Lx, Ly, Lz = trace_fieldline_ODE([-L,0,0], 0, '0', 1, ray_datenum)
    plt.plot(Lx, Lz, color='b', linewidth=1, linestyle='dashed')
    Lx, Ly, Lz = trace_fieldline_ODE([-L,0,0], 0, '0', -1, ray_datenum)
    plt.plot(Lx, Lz, color='b', linewidth=1, linestyle='dashed')
print('finished plotting field lines')

# plot field line from orbital position
Blines = []

Blines.append(trace_fieldline_ODE(Bstart, 0, '0', 1, ray_datenum))
Blines.append(trace_fieldline_ODE(Bstart, 0, '0', -1, ray_datenum))

for blinex, bliney, blinez in Blines:
    raytime = []

    # create list of the same times
    for m in range(int(len(blinex))):
        raytime.append(ray_datenum)
    
    # convert to MAG sph coords
    bpos = np.column_stack((blinex, bliney, blinez))
    bpos = coord.Coords(bpos, 'GEO', 'car', units=['Re', 'Re', 'Re'])
    bpos.ticks = Ticktock(raytime, 'UTC') # add ticks
    MAGsph_bline = bpos.convert('MAG', 'sph')
    
    btemp_rad = []
    btemp_lat = []
    btemp_lon = []
    btemp_rad.append(MAGsph_bline.radi)
    btemp_lat.append(MAGsph_bline.lati)
    btemp_lon.append(MAGsph_bline.long)

    brlon = [bl + th for bl in btemp_lon]
    bcoords = [np.column_stack([br, bl, brl]) for br, bl, brl in zip(btemp_rad, btemp_lat, brlon)]

    Rot_bline = coord.Coords(bcoords[0], 'MAG', 'sph', units=['Re', 'deg', 'deg'])
    MAGcar_bline = Rot_bline.convert('MAG', 'car')
    plt.plot(MAGcar_bline.x, MAGcar_bline.z, color='r', linewidth=1, linestyle='dashed')

# -------------------------------- GET FOOTPRINT --------------------------------
# also in GEO car, so need to use bstart 
GDZsph_foot = findFootprints(ray_datenum, Bstart, 'same')
GDZsph_foot.units = ['km', 'deg', 'deg']
GDZsph_foot.ticks = Ticktock(ray_datenum, 'UTC')
MAGsph_foot = GDZsph_foot.convert('MAG', 'sph')

# rotate and plot
Rot_foot = coord.Coords([float(MAGsph_foot.radi), float(MAGsph_foot.lati), float(MAGsph_foot.long + th)], 'MAG', 'sph', units=['Re', 'deg', 'deg'])
MAGcar_foot = Rot_foot.convert('MAG', 'car')
plt.plot(MAGcar_foot.x, MAGcar_foot.z, '-ro', label='Bfield footpoint')

# -------------------------------- FORMATTING --------------------------------
ax.set_aspect('equal')
max_lim = 4

plt.xticks(np.arange(-max_lim, max_lim, step=1))
plt.yticks(np.arange(-max_lim, max_lim, step=1))
plt.xlabel('L (R$_E$)')
plt.ylabel('L (R$_E$)')
plt.xlim([-max_lim, max_lim])
plt.ylim([-2.5, 2.5])

mytitle = str(freq[0]/1e3) + 'kHz rays at ' + str(ray_datenum)
plt.title(mytitle)
ax.legend(loc = 'lower center', fontsize =  'x-small')

savename = '/Users/rileyannereid/Desktop/' + str(freq[0]/1e3) + 'kray' + str(ray_datenum.year) + str(ray_datenum.month) + str(ray_datenum.day) + str(ray_datenum.hour) + str(ray_datenum.minute) + '.svg'
plt.savefig(savename, format='svg')
#plt.close()
plt.show()

# -------------------------------- END --------------------------------