
"""
here is a script that will call run_rays and plot the trajectory with
normalized wave power as a color scale
ONLY WORKS FOR A SINGLE POSITION - CHANGE THE TIME FOR DIFF POSITIONS
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
from run_model_dump import modeldump
# coordinate mania!
# TLES give us GEI, raytracer needs SM, and IGRF funcs needs GEO

# -------------------------------- SET TIME and OTHER SETTINGS --------------------------------

year = 2020
month = 5
day = 14
hours = 6
minutes = 59
seconds = 0

# change time information here - use UTC -
ray_datenum = dt.datetime(2020, 6, 14, 6, 59)

freq = [26.1e3] # Hz
thetalist = [0] # in deg -- what angles to launch at? 

checkdir = 0
crs_out = 'MAG'  # theres a bug with MAG coords -- maybe its the fieldlines? 
datadir = '/home/rileyannereid/workspace/SR-output/WSMR/'

# -------------------------------- GET POSITIONS --------------------------------

# define a WSMR at 1000 km up
GEOsph_wsmrpt = coord.Coords([1000e3 + R_E, 32.943896, -106.41965], 'GEO', 'sph', units=['m', 'deg', 'deg'])
GEOsph_wsmrpt.ticks = Ticktock(ray_datenum, 'UTC')
SMcar_wsmrpt = GEOsph_wsmrpt.convert('SM', 'car')
GEOcar_wsmrpt = GEOsph_wsmrpt.convert('GEO', 'car')
MAGsph_wsmrpt = GEOsph_wsmrpt.convert('MAG', 'sph')
Bstart = [float(GEOcar_wsmrpt.x) / R_E, float(GEOcar_wsmrpt.y) / R_E, float(GEOcar_wsmrpt.z) / R_E]

# -------------------------------- DEFINE RAY DIRECTIONS --------------------------------
# start position of raytracer
K = np.array([[1,0], [0,1]])
m = np.array([0, 0]).reshape(2, 1)
d = 2
n = int(100)
z = np.random.multivariate_normal(mean=m.reshape(d,), cov=K, size=n)
y = np.transpose(z)
dx = y[0] * 1000e3
dy = y[1] * 1000e3

newc = [(float(GEOcar_wsmrpt.x + sx), float(GEOcar_wsmrpt.y + sy), float(GEOcar_wsmrpt.z)) for sx, sy in zip(dx, dy)]
sc = coord.Coords(newc, 'GEO', 'car', units = ['m', 'm', 'm'])
sc.ticks = Ticktock([ray_datenum for i in np.ones(len(dx))], 'UTC')
ndone = sc.convert('GEO', 'sph')
SMpos = sc.convert('SM', 'car')

rayt = ray_datenum
modeldump(rayt.year, rayt.month, rayt.day, rayt.hour, rayt.minute, rayt.second)

dir = -1  # south
dirstr = 'south'

# fill for raytracer call
cleanpos = [(float(pos.x), float(pos.y), float(pos.z)) for pos in SMpos]
positions = cleanpos # fill from pdf
directions = [np.zeros(3) for i in range(n-1)]

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
    #print(len(tvec_datetime))
    tmp_coords.ticks = Ticktock(tvec_datetime, 'UTC')  # add ticks
    tmp_coords.sim_time = r['time']
    new_coords = tmp_coords.convert(crs_out, 'sph') # needs to be sph for rotation
    rays.append(new_coords)

# -------------------------------- PLOTTING --------------------------------
fig, ax = plt.subplots(1,1, sharex=True, sharey=True)
lw = 2  # linewidth

th = -MAGsph_wsmrpt.long

# rotate rays to be at prime merid also and plot
#for r, d in zip(rays, dlist):
for r in rays:
    if len(r.radi) > 2:
        rrad = []
        rlat = []
        rlon = []
        rrad.append(r.radi)
        rlat.append(r.lati)
        rlon.append(r.long)

        rrlon = [rl + th for rl in rlon]
        
        rcoords = [np.column_stack([rr, rl, rrl]) for rr, rl, rrl in zip(rrad, rlat, rrlon)]
       
        tvec_datetime = [ray_datenum for s in range(len(rrlon[0]))]
        
        Rot_ray = coord.Coords(rcoords[0], crs_out, 'sph', units=['m', 'deg', 'deg'])
        Rot_ray.ticks = Ticktock(tvec_datetime, 'UTC')
        outcar_ray = Rot_ray.convert(crs_out, 'car')

    if len(outcar_ray.x) > 1:
        # plotp = ax.scatter(MAGcar_ray.x / R_E, MAGcar_ray.z / R_E, c=d, s = 1, cmap = 'Reds', vmin = 0, vmax = 1.5, zorder = 103)
        plotp = ax.scatter(outcar_ray.x / R_E, outcar_ray.z / R_E, c = 'Red', s = 1, zorder = 103)

# add in color bar - will be just for the last ray, but bounds are set
# plt.colorbar(plotp, label = 'Normalized wave power')

# -------------------------------- EARTH AND IONO --------------------------------
earth = plt.Circle((0, 0), 1, color='b', alpha=0.5, zorder=100)
iono = plt.Circle((0, 0), (R_E + H_IONO) / R_E, color='g', alpha=0.5, zorder=99)
ax.add_artist(earth)
ax.add_artist(iono)

# -------------------------------- PLASMASPHERE --------------------------------
# bad pls change this @ riley!
# NOT IN CORRECT COORDS
path2plasma = '/home/rileyannereid/workspace/Stanford_Raytracer/example_scripts/modeldumps/'
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
L_shells = [2, 3, 4, -2, -3, -4]  # Field lines to draw

for L in L_shells:
    Blines = []
    if L < 0:
        rot = th
    else:
        rot = -th

    Lstart = [L, 0, rot]
    Lcoords = coord.Coords(Lstart, crs_out, 'sph', units=['Re', 'deg', 'deg'])
    Lcoords.ticks = Ticktock(ray_datenum, 'UTC')
    GEO_Lcoords = Lcoords.convert('GEO', 'car')
    Blines.append(trace_fieldline_ODE([float(GEO_Lcoords.x),float(GEO_Lcoords.y),float(GEO_Lcoords.z)], 0, '0', 1, ray_datenum))
    Blines.append(trace_fieldline_ODE([float(GEO_Lcoords.x),float(GEO_Lcoords.y),float(GEO_Lcoords.z)], 0, '0', -1, ray_datenum))

    for blinex, bliney, blinez in Blines:
        raytime = []
        # create list of the same times
        for m in range(int(len(blinex))):
            raytime.append(ray_datenum)
        
        # convert coords
        bpos = np.column_stack((blinex, bliney, blinez))
        bpos = coord.Coords(bpos, 'GEO', 'car', units=['Re', 'Re', 'Re'])
        bpos.ticks = Ticktock(raytime, 'UTC') # add ticks
        outsph_bline = bpos.convert(crs_out, 'sph')

        btemp_rad = []
        btemp_lat = []
        btemp_lon = []
        btemp_rad.append(outsph_bline.radi)
        btemp_lat.append(outsph_bline.lati)
        btemp_lon.append(outsph_bline.long)

        brlon = [bl * 0 for bl in btemp_lon]
        bcoords = [np.column_stack([br, bl, brl]) for br, bl, brl in zip(btemp_rad, btemp_lat, brlon)]
        Rot_bline = coord.Coords(bcoords[0], crs_out, 'sph', units=['Re', 'deg', 'deg'])
        Rot_bline.ticks = Ticktock(raytime, 'UTC') # add ticks
        outcar_bline = Rot_bline.convert(crs_out, 'car')
        plt.plot(np.sign(rot) * outcar_bline.x, np.sign(rot) * outcar_bline.z, color='b', linewidth=1, linestyle='dashed')

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
    
    # convert
    bpos = np.column_stack((blinex, bliney, blinez))
    bpos = coord.Coords(bpos, 'GEO', 'car', units=['Re', 'Re', 'Re'])
    bpos.ticks = Ticktock(raytime, 'UTC') # add ticks
    outsph_bline = bpos.convert(crs_out, 'sph')

    btemp_rad = []
    btemp_lat = []
    btemp_lon = []
    btemp_rad.append(outsph_bline.radi)
    btemp_lat.append(outsph_bline.lati)
    btemp_lon.append(outsph_bline.long)

    brlon = [bl * 0 for bl in btemp_lon]
    bcoords = [np.column_stack([br, bl, brl]) for br, bl, brl in zip(btemp_rad, btemp_lat, brlon)]
    Rot_bline = coord.Coords(bcoords[0], crs_out, 'sph', units=['Re', 'deg', 'deg'])
    Rot_bline.ticks = Ticktock(raytime, 'UTC') # add ticks
    outcar_bline = Rot_bline.convert(crs_out, 'car')
    plt.plot(outcar_bline.x, outcar_bline.z, color='r', linewidth=1, linestyle='dashed')

# -------------------------------- GET FOOTPRINT --------------------------------
# also in GEO car, so need to use bstart 
GDZsph_foot = findFootprints(ray_datenum, Bstart, dirstr)
GDZsph_foot.units = ['km', 'deg', 'deg']
GDZsph_foot.ticks = Ticktock(ray_datenum, 'UTC')
outsph_foot = GDZsph_foot.convert(crs_out, 'sph')

# rotate and plot
Rot_foot = coord.Coords([float(outsph_foot.radi), float(outsph_foot.lati), float(outsph_foot.long + th)], crs_out, 'sph', units=['Re', 'deg', 'deg'])
Rot_foot.ticks = Ticktock(ray_datenum, 'UTC')
outcar_foot = Rot_foot.convert(crs_out, 'car')
plt.plot(outcar_foot.x, outcar_foot.z, '-ro', label='Bfield footpoint')

# -------------------------------- FORMATTING --------------------------------
ax.set_aspect('equal')
max_lim = 4

plt.xticks(np.arange(-max_lim, max_lim, step=1))
plt.yticks(np.arange(-max_lim, max_lim, step=1))
plt.xlabel('L (R$_E$)')
plt.ylabel('L (R$_E$)')
plt.xlim([-max_lim, max_lim])
plt.ylim([-2.5, 2.5])

mytitle = str(freq[0]/1e3) + 'kHz rays at ' + str(ray_datenum) + ' in ' + crs_out + ' coords'
plt.title(mytitle)
ax.legend(loc = 'lower center', fontsize =  'x-small')

savename = datadir + str(freq[0]/1e3) + 'kHz' + str(ray_datenum.year) + str(ray_datenum.month) + str(ray_datenum.day) + str(ray_datenum.hour) + str(ray_datenum.minute) + '2Dview.png'
#plt.savefig(savename, format='svg')
plt.savefig(savename, format='png')
#plt.close()
plt.show()

# ------------------------------------------- END --------------------------------------------