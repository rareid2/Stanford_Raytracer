"""

here is a script that will call run_rays and plot the trajectory with
normalized wave power as a color scale

"""

# import needed packages
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import os
import sys
import datetime as dt
# import functions and settings from example scripts directory
from example_scripts.raytracer_utils import readdump, read_rayfile, read_rayfiles, read_damp
from example_scripts.run_rays import run_rays
from example_scripts.raytracer_settings import *
from example_scripts.IGRF_funcs import B_dir, trace_fieldline_ODE
# Spacepy (for coordinate transforms)
from spacepy import coordinates as coord
from spacepy.time import Ticktock
# for color bar plotting
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm, ListedColormap, BoundaryNorm
# for satellite orbits
from get_TLE import get_TLE

# TODO: coordinates?
# --------------------------- Orbits -------------------------------
# DSX TLE:
line1 = '1 44344U 19036F   20113.30349832 -.00000013 +00000-0 +00000-0 0  9992'
line2 = '2 44344 042.2517 093.1041 1975016 129.9693 249.1586 04.54371641013721'
DSX_pos, DSX_t = get_TLE(line1, line2, 'DSX')
x_DSX = DSX_pos[0]
y_DSX = DSX_pos[1]
z_DSX = DSX_pos[2]

# VPM TLE:
line1 = '1 45120U 19071K   20113.86456076  .00004088  00000-0  13479-3 0  9996'
line2 = '2 45120  51.6417 271.7875 0011770 209.5167 150.5152 15.33686764012484'
VPM_pos, VPM_t = get_TLE(line1, line2, 'VPM')
x_VPM = VPM_pos[0]
y_VPM = VPM_pos[1]
z_VPM = VPM_pos[2]

# --------------------------- Ray Tracing --------------------------
# define lists here - must be lists even if only one arg

freq = [26e3]
#n_pos = np.linspace(0,2,1)
n_pos = [10]
positions = [np.array([x_DSX[int(npos)], y_DSX[int(npos)], z_DSX[int(npos)]]) for npos in n_pos]

"""
thetalist is angle from local magnetic field direction in XZ plane
for example: 0 deg = parallel Bfield 
15 = 15 deg counterclockwise to Bfield
"""

thetalist = [0, 45] # in deg

# initialize - leave empty
directions = []
Bxlines = []
Bzlines =[]

for position in positions:
    # grab XZ position
    startpoint = [position[0]/R_E, position[2]/R_E]
    path = 1
    # get bfield direction
    Bx, Bz = B_dir(0, startpoint, 0, '0', path)

    # grab full magnetic field line for plotting later
    Bxline, Bzline = trace_fieldline_ODE(startpoint, 0, '0', path)
    Bxlines.append(np.array(Bxline))
    Bzlines.append(np.array(Bzline))
    Bxline, Bzline = trace_fieldline_ODE(startpoint, 0, '0', -1*path)
    Bxlines.append(np.array(Bxline))
    Bzlines.append(np.array(Bzline))

    # rotate around direction of field line
    for theta in thetalist:
        # rotate using theta
        Rv = np.array([Bx * np.cos(D2R * theta) - Bz * np.sin(D2R * theta),
                       Bx * np.sin(D2R * theta) + Bz * np.cos(D2R * theta)])
        dir_vec = np.array([Rv[0], 0, Rv[1]])
        direction = dir_vec/np.linalg.norm(dir_vec)
        # add that direction
        directions.append(direction)

# number of rays
n_rays = len(freq) * len(positions) * len(directions)
print('about to run: ', n_rays, ' rays')
# run!
run_rays(freq, positions, directions)

# ---------------------- Load output directory -------------------------

# Load all the rayfiles in the output directory
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

# quick check: did the rays propagate?
raylist = [checkray for checkray in raylist if not len(checkray["time"]) < 2]

# abandon if not
if raylist == []:
    sys.exit(0)

# ------------------------ Coordinate Conversion --------------------------

# convert to desired coordinate system into vector list rays
rays = []
for r in raylist:
    tmp_coords = coord.Coords(list(zip(r['pos'].x, r['pos'].y, r['pos'].z)), 'SM', 'car', units=['m', 'm', 'm'])
    tvec_datetime = [ray_datenum + dt.timedelta(seconds=s) for s in r['time']]
    tmp_coords.ticks = Ticktock(tvec_datetime)  # add ticks
    tmp_coords.sim_time = r['time']
    rays.append(tmp_coords)

#----------------------------- Plot rays ----------------------------------
#fig, ax = plt.subplots(1,1, sharex=True, sharey=True)
lw = 2  # linewidth
fig = plt.figure()
ax = plt.axes(projection='3d')

#initialize
r_length = []
rx = []
ry = []
rz = []

for r in rays:
    rx.append(r.x / R_E)
    ry.append(r.y / R_E)
    rz.append(r.z / R_E)
    r_length.append(len(r))

dlist = []
for d in damplist:
    damp = d["damping"]
    damp = np.squeeze(np.array(damp))
    #if len(damp) < max(r_length):
    #    leftover = max(r_length) - len(damp)
    #    damp = np.concatenate((damp, np.zeros(int(leftover))), axis=0)
    dlist.append(damp)

#def myplot(ax, xs, ys, zs, cs, cmap):
#    for x, y, z, c in zip(xs, ys, zs, cs):
#        plot = LineCollection([np.column_stack((x, y, z))], cmap=cmap, zorder=102)
#        plot.set_array(c)
#        ax.add_collection(plot)
#    return plot

# compress for plotting in 3D
darray = np.ndarray.flatten(np.array(dlist))
rx = np.ndarray.flatten(np.array(rx))
ry = np.ndarray.flatten(np.array(ry))
rz = np.ndarray.flatten(np.array(rz))

# plot it!
line = ax.scatter3D(rx, ry, rz, c=darray, cmap='Reds');
fig.colorbar(line, ax=ax, label = 'Normalized wave power')

# --------------------------- Earth and Ionosphere ---------------------------
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
earthx = np.outer(np.cos(u), np.sin(v))
earthy = np.outer(np.sin(u), np.sin(v))
earthz = np.outer(np.ones(np.size(u)), np.cos(v))
ax.plot_surface(earthx, earthy, earthz, color='b', zorder = 100)

# add in the ionosphere
iono_r = (R_E + H_IONO) / R_E
ionox = iono_r * np.outer(np.cos(u), np.sin(v))
ionoy = iono_r * np.outer(np.sin(u), np.sin(v))
ionoz = iono_r * np.outer(np.ones(np.size(u)), np.cos(v))
ax.plot_surface(ionox, ionoy, ionoz, color='g', alpha = 0.2, zorder = 99)

# ------------------------ Field Lines --------------------------
L_shells = [2, 3, 4]  # Field lines to draw
# (from IGRF13 model)
for L in L_shells:
    Lx, Lz = trace_fieldline_ODE([L,0], 0, '0', 1)
    plt.plot(Lx, np.zeros(len(Lx)), Lz, color='b', linewidth=1, linestyle='dashed')
    Lx, Lz = trace_fieldline_ODE([L,0], 0, '0', -1)
    plt.plot(Lx, np.zeros(len(Lx)), Lz, color='b', linewidth=1, linestyle='dashed')
    Lx, Lz = trace_fieldline_ODE([-L,0], 0, '0', 1)
    plt.plot(Lx, np.zeros(len(Lx)), Lz, color='b', linewidth=1, linestyle='dashed')
    Lx, Lz = trace_fieldline_ODE([-L,0], 0, '0', -1)
    plt.plot(Lx, np.zeros(len(Lx)), Lz, color='b', linewidth=1, linestyle='dashed')

# plot field line from orbital position
for linex, linez in zip(Bxlines, Bzlines):
    plt.plot(linex, linez, color='r', linewidth=1, linestyle='dashed')
fig.show()
"""
# ---------------------- Satellite Orbits   -----------------------------
plt.plot(x_DSX/R_E, z_DSX/R_E, c='y', zorder = 105, label = 'DSX')
plt.plot(x_VPM/R_E, z_VPM/R_E, c='y', zorder = 105, label = 'VPM')
for npos in n_pos:
    plt.plot(x_DSX[int(npos)]/R_E, z_DSX[int(npos)]/R_E, '-bo', zorder = 103)

# ------------------------- Plasmasphere ------------------------------
plasma_model_dump = os.path.join(ray_out_dir, 'model_dump_mode_1_XZ.dat')
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

# ----------------------------- More Formatting ----------------------------
ax.set_aspect('equal')
max_lim = max(L_shells)+1

plt.xticks(np.arange(-max_lim, max_lim, step=1))
plt.yticks(np.arange(-max_lim, max_lim, step=1))
plt.xlabel('L (R$_E$)')
plt.ylabel('L (R$_E$)')
plt.xlim([-max_lim, max_lim])
plt.ylim([-2.5, 2.5])

#ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.25), ncol=2)
str_freq = str(int(freq[0] / 1e3))
str_orbital_pos = str(n_pos)
fig_title = str_freq + ' kHz rays in XZ plane'
plt.title(fig_title)

# ------------------------------- Saving ---------------------------------------
#savename = 'plots/XZ_' + str_freq + 'kHz_%03d.png' %p
savename = 'plots/XZ_' + str_freq + 'kHz_SVGTEST.svg'
fig.savefig(savename, format='svg')

plt.show()
plt.close()
"""