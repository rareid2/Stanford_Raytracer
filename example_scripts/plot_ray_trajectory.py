"""

here is a script that will call run_rays and plot the trajectory with
normalized wave power as a color scale
this is currently set for XZ coordinates

"""

# import needed packages
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import datetime as dt

# import functions and settings from this example scripts directory
from raytracer_utils import readdump, read_rayfile, read_rayfiles, read_damp
from run_rays import run_rays
from raytracer_settings import *

# for IGRF
from IGRF_funcs import IGRFline, IGRFdirection

# Spacepy (for coordinate transforms)
from spacepy import coordinates as coord
from spacepy.time import Ticktock

# for color bar plotting
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm, ListedColormap, BoundaryNorm

# for satellites
from get_TLE import get_TLE


# ------------------ Orbits -------------------------------
# DSX TLE:
line1 = '1 44344U 19036F   20099.44261897 -.00000008 +00000-0 +00000-0 0  9998'
line2 = '2 44344 042.2458 098.1824 1975230 124.0282 256.3811 04.54371606013099'

x_DSX, y_DSX, z_DSX = get_TLE(line1, line2, 'DSX')
x_VPM, y_VPM, z_VPM = get_TLE(line1, line2, 'VPM')

# ------------------ Ray Tracing --------------------------
# define lists here - must be lists even if only one arg
freq = [26e3]
n_pos = 40
positions = [np.array([x_DSX[n_pos], y_DSX[n_pos], z_DSX[n_pos]])]
thetalist = [0] # in deg
directions = []

# theta is counterclockwise direction, use 0 for field-aligned

# find the direction of the local magnetic field line using IGRFdirection

for position in positions:

    # convert to lla real quick
    cart_coords = [position[0], position[1], position[2]]
    cvals = coord.Coords(cart_coords, 'GEO', 'car')

    # add time vector - just one input
    cvals.ticks = Ticktock(ray_datenum)  # add ticks
    # convert to lla
    llacoord = cvals.convert('GEO', 'sph')

    lat2, lon2, alt2 = IGRFline(llacoord.lati, llacoord.long, llacoord.radi)
    sph_coords = list(zip(alt2, lat2, lon2))
    cvals = coord.Coords(sph_coords, 'GEO', 'sph')
    tvec_datetime = [ray_datenum + dt.timedelta(seconds=s) for s in range(len(sph_coords))]
    cvals.ticks = Ticktock(tvec_datetime)  # add ticks
    field_line = cvals.convert('GEO', 'car')

    # find direction of local magnetic field line
    dr, dp, dth = IGRFdirection(llacoord.lati, llacoord.long, llacoord.radi)
    dx = dr * np.sin(dth) * np.cos(dp)
    dy = dr * np.sin(dth) * np.sin(dp)
    dz = dr * np.cos(dth)

    # desired angles:
    for theta in thetalist:

        dline = np.array([dx, dy, dz])
        unit_d = dline/np.linalg.norm(dline)

        # rotation matrix around y-axis
        Ry = np.array([[np.cos(theta*D2R), 0, np.sin(theta*D2R)],
                      [0, 1, 0,], [-np.sin(theta*D2R), 0, np.cos(theta*D2R)]])
        direction = np.matmul(Ry, unit_d)
        if direction.all() == unit_d.all():
            direction = np.zeros(3)
        directions.append(direction)

# need to be unit vectors
#theta = np.array([45, 315])
#thetax = np.cos(np.deg2rad(theta))
#thetaz = np.sin(np.deg2rad(theta))
#directions = []
#for angx, angz in zip(thetax, thetaz):
#    directions.append(np.array([angx, 0, angz]))

# number of rays
n_rays = len(freq) * len(positions) * len(directions)

# run!
run_rays(freq, positions, directions)

# -------------- Load output directory ----------

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

# ------------------ Coordinate Conversion --------------------------

# convert to desired coordinate system into vector list rays
rays = []
for r in raylist:
    tmp_coords = coord.Coords(list(zip(r['pos'].x, r['pos'].y, r['pos'].z)), 'SM', 'car', units=['m', 'm', 'm'])
    tvec_datetime = [ray_datenum + dt.timedelta(seconds=s) for s in r['time']]
    tmp_coords.ticks = Ticktock(tvec_datetime)  # add ticks
    tmp_coords.sim_time = r['time']
    rays.append(tmp_coords)

#-------------------------- Plot rays ----------------------------------
fig, ax = plt.subplots(1,1, sharex=True, sharey=True)
lw = 2  # linewidth

#initialize
r_length = []
rx = []
rz = []

for r in rays:
    rx.append(r.x / R_E)
    rz.append(r.z / R_E)
    r_length.append(len(r))

dlist = []
for d in damplist:
    damp = d["damping"]
    damp = np.squeeze(np.array(damp))
    if len(damp) < max(r_length):
        leftover = max(r_length) - len(damp)
        damp = np.concatenate((damp, np.zeros(int(leftover))), axis=0)
    dlist.append(damp)

def myplot(ax, xs, ys, zs, cmap):
    for x, y, z in zip(xs, ys, zs):
        plot = LineCollection([np.column_stack((x, y))], cmap=cmap, zorder=102)
        plot.set_array(z)
        ax.add_collection(plot)
    return plot

#line = myplot(ax, rx, rz, dlist, 'Reds')

#fig.colorbar(line, ax=ax, label = 'Normalized wave power')

# -------------------------- Figure formatting ---------------------------
L_shells = [2, 3, 4, 5]  # Field lines to draw

# -------- Earth and Iono --------
earth = plt.Circle((0, 0), 1, color='b', alpha=1, zorder=100)
iono = plt.Circle((0, 0), (R_E + H_IONO) / R_E, color='g', alpha=0.5, zorder=99)

ax.add_artist(earth)
ax.add_artist(iono)

# -------- fieldlines -------- (from IGRF13 model)

for L in L_shells:
#   Plot IGRF field lines
    lat, lon, alt = IGRFline(0, 0, L * R_E)
    sph_coords = list(zip(alt, lat, lon))
    cvals = coord.Coords(sph_coords, 'GEO', 'sph')
    tvec_datetime = [ray_datenum + dt.timedelta(seconds=s) for s in range(len(sph_coords))]
    cvals.ticks = Ticktock(tvec_datetime)  # add ticks
    newcoord = cvals.convert('GEO', 'car')

    plt.plot(newcoord.x/R_E, newcoord.z/R_E, color='b', linewidth=1, linestyle='dashed')
    plt.plot(-newcoord.x/R_E, newcoord.z/R_E, color='b', linewidth=1, linestyle='dashed')
    plt.plot(newcoord.x/R_E, -newcoord.z/R_E, color='b', linewidth=1, linestyle='dashed')
    plt.plot(-newcoord.x/R_E, -newcoord.z/R_E, color='b', linewidth=1, linestyle='dashed')

# plot field line from orbital position
plt.plot(field_line.x/R_E, field_line.z/R_E, color='r', linewidth=1, linestyle='dashed')

for L in L_shells:
    # Plot dipole field lines for both profile views
    lam = np.linspace(-80, 80, 181)
    L_r = L * pow(np.cos(lam * D2R), 2)
    Lx = L_r * np.cos(lam * D2R)
    Lz = L_r * np.sin(lam * D2R)
    ax.plot(Lx, Lz, color='g', linewidth=1, linestyle='dashed')  # Field line
    ax.plot(-Lx, Lz, color='g', linewidth=1, linestyle='dashed')  # Field line (other side)

# -------- sat orbits   --------

#plt.plot(x_DSX/R_E, z_DSX/R_E, c='y', zorder = 101, label = 'DSX')
#plt.plot(x_VPM/R_E, z_VPM/R_E, c='y', zorder = 102, label = 'VPM')
plt.plot(x_DSX[n_pos]/R_E, z_DSX[n_pos]/R_E, '-bo', zorder = 103)

# -------- plasmapause --------
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


# -------- figure formatting --------
ax.set_aspect('equal')

max_lim = max(L_shells)+1

plt.xticks(np.arange(-max_lim, max_lim, step=1))
plt.yticks(np.arange(-max_lim, max_lim, step=1))
plt.xlabel('L (R$_E$)')
plt.ylabel('L (R$_E$)')
plt.xlim([-max_lim, max_lim])
plt.ylim([-2.5, 2.5])

ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.25), ncol=2)
str_freq = str(int(freq[0] / 1e3))
str_orbital_pos = str(n_pos)
fig_title = str_freq + ' kHz rays \n in SM coordinates in XZ plane'
plt.title(fig_title)

# -------- saving --------
#savename = 'plots/XZ_' + str_freq + 'kHz_%03d.png' %p
#savename = 'plots/XZ_' + str_freq + 'kHz_test2.png'
#fig.savefig(savename)

plt.show()
#plt.close()
