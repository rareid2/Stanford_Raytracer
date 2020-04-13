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

# Spacepy (for coordinate transforms)
from spacepy import coordinates as coord
from spacepy.time import Ticktock

# for color bar plotting
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm, ListedColormap, BoundaryNorm

# ------------------ Orbits -------------------------------
from get_TLE import x_sat, y_sat, z_sat

# grab the satellite orbits
x_DSX = x_sat[0]
y_DSX = y_sat[0]
z_DSX = z_sat[0]

x_VPM = x_sat[1]
y_VPM = y_sat[1]
z_VPM = z_sat[1]

# ------------------ Ray Tracing --------------------------
# create lists - must be lists even if only one arg
freq = [26e3]
n_pos = 300
positions = [np.array([x_DSX[n_pos], y_DSX[n_pos], z_DSX[n_pos]])]
#positions = [np.array([-2*R_E,y_DSX[n_pos], 0])]

# need to be unit vectors
theta = np.array([45, 60, 75, 90, 105, 120, 135, 225, 240, 255, 270, 285, 300, 315])
theta = np.array([315, 300])
thetax = np.cos(np.deg2rad(theta))
thetaz = np.sin(np.deg2rad(theta))
directions = []
for angx, angz in zip(thetax, thetaz):
    directions.append(np.array([angx, 0, angz]))

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
        plot = LineCollection([np.column_stack((x, y))], cmap=cmap)
        plot.set_array(z)
        ax.add_collection(plot)
    return plot

line = myplot(ax, rx, rz, dlist, 'Reds')

fig.colorbar(line, ax=ax, label = 'Normalized wave power')

# -------------------------- Figure formatting ---------------------------
L_shells = [2, 3, 4, 5]  # Field lines to draw

# -------- Earth and Iono --------
earth = plt.Circle((0, 0), 1, color='b', alpha=1, zorder=100)
iono = plt.Circle((0, 0), (R_E + H_IONO) / R_E, color='g', alpha=0.5, zorder=99)

ax.add_artist(earth)
ax.add_artist(iono)

# -------- fieldlines -------- (dipole model)
for L in L_shells:
    # Plot dipole field lines for both profile views
    lam = np.linspace(-80, 80, 181)
    L_r = L * pow(np.cos(lam * D2R), 2)
    Lx = L_r * np.cos(lam * D2R)
    Lz = L_r * np.sin(lam * D2R)
    ax.plot(Lx, Lz, color='b', linewidth=1, linestyle='dashed')  # Field line
    ax.plot(-Lx, Lz, color='b', linewidth=1, linestyle='dashed')  # Field line (other side)

# -------- sat orbits   --------

plt.plot(x_DSX/R_E, z_DSX/R_E, c='y', zorder = 101, label = 'DSX')
plt.plot(x_VPM/R_E, z_VPM/R_E, c='y', zorder = 102, label = 'VPM')
plt.plot(x_DSX[300]/R_E, z_DSX[300]/R_E, '-bo', zorder = 103)

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
#

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
