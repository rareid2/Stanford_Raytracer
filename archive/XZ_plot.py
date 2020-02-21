# import packages
import numpy as np
import matplotlib.pyplot as plt
import os
import datetime as dt

# import functions from this example scripts directory
from run_them_rays import run_them_rays
from orbits import orbits
from raytracer_utils import readdump, read_rayfile, read_rayfiles, read_damp

# Spacepy (for coordinate transforms)
from spacepy import coordinates as coord
from spacepy.time import Ticktock

# -------------- constants -------------
D2R = (np.pi / 180.0)
R2D = 1.0 / D2R
R_E = 6371e3  # m
H_IONO = 1000e3

# ------------- input parameters -------
# get orbits with set time points
n = 200
dsxorbit, vpmorbit = orbits(n)

freq = 1e3  # Hz
orbital_pos = 0  # from 0 to n
directions = [np.array([1, 0, 0]), np.array([0, 0, 1])]
"""
directions = list()
theta = 40  # deg
for i in range(1, 11):
    dir = np.array([np.cos(np.deg2rad(theta+(5*i))), 0, np.sin(np.deg2rad(theta+(5*i)))])
    directions.append(dir)

theta = 220  # deg
for i in range(1, 11):
    dir = np.array([np.cos(np.deg2rad(theta+(5*i))), 0, np.sin(np.deg2rad(theta+(5*i)))])
    directions.append(dir)
"""
"""
d1 = np.array([1, 0, 1])/np.linalg.norm(np.array([1, 0, 1]))
d2 = np.array([0, 0, 1])/np.linalg.norm(np.array([0, 0, 1]))
d3 = np.array([1, 0, 0])/np.linalg.norm(np.array([1, 0, 0]))
d4 = np.array([0, 0, -1])/np.linalg.norm(np.array([0, 0, -1]))
d5 = np.array([-1, 0, -1])/np.linalg.norm(np.array([-1, 0, -1]))
d6 = np.array([-1, 0, 0])/np.linalg.norm(np.array([-1, 0, 0]))
d7 = np.array([-1, 0, 1])/np.linalg.norm(np.array([-1, 0, 1]))
directions = [np.zeros(3), d1, d2, d3, d4, d5, d6, d7]
"""

fig, ax = plt.subplots()

for direction in directions:
    run_them_rays(orbital_pos, freq, direction)

    fig.gca()

    dir_str  = int(np.round(np.rad2deg(np.arcsin(direction[0]))))

    # Load n' plot the Ngo files:
    project_root = os.getcwd()
    ray_out_dir = os.path.join(project_root, "test_outputs")

    # Load all the rayfiles in the output directory
    d = os.listdir(ray_out_dir)
    file_titles = ['example_ray_mode1']

    raylist = []

    for r in file_titles:
        raylist += read_rayfile(os.path.join(ray_out_dir, r + '.ray'))

    flashtime = dt.datetime(2010, 1, 1, 0, 0, 0)

    # Put the rays into a friendly system (and so we can use coordinate transforms)
    rays = []
    for r in raylist:
        tmp_coords = coord.Coords(list(zip(r['pos'].x, r['pos'].y, r['pos'].z)), 'SM', 'car', units=['m', 'm', 'm'])
        tvec_datetime = [flashtime + dt.timedelta(seconds=s) for s in r['time']]
        tmp_coords.ticks = Ticktock(tvec_datetime)  # add ticks
        tmp_coords.sim_time = r['time']
        rays.append(tmp_coords)

    # Plot rays:
    lw = 1  # linewidth
    for name, r in zip(file_titles, rays):
        ax.plot(r.x / R_E, r.z / R_E, linewidth=lw, label=dir_str)

# -------- figure formatting --------
L_shells = [2, 3]  # Field lines to draw

# Plot the earth
earth = plt.Circle((0, 0), 1, color='0.5', alpha=1, zorder=100)
iono = plt.Circle((0, 0), (R_E + H_IONO) / R_E, color='c', alpha=0.5, zorder=99)

ax.add_artist(earth)
ax.add_artist(iono)

# Plot the fieldlines (dipole model; could use something more complex)
for L in L_shells:
    # Plot dipole field lines for both profile views
    lam = np.linspace(-80, 80, 181)
    L_r = L * pow(np.cos(lam * D2R), 2)
    Lx = L_r * np.cos(lam * D2R)
    Ly = L_r * np.sin(lam * D2R)
    ax.plot(Lx, Ly, color='r', linewidth=1, linestyle='dashed')  # Field line
    ax.plot(-Lx, Ly, color='r', linewidth=1, linestyle='dashed')  # Field line (other side)

# add in DSX and VPM
ax.plot(dsxorbit[:, 0] / R_E, dsxorbit[:, 2] / R_E, label='DSX')
ax.plot(vpmorbit[:, 0] / R_E, vpmorbit[:, 2] / R_E, label='VPM')

# tidy up axes
plt.xticks(np.arange(-3, 3, step=1))
plt.yticks(np.arange(-3, 3, step=1))
plt.xlabel('L (R$_E$)')
plt.ylabel('L (R$_E$)')

#plt.legend()

str_freq = str(int(freq / 1e3))
str_orbital_pos = str(orbital_pos)
fig_title = str_freq + ' kHz rays' + ' at position ' + str_orbital_pos + ' orbit in SM coordinates in XZ plane'
plt.title(fig_title)

# save it!
savename = 'plots/' + 'XZ_' + str_freq + 'kHz_pos' + str_orbital_pos + '_mode1.png'
plt.savefig(savename)

plt.show()
plt.close()