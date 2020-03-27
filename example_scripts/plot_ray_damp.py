# script to plot damping of rays

import numpy as np
import matplotlib.pyplot as plt
import os
import datetime as dt

# import functions from this example scripts directory
from raytracer_utils import readdump, read_rayfile, read_rayfiles, read_damp
from run_rays import run_rays

# Spacepy (for coordinate transforms)
from spacepy import coordinates as coord
from spacepy.time import Ticktock

from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm

# ------------------Constants --------------------------
R_E = 6371e3  # m

# ------------------ create Orbits --------------------------
from DSX_TLE import x_DSX, y_DSX, z_DSX

# create lists
freq = 10e3
n_pos = 0
positions = [np.array([x_DSX[n_pos],y_DSX[n_pos],z_DSX[n_pos]])]
directions = [np.array([0,0,0])]
run_rays(freq, positions, directions)

# ----------- PLOTTING -----------------
fig, ax = plt.subplots()

# -------------- Load output directory ----------
project_root = os.getcwd();
ray_out_dir = os.path.join(project_root, "test_outputs");

yearday = '2010001'  # YYYYDDD
milliseconds_day = 0  # milliseconds into the day
ray_datenum = dt.datetime(2010, 1, 1, 0, 0, 0)

# Load all the rayfiles in the output directory
d = os.listdir(ray_out_dir)
file_titles = ['example_ray_mode1']

damplist = []
for r in file_titles:
    damplist += read_damp(os.path.join(ray_out_dir, r + '.damp'))

# plot damping
for name, d in zip(file_titles, damplist):
    ax.plot(d["time"], d["damping"])

plt.xlabel('Time [s]')
plt.ylabel('Normalized wave power')
plt.title('Damping')
plt.legend()
plt.legend.linespaceing = 0.1
plt.legend.loc = 'upper right'

# save it!
savename = 'plots/1kHz_damp_mode1.png'
plt.savefig(savename)

plt.show()
plt.close()
