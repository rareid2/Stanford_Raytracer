
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

<<<<<<< HEAD
import xlsxwriter


workbook = xlsxwriter.Workbook('more_conjunctions.xlsx')
worksheet = workbook.add_worksheet()

worksheet.write(0, 0, 'date')
worksheet.write(1, 0, 'freq')
worksheet.write(2, 0, 'angle')
worksheet.write(3, 0, 'DSXpos')
worksheet.write(4, 0, 'VPMpos')
worksheet.write(6, 0, 'rx')
worksheet.write(6, 1, 'ry')
worksheet.write(6, 2, 'rz')
worksheet.write(6, 3, 'amp')
worksheet.write(6, 4, 'Lx')
worksheet.write(6, 5, 'Ly')
worksheet.write(6, 6, 'Lz')
worksheet.write(6, 7, '-Lx')
worksheet.write(6, 8, '-Ly')
worksheet.write(6, 9, '-Lz')


# TODO: coordinates?
# --------------------------- Orbits -------------------------------
# DSX TLE:
line1 = '1 44344U 19036F   20116.38366941 -.00000011 +00000-0 +00000+0 0  9990'
line2 = '2 44344 042.2529 091.9758 1974961 131.2888 247.5136 04.54371596013862'

=======
# TODO: coordinates?
# --------------------------- Orbits -------------------------------
# DSX TLE:
line1 = '1 44344U 19036F   20113.30349832 -.00000013 +00000-0 +00000-0 0  9992'
line2 = '2 44344 042.2517 093.1041 1975016 129.9693 249.1586 04.54371641013721'
>>>>>>> a40385af70d8379ef79a0df25c5c9dc386b09bc8
DSX_pos, DSX_t = get_TLE(line1, line2, 'DSX')
x_DSX, y_DSX, z_DSX = DSX_pos[0], DSX_pos[1], DSX_pos[2]

# VPM TLE:
<<<<<<< HEAD
line1 = '1 45120U 19071K   20116.86188762  .00004015  00000-0  13248-3 0  9997'
line2 = '2 45120  51.6428 257.3005 0011715 220.9363 139.0740 15.33712868012946'
=======
line1 = '1 45120U 19071K   20113.86456076  .00004088  00000-0  13479-3 0  9996'
line2 = '2 45120  51.6417 271.7875 0011770 209.5167 150.5152 15.33686764012484'
>>>>>>> a40385af70d8379ef79a0df25c5c9dc386b09bc8
VPM_pos, VPM_t = get_TLE(line1, line2, 'VPM')
x_VPM, y_VPM, z_VPM = VPM_pos[0], VPM_pos[1], VPM_pos[2]

# --------------------------- Ray Tracing --------------------------
# define lists here - must be lists even if only one arg
<<<<<<< HEAD
# started at midnight on the 29th
freq = [8.2e3]
#n_pos = np.linspace(0,49,1)

positions = []

orbitdata = np.genfromtxt('orbit_data.txt')

apa_new = np.genfromtxt('apa.txt')
apa_new = [int(apa) for apa in apa_new]
apa_new = apa_new[0::100]
mytime = []
vpmspots = []
#apa_new = [1,2,3,4]
for ap in apa_new:
    pos = orbitdata[ap]
    vpmspot = pos[3:6]
    pos = pos[0:3]
    positions.append(pos)
    vpmspots.append(vpmspot)
    mytime.append(ray_datenum + dt.timedelta(minutes=ap))

print('got all positions')

mypositions = positions

#positions = [np.array([x_DSX[int(npos)], y_DSX[int(npos)], z_DSX[int(npos)]]) for npos in n_pos]
looplen = len(mypositions)
looplen = 2
for numcount in range(looplen):
    positions = [mypositions[numcount]]
    ray_datenum = mytime[numcount]
    """
    thetalist is angle from local magnetic field direction in XZ plane
    for example: 0 deg = parallel Bfield 
    15 = 15 deg counterclockwise to Bfield
    """

    thetalist = [] # in deg

    # initialize - leave empty
    directions = []
    Blines = []

    for position in positions:
        # grab position
        startpoint = [position[0]/R_E, position[1]/R_E, position[2]/R_E]
        # get bfield direction
        Bx, By, Bz = B_dir(0, startpoint, 0, '0', 1)
        dirB = np.reshape(np.array([Bx, By, Bz]), (1,3))

        # grab full magnetic field line for plotting later
        Blines.append(trace_fieldline_ODE(startpoint, 0, '0', 1))
        Blines.append(trace_fieldline_ODE(startpoint, 0, '0', -1))

        # rotate around direction of field line around x axis
        for theta in thetalist:
            R = [ [1, 0, 0], [0, np.cos(D2R * theta), - np.sin(D2R * theta)],
                  [0, np.sin(D2R * theta), np.cos(D2R * theta)] ]
            direction = np.matmul(dirB, np.reshape(np.array(R), (3, 3)))
            direction = direction/np.linalg.norm(direction)
            # add that normalized direction
            directions.append(np.squeeze(direction))

    # number of rays
    # n_rays = len(freq) * len(positions) * len(directions)
    # print('about to run: ', n_rays, ' rays')

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
    fig, ax = plt.subplots(1,1, sharex=True, sharey=True)
    lw = 2  # linewidth

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

    def myplot(ax, xs, ys, zs, cmap):
        for x, y, z in zip(xs, ys, zs):
            plot = LineCollection([np.column_stack((x, y))], cmap=cmap, zorder=102)
            plot.set_array(z)
            ax.add_collection(plot)
        return plot

    line = myplot(ax, rx, rz, dlist, 'Reds')
    fig.colorbar(line, ax=ax, label = 'Normalized wave power')

    # --------------------------- Figure formatting ---------------------------
    L_shells = [2, 3, 4]  # Field lines to draw

    # Earth and Iono
    earth = plt.Circle((0, 0), 1, color='b', alpha=1, zorder=100)
    iono = plt.Circle((0, 0), (R_E + H_IONO) / R_E, color='g', alpha=0.5, zorder=99)
    ax.add_artist(earth)
    ax.add_artist(iono)

    # ---------------------- Field Lines -----------------------
    # (from IGRF13 model)
    for L in L_shells:
        Lx, Ly, Lz = trace_fieldline_ODE([L,0,0], 0, '0', 1)
        plt.plot(Lx, Lz, color='b', linewidth=1, linestyle='dashed')
        Lx, Ly, Lz = trace_fieldline_ODE([L,0,0], 0, '0', -1)
        plt.plot(Lx, Lz, color='b', linewidth=1, linestyle='dashed')
        Lx, Ly, Lz = trace_fieldline_ODE([-L,0,0], 0, '0', 1)
        plt.plot(Lx, Lz, color='b', linewidth=1, linestyle='dashed')
        Lx, Ly, Lz = trace_fieldline_ODE([-L,0,0], 0, '0', -1)
        plt.plot(Lx, Lz, color='b', linewidth=1, linestyle='dashed')
    print('finished plotting field lines')

    # plot field line from orbital position
    for blinex, bliney, blinez in Blines:
        plt.plot(blinex, blinez, color='r', linewidth=1, linestyle='dashed')

    # ---------------------- Satellite Orbits   -----------------------------
    plt.plot(position[0] / R_E, position[2] / R_E, 'r', zorder = 105, label = 'DSX')
    plt.plot(vpmspots[numcount][0] / R_E, vpmspots[numcount][2] / R_E, 'r', zorder = 105, label = 'VPM')
    #for npos in n_pos:
    plt.plot(x_DSX / R_E, z_DSX / R_E, 'y', zorder = 103)
    plt.plot(x_VPM / R_E, z_VPM / R_E, 'y', zorder = 103)

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
    fig_title = str(freq[0]/1e3) + ' kHz rays in XZ plane'
    plt.title(fig_title)

    # ------------------------------- Saving ---------------------------------------
    #savename = 'plots/XZ_' + str_freq + 'kHz_%03d.png' %p
    savename = 'plots/ray_' + str(ray_datenum) + '.png'
    fig.savefig(savename, format='png')

    plt.show()
    plt.close()

    #vpmspot = VPMdata[apa_new[numcount],3:6]

    worksheet.write(0, (11 * numcount) + 1, str(ray_datenum))
    worksheet.write(1, (11 * numcount) + 1, freq[0])
    worksheet.write(2, (11 * numcount) + 1, thetalist[0])
    worksheet.write(3, (11 * numcount) + 1, positions[0][0])
    worksheet.write(3, (11 * numcount) + 2, positions[0][1])
    worksheet.write(3, (11 * numcount) + 3, positions[0][2])
    worksheet.write(4, (11 * numcount) + 1, vpmspots[0][0])
    worksheet.write(4, (11 * numcount) + 2, vpmspots[0][1])
    worksheet.write(4, (11 * numcount) + 3, vpmspots[0][2])

    # Iterate over the data and write it out row by row.
    for i in range(len(rx[0])):
        worksheet.write(i + 7, (11 * numcount), rx[0][i])
        worksheet.write(i + 7, (11 * numcount) + 1, ry[0][i])
        worksheet.write(i + 7, (11 * numcount) + 2, rz[0][i])
        worksheet.write(i + 7, (11 * numcount) + 3, dlist[0][i])
        for m in range(len(Blines[0])):
            worksheet.write(i + 7, (11 * numcount) + 4, Blines[0][0][m])
            worksheet.write(i + 7, (11 * numcount) + 5, Blines[0][1][m])
            worksheet.write(i + 7, (11 * numcount) + 6, Blines[0][2][m])
        for n in range(len(Blines[1])):
            worksheet.write(i + 7, (11 * numcount) + 7, Blines[1][0][n])
            worksheet.write(i + 7, (11 * numcount) + 8, Blines[1][1][n])
            worksheet.write(i + 7, (11 * numcount) + 9, Blines[1][2][n])
    print('wrote to xls')

workbook.close()  # mistake here
=======

freq = [8.2e3]
n_pos = np.linspace(0,49,2)
positions = [np.array([x_DSX[int(npos)], y_DSX[int(npos)], z_DSX[int(npos)]]) for npos in n_pos]

"""
thetalist is angle from local magnetic field direction in XZ plane
for example: 0 deg = parallel Bfield 
15 = 15 deg counterclockwise to Bfield
"""

thetalist = [0] # in deg

# initialize - leave empty
directions = []
Blines = []

for position in positions:
    # grab position
    startpoint = [position[0]/R_E, position[1]/R_E, position[2]/R_E]
    # get bfield direction
    Bx, By, Bz = B_dir(0, startpoint, 0, '0', 1)
    dirB = np.reshape(np.array([Bx, By, Bz]), (1,3))

    # grab full magnetic field line for plotting later
    Blines.append(trace_fieldline_ODE(startpoint, 0, '0', 1))
    Blines.append(trace_fieldline_ODE(startpoint, 0, '0', -1))

    # rotate around direction of field line around x axis
    for theta in thetalist:
        R = [ [1, 0, 0], [0, np.cos(D2R * theta), - np.sin(D2R * theta)],
              [0, np.sin(D2R * theta), np.cos(D2R * theta)] ]
        direction = np.matmul(dirB, np.reshape(np.array(R), (3, 3)))
        direction = direction/np.linalg.norm(direction)
        # add that normalized direction
        directions.append(np.squeeze(direction))

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
    #if len(damp) < max(r_length):
    #    leftover = max(r_length) - len(damp)
    #    damp = np.concatenate((damp, np.zeros(int(leftover))), axis=0)
    dlist.append(damp)

# normalize
normdlist = [float(i)/sum(dlist[0]) for i in dlist[0]]

def myplot(ax, xs, ys, zs, cmap):
    for x, y, z in zip(xs, ys, zs):
        plot = LineCollection([np.column_stack((x, y))], cmap=cmap, zorder=102)
        plot.set_array(z)
        ax.add_collection(plot)
    return plot

line = myplot(ax, rx, rz, normdlist, 'Reds')
fig.colorbar(line, ax=ax, label = 'Normalized wave power')

# --------------------------- Figure formatting ---------------------------
L_shells = [2, 3, 4]  # Field lines to draw

# Earth and Iono
earth = plt.Circle((0, 0), 1, color='b', alpha=1, zorder=100)
iono = plt.Circle((0, 0), (R_E + H_IONO) / R_E, color='g', alpha=0.5, zorder=99)
ax.add_artist(earth)
ax.add_artist(iono)

# ---------------------- Field Lines -----------------------
# (from IGRF13 model)
for L in L_shells:
    Lx, Ly, Lz = trace_fieldline_ODE([L,0,0], 0, '0', 1)
    plt.plot(Lx, Lz, color='b', linewidth=1, linestyle='dashed')
    Lx, Ly, Lz = trace_fieldline_ODE([L,0,0], 0, '0', -1)
    plt.plot(Lx, Lz, color='b', linewidth=1, linestyle='dashed')
    Lx, Ly, Lz = trace_fieldline_ODE([-L,0,0], 0, '0', 1)
    plt.plot(Lx, Lz, color='b', linewidth=1, linestyle='dashed')
    Lx, Ly, Lz = trace_fieldline_ODE([-L,0,0], 0, '0', -1)
    plt.plot(Lx, Lz, color='b', linewidth=1, linestyle='dashed')
print('finished plotting field lines')

# plot field line from orbital position
for blinex, bliney, blinez in Blines:
    plt.plot(blinex, blinez, color='r', linewidth=1, linestyle='dashed')

# ---------------------- Satellite Orbits   -----------------------------
plt.plot(x_DSX/R_E, z_DSX/R_E, c='y', zorder = 105, label = 'DSX')
plt.plot(x_VPM/R_E, z_VPM/R_E, c='y', zorder = 105, label = 'VPM')
for npos in n_pos:
    plt.plot(x_DSX[int(npos)]/R_E, z_DSX[int(npos)]/R_E, '-bo', zorder = 103)
    plt.plot(x_VPM[int(npos)] / R_E, z_VPM[int(npos)] / R_E, '-ro', zorder=103)
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
fig_title = str(freq[0]/1e3) + ' kHz rays in XZ plane'
plt.title(fig_title)

# ------------------------------- Saving ---------------------------------------
#savename = 'plots/XZ_' + str_freq + 'kHz_%03d.png' %p
savename = 'plots/ray_' + str(ray_datenum) + '.png'
fig.savefig(savename, format='png')

fig.show()
plt.close()
>>>>>>> a40385af70d8379ef79a0df25c5c9dc386b09bc8
