# script to plot trajectory of rays in XZ coordinates with attenuation colorbar

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

# for color bar plotting
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm

# ------------------ Constants --------------------------
D2R = (np.pi / 180.0)
R2D = 1.0 / D2R
R_E = 6371e3  # m
H_IONO = 1000e3

# ------------------ Orbits --------------------------
from DSX_TLE import x_DSX, y_DSX, z_DSX, slope

# ------------------ Ray Tracing --------------------------
# create lists
freq = 26e3
p = 63
iterate_pos = np.linspace(300, 600, 6)

for n_pos in iterate_pos:
    n_pos = int(n_pos)

    positions = [np.array([x_DSX[n_pos], y_DSX[n_pos], z_DSX[n_pos]])]

    """
    nearx = int(x_DSX[n_pos] / R_E)
    nearz = int(z_DSX[n_pos] / R_E)
    # figure out angle between ray start and field line
    L_shells = [nearx]  # Field lines to draw
    for L in L_shells:
        # dipole field lines for both profile views
        lam = np.linspace(-80, 80, 181)
        L_r = L * pow(np.cos(lam * D2R), 2)
        Lx = L_r * np.cos(lam * D2R)
        Lz = L_r * np.sin(lam * D2R)

    idx = (np.abs(Lx -(x_DSX[n_pos] / R_E))).argmin()
    idz = (np.abs(Lz - (z_DSX[n_pos] / R_E))).argmin()

    btheta = np.arctan2((Lz[idz] - Lz[idz + 1]), (Lx[idx] - Lz[idx + 1]))
    btheta = np.rad2deg(btheta)
    

    if 0 < n_pos < 10:
        btheta = 0
    if 10 < n_pos < 175:
        btheta = 180
    directionlist = []
    thetastep = 10

    for theta in range(45, 315, thetastep):
        dir = np.array([np.cos(np.deg2rad(np.abs(theta - btheta))), 0, np.sin(np.deg2rad(np.abs(theta - btheta)))])
        if theta <= 135:
            directionlist.append(dir)
        if theta >= 225:
            directionlist.append(dir)

    #directionlist[int(45 / thetastep)] = np.array([0, 0, 0])
    """

    btheta = np.rad2deg(np.arctan(slope[n_pos]))
    directionlist = []
    thetastep = 15


    for theta in range(45, 315, thetastep):
        dir = np.array([np.cos(np.deg2rad(theta - btheta)), 0, np.sin(np.deg2rad(theta - btheta))])
        if theta <= 135:
            directionlist.append(dir)
        if theta >= 225:
            directionlist.append(dir)

    for direction in directionlist:
        directions = [direction]
        # run!
        run_rays(freq, positions, directions)

        # -------------- Load output directory ----------
        project_root = os.getcwd();
        ray_out_dir = os.path.join(project_root, "test_outputs");

        yearday = '2010001'  # YYYYDDD
        milliseconds_day = 0  # milliseconds into the day
        ray_datenum = dt.datetime(2010, 1, 1, 0, 0, 0)

        # Load all the rayfiles in the output directory
        d = os.listdir(ray_out_dir)
        file_titles = ['example_ray_mode1']

        raylist = []

        for r in file_titles:
            raylist += read_rayfile(os.path.join(ray_out_dir, r + '.ray'))

        damplist = []

        for r in file_titles:
            damplist += read_damp(os.path.join(ray_out_dir, r + '.damp'))

        # ------------------ Coordinate Conversion --------------------------
        rays = []
        for r in raylist:
            tmp_coords = coord.Coords(list(zip(r['pos'].x, r['pos'].y, r['pos'].z)), 'SM', 'car', units=['m', 'm', 'm'])
            tvec_datetime = [ray_datenum + dt.timedelta(seconds=s) for s in r['time']]
            tmp_coords.ticks = Ticktock(tvec_datetime)  # add ticks
            tmp_coords.sim_time = r['time']
            rays.append(tmp_coords)

        #------------------ Plot rays --------------------------
        fig, ax = plt.subplots(1,1, sharex=True, sharey=True)
        lw = 2  # linewidth

        for r, dir0 in zip(rays, directions):
            rx = []
            rz = []

            if np.count_nonzero(dir0) == 0:
                dir_str = 'field-aligned'
            else:
                angdiff = np.abs(btheta - np.rad2deg(np.arctan(dir0[2] / dir0[0])))
                if angdiff >= 180:
                    dir_str = int(360 - angdiff)
                else:
                    dir_str = int(angdiff)
            dir_str = 'ray'
            rx.append(r.x / R_E)
            rz.append(r.z / R_E)

            damp = []
            for d in damplist:
                damp.append(d["damping"])

            if np.size(damp) - np.size(rx) > 0:
                damp = damp[0:np.size(rx)]
            elif np.size(damp) - np.size(rx) < 0:
                damplast = damp[-1]
                for n in range(np.size(damp), np.size(rx)):
                    damp.append(damplast)

            rx = np.array(rx)
            rz = np.array(rz)

            damp = np.squeeze(np.array(damp))

            points = np.array([rx, rz]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)

            lc = LineCollection(segments, cmap='Spectral', label = dir_str)
            lc.set_array(damp)
            lc.set_linewidth(lw)
            line = ax.add_collection(lc)
            fig.colorbar(line, ax=ax, label = 'Normalized wave power')

        # -------- figure formatting --------
        L_shells = [2, 3, 4, 5]  # Field lines to draw

        # -------- Earth and Iono --------
        earth = plt.Circle((0, 0), 1, color='0.5', alpha=1, zorder=100)
        iono = plt.Circle((0, 0), (R_E + H_IONO) / R_E, color='c', alpha=0.5, zorder=99)

        ax.add_artist(earth)
        ax.add_artist(iono)

        # -------- fieldlines -------- (dipole model; could use something more complex)
        for L in L_shells:
            # Plot dipole field lines for both profile views
            lam = np.linspace(-80, 80, 181)
            L_r = L * pow(np.cos(lam * D2R), 2)
            Lx = L_r * np.cos(lam * D2R)
            Lz = L_r * np.sin(lam * D2R)
            ax.plot(Lx, Lz, color='r', linewidth=1, linestyle='dashed')  # Field line
            ax.plot(-Lx, Lz, color='r', linewidth=1, linestyle='dashed')  # Field line (other side)

        # -------- DSX --------
        plt.plot(x_DSX/R_E, z_DSX/R_E, c='g', label = 'DSX')

        # -------- plasmapause --------
        lam = np.linspace(-80, 80, 181)
        L_r = 5 * pow(np.cos(lam * D2R), 2)
        Lx = L_r * np.cos(lam * D2R)
        Ly = L_r * np.sin(lam * D2R)
        ax.plot(Lx, Ly, color='b', linewidth=3, linestyle='dashed')  # Field line
        ax.plot(-Lx, Ly, color='b', linewidth=3, linestyle='dashed')  # Field line (other side)

        """
        # add in plasmasphere
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
        plt.pcolormesh(px, py, np.log10(Ne_xz))
        """

        # -------- figure formatting --------
        ax.set_aspect('equal')

        plt.xticks(np.arange(-6, 6, step=1))
        plt.yticks(np.arange(-6, 6, step=1))
        plt.xlabel('L (R$_E$)')
        plt.ylabel('L (R$_E$)')
        plt.xlim([-6, 6])
        plt.ylim([-2.5, 2.5])

        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.25), ncol=2)

        str_freq = str(int(freq / 1e3))
        str_orbital_pos = str(n_pos)
        fig_title = str_freq + ' kHz rays \n in SM coordinates in XZ plane'
        plt.title(fig_title)

        # -------- saving --------
        savename = 'plots/XZ_' + str_freq + 'kHz_%03d.png' %p
        fig.savefig(savename)

        fig.show()
        plt.close()

        # iterate
        p += 1