def run_them_rays(freq, orbital_pos):
    # import packages
    import numpy as np
    import matplotlib.pyplot as plt
    import os
    import datetime as dt

    # import functions from this example scripts directory
    #from orbits import orbits
    from TLE import x_DSX, y_DSX, z_DSX
    from raytracer_utils import readdump, read_rayfile, read_rayfiles, read_damp

    # Spacepy (for coordinate transforms)
    from spacepy import coordinates as coord
    from spacepy.time import Ticktock

    # ------------------Constants --------------------------
    D2R = (np.pi / 180.0)
    R2D = 1.0 / D2R
    R_E = 6371e3  # m
    H_IONO = 1000e3

    # -------------- Simulation parameters -----------------
    t_max = 5   # Maximum duration in seconds
    dt0 = 1e-3  # Initial timestep in seconds
    dtmax = 0.1  # Maximum allowable timestep in seconds
    root = 2  # Which root of the Appleton-Hartree equation
    # (1 = negative, 2 = positive)
    # (2=whistler in magnetosphere)
    fixedstep = 0  # Don't use fixed step sizes, that's a bad idea.
    maxerr = 5.0e-4  # Error bound for adaptive timestepping
    maxsteps = 2e5  # Max number of timesteps (abort if reached)

    use_IGRF = 1  # Magnetic field model (1 for IGRF, 0 for dipole)
    use_tsyg = 1  # Use the Tsyganenko magnetic field model corrections
    minalt = R_E + 300e3  # cutoff altitude in meters
    # TODO: make a max cutoff alt

    pos0 = np.array([x_DSX[orbital_pos],y_DSX[orbital_pos],z_DSX[orbital_pos]])

    # -------------- Starting ray parameters -------------
    # freq = 1e3  # Hz
    w0 = freq * 2.0 * np.pi
    # only pass in unit vectors
    # n = 200
    # dsxorbit, vpmorbit = orbits(n)
    # orbital_pos = 0
    # pos0 = dsxorbit[orbital_pos, 0:3]
    #pos0 = [12677107.852887288, 1094813.8049563468, 671478.3829286593]

    directions = []

    theta = 40  # deg
    for i in range(1, 22):
        dir = np.array([np.cos(np.deg2rad(theta + (5 * i))), 0, np.sin(np.deg2rad(theta + (5 * i)))])
        directions.append(dir)

    directions_neg = list()
    theta = 220  # deg
    for i in range(1, 22):
        dir = np.array([np.cos(np.deg2rad(theta + (5 * i))), 0, np.sin(np.deg2rad(theta + (5 * i)))])
        directions.append(dir)


    # -------------- Environmental parameters -------------
    yearday = '2010001'  # YYYYDDD
    milliseconds_day = 0  # milliseconds into the day
    ray_datenum = dt.datetime(2010, 1, 1, 0, 0, 0)

    Kp = 2
    AE = 1.6
    Pdyn = 4
    Dst = 1.0
    ByIMF = 0.0
    BzIMF = -5
    # Tsyganenko correction parameters
    W = [0.132, 0.303, 0.083, 0.070, 0.211, 0.308]  # Doesn't matter if we're not using Tsyg

    # Which plasmasphere models should we run?
    #   1 - Legacy (Ngo) model
    #   2 - GCPM (Accurate, but * s l o w * )
    #   3 - Uniformly interpolated precomputed model
    #   4 - Randomly interpolated precomputed model
    #   5 - (not real)
    #   6 - Simplified GCPM from Austin Sousa's thesis

    modes_to_do = [1]

    # Should we include a geometric focusing term in the damping?
    include_geom_factor = 0  # 1 for yes

    # -------------- Set up the output directory ----------
    project_root = os.getcwd()
    ray_out_dir = os.path.join(project_root, "test_outputs")

    print("output directory:", ray_out_dir)
    if not os.path.exists(ray_out_dir):
        os.mkdir(ray_out_dir)

    # ------ Write the ray input file ---
    ray_inpfile = os.path.join(project_root, "ray_inpfile.txt")

    f = open(ray_inpfile, 'w')

    for dir0 in directions:
        f.write('%1.15e %1.15e %1.15e %1.15e %1.15e %1.15e %1.15e\n' % (
            pos0[0], pos0[1], pos0[2], dir0[0], dir0[1], dir0[2], w0))
    f.close()

    # GCPM model and damping code needs to be run in the same directory
    # as the binary file (and all the misc data files)
    cwd = os.getcwd()
    os.chdir('../bin')

    for mode in modes_to_do:

        # The output file paths
        ray_outfile = os.path.join(ray_out_dir, 'example_ray_mode%d.ray' % mode)
        damp_outfile = os.path.join(ray_out_dir, 'example_ray_mode%d.damp' % mode)

        # The base command -- with parameters common for all modes
        base_cmd = './raytracer --outputper=%d --dt0=%g --dtmax=%g' % (1, dt0, dtmax) + \
                   ' --tmax=%g --root=%d --fixedstep=%d --maxerr=%g' % (t_max, root, fixedstep, maxerr) + \
                   ' --modelnum=%d --maxsteps=%d --minalt=%d --inputraysfile="%s"' % (
                       mode, maxsteps, minalt, ray_inpfile) + \
                   ' --outputfile="%s" --yearday=%s --milliseconds_day=%d' % (ray_outfile, yearday, milliseconds_day) + \
                   ' --use_tsyganenko=%d --use_igrf=%d --tsyganenko_Pdyn=%g' % (use_tsyg, use_IGRF, Pdyn) + \
                   ' --tsyganenko_Dst=%g --tsyganenko_ByIMF=%g --tsyganenko_BzIMF=%g' % (Dst, ByIMF, BzIMF) + \
                   ' --tsyganenko_W1=%g --tsyganenko_W2=%g --tsyganenko_W3=%g' % (W[0], W[1], W[2]) + \
                   ' --tsyganenko_W4=%g --tsyganenko_W5=%g --tsyganenko_W6=%g' % (W[3], W[4], W[5])

        base_damp_cmd = './damping --inp_file "%s" --out_file "%s" ' % (ray_outfile, damp_outfile) + \
                        ' --Kp %g --AE %g' % (Kp, AE) + \
                        ' --yearday %s --msec %d' % (yearday, milliseconds_day) + \
                        ' --geom_factor=%d' % include_geom_factor

        if mode == 1:
            # -------------- Test the Ngo model ------------------
            configfile = os.path.join(project_root, "newray_default.in")
            damp_mode = 0

            ray_cmd = base_cmd + ' --ngo_configfile="%s"' % (configfile)
            damp_cmd = base_damp_cmd + ' --mode %d' % damp_mode

        if mode == 2:
            # -------------- Test the full GCPM model ------------------
            damp_mode = 1

            ray_cmd = base_cmd + ' --gcpm_kp=%g' % (Kp)
            damp_cmd = base_damp_cmd + ' --mode %d' % damp_mode

        if mode == 3:
            # -------------- Test the uniformly-sampled GCPM model ------------------
            mode3_interpfile = os.path.join(project_root, 'precomputed_grids',
                                            'gcpm_kp4_2001001_L10_80x80x80_noderiv.txt')
            damp_mode = 1

            ray_cmd = base_cmd + ' --interp_interpfile="%s"' % (mode3_interpfile)
            damp_cmd = base_damp_cmd + ' --mode %d' % damp_mode

        if mode == 4:
            # -------------- Test the randomly-sampled GCPM model ------------------
            mode4_modelfile = os.path.join(project_root,
                                           'precomputed_grids', 'precomputed_model_gcpm_2010001_0_kp2_L10_random.dat')

            # Mode4 interpolation parameters:
            scattered_interp_window_scale = 1.2
            scattered_interp_order = 2
            scattered_interp_exact = 0  # Try 0 if there's weirdness at discontinuities
            scattered_interp_local_window_scale = 5

            damp_mode = 1

            ray_cmd = base_cmd + ' --interp_interpfile=%s' % (mode4_modelfile) + \
                      ' --scattered_interp_window_scale=%d' % (scattered_interp_window_scale) + \
                      ' --scattered_interp_order=%d' % (scattered_interp_order) + \
                      ' --scattered_interp_exact=%d' % (scattered_interp_exact) + \
                      ' --scattered_interp_local_window_scale=%d' % (scattered_interp_local_window_scale)
            damp_cmd = base_damp_cmd + ' --mode %d' % damp_mode

        if mode == 6:
            # -------------- Test the Simplified GCPM model ------------------
            MLT = 0
            fixed_MLT = 1  # Force the raytracer to stay in the meridonal plane?
            damp_mode = 0

            ray_cmd = base_cmd + ' --MLT="%g" --fixed_MLT=%g --kp=%g' % (MLT, fixed_MLT, Kp)
            damp_cmd = base_damp_cmd + ' --mode %d' % damp_mode

        # Run it!

        print("------- Running mode %d -------" % mode)
        print("Command is:")
        print(ray_cmd)
        print()

        os.system(ray_cmd)

        print("------- Running damping, mode %d -------" % damp_mode)

        print(damp_cmd)
        os.system(damp_cmd)

    # Move back to the working directory
    os.chdir(cwd)

    # ----------- PLOTTING -----------------
    fig, ax = plt.subplots()

    # Load all the rayfiles in the output directory
    d = os.listdir(ray_out_dir)
    file_titles = ['example_ray_mode1']

    raylist = []

    for r in file_titles:
        raylist += read_rayfile(os.path.join(ray_out_dir, r + '.ray'))

    rays = []
    for r in raylist:
        tmp_coords = coord.Coords(list(zip(r['pos'].x, r['pos'].y, r['pos'].z)), 'SM', 'car', units=['m', 'm', 'm'])
        tvec_datetime = [ray_datenum + dt.timedelta(seconds=s) for s in r['time']]
        tmp_coords.ticks = Ticktock(tvec_datetime)  # add ticks
        tmp_coords.sim_time = r['time']
        rays.append(tmp_coords)

    # TODO: I took out a line here that was parsing throuhg multiple files, may want to add it back
    # Plot rays:
    lw = 1  # linewidth
    for r, dir0 in zip(rays, directions):
        dir_str = int(np.round(np.rad2deg(np.arcsin(dir0[0]))))
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
    plt.plot(x_DSX/(R_E),z_DSX/(R_E), c='g')
    #ax.plot(vpmorbit[:, 0] / R_E, vpmorbit[:, 2] / R_E, label='VPM')

    # tidy up axes
    plt.xticks(np.arange(-3, 3, step=1))
    plt.yticks(np.arange(-3, 3, step=1))
    plt.xlabel('L (R$_E$)')
    plt.ylabel('L (R$_E$)')

    plt.legend()

    str_freq = str(int(freq / 1e3))
    str_orbital_pos = str(orbital_pos)
    fig_title = str_freq + ' kHz rays' + ' at position ' + str_orbital_pos + ' orbit in SM coordinates in XZ plane'
    plt.title(fig_title)

    # save it!
    savename = 'plots/' + 'XZ_' + str_freq + 'kHz_pos' + str_orbital_pos + '_mode1.png'
    plt.savefig(savename)

    plt.show()
    plt.close()

    """
    # check if conjunction
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='polar')

    r_dist = []
    c_dir = []

    condition = 500e3  # m

    
    print('For DSX at ', dsxorbit[orbital_pos, 0], dsxorbit[orbital_pos, 2], '[m]:')
    for r, dir0 in zip(rays, directions):
        dir_str = int(np.round(np.rad2deg(np.arcsin(dir0[0]))))
        for m in range(n):
            dist = np.sqrt((r.x[-1] - vpmorbit[m, 0]) ** 2 + (r.z[-1] - vpmorbit[m, 2]) ** 2)
            if dist <= condition:
                r_dist.append(dist / condition)
                c_dir.append(np.deg2rad(dir_str))
                print(dir_str, 'degree ray conjunction at ', vpmorbit[m, 0], vpmorbit[m, 2], '[m],  dist= ', dist,
                      ' [m]')
    
    ax.scatter(c_dir, r_dist)
    fig_title = str_freq + ' kHz rays' + ' at position ' + str_orbital_pos + ' orbit in SM coordinates in XZ plane'
    plt.title(fig_title)

    # save it!
    savename = 'plots/' + 'polar_' + str_freq + 'kHz_pos' + str_orbital_pos + '_mode1.png'
    plt.savefig(savename)

    plt.show()
    plt.close()
    """


# run it for freq range of interest at apogee and perigee of DSX
#for f in range(30, 31):
for o in [400, 425, 450, 475, 500, 525, 550, 575, 0, 25, 50, 75, 175, 200, 225, 250, 300, 325]:
    run_them_rays(26e3, o)

""" damping
# import packages
import numpy as np
import matplotlib.pyplot as plt
import os
import datetime as dt

# import functions from this example scripts directory
from orbits import orbits
from raytracer_utils import readdump, read_rayfile, read_rayfiles, read_damp

# Spacepy (for coordinate transforms)
from spacepy import coordinates as coord
from spacepy.time import Ticktock

fig2, ax2 = plt.subplots(1, 1)

# Load n' plot the Ngo files:
project_root = os.getcwd()
ray_out_dir = os.path.join(project_root, "test_outputs")

# Load all the rayfiles in the output directory
d = os.listdir(ray_out_dir)
file_titles = ['example_ray_mode1']

damplist = []
for r in file_titles:
    damplist += read_damp(os.path.join(ray_out_dir, r + '.damp'))

# plot damping
for name, d in zip(file_titles, damplist):
    ax2.plot(d["time"], d["damping"])

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
"""
