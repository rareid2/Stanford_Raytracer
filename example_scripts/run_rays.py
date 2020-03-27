# function to run the raytracer!

# inputs: freq, positions, directions
# freq: python list of ints in Hz
# positions: python list of arrays in SM cartesian coordinates
# directions: python list of arrays as unit vectors in cartesian coordinates

# for every position, every combination of frequency and direction will be launched
# total launched rays = (np.size(freq) * np.size(directions))* np.size(positions)

# import packages
import numpy as np  # for math
import matplotlib.pyplot as plt  # for plotting
import os  # for running commands
import datetime as dt  # for coordinate transforms

# Spacepy (for coordinate transforms)
from spacepy import coordinates as coord
from spacepy.time import Ticktock

def run_rays(freq, positions, directions):

    # Constants
    D2R = (np.pi / 180.0)
    R2D = 1.0 / D2R
    R_E = 6371e3  # m
    H_IONO = 1000e3

    # Simulation parameters
    t_max = 20       # Maximum duration in seconds
    dt0 = 1e-3       # Initial timestep in seconds
    dtmax = 0.1      # Maximum allowable timestep in seconds
    root = 2         # Which root of the Appleton-Hartree equation
                     # (1 = negative, 2 = positive)
                     # (2=whistler in magnetosphere)
    fixedstep = 0    # Don't use fixed step sizes, that's a bad idea.
    maxerr = 5.0e-4  # Error bound for adaptive timestepping
    maxsteps = 2e5   # Max number of timesteps (abort if reached)
    use_IGRF = 1     # Magnetic field model (1 for IGRF, 0 for dipole)
    use_tsyg = 1     # Use the Tsyganenko magnetic field model corrections
    minalt = R_E     # cutoff altitude in meters
    # TODO: make a max cutoff alt or check if inside the plasmasphere

    # Environmental parameters
    yearday = '2020001'   # YYYYDDD
    milliseconds_day = 0  # milliseconds into the day
    ray_datenum = dt.datetime(2020, 1, 1, 0, 0, 0)

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

    # Set up the output directory
    project_root = os.getcwd()  # grabs current full path
    ray_out_dir = os.path.join(project_root, "test_outputs")

    # Create directory for outputs if doesn't already exist
    print("output directory:", ray_out_dir)
    if not os.path.exists(ray_out_dir):
        os.mkdir(ray_out_dir)

    # Write the ray input file
    ray_inpfile = os.path.join(project_root, "ray_inpfile.txt")
    f = open(ray_inpfile, 'w')

    # Go through list of positions, write a new ray for every direction and freq at each pos
    for pos0 in positions:
        for dir0 in directions:
            for fr in freq:
                w0 = fr * 2.0 * np.pi
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
            # Test the Ngo model
            configfile = os.path.join(project_root, "newray_default.in")
            damp_mode = 0

            ray_cmd = base_cmd + ' --ngo_configfile="%s"' % (configfile)
            damp_cmd = base_damp_cmd + ' --mode %d' % damp_mode

        if mode == 2:
            # Test the full GCPM model
            damp_mode = 1

            ray_cmd = base_cmd + ' --gcpm_kp=%g' % (Kp)
            damp_cmd = base_damp_cmd + ' --mode %d' % damp_mode

        if mode == 3:
            # Test the uniformly-sampled GCPM model
            mode3_interpfile = os.path.join(project_root, 'precomputed_grids',
                                            'gcpm_kp4_2001001_L10_80x80x80_noderiv.txt')
            damp_mode = 1

            ray_cmd = base_cmd + ' --interp_interpfile="%s"' % (mode3_interpfile)
            damp_cmd = base_damp_cmd + ' --mode %d' % damp_mode

        if mode == 4:
            # Test the randomly-sampled GCPM model
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
            # Test the Simplified GCPM model
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

run_rays([1e3], [np.array([0,0,0])], [np.array([0,0,0])])