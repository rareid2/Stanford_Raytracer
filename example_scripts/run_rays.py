"""

function to run the raytracer!

inputs: freq, positions, directions
freq: python list of ints in Hz
positions: python list of arrays in SM cartesian coordinates
directions: python list of arrays as unit vectors in cartesian coordinates

for every position, every combination of frequency and direction will be launched
total launched rays = (np.size(freq) * np.size(directions))* np.size(positions)

"""

# import packages
import numpy as np        
import os                        

# get settings
from raytracer_settings import *

def run_rays(freq, positions, directions, yearday, milliseconds_day, tmpdir):

    #  ---------------------------- SET UP FILE SETTINGS  ------------------------------
    project_root = os.getcwd()  # grabs current full path
    tempdirectory = tmpdir
    # tempdirectory = os.getcwd() 
    radstr = ''.join(random.choices(string.ascii_lowercase, k=16))
    radstr = 'ray_inpfile.txt'

    # Set input file path
    ray_inpfile = os.path.join(tempdirectory, radstr)
    # uncomment this if you want to be able to traceback the input and output files
    print(f'tempfile is {ray_inpfile}')

    # Set config file for Ngo Model (Mode 1)
    configfile = os.path.join(project_root, "newray_default.in")

    # Set config file for mode 3
    mode3_interpfile = os.path.join(project_root, 'precomputed_grids',
                                                'gcpm_kp4_2001001_L10_80x80x80_noderiv.txt')
    # Set config file for mode 4
    mode4_modelfile = os.path.join(project_root,
                                'precomputed_grids', 'precomputed_model_gcpm_2010001_0_kp2_L10_random.dat')

    # Set up the output directory
    ray_out_dir = tempdirectory

    # Create directory for outputs if doesn't already exist
    # print("output directory:", ray_out_dir)
    #if not os.path.exists(ray_out_dir):
    #    os.mkdir(ray_out_dir)

    # ---------------------------- END SET UP FILE SETTINGS  ------------------------------

    # ------------------------------- START THE RAYTRACER  --------------------------------

    # Write the ray input file
    f = open(ray_inpfile, 'w')

    # Go through list of positions, write a new ray for every direction and freq at each pos
    for pos0, dir0 in zip(positions, directions):
        for fr in freq:
            w0 = fr * 2.0 * np.pi
            f.write('%1.15e %1.15e %1.15e %1.15e %1.15e %1.15e %1.15e\n' % (
                pos0[0], pos0[1], pos0[2], dir0[0], dir0[1], dir0[2], w0))
    f.close()

    # GCPM model and damping code needs to be run in the same directory
    # as the binary file (and all the misc data files)
    cwd = os.getcwd()
    # os.chdir('example_scripts')
    os.chdir('../bin')

    for mode in modes_to_do:

        # Set output file path
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
            damp_mode = 1

            ray_cmd = base_cmd + ' --interp_interpfile="%s"' % (mode3_interpfile)
            damp_cmd = base_damp_cmd + ' --mode %d' % damp_mode

        if mode == 4:
            # Test the randomly-sampled GCPM model
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

        #print("------- Running mode %d -------" % mode)
        #print("Command is:")
        print(ray_cmd)
        #print()

        os.system(ray_cmd)

        #print("------- Running damping, mode %d -------" % damp_mode)

        print(damp_cmd)
        os.system(damp_cmd)

    # Move back to the working directory
    os.chdir(cwd)

    # print('raytracer done')

#  ------------------------------- END THE RAYTRACER  --------------------------------