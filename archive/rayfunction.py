## this is a script to run_example_rays as a function

# ------------------ Inputs --------------------------
#t_max Maximum duration in seconds
#inp_lat = 21;             # launch latitude (geomagnetic)
#inp_lon = 104;            # launch longitude (geomagnetic)
#inp_alt = (900+R_E)*1e3   # launch altitude in meters
#freq = 440;               # ray frequency in Hz. How about concert A?
#launch_direction = "field-aligned"

#example call run_rays_func(20,21,100,(900+R_E)*1e3,440,[0,1,0])
#output = build the required files to be parsed and plotted

def run_rays_func(t_max,inp_lat,inp_lon,inp_alt,freq,launch_direction):

    import numpy as np
    import subprocess
    import os
    import datetime

    # Spacepy (for coordinate transforms)

    from spacepy import coordinates as coord
    from spacepy.time import Ticktock

    # ------------------Constants --------------------------
    R_E = 6371 # in km

    # -------------- Simulation parameters -----------------

    dt0 = 1e-3      # Initial timestep in seconds
    dtmax = 0.1     # Maximum allowable timestep in seconds
    root = 2        # Which root of the Appleton-Hartree equation
                    # (1 = negative, 2 = positive)
                    # (2=whistler in magnetosphere)
    fixedstep = 0   # Don't use fixed step sizes, that's a bad idea.
    maxerr = 5.0e-4 # Error bound for adaptive timestepping
    maxsteps = 2e5  # Max number of timesteps (abort if reached)
    use_IGRF = 0    # Magnetic field model (1 for IGRF, 0 for dipole)
    use_tsyg = 0    # Use the Tsyganenko magnetic field model corrections
    minalt   = (R_E + 400)*1e3 # cutoff altitude in meters

    # -------------- Environmental parameters -------------
    # setting this for the first day of 2010, env parameters should be accurate

    yearday = '2010001'		#YYYYDDD
    milliseconds_day = 0;	# milliseconds into the day
    ray_datenum = datetime.datetime(2010, 1, 1, 0, 0, 0);

    Kp    = 5
    AE    = 1.6
    Pdyn  = 4
    Dst   = 1.0
    ByIMF = 0.0
    BzIMF = -5
    # Tsyganenko correction parameters
    W = [0.132,    0.303,    0.083,    0.070,    0.211,    0.308 ]  # Doesn't matter if we're not using Tsyg

    # Which plasmasphere models should we run?
    #   1 - Legacy (Ngo) model
    #   2 - GCPM (Accurate, but * s l o w * )
    #   3 - Uniformly interpolated precomputed model
    #   4 - Randomly interpolated precomputed model
    #   5 - (not real)
    #   6 - Simplified GCPM from Austin Sousa's thesis

    modes_to_do = [6] #try mode 6 for Kp and MLT

    # Should we include a geometric focusing term in the damping?
    include_geom_factor = 1 # 1 for yes

    # -------------- Set up the output directory ----------
    project_root = os.getcwd();
    ray_out_dir = os.path.join(project_root, "test_outputs");

    #print("output directory:",ray_out_dir)
    if not os.path.exists(ray_out_dir):
    	os.mkdir(ray_out_dir)


    # ------ Write the ray input file ---

    ray_inpfile  = os.path.join(project_root,"ray_inpfile.txt")

    f = open(ray_inpfile,'w')

    # --- Rotate to SM using spacepy coordinate transform ---
    tmp_coords = coord.Coords((inp_alt, inp_lat, inp_lon),'MAG','sph',units=['m','deg','deg'])
    tmp_coords.ticks = Ticktock(ray_datenum) # add ticks
    pos0 = tmp_coords.convert('SM','car')

    nonzerocounter = np.count_nonzero(launch_direction)
    if nonzerocounter!=0:
        dir0 = launch_direction/(np.linalg.norm(launch_direction))
    else:
        dir0 = np.zeros(3)

    w0   = freq*2.0*np.pi

    f.write('%1.15e %1.15e %1.15e %1.15e %1.15e %1.15e %1.15e\n'%(pos0.x, pos0.y, pos0.z, dir0[0], dir0[1], dir0[2], w0))
    f.close()

    # GCPM model and damping code needs to be run in the same directory
    # as the binary file (and all the misc data files)
    cwd = os.getcwd();
    os.chdir('../bin')

    for mode in modes_to_do:

        # The output file paths
        ray_outfile  = os.path.join(ray_out_dir,'example_ray_mode%d.ray'%mode)
        damp_outfile = os.path.join(ray_out_dir,'example_ray_mode%d.damp'%mode)


        # The base command -- with parameters common for all modes
        base_cmd= './raytracer --outputper=%d --dt0=%g --dtmax=%g'%(1, dt0, dtmax) + \
                 ' --tmax=%g --root=%d --fixedstep=%d --maxerr=%g'%(t_max, root, fixedstep, maxerr) + \
                 ' --modelnum=%d --maxsteps=%d --minalt=%d --inputraysfile="%s"'%(mode,maxsteps, minalt, ray_inpfile) + \
                 ' --outputfile="%s" --yearday=%s --milliseconds_day=%d'%(ray_outfile, yearday, milliseconds_day) + \
                 ' --use_tsyganenko=%d --use_igrf=%d --tsyganenko_Pdyn=%g'%(use_tsyg, use_IGRF, Pdyn) + \
                 ' --tsyganenko_Dst=%g --tsyganenko_ByIMF=%g --tsyganenko_BzIMF=%g'%( Dst, ByIMF, BzIMF ) + \
                 ' --tsyganenko_W1=%g --tsyganenko_W2=%g --tsyganenko_W3=%g'%(W[0], W[1], W[2]) + \
                 ' --tsyganenko_W4=%g --tsyganenko_W5=%g --tsyganenko_W6=%g'%(W[3], W[4], W[5])

        base_damp_cmd =  './damping --inp_file "%s" --out_file "%s" '%(ray_outfile, damp_outfile) + \
                    ' --Kp %g --AE %g'%(Kp, AE) + \
                    ' --yearday %s --msec %d'%(yearday, milliseconds_day) +\
                    ' --geom_factor=%d'%include_geom_factor

        if mode == 1:
            # -------------- Test the Ngo model ------------------
            configfile = os.path.join(project_root,"newray_default.in")
            damp_mode = 0

            ray_cmd = base_cmd + ' --ngo_configfile="%s"'%(configfile)
            damp_cmd =  base_damp_cmd + ' --mode %d'%damp_mode

        if mode == 2:
            # -------------- Test the full GCPM model ------------------
            damp_mode = 1

            ray_cmd = base_cmd + ' --gcpm_kp=%g'%(Kp)
            damp_cmd =  base_damp_cmd + ' --mode %d'%damp_mode

        if mode == 3:
            # -------------- Test the uniformly-sampled GCPM model ------------------
            mode3_interpfile = os.path.join(project_root,'precomputed_grids',
                'gcpm_kp4_2001001_L10_80x80x80_noderiv.txt')
            damp_mode = 1

            ray_cmd = base_cmd +' --interp_interpfile="%s"'%(mode3_interpfile)
            damp_cmd =  base_damp_cmd + ' --mode %d'%damp_mode

        if mode == 4:
            # -------------- Test the randomly-sampled GCPM model ------------------
            mode4_modelfile = os.path.join(project_root,
                'precomputed_grids','precomputed_model_gcpm_2010001_0_kp2_L10_random.dat')

            # Mode4 interpolation parameters:
            scattered_interp_window_scale = 1.2
            scattered_interp_order = 2
            scattered_interp_exact = 0   # Try 0 if there's weirdness at discontinuities
            scattered_interp_local_window_scale = 5

            damp_mode = 1

            ray_cmd = base_cmd + ' --interp_interpfile=%s'%(mode4_modelfile) +\
                    ' --scattered_interp_window_scale=%d'%(scattered_interp_window_scale) +\
                    ' --scattered_interp_order=%d'%(scattered_interp_order) +\
                    ' --scattered_interp_exact=%d'%(scattered_interp_exact) +\
                    ' --scattered_interp_local_window_scale=%d'%(scattered_interp_local_window_scale)
            damp_cmd =  base_damp_cmd + ' --mode %d'%damp_mode

        if mode == 6:
            # -------------- Test the Simplified GCPM model ------------------
            MLT = 0
            fixed_MLT = 1 # Force the raytracer to stay in the meridonal plane?
            damp_mode = 0

            ray_cmd = base_cmd + ' --MLT="%g" --fixed_MLT=%g --kp=%g'%(MLT, fixed_MLT, Kp)
            damp_cmd =  base_damp_cmd + ' --mode %d'%damp_mode

        # Run it!

        #print("------- Running mode %d -------"%mode)
        #print("Command is:")
        #print(ray_cmd)
        #print();

        os.system(ray_cmd)

        #print("------- Running damping, mode %d -------"%damp_mode)

        #print(damp_cmd)
        os.system(damp_cmd)

    # Move back to the working directory
    os.chdir(cwd)

    print('ray done!')

    return
