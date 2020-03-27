import numpy as np
import subprocess
import os
import datetime

# Dump plasmasphere models
print("Dumping plasmasphere models")

project_root = os.getcwd();
ray_out_dir = os.path.join(project_root, "test_outputs");
configfile = os.path.join(project_root, "newray_default.in")

R_E = 6371 # km

yearday = '2010001'  # YYYYDDD
milliseconds_day = 0;  # milliseconds into the day

Kp = 2
AE = 1.6
Pdyn = 4
Dst = 1.0
ByIMF = 0.0
BzIMF = -5
# Tsyganenko correction parameters
W = [0.132, 0.303, 0.083, 0.070, 0.211, 0.308]  # Doesn't matter if we're not using Tsyg

modes_to_do = [1]

mag_dump = False  # True for mag dipole coords, false for SM
use_IGRF = 1  # Magnetic field model (1 for IGRF, 0 for dipole)
use_tsyg = 1  # Use the Tsyganenko magnetic field model corrections

mode3_modelfile = os.path.join(project_root,
                               'precomputed_grids', 'gcpm_kp4_2001001_L10_80x80x80_noderiv.txt')
mode4_modelfile = os.path.join(project_root,
                               'precomputed_grids', 'precomputed_model_gcpm_2010001_0_kp2_L10_random.dat')

# Mode4 interpolation parameters:
scattered_interp_window_scale = 1.5
scattered_interp_order = 2
scattered_interp_exact = 1  # Try 0 if there's weirdness at discontinuities
scattered_interp_local_window_scale = 5

# GCPM model and damping code needs to be run in the same directory
# as the binary file (and all the misc data files)
cwd = os.getcwd();
os.chdir('../bin')

for mode in modes_to_do:
    for plane in ['XZ', 'XY', 'YZ']:

        print("doing model %d, %s plane" % (mode, plane))
        print(str(Kp))
        maxD = 10.0 * R_E * 1e3
        if plane == 'XZ':
            minx = -maxD
            maxx = maxD
            miny = 0
            maxy = 0
            minz = -maxD
            maxz = maxD
            nx = 200
            ny = 1
            nz = 200
        if plane == 'XY':
            minx = -maxD
            maxx = maxD
            miny = -maxD
            maxy = maxD
            minz = 0
            maxz = 0
            nx = 200
            ny = 200
            nz = 1
        if plane == 'YZ':
            minx = 0
            maxx = 0
            miny = -maxD
            maxy = maxD
            minz = -maxD
            maxz = maxD
            nx = 1
            ny = 200
            nz = 200

        # model_outfile='model_dump_mode%d_%d_%s.dat'%(modelnum, use_IGRF, plane)
        if mag_dump:
            model_outfile = os.path.join(ray_out_dir, 'model_dump_mode_%d_%s_MAG.dat' % (mode, plane))
        else:
            model_outfile = os.path.join(ray_out_dir, 'model_dump_mode_%d_%s.dat' % (mode, plane))

        cmd = './dumpmodel' + \
              ' --modelnum=%d --yearday=%s --milliseconds_day=%d ' % (mode, yearday, milliseconds_day) + \
              '--minx=%g --maxx=%g ' % (minx, maxx) + \
              '--miny=%g --maxy=%g ' % (miny, maxy) + \
              '--minz=%g --maxz=%g ' % (minz, maxz) + \
              '--nx=%g --ny=%g --nz=%g ' % (nx, ny, nz) + \
              '--filename="%s" ' % (model_outfile) + \
              '--use_igrf=%g --use_tsyganenko=%g ' % (use_IGRF, 0) + \
              '--tsyganenko_Pdyn=%g ' % (Pdyn) + \
              '--tsyganenko_Dst=%g ' % (Dst) + \
              '--tsyganenko_ByIMF=%g ' % (ByIMF) + \
              '--tsyganenko_BzIMF=%g ' % (BzIMF) + \
              '--tsyganenko_W1=%g ' % (W[0]) + \
              '--tsyganenko_W2=%g ' % (W[1]) + \
              '--tsyganenko_W3=%g ' % (W[2]) + \
              '--tsyganenko_W4=%g ' % (W[3]) + \
              '--tsyganenko_W5=%g ' % (W[4]) + \
              '--tsyganenko_W6=%g ' % (W[5]) + \
              '--gcpm_kp=%g ' % (Kp) + \
              '--ngo_configfile="%s" ' % configfile
        # '--ngo_configfile=%s '%(os.path.join(working_path,'newray.in')) +\

        if mode == 3:
            cmd += ' --interp_interpfile="%s"' % mode3_modelfile

        if mode == 4:
            cmd += '--interp_interpfile=%s ' % (mode4_modelfile) + \
                   '--scattered_interp_window_scale=%d ' % (scattered_interp_window_scale) + \
                   '--scattered_interp_order=%d ' % (scattered_interp_order) + \
                   '--scattered_interp_exact=%d ' % (scattered_interp_exact) + \
                   '--scattered_interp_local_window_scale=%d ' % (scattered_interp_local_window_scale)

        if mode == 6:
            cmd += ' --kp=%g ' % (Kp)
            if mag_dump:
                cmd += ' --mag_coords=1 '
        print(cmd)

        os.system(cmd)
        # os.system('mv %s %s'%(model_outfile, os.path.join(ray_out_dir, model_outfile)))

# Move back to the working directory
os.chdir(cwd)