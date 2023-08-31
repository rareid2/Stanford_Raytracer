import numpy as np
import datetime as dt
import os
import json
from multiprocessing import Pool, cpu_count
from constants import *
from convert_coords import convert2
# ------------------- helper func for the datetimes ---------------------
def get_yearmiliday(ray_datenum):
    # convert for raytracer settings for convenience
    days_in_the_year = ray_datenum.timetuple().tm_yday
    days_in_the_year = format(days_in_the_year, '03d')

    # yearday and miliseconds day are used by raytracer
    yearday = str(ray_datenum.year)+ str(days_in_the_year)   # YYYYDDD
    milliseconds_day = ray_datenum.hour*3.6e6 + ray_datenum.minute*6.0e4 + ray_datenum.second*1.0e3

    return yearday, milliseconds_day
# -------------------------------------- getting the bdir! ---------------------------------
def getBdir_src(ray_datenum, positions, thetas, phis, rayfile_directory):

    yearday, milliseconds_day = get_yearmiliday(ray_datenum)

    # get over to where the binary is
    cwd = os.getcwd()
    os.chdir(raytracer_dir)

    model_outfile = os.path.join(rayfile_directory, 'bmodel_dump.dat')
    model_inpfile = os.path.join(rayfile_directory, 'bmodel_input.in')

    # Write input positions file
    f = open(model_inpfile, 'w')

    # Go through list of positions, write a new line for each
    for pos0 in positions:
        f.write('%1.15e %1.15e %1.15e\n' % (
            pos0[0], pos0[1], pos0[2]))
    f.close()

    # get number of rays for memory allocation
    nrays = len(positions)

    cmd = './getBfield ' + \
        '--filename="%s" ' % (model_outfile) + \
        '--infilename="%s" ' % (model_inpfile) + \
        ' --yearday=%s --milliseconds_day=%d ' % (yearday, milliseconds_day) + \
        '--use_igrf=%g --use_tsyganenko=%g ' % (use_IGRF, use_tsyg) + \
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
        '--nrays=%g' % nrays

    # execute the command
    os.system(cmd)

    # save the outputs - converted directions
    converted_dirs = []

    # intermediate output
    Bunits = []

    # open the file and save stuff
    f = open(model_outfile,'r')
    three_count = 0
    # read output file with all the directions
    for line in f.readlines():
        # comes out line by line, need to read in directions one at a time
        if three_count==0:
            B0 = []
        if three_count<2:
            B0.append(float(line))
            three_count+=1
        else:
            B0.append(float(line))
            three_count=0

            # reset, that means we have the direc
            Bmag = np.sqrt(B0[0]**2 + B0[1]**2 + B0[2]**2)
            Bunit = [B0[0]/Bmag, B0[1]/Bmag, B0[2]/Bmag]
            Bunits.append(Bunit)

    f.close()

    # now change to new directions with thetas and phis
    # theta is alpha
    for Bunit, theta, phi in zip(Bunits,thetas,phis):

        if Bunit[0] == 0 or Bunit[2] == 0:
            r1 = [1,(-1*Bunit[0]-Bunit[2])/Bunit[1],1]
        else:
            r1 = [1,1,(-1*Bunit[1]-Bunit[0])/Bunit[2]]

        r1 = np.array(r1)/np.linalg.norm(np.array(r1))
        r2  = np.cross(r1,Bunit)
        T_rotate = np.column_stack((r1,r2,Bunit))

        r = 1/(np.cos(theta*D2R))
        cone_vec = np.array([r*np.sin(theta*D2R)*np.cos(phi*D2R),r*np.sin(theta*D2R)*np.sin(phi*D2R),r*np.cos(theta*D2R)])
        cone_vec = np.matmul(T_rotate,np.transpose(cone_vec))

        cone_vec = cone_vec/np.linalg.norm(cone_vec)
        if theta > 90 and theta < 270:
            zsign = -1
        elif theta < 90 and theta > -90:
            zsign = 1
        else:
            print('ERROR INPUT THETA IS NOT BETWEEN -90 and 90 or 90 and 270 deg')
        converted_dirs.append(zsign*cone_vec)

    # for plotting

    #ax = plt.axes(projection='3d')
    #ax.quiver(0,0,0,converted_dirs[0][0],converted_dirs[0][1],converted_dirs[0][2],length=1)
    #ax.quiver(0,0,0,Bunit[0],Bunit[1],Bunit[2],length=1,color='g')
    #ax.axes.set_xlim3d(left=0, right=1.5)
    #ax.axes.set_ylim3d(bottom=0, top=1.5)
    #ax.axes.set_zlim3d(bottom=0, top=1.5)
    #plt.show()
    #plt.close()
    os.chdir(cwd)
    return converted_dirs
# ----------------------------- run -----------------------------------------------
def single_run_rays(ray_datenum, positions, directions, freqs, rayfile_directory, md, wid):

    # this is just a silly way to get the names of the workers to keep everyone separate
    fp = 'worker_'+str(wid)
    yearday, milliseconds_day = get_yearmiliday(ray_datenum)

    # Set input file path
    ray_inpfile = os.path.join(rayfile_directory, 'ray_inpfile_'+fp+'.txt')

    # Write the ray input file
    f = open(ray_inpfile, 'w')

    # Go through list of positions, write a new ray for every direction and freq at each pos
    for pos0, dir0, fr in zip(positions, directions, freqs):
        w0 = fr * 2.0 * np.pi
        f.write('%1.15e %1.15e %1.15e %1.15e %1.15e %1.15e %1.15e\n' % (
            pos0[0], pos0[1], pos0[2], dir0[0], dir0[1], dir0[2], w0))
    f.close()

    # GCPM model and damping code needs to be run in the same directory
    # as the binary file (and all the misc data files)
    cwd = os.getcwd()
    os.chdir(raytracer_dir)

    # Set output file path
    ray_outfile = os.path.join(rayfile_directory, 'ray_out_mode%d_%s.ray' % (md,fp))

    # The base command -- with parameters common for all modes
    base_cmd = './raytracer --outputper=%d --dt0=%g --dtmax=%g' % (outputper, dt0, dtmax) + \
                ' --tmax=%g --root=%d --fixedstep=%d --maxerr=%g' % (t_max, root, fixedstep, maxerr) + \
                ' --modelnum=%d --maxsteps=%d --minalt=%d --inputraysfile="%s"' % (
                    md, maxsteps, minalt, ray_inpfile) + \
                ' --outputfile="%s" --yearday=%s --milliseconds_day=%d' % (ray_outfile, yearday, milliseconds_day) + \
                ' --use_tsyganenko=%d --use_igrf=%d --tsyganenko_Pdyn=%g' % (use_tsyg, use_IGRF, Pdyn) + \
                ' --tsyganenko_Dst=%g --tsyganenko_ByIMF=%g --tsyganenko_BzIMF=%g' % (Dst, ByIMF, BzIMF) + \
                ' --tsyganenko_W1=%g --tsyganenko_W2=%g --tsyganenko_W3=%g' % (W[0], W[1], W[2]) + \
                ' --tsyganenko_W4=%g --tsyganenko_W5=%g --tsyganenko_W6=%g' % (W[3], W[4], W[5])

    if md == 6:
        # Test the Simplified GCPM model
        MLT = 0
        fixed_MLT = 0 # Force the raytracer to stay in the meridonal plane?

        ray_cmd = base_cmd + ' --MLT="%g" --fixed_MLT=%g --kp=%g' % (MLT, fixed_MLT, Kp)

    if md == 7:
        ray_cmd = base_cmd + ' --gcpm_kp=%g' % (Kp)

    os.system(ray_cmd)

    # Move back to the working directory
    os.chdir(cwd)

    print('raytracer done')

    return fp
# -------------------------------- PARALLELIZE --------------------------------
def GetArgs(nmbrcores, nrays, ray_datenum, positions, directions, freqs, output_directory, md):

    # lists the same length as the number of cores
    time_list = [ray_datenum for i in range(nmbrcores)]
    output_directory_list = [output_directory for i in range(nmbrcores)]
    mds = [md for i in range(nmbrcores)]

    # break up all the arguments
    direction_list = [directions[int(i * (nrays/nmbrcores)):int((i+1)*nrays/nmbrcores)] for i in range(nmbrcores)]
    position_list = [positions[int(i * (nrays/nmbrcores)):int((i+1)*nrays/nmbrcores)] for i in range(nmbrcores)]
    freq_list = [freqs[int(i * (nrays/nmbrcores)):int((i+1)*nrays/nmbrcores)] for i in range(nmbrcores)]

    # unique name
    worker_id = [i for i in range(1, nmbrcores+1)]

    # send in all the arguments
    lstarg = zip(time_list, position_list, direction_list, freq_list, output_directory_list, mds, worker_id)

    return lstarg
# -------------------------------- MAIN --------------------------------
if __name__ == '__main__':

    #  -------------------------- import settings from config file -----------------------------
    rayfile_directory = settings['rayfile_directory']
    mode = settings['mode']

    if mode == 'GCPM':
        md = 6
    if mode == 'DE':
        md = 7

    # -------------------- process datetime info -----------------------------------
    ray_datenum = dt.datetime.strptime(settings['date_time_iso_utc'],"%Y-%m-%dT%H:%M:%S.%fZ")
    # can't be timezone naive
    ray_datenum.replace(tzinfo=dt.timezone.utc)

    # only need to run this once! not in parallel!
    # Create directory for inputs/outputs if doesn't already exist
    # OLD: ray_out_dir = os.path.join(rayfile_directory,dt.datetime.strftime(ray_datenum, '%Y-%m-%d_%H_%M_%S'))
    ray_out_dir = os.path.join(rayfile_directory, runName)
    if not os.path.exists(ray_out_dir):
        os.makedirs(ray_out_dir)

    # --------------- get each ray position, freq, direcion --------------------------
    positions = []
    freqs = []
    directions = []
    thetas = []
    phis = []

    for ri,r in enumerate(settings['rays']):
        alt = r['alt']
        lat = r['lat']
        lon = r['lon']
        ray_start = convert2([[alt+R_E,lat,lon]], [ray_datenum],'GEO','sph',['m','deg','deg'], 'SM', 'car', ['m','m','m'])
        positions.append(ray_start[0])

        freq = r['freq']
        freqs.append(freq)

        alpha = r['alpha']
        phi = r['phi']

        thetas.append(alpha)
        phis.append(phi)

    # get directions as a batch
    directions = getBdir_src(ray_datenum, positions, thetas, phis, rayfile_directory)

    print('initialized all rays')

    nmbrcores = os.cpu_count() - 1
    #nmbrcores = 7
    nrays = len(positions) # same as length of directions and freqs
    print(nmbrcores, ' cores running ', nrays, ' total')

    # Zips everything for input into parallel processing
    lstarg = GetArgs(nmbrcores, nrays, ray_datenum, positions, directions,
                     freqs, ray_out_dir, md)

    with Pool(nmbrcores) as p:
        results = p.starmap(single_run_rays, lstarg)
