import numpy as np
import datetime as dt
import os
import json
from multiprocessing import Pool, cpu_count
import multiprocessing
import pandas as pd
from constants import *
from convert_coords import convert2

# HERE ARE HELPER FUNCS TO RUN rayinput.py----------------------------------------------------------

# -------------------------------------- getting the bdir! ---------------------------------
def getBdir(ray_start, ray_datenum, rayfile_directory, thetas, phis, md, freq):
    positions = ray_start
    directions = [(0,0,0)]
    freqs = [freq]
    # going to run a ray real quick
    single_run_rays(ray_datenum, positions, directions, freqs, rayfile_directory, md)

    # Load all the rayfiles in the output directory
    ray_out_dir = rayfile_directory + '/'+dt.datetime.strftime(ray_datenum, '%Y-%m-%d_%H_%M_%S')
    file_titles = os.listdir(ray_out_dir)

    # create empty lists to fill with ray files and damp files
    raylist = []
    for filename in file_titles:
        if '.ray' in filename and str(md) in filename and str('Main') in filename:
            print(filename)
            raylist += read_rayfile(os.path.join(ray_out_dir, filename))

    # get b direction for this ray
    for r in raylist:
        B0 = [r['B0'].x[0], r['B0'].y[0], r['B0'].z[0]]
        # vec in T in SM car coordinates
        # create unit vector 
        Bmag = np.sqrt(r['B0'].x[0]**2 + r['B0'].y[0]**2 + r['B0'].z[0]**2)
        Bunit = [r['B0'].x[0]/Bmag, r['B0'].y[0]/Bmag, r['B0'].z[0]/Bmag]

    # now we have Bunit in SM car
    # let's put it in spherical (easier for changing the wavenormal)
    sph_dir = convert2([Bunit], ray_datenum, 'SM', 'car', ['Re','Re','Re'], 'SM', 'sph', ['Re','deg','deg'])     
    converted_dirs = []
    for theta, phi in zip(thetas,phis):
        new_dir = [sph_dir[0][0],sph_dir[0][1]+theta, sph_dir[0][2]+phi]
        converted_dir = convert2([new_dir], ray_datenum, 'SM', 'sph', ['Re','deg','deg'], 'SM', 'car', ['Re','Re','Re']) 
        
        converted_dirs.append(converted_dir[0])

    return converted_dirs

# ------------------- helper func for the datetimes ---------------------
def get_yearmiliday(ray_datenum):
    # convert for raytracer settings for convenience
    days_in_the_year = ray_datenum.timetuple().tm_yday
    days_in_the_year = format(days_in_the_year, '03d')

    # yearday and miliseconds day are used by raytracer
    yearday = str(ray_datenum.year)+ str(days_in_the_year)   # YYYYDDD
    milliseconds_day = ray_datenum.hour*3.6e6 + ray_datenum.minute*6.0e4 + ray_datenum.second*1.0e3

    return yearday, milliseconds_day


# ----------------------------- run -----------------------------------------------
def single_run_rays(ray_datenum, positions, directions, freqs, rayfile_directory, md):

    # this is just a silly way to get the names of the workers to keep everyone separate
    fp = str(multiprocessing.current_process())[19:36]
    fp = fp.replace("'",'')
    if len(fp) < 17:
        fp = fp[:-1] + '0' + fp[-1]
    # make numbers double digit

    yearday, milliseconds_day = get_yearmiliday(ray_datenum)

    project_root = os.getcwd()  # grabs current full path

    # Create directory for inputs/outputs if doesn't already exist
    ray_out_dir = rayfile_directory + '/'+ dt.datetime.strftime(ray_datenum, '%Y-%m-%d_%H_%M_%S')

    if not os.path.exists(ray_out_dir):
        os.makedirs(ray_out_dir)

    # Set input file path
    ray_inpfile = os.path.join(ray_out_dir, 'ray_inpfile'+fp+'.txt')

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
    ray_outfile = os.path.join(ray_out_dir, 'ray_out_mode%d_%s.ray' % (md,fp))

    # The base command -- with parameters common for all modes
    base_cmd = './raytracer --outputper=%d --dt0=%g --dtmax=%g' % (1, dt0, dtmax) + \
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

    return 

# -------------------------------- PARALLELIZE --------------------------------
def parallel_run_rays(ray_datenum, positions, directions, freqs, output_directory, md):

    # parallel
    nmbrcores = cpu_count()
    nrays = len(positions) # same as length of directions and freqs
    print(nmbrcores, ' cores running ', nrays, ' total')

    # lists the same length as the number of cores
    time_list = [ray_datenum for i in range(nmbrcores)]
    output_directory_list = [output_directory for i in range(nmbrcores)]
    mds = [md for i in range(nmbrcores)]

    # break up all the arguments 
    direction_list = [directions[int(i * (nrays/nmbrcores)):int((i+1)*nrays/nmbrcores)] for i in range(nmbrcores)]
    position_list = [positions[int(i * (nrays/nmbrcores)):int((i+1)*nrays/nmbrcores)] for i in range(nmbrcores)]
    freq_list = [freqs[int(i * (nrays/nmbrcores)):int((i+1)*nrays/nmbrcores)] for i in range(nmbrcores)]

    # send in all the arguments
    lstarg = zip(time_list, position_list, direction_list, freq_list, output_directory_list, mds)
    
    with Pool(nmbrcores) as p:
        results = p.starmap(single_run_rays, lstarg)

# ---------------------------- read output using pandas -------------------------
def read_rayfile(rayfile):
    ''' Load output from Forest's raytracer'''
    x = pd.read_csv(rayfile,delim_whitespace=True, header=None)
    # % fields: raynum, stopcond, time, pos(3), vprel(3), vgrel(3), n(3),
    # % B0(3), Nspec, qs(Nspec), ms(Nspec), Ns(Nspec), nus(Nspec)

    raynums = np.unique(x[0])
    numrays = len(raynums)

    out = []
    for ii in range(numrays):
        tmp = x[x[0]==raynums[ii]]
        tmp.reset_index(drop=True)
        #data = pd.DataFrame(columns=['raynum','time','pos','vprel','vgrel','n','B0','qs','ms','Ns','nus'])
        data = dict()
        data['time'] = tmp[2]
        data['pos'] = tmp.loc[:,3:5]
        data['pos'].columns=['x','y','z']
        data['vprel'] = tmp.loc[:,6:8]
        data['vprel'].columns=['x','y','z']
        data['vgrel'] = tmp.loc[:,9:11]
        data['vgrel'].columns=['x','y','z']
        data['n'] = tmp.loc[:,12:14]
        data['n'].columns=['x','y','z']
        data['B0'] = tmp.loc[:,15:17]
        data['B0'].columns=['x','y','z']
        
        # #data = pd.DataFrame(columns=['raynum','time','pos','vprel','vgrel','n','B0','qs','ms','Ns','nus'])
        # data = dict()
        # data['time'] = tmp[2]

        # data['pos']   = coord.Coords(zip(tmp.loc[:,3], tmp.loc[:,4],  tmp.loc[:,5]), 'SM','car')

        # data['pos']   = data['pos'].
        # data['vprel'] = tmp.loc[:,6:8]
        # data['vprel'].columns=['x','y','z']
        # data['vgrel'] = tmp.loc[:,9:11]
        # data['vgrel'].columns=['x','y','z']
        # data['n'] = tmp.loc[:,12:14]
        # data['n'].columns=['x','y','z']
        # data['B0'] = tmp.loc[:,15:17]
        # data['B0'].columns=['x','y','z']

        data['w'] = tmp.iloc[0,18]       # Frequency
        data['Nspec'] = tmp.iloc[0,19]   # Number of species
        data['stopcond'] = tmp.iloc[0,1] # Stop condition
        data['qs'] =  tmp.loc[:,20 + 0*data['Nspec']:20 + 1*data['Nspec'] - 1]
        data['ms'] =  tmp.loc[:,20 + 1*data['Nspec']:20 + 2*data['Nspec'] - 1]
        data['Ns'] =  tmp.loc[:,20 + 2*data['Nspec']:20 + 3*data['Nspec'] - 1]
        data['nus']=  tmp.loc[:,20 + 3*data['Nspec']:20 + 4*data['Nspec'] - 1]

        # read damping data too, if we have it:
        if (tmp.shape[1] > 20+4*data['Nspec']):
            print('we have it')
            data['damping'] = tmp.loc[:,20+4*data['Nspec']]

            
        out.append(data)
    return out
