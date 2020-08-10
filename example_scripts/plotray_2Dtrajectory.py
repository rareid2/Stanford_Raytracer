
"""
here is a script that will call run_rays and plot the trajectory with
normalized wave power as a color scale
ONLY WORKS FOR A SINGLE POSITION - CHANGE THE TIME FOR DIFF POSITIONS
plotting in MAG coordinates
"""

# import needed packages
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
import sys
import datetime as dt
from raytracer_utils import readdump, read_rayfile, read_rayfiles, read_damp
from run_rays import run_rays
from raytracer_settings import *
from IGRF_funcs import B_dir, trace_fieldline_ODE, findFootprints, B_direasy
from spacepy import coordinates as coord
from spacepy import irbempy
from spacepy.time import Ticktock
from TLE_funcs import TLE2pos
import tempfile
from run_model_dump import modeldump

# coordinate mania!
# TLES give us GEI, raytracer needs SM, and IGRF funcs needs GEO

# -------------------------------- SET TIME and OTHER SETTINGS --------------------------------
# change time information here - use UTC -

#dates = [dt.datetime(2020, 7, 25, 0, 6, 0)]
#fs = [8.2e3]
#bs = ['burst']


def readconjlog(fnameconj):
    dates = []

    infile = open(fnameconj,'r')
    chdatel = 0

    # goes thru line by line
    for line in infile:
        out = line
        # skip if empty
        if not out.strip():
            continue

        year = 2020
            
        month = int(out[5:7])
        day = int(out[8:10])
        hour = int(out[11:13])
        minute = int(out[14:16])
        second = int(out[17:19])
        
        # convert to fractions of a day
        minphr = 60
        
        chmin = minute/minphr
        chdate = hour + chmin

        if np.abs(chdate - chdatel) > 30/minphr: # if more than 30 minutes apart, this is a new conjunction!
            dates.append(dt.datetime(year, month, day, hour, minute))
        chdatel = chdate

    infile.close()
    return dates
weekstart = 17
stdates = [dt.datetime(2020, 8, d, 0, 0) for d in range(weekstart, weekstart+8)]

# conj list? 
clist = 1

condtime = []
fs = []
bs = []

if clist == 1: 
    for cdate in stdates:
        
        year = cdate.year
        month = cdate.month
        day = cdate.day
        hours = cdate.hour
        minutes = cdate.minute
        seconds = cdate.second
        ray_datenum = cdate

        datadir = '/home/rileyannereid/workspace/SR-output/' + 'fullday' + '/'
        datadir = datadir + str(ray_datenum.month) + str(ray_datenum.day) + str(ray_datenum.year) + '/'
        fname = datadir + str(ray_datenum.month) + str(ray_datenum.day) + str(ray_datenum.year) + 'conjlist.txt'
        newdates = readconjlog(fname)

        condtime.extend(newdates)
        for ii in condtime:
            fs.append(28e3)
            bs.append('fullday')
        for ii in condtime:
            fs.append(8.2e3)
            bs.append('fullday')
        condtime.extend(newdates)



    dates = condtime


else:
    weekstart = 17
    dates = [dt.datetime(2020, 8, d, 0, 0) for d in range(weekstart, weekstart+8)]
    fs = [28e3 for d in range(weekstart, weekstart+8)] # arbitrary
    bs = ['fullday' for d in range(weekstart, weekstart+8)]



for cdate, cf, bsstr in zip(dates, fs, bs):
    year = cdate.year
    month = cdate.month
    day = cdate.day
    hours = cdate.hour
    minutes = cdate.minute
    seconds = cdate.second

    freq = [cf]

    ray_datenum = dt.datetime(year, month, day, hours, minutes, seconds)
    modeldump(year, month, day, hours, minutes, seconds) # run model dump to update plasmasphere
    thetalist = [45] # in deg -- what angles to launch at? 

    checkdir = 0
    crs_out = 'MAG'  # theres a bug with GEO coords -- maybe its the fieldlines? 
    datadir = '/home/rileyannereid/workspace/SR-output/' + bsstr + '/'
    rename = str(ray_datenum.month) + str(ray_datenum.day) + str(ray_datenum.year)

    datadir = datadir + rename + '/' + str(freq[0]/1e3) + 'kHz' + rename + str(ray_datenum.hour) + str(ray_datenum.minute) + '/'
    try:
        os.mkdir(datadir)
    except OSError:
        print ("Creation of the directory %s failed" % datadir)
    else:
        print ("Successfully created the directory %s" % datadir)
    # -------------------------------- GET POSITIONS --------------------------------
    # get DSX and VPM positions for... 
    r, tvec = TLE2pos(1, ray_datenum, 1)

    # redefine time here -- more accurate
    ray_datenum = tvec[0]

    # convert to meters
    dsx = [rpos*1e3 for rpos in r[0]]
    vpm = [rpos*1e3 for rpos in r[1]]

    # only grab first one - weird bug fix with JD dates
    dsx = [dsx[0]]
    vpm = [vpm[0]]

    # convert startpoint to SM car for raytracer
    GEIcar_dsx = coord.Coords(dsx, 'GEI', 'car', units=['m', 'm', 'm'])
    GEIcar_dsx.ticks = Ticktock(ray_datenum, 'UTC') # add ticks

    SMcar_dsx = GEIcar_dsx.convert('SM', 'car') # needed for raytracer
    GEOcar_dsx = GEIcar_dsx.convert('GEO', 'car') # needed for Bfield calcs

    # convert vpm -- to check which hemi and plot later
    GEIcar_vpm = coord.Coords(vpm, 'GEI', 'car', units=['m', 'm', 'm'])
    GEIcar_vpm.ticks = Ticktock(ray_datenum, 'UTC') # add ticks
    outsph_vpm = GEIcar_vpm.convert(crs_out, 'sph')

    # -------------------------------- DEFINE RAY DIRECTIONS --------------------------------
    # start position of raytracer
    position = [float(SMcar_dsx.x), float(SMcar_dsx.y), float(SMcar_dsx.z)]

    # check with hemi we are in
    if outsph_vpm.lati > 0:
        dir = 1   # north
        dirstr = 'north'
    else:
        dir = -1  # south
        dirstr = 'south'

    Bstart = [float(GEOcar_dsx.x)/R_E, float(GEOcar_dsx.y)/R_E, float(GEOcar_dsx.z)/R_E]
    Bx, By, Bz = B_direasy(ray_datenum, Bstart, dir)

    # convert direction to SM coordinates for raytracer
    dirB = np.reshape(np.array([Bx, By, Bz]), (1, 3))
    dirB = coord.Coords(dirB[0], 'GEO', 'car', units=['Re', 'Re', 'Re'])
    dirB.ticks = Ticktock(ray_datenum, 'UTC') # add ticks
    SMsph_dirB = dirB.convert('SM', 'sph')

    # fill for raytracer call
    positions = []
    directions = []

    # rotate directions
    for theta in thetalist:
        # increase (or decrease) polar angle
        newth = float(SMsph_dirB.lati) + theta
        Rot_dirB = [float(SMsph_dirB.radi), newth, float(SMsph_dirB.long)] 
        Rot_dirB = coord.Coords(Rot_dirB, 'SM', 'sph')
        Rot_dirB.ticks = Ticktock(ray_datenum, 'UTC') # add ticks
        SMcar_dirB = Rot_dirB.convert('SM', 'car')

        if checkdir == 1:
            outcar_dirB = dirB.convert(crs_out, 'car') # before theta addition
            outcar_dirBrot = Rot_dirB.convert(crs_out, 'car')
            fig, ax = plt.subplots(1,1, sharex=True, sharey=True)
            ax.quiver(0, 0, outcar_dirB.x, outcar_dirB.z, label='th', color = 'r')
            ax.quiver(0, 0, outcar_dirBrot.x, outcar_dirBrot.z, label='B0', color = 'b')
            plt.legend()
            plt.title(str(theta) + ' deg from B0 in ' + crs_out + ' coords')
            plt.show()
            plt.close()

        direction = [float(SMcar_dirB.x), float(SMcar_dirB.y), float(SMcar_dirB.z)]
        
        # add the normalized direction (or zeros)
        directions.append(np.squeeze(direction))

        # make sure position list matches direction list
        positions.append(position)

    # -------------------------------- RUN RAYS --------------------------------
    # convert for raytracer settings
    days_in_the_year = ray_datenum.timetuple().tm_yday
    days_in_the_year = format(days_in_the_year, '03d')

    # yearday and miliseconds day are used by raytracer
    yearday = str(year)+ str(days_in_the_year)   # YYYYDDD
    milliseconds_day = hours*3.6e6 + minutes*6e4 + seconds*1e3

    # run it!
    tmpdir = tempfile.mkdtemp() 
    run_rays(freq, positions, directions, yearday, milliseconds_day, tmpdir)

    # -------------------------------- LOAD OUTPUT --------------------------------
    # Load all the rayfiles in the output directory
    ray_out_dir = tmpdir
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

    # -------------------------------- CONVERT COORDINATES --------------------------------
    # convert to desired coordinate system into vector list rays
    rays = []
    for r in raylist:
        tmp_coords = coord.Coords(list(zip(r['pos'].x, r['pos'].y, r['pos'].z)), 'SM', 'car', units=['m', 'm', 'm'])
        tvec_datetime = [ray_datenum + dt.timedelta(seconds=s) for s in r['time']]
        tmp_coords.ticks = Ticktock(tvec_datetime, 'UTC')  # add ticks
        tmp_coords.sim_time = r['time']
        new_coords = tmp_coords.convert(crs_out, 'sph') # needs to be sph for rotation
        rays.append(new_coords)

    dlist = []
    for d in damplist:
        damp = d["damping"]
        damp = np.squeeze(np.array(damp))
        dlist.append(damp)
        #print(damp)
        
    # -------------------------------- PLOTTING --------------------------------
    fig, ax = plt.subplots(1,1, sharex=True, sharey=True)
    lw = 2  # linewidth

    # rotate plot to be in plane of view
    outsph_dsx = GEIcar_dsx.convert(crs_out, 'sph') # add ticks
    th = -outsph_dsx.long

    # rotate long to be at prime merid
    Rot_dsx = coord.Coords([float(outsph_dsx.radi), float(outsph_dsx.lati), float(outsph_dsx.long + th)], crs_out, 'sph', units=['m', 'deg', 'deg'])
    Rot_vpm = coord.Coords([float(outsph_vpm.radi), float(outsph_vpm.lati), float(outsph_vpm.long + th)], crs_out, 'sph', units=['m', 'deg', 'deg'])
    Rot_dsx.ticks = Ticktock(ray_datenum, 'UTC')
    Rot_vpm.ticks = Ticktock(ray_datenum, 'UTC')
    outcar_dsx = Rot_dsx.convert(crs_out, 'car')
    outcar_vpm = Rot_vpm.convert(crs_out, 'car')

    # plot sat locations
    plt.plot(outcar_dsx.x / R_E, outcar_dsx.z / R_E, '-go', zorder=105, label='DSX')
    plt.plot(outcar_vpm.x / R_E, outcar_vpm.z / R_E, '-yo', zorder=104, label='VPM')

    # rotate rays to be at prime merid also and plot
    #for r, d in zip(rays, dlist):
    for r in rays:
        rrad = []
        rlat = []
        rlon = []
        rrad.append(r.radi)
        rlat.append(r.lati)
        rlon.append(r.long)

        rrlon = [rl + th for rl in rlon]
        rcoords = [np.column_stack([rr, rl, rrl]) for rr, rl, rrl in zip(rrad, rlat, rrlon)]

        Rot_ray = coord.Coords(rcoords[0], crs_out, 'sph', units=['m', 'deg', 'deg'])
        Rot_ray.ticks = Ticktock(tvec_datetime, 'UTC')
        outcar_ray = Rot_ray.convert(crs_out, 'car')

        if len(outcar_ray.x) > 1:
            # plotp = ax.scatter(MAGcar_ray.x / R_E, MAGcar_ray.z / R_E, c=d, s = 1, cmap = 'Reds', vmin = 0, vmax = 1.5, zorder = 103)
            plotp = ax.scatter(outcar_ray.x / R_E, outcar_ray.z / R_E, c = 'Red', s = 1, zorder = 103)

    # add in color bar - will be just for the last ray, but bounds are set
    # plt.colorbar(plotp, label = 'Normalized wave power')

    # -------------------------------- EARTH AND IONO --------------------------------
    earth = plt.Circle((0, 0), 1, color='b', alpha=0.5, zorder=100)
    iono = plt.Circle((0, 0), (R_E + H_IONO) / R_E, color='g', alpha=0.5, zorder=99)
    ax.add_artist(earth)
    ax.add_artist(iono)

    # -------------------------------- PLASMASPHERE --------------------------------
    # bad pls change this @ riley!
    # NOT IN CORRECT COORDS
    path2plasma = '/home/rileyannereid/workspace/Stanford_Raytracer/example_scripts/modeldumps/'
    plasma_model_dump = os.path.join(path2plasma, 'model_dump_mode_1_XZ.dat')
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

    # ---------------------------------- BFIELD -----------------------------------
    # (from IGRF13 model)
    L_shells = [2, 3, 4, -2, -3, -4]  # Field lines to draw

    for L in L_shells:
        Blines = []
        if L < 0:
            rot = th
        else:
            rot = -th

        Lstart = [L, 0, rot]
        Lcoords = coord.Coords(Lstart, crs_out, 'sph', units=['Re', 'deg', 'deg'])
        Lcoords.ticks = Ticktock(ray_datenum, 'UTC')
        GEO_Lcoords = Lcoords.convert('GEO', 'car')
        Blines.append(trace_fieldline_ODE([float(GEO_Lcoords.x),float(GEO_Lcoords.y),float(GEO_Lcoords.z)], 0, '0', 1, ray_datenum))
        Blines.append(trace_fieldline_ODE([float(GEO_Lcoords.x),float(GEO_Lcoords.y),float(GEO_Lcoords.z)], 0, '0', -1, ray_datenum))

        for blinex, bliney, blinez in Blines:
            raytime = []
            # create list of the same times
            for m in range(int(len(blinex))):
                raytime.append(ray_datenum)
            
            # convert coords
            bpos = np.column_stack((blinex, bliney, blinez))
            bpos = coord.Coords(bpos, 'GEO', 'car', units=['Re', 'Re', 'Re'])
            bpos.ticks = Ticktock(raytime, 'UTC') # add ticks
            outsph_bline = bpos.convert(crs_out, 'sph')

            btemp_rad = []
            btemp_lat = []
            btemp_lon = []
            btemp_rad.append(outsph_bline.radi)
            btemp_lat.append(outsph_bline.lati)
            btemp_lon.append(outsph_bline.long)

            brlon = [bl * 0 for bl in btemp_lon]
            bcoords = [np.column_stack([br, bl, brl]) for br, bl, brl in zip(btemp_rad, btemp_lat, brlon)]
            Rot_bline = coord.Coords(bcoords[0], crs_out, 'sph', units=['Re', 'deg', 'deg'])
            Rot_bline.ticks = Ticktock(raytime, 'UTC') # add ticks
            outcar_bline = Rot_bline.convert(crs_out, 'car')
            plt.plot(np.sign(rot) * outcar_bline.x, np.sign(rot) * outcar_bline.z, color='b', linewidth=1, linestyle='dashed')

    print('finished plotting field lines')

    # plot field line from orbital position
    Blines = []

    Blines.append(trace_fieldline_ODE(Bstart, 0, '0', 1, ray_datenum))
    Blines.append(trace_fieldline_ODE(Bstart, 0, '0', -1, ray_datenum))

    for blinex, bliney, blinez in Blines:
        raytime = []

        # create list of the same times
        for m in range(int(len(blinex))):
            raytime.append(ray_datenum)
        
        # convert
        bpos = np.column_stack((blinex, bliney, blinez))
        bpos = coord.Coords(bpos, 'GEO', 'car', units=['Re', 'Re', 'Re'])
        bpos.ticks = Ticktock(raytime, 'UTC') # add ticks
        outsph_bline = bpos.convert(crs_out, 'sph')

        btemp_rad = []
        btemp_lat = []
        btemp_lon = []
        btemp_rad.append(outsph_bline.radi)
        btemp_lat.append(outsph_bline.lati)
        btemp_lon.append(outsph_bline.long)

        brlon = [bl * 0 for bl in btemp_lon]
        bcoords = [np.column_stack([br, bl, brl]) for br, bl, brl in zip(btemp_rad, btemp_lat, brlon)]
        Rot_bline = coord.Coords(bcoords[0], crs_out, 'sph', units=['Re', 'deg', 'deg'])
        Rot_bline.ticks = Ticktock(raytime, 'UTC') # add ticks
        outcar_bline = Rot_bline.convert(crs_out, 'car')
        plt.plot(outcar_bline.x, outcar_bline.z, color='r', linewidth=1, linestyle='dashed')

    # -------------------------------- GET FOOTPRINT --------------------------------
    # also in GEO car, so need to use bstart 
    GDZsph_foot = findFootprints(ray_datenum, Bstart, dirstr)
    GDZsph_foot.units = ['km', 'deg', 'deg']
    GDZsph_foot.ticks = Ticktock(ray_datenum, 'UTC')
    outsph_foot = GDZsph_foot.convert(crs_out, 'sph')

    # rotate and plot
    Rot_foot = coord.Coords([float(outsph_foot.radi), float(outsph_foot.lati), float(outsph_foot.long + th)], crs_out, 'sph', units=['Re', 'deg', 'deg'])
    Rot_foot.ticks = Ticktock(ray_datenum, 'UTC')
    outcar_foot = Rot_foot.convert(crs_out, 'car')
    plt.plot(outcar_foot.x, outcar_foot.z, '-ro', label='Bfield footpoint')

    # -------------------------------- FORMATTING --------------------------------
    ax.set_aspect('equal')
    max_lim = 4

    plt.xticks(np.arange(-max_lim, max_lim, step=1))
    plt.yticks(np.arange(-max_lim, max_lim, step=1))
    plt.xlabel('L (R$_E$)')
    plt.ylabel('L (R$_E$)')
    plt.xlim([-max_lim, max_lim])
    plt.ylim([-2.5, 2.5])

    mytitle = str(freq[0]/1e3) + 'kHz rays at ' + str(ray_datenum) + ' in ' + crs_out + ' coords'
    plt.title(mytitle)
    ax.legend(loc = 'lower center', fontsize =  'x-small')

    savename = datadir + str(freq[0]/1e3) + 'kHz' + str(ray_datenum.year) + str(ray_datenum.month) + str(ray_datenum.day) + str(ray_datenum.hour) + str(ray_datenum.minute) + '2Dview.png'
    #plt.savefig(savename, format='svg')
    plt.savefig(savename, format='png')
    #plt.close()
    #plt.show()

# final plot
"""
if checkdir==1:
    fig, ax = plt.subplots(1,1, sharex=True, sharey=True)

    outcarray = [pos for pos in outcar_ray]
    loopstop = len(outcarray) - 2

    for posind, pos in enumerate(outcarray):
        if posind > loopstop:
            break
        plt.quiver(pos.x/R_E, pos.z/R_E, (outcarray[posind+1].x - pos.x) / R_E, (outcarray[posind+1].z - pos.z) / R_E)
    
    plt.plot(outcar_bline.x, outcar_bline.z, color='r', linewidth=1, linestyle='dashed', label='bfieldline')
    earth = plt.Circle((0, 0), 1, color='b')
    ax.add_artist(earth)
    plt.title('ray launched at ' + str(ray_datenum) + ' at ' + str(thetalist[0]) + ' deg from B0 in ' + crs_out + ' coords')
    plt.show()
"""
# ------------------------------------------- END --------------------------------------------