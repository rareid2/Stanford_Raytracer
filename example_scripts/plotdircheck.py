# import needed packages
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
import math
import sys
import aacgmv2
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
from plotrefractivesurface import getLshell, stix_parameters

R2D = 180./np.pi
D2R = np.pi/180.
Hz2Rad = 2.*np.pi
Rad2Hz = 1./Hz2Rad
eo   = 8.854e-12   # C^2/Nm^2 
c    = 2.998e8     # m/s
Q_EL = 1.602e-19   # C
M_EL = 9.1e-31     # kg
M_P = 1.67e-27     # kg
R_E = 6371e3  # m

sys.path.insert(1, '/home/rileyannereid/workspace/scratches/')
from simplegifs import simplegifs

#-------------------------------- GET SAT POS --------------------------------
def getDSXVPMpos(ray_datenum, duration):
    # get DSX and VPM positions for... 
    r, tvec = TLE2pos(duration, ray_datenum, 1)

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
    MAGsph_vpm = GEIcar_vpm.convert('MAG', 'sph')

    # check with hemi we are in
    if MAGsph_vpm.lati > 0:
        dir = 1   # north
        dirstr = 'north'
    else:
        dir = -1  # south
        dirstr = 'south'

    return SMcar_dsx, GEOcar_dsx, MAGsph_vpm, dir, dirstr
#---------------------------------------------------------------------------


#-----------------------------------DEF DIRECS ------------------------------
# input is a GEO cartesian position in m
# theta is deg from fieldline
# dir is direction along fieldline

def getDIR(dir, ray_datenum, GEOposition, theta):
    Bstart =  [float(GEOposition.x)/R_E, float(GEOposition.y)/R_E, float(GEOposition.z)/R_E]

    Bx, By, Bz = B_direasy(ray_datenum, Bstart, dir)

    # convert direction to SM coordinates for raytracer
    dirB = np.reshape(np.array([Bx, By, Bz]), (1, 3))
    dirB = coord.Coords(dirB[0], 'GEO', 'car', units=['Re', 'Re', 'Re'])
    dirB.ticks = Ticktock(ray_datenum, 'UTC') # add ticks
    SMsph_dirB = dirB.convert('SM', 'sph')

    # increase (or decrease) polar angle
    newth = float(SMsph_dirB.lati) + theta
    Rot_dirB = [float(SMsph_dirB.radi), newth, float(SMsph_dirB.long)] 
    Rot_dirB = coord.Coords(Rot_dirB, 'SM', 'sph')
    Rot_dirB.ticks = Ticktock(ray_datenum, 'UTC') # add ticks
    SMcar_dirB = Rot_dirB.convert('SM', 'car')

    direction = [float(SMcar_dirB.x), float(SMcar_dirB.y), float(SMcar_dirB.z)]
    
    return direction
#---------------------------------------------------------------------------


#---------------------------------------------------------------------------
def rotateplane(th, MAGsphpos, ray_datenum):
    # rotate long to be at prime merid
    Rotpos = coord.Coords([float(MAGsphpos.radi), float(MAGsphpos.lati), float(MAGsphpos.long + th)], 'MAG', 'sph', units=['m', 'deg', 'deg'])
    Rotpos.ticks = Ticktock(ray_datenum, 'UTC')
    MAGcarpos = Rotpos.convert('MAG', 'car')
    return MAGcarpos
#---------------------------------------------------------------------------


#---------------------------------------------------------------------------
# plot field line from orbital position
def getDSXfieldline(GEOposition, ray_datenum, th):
    Bstart =  [float(GEOposition.x)/R_E, float(GEOposition.y)/R_E, float(GEOposition.z)/R_E]
    bline_dsx = []
    Blines = []

    Blines.append(trace_fieldline_ODE(Bstart, 0, '0', 1, ray_datenum))
    Blines.append(trace_fieldline_ODE(Bstart, 0, '0', -1, ray_datenum))

    for blinex, bliney, blinez in Blines:
        raytimes = [ray_datenum for m in range(int(len(blinex)))]

        # convert
        bpos = np.column_stack((blinex, bliney, blinez))
        bpos = coord.Coords(bpos, 'GEO', 'car', units=['Re', 'Re', 'Re'])
        bpos.ticks = Ticktock(raytimes, 'UTC') # add ticks
        MAGsph_bline = bpos.convert('MAG', 'sph')
        mc_bline = []
        for msph_bline in MAGsph_bline:
            MAGcar_bline_dsx = rotateplane(th, msph_bline, ray_datenum)
            mc_bline.append(MAGcar_bline_dsx)
        bline_dsx.append(mc_bline)

    return bline_dsx
#---------------------------------------------------------------------------

#---------------------------------------------------------------------------
def plotray2Ddir(ray_datenum, fs, intcheck):

    year = ray_datenum.year
    month = ray_datenum.month
    day = ray_datenum.day
    hours = ray_datenum.hour
    minutes = ray_datenum.minute
    seconds = ray_datenum.second

    freq = [fs]

    modeldump(year, month, day, hours, minutes, seconds) # run model dump to update plasmasphe

    datadir = '/home/rileyannereid/workspace/SR-output/' + 'kvecs2/'

    # -------------------------------- DEFINE RAY PARAMS --------------------------------
    # start position of raytracer
    SMcar_dsx, GEOcar_dsx, MAGsph_vpm, dir, dirstr = getDSXVPMpos(ray_datenum, 1)
    position = [float(SMcar_dsx.x), float(SMcar_dsx.y), float(SMcar_dsx.z)]

    theta = 45
    direction = getDIR(dir, ray_datenum, GEOcar_dsx, theta)

    # testing for fig in J. Bortnik thesis
    #jbpos = coord.Coords([1000e3+R_E, 30, 0], 'MAG', 'sph', units=['m', 'deg', 'deg'])
    #jbpos.ticks = Ticktock(ray_datenum, 'UTC')
    #newjbos = jbpos.convert('SM', 'car')
    #position = [float(newjbos.x), float(newjbos.y), float(newjbos.z)]
    #direction = coord.Coords([1,0,0], 'MAG', 'sph', units = ['m', 'm', 'm'])
    #direction.ticks = Ticktock(ray_datenum, 'UTC')
    #newdir = direction.convert('SM', 'car')
    #direction = [float(newdir.x), float(newdir.y), float(newdir.z)]
    
    # -------------------------------- RUN RAYS --------------------------------
    # convert for raytracer settings
    days_in_the_year = ray_datenum.timetuple().tm_yday
    days_in_the_year = format(days_in_the_year, '03d')

    # yearday and miliseconds day are used by raytracer
    yearday = str(year)+ str(days_in_the_year)   # YYYYDDD
    
    minutes = minutes + 2
    milliseconds_day = hours*3.6e6 + minutes*6e4 + seconds*1e3

    # run it!
    tmpdir = tempfile.mkdtemp() 
    
    run_rays(freq, [position], [direction], yearday, milliseconds_day, tmpdir)

    # -------------------------------- LOAD OUTPUT --------------------------------
    # Load all the rayfiles in the output directory
    ray_out_dir = tmpdir
    file_titles = os.listdir(ray_out_dir)

    # create empty lists to fill with ray files and damp files
    raylist = []

    for filename in file_titles:
        if '.ray' in filename:
            raylist += read_rayfile(os.path.join(ray_out_dir, filename))

    # -------------------------------- CONVERT COORDINATES --------------------------------
    # convert to desired coordinate system into vector list rays
    for r in raylist:
        ray = r # for output
        w = r['w']

        tmp_coords = coord.Coords(list(zip(r['pos'].x, r['pos'].y, r['pos'].z)), 'SM', 'car', units=['m', 'm', 'm'])
        tmp_kcoords = coord.Coords(list(zip((w/c) * r['n'].x, (w/c) * r['n'].y, (w/c) * r['n'].z)), 'SM', 'car', units=['m', 'm', 'm'])
        
        # convert to a unit vector first
        unitk = [(float(tmp_kcoords[s].x), float(tmp_kcoords[s].y), float(tmp_kcoords[s].z)) / np.sqrt(tmp_kcoords[s].x**2 + tmp_kcoords[s].y**2 + tmp_kcoords[s].z**2) for s in range(len(tmp_kcoords))]
        unitk_coords = coord.Coords(unitk, 'SM', 'car', units=['m', 'm', 'm'])
        tmp_kcoords = unitk_coords

        tvec_datetime = [ray_datenum + dt.timedelta(seconds=s) for s in r['time']]
        tmp_coords.ticks = Ticktock(tvec_datetime, 'UTC')  # add ticks
        tmp_kcoords.ticks = Ticktock(tvec_datetime, 'UTC')  # add ticks
        new_kcoords = tmp_kcoords.convert('MAG', 'car')
        MAGcar_k = new_kcoords

        tmp_coords.sim_time = r['time']
        new_coords = tmp_coords.convert('MAG', 'car')
        MAGcar_ray = new_coords
        
    # -------------------------------- PLOTTING --------------------------------
    fig, ax = plt.subplots(1,1,  figsize=(5, 10), sharex=True, sharey=True)

    lw = 2  # linewidth

    # rotate plot to be in plane of view
    MAGsph_dsx = GEOcar_dsx.convert('MAG', 'sph') # add ticks
    th = -MAGsph_dsx.long

    MAGcar_dsx = rotateplane(th, MAGsph_dsx, ray_datenum)
    MAGcar_vpm = rotateplane(th, MAGsph_vpm, ray_datenum)

    # plot sat locations
    plt.plot(MAGcar_dsx.x / R_E, MAGcar_dsx.z / R_E, '-go', zorder=105, label='DSX')
    plt.plot(MAGcar_vpm.x / R_E, MAGcar_vpm.z / R_E, '-yo', zorder=104, label='VPM')

    # rotate rays
    MAGsph_ray = MAGcar_ray.convert('MAG', 'sph')
    MAGsph_k = MAGcar_k.convert('MAG', 'sph')
    
    # have to loop to use rotate func
    mc_rayx = []
    mc_kx = []
    mc_rayz = []
    mc_kz = []

    for msph_ray, msph_k in zip(MAGsph_ray, MAGsph_k):
        MAGcar_ray = rotateplane(th, msph_ray, ray_datenum)
        mc_rayx.append(MAGcar_ray.x / R_E)
        mc_rayz.append(MAGcar_ray.z / R_E)
        MAGcar_k = rotateplane(th, msph_k, ray_datenum)
        mc_kx.append(MAGcar_k.x / R_E)
        mc_kz.append(MAGcar_k.z / R_E)

    ax.scatter(mc_rayx, mc_rayz, c = 'Black', s = 1, zorder = 103)
    #ax.quiver(MAGcar_ray.x[::intcheck] / R_E, MAGcar_ray.z[::intcheck] / R_E, MAGcar_k.x[::intcheck], MAGcar_k.z[::intcheck], color='black', zorder=104)
    
    for tti, tt in enumerate(mc_rayx): # just a nice way to get # of steps 
        if tti % intcheck == 0:   
            ax.text(mc_rayx[tti] - 0.1, mc_rayz[tti] - 0.05, str(int(tti/intcheck)), fontsize=10)

    # -------------------------------- EARTH AND IONO --------------------------------
    earth = plt.Circle((0, 0), 1, color='b', alpha=0.5, zorder=100)
    ax.add_artist(earth)

    # ---------------------------------- BFIELD -----------------------------------
    # (from IGRF13 model)
    L_shells = [2, 3]  # Field lines to draw

    for L in L_shells:
        Blines = []
        Lstart = [L, 0, 0]
        
        Lcoords = coord.Coords(Lstart, 'MAG', 'sph', units=['Re', 'deg', 'deg'])
        Lcoords.ticks = Ticktock(ray_datenum, 'UTC')
        GEO_Lcoords = Lcoords.convert('GEO', 'car')

        Blines.append(trace_fieldline_ODE([float(GEO_Lcoords.x),float(GEO_Lcoords.y),float(GEO_Lcoords.z)], 0, '0', 1, ray_datenum))
        Blines.append(trace_fieldline_ODE([float(GEO_Lcoords.x),float(GEO_Lcoords.y),float(GEO_Lcoords.z)], 0, '0', -1, ray_datenum))
        
        for blinex, bliney, blinez in Blines:
            raytimes = [ray_datenum for m in range(int(len(blinex)))]
            
            # convert coords
            bpos = np.column_stack((blinex, bliney, blinez))
            bpos = coord.Coords(bpos, 'GEO', 'car', units=['Re', 'Re', 'Re'])
            bpos.ticks = Ticktock(raytimes, 'UTC') # add ticks
            MAGcar_bline = bpos.convert('MAG', 'car')
            
            ax.plot(MAGcar_bline.x, MAGcar_bline.z, color='b', linewidth=1, linestyle='dashed')

    print('finished plotting field lines')

    bline_dsx = getDSXfieldline(GEOcar_dsx, ray_datenum, th)
    bline_dsx = bline_dsx[0]
    for bb in bline_dsx:
        ax.plot(bb.x, bb.z, color='r', linewidth=1, linestyle='dashed')

    # -------------------------------- FORMATTING --------------------------------
    ax.set_aspect('equal')
    max_lim = 3

    plt.xticks(np.arange(-max_lim, max_lim, step=1))
    plt.yticks(np.arange(-max_lim, max_lim, step=1))
    plt.xlabel('L (R$_E$)')
    plt.ylabel('L (R$_E$)')
    plt.xlim([0, max_lim])
    plt.ylim([-2, 2])

    mytitle = str(freq[0]/1e3) + 'kHz rays at ' + str(ray_datenum.month) + '-' + str(ray_datenum.day) + '-' + str(ray_datenum.hour) + ':' + str(minutes) + '\n' + str(theta) + ' intial angle'
    plt.title(mytitle)
    ax.legend(loc = 'lower center', fontsize =  'x-small')

    savename = datadir + str(freq[0]/1e3) + 'kHz_' + str(ray_datenum.month) + str(ray_datenum.day) + str(ray_datenum.year) + str(ray_datenum.hour) + str(minutes) + '_' + str(theta) + 'initialangle' + '.svg'
    plt.savefig(savename, format='svg')
    #plt.show()
    plt.close()

    return ray
# --------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------
def plotrefractive(ray, ray_datenum, intcheck):
    # set up phi vec
    phi_vec = np.linspace(0,360,int(1e5))*D2R
    w= ray['w']
    
    savenames = []
    for tti,tt in enumerate(ray['time']):
        
        if tti % intcheck == 0:
            t = tti

            bb = ray['B0'].iloc[t]
            bbdir = [bb.x / np.sqrt(bb.x**2 + bb.y**2 + bb.z**2), bb.y / np.sqrt(bb.x**2 + bb.y**2 + bb.z**2), bb.z /  np.sqrt(bb.x**2 + bb.y**2 + bb.z**2)]
            SMcar_bb = coord.Coords(bbdir, 'SM', 'car')
            SMcar_bb.ticks = Ticktock([ray_datenum for i in range(len(SMcar_bb.x))], 'UTC')
            MAGcar_bb = SMcar_bb.convert('MAG', 'car') 
            
            kk = (w/c) * ray['n'].iloc[t]
            kkdir = [kk.x / np.sqrt(kk.x**2 + kk.y**2 + kk.z**2), kk.y / np.sqrt(kk.x**2 + kk.y**2 + kk.z**2), kk.z /  np.sqrt(kk.x**2 + kk.y**2 + kk.z**2)]
            SMcar_kk = coord.Coords(kkdir, 'SM', 'car')
            SMcar_kk.ticks = Ticktock([ray_datenum for i in range(len(SMcar_bb.x))], 'UTC')
            MAGcar_kk = SMcar_kk.convert('MAG', 'car') 

            # anooying way of getting dot product but oh well need to check everything am tired
            vector_1 = [float(MAGcar_bb.x), float(MAGcar_bb.z)]
            vector_2 = [float(MAGcar_kk.x), float(MAGcar_kk.z)]
            unit_vector_1 = vector_1 / np.linalg.norm(vector_1)
            unit_vector_2 = vector_2 / np.linalg.norm(vector_2)
            dot_product = np.dot(unit_vector_1, unit_vector_2)
            ang = np.arccos(dot_product)

            # get stix param
            R, L, P, S, D = stix_parameters(ray, t, w)
            root = -1 # why ?? CAUSE WHISTLER

            k_vec = np.zeros_like(phi_vec)
            eta_vec = np.zeros_like(phi_vec)

            # solution from antenna white paper!
            resangle = np.arctan(np.sqrt(-P/S))

            for phi_ind, phi  in enumerate(phi_vec):

                # Solve the cold plasma dispersion relation
                cos2phi = pow(np.cos(phi),2)
                sin2phi = pow(np.sin(phi),2)

                A = S*sin2phi + P*cos2phi
                B = R*L*sin2phi + P*S*(1.0+cos2phi)

                discriminant = B*B - 4.0*A*R*L*P

                n1sq = (B + np.sqrt(discriminant))/(2.0*A)
                n2sq = (B - np.sqrt(discriminant))/(2.0*A)

                # negative refers to the fact that ^^^ B - sqrt
                n1 = np.sqrt(n1sq)
                n2 = np.sqrt(n2sq)
                
                k_vec[phi_ind] = w*n2/c
                eta_vec[phi_ind] = n2

            # plot it --------------------------
            fig, ax = plt.subplots(1,1)

            # rotate with bfield
            bang = np.arccos(MAGcar_bb.z / np.sqrt(MAGcar_bb.x**2 + MAGcar_bb.z**2))
            if MAGcar_bb.x < 0:
                bang = -bang
            else:
                bang = bang

            # plot the surface
            ax.plot(eta_vec*np.sin(phi_vec + bang), eta_vec*np.cos(phi_vec + bang), 'gray', LineWidth = 1, label = 'e + ions')

            # find eta at kvec
            etaind = min(range(len(phi_vec)), key=lambda i: abs(phi_vec[i]-ang))
            etaang = eta_vec[etaind]

            # note to self -- left off correctly plotting the kvector (yay! ) but still  not getting the line crossing
            # correct - maybe try a different method or use an intersection method? 
            # next, get the wavenormal figured out, but this is looking CORRECT!

            ax.plot([-1e5*MAGcar_bb.x, 1e5*MAGcar_bb.x], [-1e5*MAGcar_bb.z, 1e5*MAGcar_bb.z], 'b', linestyle='--', label = 'B0')
            ax.plot([0, etaang*np.sin(ang + bang)], [0, etaang*np.cos(ang + bang)], 'r', label = 'kvec')
            #ax.quiver(0, 0, MAGcar_kk.x, MAGcar_kk.z)

            # find normal at that point
            f1 = eta_vec[etaind-1]*np.cos(phi_vec[etaind-1]+bang)
            f2 = eta_vec[etaind+1]*np.cos(phi_vec[etaind+1]+bang)
            x1 = eta_vec[etaind-1]*np.sin(phi_vec[etaind-1]+bang)
            x2 = eta_vec[etaind+1]*np.sin(phi_vec[etaind+1]+bang)
            
            # this was all Sam
            taneta = (f2-f1)/(x2-x1)
            normeta = -1/taneta
            ntheta = np.arctan2(x1-x2, f2-f1)
            if ntheta < 0: 
                ntheta = ntheta + np.pi

            pp = ntheta
            # bottom surface
            #if ang > np.pi/2:
            #    ntheta = ntheta + np.pi

            if ntheta > np.pi/2 and ang < np.pi/2:
                ntheta = ntheta + np.pi
            
            # make a line
            # y = mx + b
            intercept = (etaang * np.cos(ang + bang)) - (normeta * etaang * np.sin(ang + bang))
            #ax.scatter(x1, f1, label='1')
            #ax.scatter(x2, f2, label='2')
            #ax.plot([etaang * np.sin(ang + bang), etaang * np.cos(ang + bang), 1, normeta * 1 + intercept) 
            ax.quiver(etaang*np.sin(ang + bang), etaang*np.cos(ang + bang), np.cos(ntheta),  np.sin(ntheta))
            #ax.quiver(etaang*np.sin(ang + bang), etaang*np.cos(ang + bang)

            # ------------------------------------ formatting ----------------------------------------
            xlim1 = -300
            xlim2 = -xlim1
            ylim1 = xlim1
            ylim2 = -ylim1
            scale = xlim2/500

            ax.set_xlim([xlim1, xlim2])
            ax.set_ylim([ylim1, ylim2])
            ax.set_xlabel('Transverse Refractive Component')

            #ax.annotate('fp = ' + str(float(round((np.sqrt(wp2)/(2*np.pi))/1e3, 1))) + ' kHz', xy=(xlim1+(scale*100),ylim2-(scale*200)))

            resonanceangle = float(R2D*resangle)
            ax.annotate('${\Theta}$ < ' + str(round(resonanceangle, 2)) + 'deg', xy=(xlim2 - (scale*200), scale*50))
            
            rename = str(ray_datenum.month) + str(ray_datenum.day) + str(ray_datenum.year)

            ax.set_title(str(round(w/(2*np.pi*1e3),1)) + ' kHz Refractive Surface at ' + str(ray_datenum) + '\n' + r"$\bf{" + str(int(tti / intcheck)) + "}$")
            plt.legend(loc='upper right')
            datadir = '/home/rileyannereid/workspace/SR-output/' + 'kvecs2/'

            imgdir = datadir + str(round(w/(2*np.pi*1e3),1)) + 'kHz' + rename + str(ray_datenum.hour) + str(ray_datenum.minute) + 'refractivesurfaces/'
            
            try:
                os.mkdir(imgdir)
            except OSError:
                pass
            else:
                pass
            
            plt.savefig(imgdir + str(round(w/(2*np.pi*1e3),1)) + 'kHz' + rename + str(ray_datenum.hour) + str(ray_datenum.minute) + 'refractivesurface' + str(tti) + '.png', format='png')
            
            # debug
            if tti / intcheck > 70:
                print(ntheta)
                print('pp=', pp)
                #plt.show()
            savenames.append(imgdir + str(round(w/(2*np.pi*1e3),1)) + 'kHz' + rename + str(ray_datenum.hour) + str(ray_datenum.minute) + 'refractivesurface' + str(tti) + '.png')
            plt.close()

    simplegifs(savenames, datadir + str(round(w/(2*np.pi*1e3),1)) + 'kHz' + rename + str(ray_datenum.hour) + str(ray_datenum.minute) + '.gif')
    
    return

# ------------------------------------------- END --------------------------------------------

intcheck = 1
ray_datenum = dt.datetime(2020,9,14,22,53)
fs = 8.2e3
ray = plotray2Ddir(ray_datenum, fs, intcheck)
plotrefractive(ray, ray_datenum, intcheck)