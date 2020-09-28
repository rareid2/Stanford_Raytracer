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

# Define tangent line
# y = m*(x - x1) + y1
def tanline(x, x1, y1, m):
    return m*(x - x1) + y1

def plotraydir(dates, fs, bs):

    for cdate, cf, bsstr in zip(dates, fs, bs):
        year = cdate.year
        month = cdate.month
        day = cdate.day
        hours = cdate.hour
        minutes = cdate.minute
        seconds = cdate.second

        freq = [cf]
        intcheck = 10

        ray_datenum = dt.datetime(year, month, day, hours, minutes, seconds)
        modeldump(year, month, day, hours, minutes, seconds) # run model dump to update plasmasphere
        #thetalist = [0, 15, 30, 45, -15, -30, -45] # in deg -- what angles to launch at? 

        checkdir = 0
        crs_out = 'MAG'  # theres a bug with GEO coords -- maybe its the fieldlines? -- DONT CAHNGE THINS
        datadir = '/home/rileyannereid/workspace/SR-output/' + 'kvecs2/'

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
        #quick check real quick

        jbpos = coord.Coords([1000e3+R_E, 30, 0], 'MAG', 'sph', units=['m', 'deg', 'deg'])
        jbpos.ticks = Ticktock(ray_datenum, 'UTC')
        newjbos = jbpos.convert('SM', 'car')
        newjbosGEO = jbpos.convert('GEO', 'car')
        position = [float(newjbos.x), float(newjbos.y), float(newjbos.z)]

        # check with hemi we are in
        if outsph_vpm.lati > 0:
            dir = 1   # north
            dirstr = 'north'
        else:
            dir = -1  # south
            dirstr = 'south'

        dir = -1
        dirstr = 'south'

        Bstart = [float(GEOcar_dsx.x)/R_E, float(GEOcar_dsx.y)/R_E, float(GEOcar_dsx.z)/R_E]
        Bstart =  [float(newjbosGEO.x)/R_E, float(newjbosGEO.y)/R_E, float(newjbosGEO.z)/R_E]
        Bx, By, Bz = B_direasy(ray_datenum, Bstart, dir)

        # convert direction to SM coordinates for raytracer
        dirB = np.reshape(np.array([Bx, By, Bz]), (1, 3))
        dirB = coord.Coords(dirB[0], 'GEO', 'car', units=['Re', 'Re', 'Re'])
        dirB.ticks = Ticktock(ray_datenum, 'UTC') # add ticks
        SMsph_dirB = dirB.convert('SM', 'sph')

        # fill for raytracer call
        positions = []
        directions = []
        thetalist = [0]

        # rotate directions
        for theta in thetalist:
            # increase (or decrease) polar angle
            newth = float(SMsph_dirB.lati) + theta
            Rot_dirB = [float(SMsph_dirB.radi), newth, float(SMsph_dirB.long)] 
            Rot_dirB = coord.Coords(Rot_dirB, 'SM', 'sph')
            Rot_dirB.ticks = Ticktock(ray_datenum, 'UTC') # add ticks
            SMcar_dirB = Rot_dirB.convert('SM', 'car')

            direction = [float(SMcar_dirB.x), float(SMcar_dirB.y), float(SMcar_dirB.z)]
            
            # add the normalized direction (or zeros)
            directions.append(np.squeeze(direction))

            # make sure position list matches direction list
            positions.append(position)
        
        #quick check
        direction = coord.Coords([1,0,0], 'MAG', 'sph', units = ['m', 'm', 'm'])
        direction.ticks = Ticktock(ray_datenum, 'UTC')
        newdir = direction.convert('SM', 'car')
        directions = [[float(newdir.x), float(newdir.y), float(newdir.z)]]

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
        kvecs = []
        for r in raylist:
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
            new_kcoords = tmp_kcoords.convert(crs_out, 'car')
            kvecs.append(new_kcoords)

            tmp_coords.sim_time = r['time']
            new_coords = tmp_coords.convert(crs_out, 'sph') # needs to be sph for rotation
            rays.append(new_coords)
            
        # -------------------------------- PLOTTING --------------------------------
        fig, ax = plt.subplots(1,1, sharex=True, sharey=True)
        lw = 2  # linewidth

        # rotate plot to be in plane of view
        outsph_dsx = GEIcar_dsx.convert(crs_out, 'sph') # add ticks
        th = -outsph_dsx.long
        th = 0

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

        # rotate rays

        for r, k in zip(rays, kvecs):
            print(k)
            rrad = []
            rlat = []
            rlon = []
            rrad.append(r.radi)
            rlat.append(r.lati)
            rlon.append(r.long)

            rrlon = [rl + th for rl in rlon]
            rcoords = [np.column_stack([rr, rl, rrl]) for rr, rl, rrl in zip(rrad, rlat, rrlon)]

            Rot_ray = coord.Coords(rcoords[0], crs_out, 'sph', units=['m', 'deg', 'deg'])

            newtvec = [ray_datenum for i in range(len(Rot_ray))] #WHYY???
            Rot_ray.ticks = Ticktock(newtvec, 'UTC')
            outcar_ray = Rot_ray.convert(crs_out, 'car')

            if len(outcar_ray.x) > 1:
                # plotp = ax.scatter(MAGcar_ray.x / R_E, MAGcar_ray.z / R_E, c=d, s = 1, cmap = 'Reds', vmin = 0, vmax = 1.5, zorder = 103)
                plotp = ax.scatter(outcar_ray.x / R_E, outcar_ray.z / R_E, c = 'Black', s = 1, zorder = 103)
                ax.quiver(outcar_ray.x[::intcheck] / R_E, outcar_ray.z[::intcheck] / R_E, k.x[::intcheck], k.z[::intcheck], zorder=104)
            for tti, tt in enumerate(r.radi): 
                if tti % intcheck == 0:   
                    ax.text(outcar_ray.x[tti] / R_E, (outcar_ray.z[tti] - 100e3) / R_E, str(int(tti/intcheck)), fontsize=10)

        # -------------------------------- EARTH AND IONO --------------------------------
        earth = plt.Circle((0, 0), 1, color='b', alpha=0.5, zorder=100)
        ax.add_artist(earth)

        # ---------------------------------- BFIELD -----------------------------------
        # (from IGRF13 model)
        L_shells = [2, 3, 4, -2, -3, -4]  # Field lines to draw

        for L in L_shells:
            Blines = []
            if L < 0:
                rot = 1
            else:
                rot = -1

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

        # get an entire ray lets just try one for now
        for r in raylist:
            ray = r
            
        w = ray['w']           

        f = w/(2*np.pi)
        print('frequency of ray is:', w/(2*np.pi), ' Hz')

        # set up phi vec
        phi_vec = np.linspace(0,360,int(1e5))*D2R

        # -------------------------------- FORMATTING --------------------------------
        ax.set_aspect('equal')
        max_lim = 3

        plt.xticks(np.arange(-max_lim, max_lim, step=1))
        plt.yticks(np.arange(-max_lim, max_lim, step=1))
        plt.xlabel('L (R$_E$)')
        plt.ylabel('L (R$_E$)')
        plt.xlim([0, max_lim])
        plt.ylim([-2, 2])

        mytitle = str(freq[0]/1e3) + 'kHz rays at ' + str(ray_datenum.month) + '-' + str(ray_datenum.day) + '-' + str(ray_datenum.hour) + ':' + str(minutes) + '\n' + str(thetalist[0]) + 'intial angle'
        plt.title(mytitle)
        ax.legend(loc = 'lower center', fontsize =  'x-small')

        savename = datadir + str(freq[0]/1e3) + 'kHz_' + str(ray_datenum.month) + str(ray_datenum.day) + str(ray_datenum.year) + str(ray_datenum.hour) + str(minutes) + '_' + str(thetalist[0]) + 'initialangle' + '.png'
        plt.savefig(savename, format='png')
        #plt.show()
        plt.close()

        # grab every xseconds seconds
        savenames = []
        for tti,tt in enumerate(ray['time']):
            
            if tti % intcheck == 0:
                # grab only t = 0
                t = tti

                bb = ray['B0'].iloc[t]
                bbdir = [bb.x / np.sqrt(bb.x**2 + bb.y**2 + bb.z**2), bb.y / np.sqrt(bb.x**2 + bb.y**2 + bb.z**2), bb.z /  np.sqrt(bb.x**2 + bb.y**2 + bb.z**2)]
                kk = (w/c) * ray['n'].iloc[t]
                kkdir = [kk.x / np.sqrt(kk.x**2 + kk.y**2 + kk.z**2), kk.y / np.sqrt(kk.x**2 + kk.y**2 + kk.z**2), kk.z /  np.sqrt(kk.x**2 + kk.y**2 + kk.z**2)]

                ang = math.acos((bbdir[0] * kkdir[0]) + (bbdir[2] * kkdir[2]) / (np.sqrt(kkdir[0]**2 + kkdir[2]**2) + np.sqrt(bbdir[0]**2 + bbdir[2]**2)))
                
                # for earlier use
                Lshell = getLshell(ray, t, ray_datenum)

                # get stix param
                R, L, P, S, D = stix_parameters(ray, t, w)

                root = -1 # why ?? CAUSE WHISTLER

                k_vec = np.zeros_like(phi_vec)
                eta_vec=np.zeros_like(phi_vec)

                # solution from antenna white paper!
                resangle = np.arctan(np.sqrt(-P/S))

                # cone = []
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
                    #if n2sq < 0:
                    #    cone.append(phi)
                    # only get whistler solution from minus root ( i think)
                    # important to call these plus and minus roots! NOT POS AND NEG
                    n2 = np.sqrt(n2sq)
                    # Order the roots -- ?
                    """
                    if abs(n1) > abs(n2):
                        k2 = w*n1/c
                        k1 = w*n2/c
                    else:
                        k1 = w*n1/c
                        k2 = w*n2/c

                    if root==1.0:
                        k = k1
                        eta = n1
                    else:
                        k = k2
                        eta = n2
                    """
                    k_vec[phi_ind] = w*n2/c
                    eta_vec[phi_ind] = n2

                # repeat for only AH solution

                # grab for ONLY electrons -- should I assume electron fill ENTIRE population...?
                Ns = float(ray['Ns'].iloc[t,0])
                Q = float(ray['qs'].iloc[t,0])
                M = float(ray['ms'].iloc[t,0])
                B   =  ray['B0'].iloc[t]
                Bmag = np.linalg.norm(B)

                # makes it easier to square here
                w2 = w*w
                wp2 = Ns*pow(Q,2)/(eo*M)
                wh = Q*Bmag/M
                wh2 = wh*wh

                root = 1

                # solve appleton-hartree eq
                numerator = wp2/w2
                denom1 = (wh2*pow(np.sin(phi_vec),2))/(2*(w2 - wp2))
                denom2 = np.sqrt(pow(denom1, 2) + wh2*pow(np.cos(phi_vec), 2)/w2)
                eta2_AH   = 1 - (numerator/(1 - denom1 + root*denom2))
                eta_AH = np.sqrt(-eta2_AH)

                # plot it 
                fig, ax = plt.subplots(1,1)

                #ax.plot(eta_AH*np.sin(phi_vec), eta_AH*np.cos(phi_vec), LineWidth = 1, label = 'e only')
                ax.plot(eta_vec*np.sin(phi_vec), eta_vec*np.cos(phi_vec), LineWidth = 1, label = 'e + ions')
                
                #figure out how to get correct eta
                etaind = min(range(len(phi_vec)), key=lambda i: abs(phi_vec[i]-ang))
                etaang = eta_vec[etaind]
                ax.plot([0, etaang*np.sin(ang)], [0, etaang*np.cos(ang)]) # plots the kvec

                # find normal at that point
                f1 = eta_vec[etaind-1]*np.cos(phi_vec[etaind-1])
                f2 = eta_vec[etaind+1]*np.cos(phi_vec[etaind+1])
                dx = phi_vec[etaind+1] - phi_vec[etaind-1]
                taneta = (f2-f1)/dx
                normeta = -1/taneta

                #xrange = np.linspace(eta_vec[etaind-1000]*np.sin(phi_vec[etaind-1000]), eta_vec[etaind+1000]*np.sin(phi_vec[etaind+1000]), 10)
                #ax.quiver(etaang*np.sin(ang), etaang*np.cos(ang), 1/np.sqrt(1+normeta**2), normeta/np.sqrt(1+normeta**2))

                # lol dont do this
                #findcone = eta_vec*np.sin(phi_vec)
                #for cone in findcone:
                    #if cone > 100:
                    #    wherecone = np.where(findcone == cone)
                    #    conetheta = phi_vec[wherecone]
                    #    conetheta = conetheta * R2D
                    #    break

                #archeight = float(eta_vec[wherecone])

                # formatting
                xlim1 = -300
                xlim2 = -xlim1
                ylim1 = xlim1
                ylim2 = -ylim1

                scale = xlim2/500

                ax.set_xlim([xlim1, xlim2])
                ax.set_ylim([ylim1, ylim2])
                ax.set_xlabel('Transverse Refractive Component')

                ax.arrow(0, ylim1, 0, 2*ylim2-1, length_includes_head=True, head_width=3, head_length=5, color = 'grey', ls = '--')
                ax.annotate('B0', xy=(scale*10,ylim2-(scale*50)))

                ax.annotate('fp = ' + str(float(round((np.sqrt(wp2)/(2*np.pi))/1e3, 1))) + ' kHz', xy=(xlim1+(scale*100),ylim2-(scale*200)))

                resonanceangle = float(R2D*resangle)

                #pac = Patch.Arc([0, 0], archeight, archeight, angle=0, theta1=0, theta2=float(resonanceangle), edgecolor = 'r')
                #ax.add_patch(pac)

                ax.annotate('${\Theta}$ < ' + str(round(resonanceangle, 2)) + 'deg', xy=(xlim2 - (scale*200), scale*50))
                rename = str(ray_datenum.month) + str(ray_datenum.day) + str(ray_datenum.year)

                ax.set_title(str(freq[0]/1e3) + ' kHz Refractive Surface at ' + str(ray_datenum) + ' ' + str(float(tti / intcheck)) + ' pos')
                plt.legend()

                imgdir = datadir + str(freq[0]/1e3) + 'kHz' + rename + str(ray_datenum.hour) + str(ray_datenum.minute) + 'refractivesurfaces/'
                try:
                    os.mkdir(imgdir)
                except OSError:
                    pass
                else:
                    pass
                
                plt.savefig(imgdir + str(freq[0]/1e3) + 'kHz' + rename + str(ray_datenum.hour) + str(ray_datenum.minute) + 'refractivesurface' + str(tti) + '.png', format='png')
                #plt.show()
                savenames.append(imgdir + str(freq[0]/1e3) + 'kHz' + rename + str(ray_datenum.hour) + str(ray_datenum.minute) + 'refractivesurface' + str(tti) + '.png')
                plt.close()

        simplegifs(savenames, datadir + str(freq[0]/1e3) + 'kHz' + rename + str(ray_datenum.hour) + str(ray_datenum.minute) + '.gif')

        # ------------------------------------------- END --------------------------------------------
        """
        fig, ax = plt.subplots(1,1, sharex=True, sharey=True)
        lw = 2  # linewidth


        # full plot!

        # repeat for each kvec @ 475
        mycolors = ['blue', 'red', 'cyan', 'brown', 'green', 'magenta', 'yellow', 'purple', 'orange', 'pink', 'maroon']

        i = 0

        for r, k in zip(rays, kvecs):
            th = 0

            rrad = []
            rlat = []
            rlon = []
            rrad.append(r.radi)
            rlat.append(r.lati)
            rlon.append(r.long)

            rrlon = [rl + th for rl in rlon]
            rcoords = [np.column_stack([rr, rl, rrl]) for rr, rl, rrl in zip(rrad, rlat, rrlon)]

            Rot_ray = coord.Coords(rcoords[0], crs_out, 'sph', units=['m', 'deg', 'deg'])
            newtvec = [ray_datenum for i in range(len(Rot_ray))] #WHYY???
            Rot_ray.ticks = Ticktock(newtvec, 'UTC')
            outcar_ray = Rot_ray.convert('GEO', 'car')
            
            sc = 1

            if i == 0: 
                B475 = [float(outcar_ray[-1].x)/R_E, float(outcar_ray[-1].y)/R_E, float(outcar_ray[-1].z)/R_E]
                Bx, By, Bz = B_direasy(ray_datenum, B475, dir)
                ax.quiver(0,0, sc*Bx, sc*Bz, label = 'Bfield@475km', angles='xy', scale_units='xy', scale=25, linewidth = 1,)

            if r.radi[-1] < 500e3 + R_E:

                # plotp = ax.scatter(MAGcar_ray.x / R_E, MAGcar_ray.z / R_E, c=d, s = 1, cmap = 'Reds', vmin = 0, vmax = 1.5, zorder = 103)
                ax.quiver(0, 0, sc*k.x[-1], sc*k.z[-1], label= thetalist[i], angles='xy', scale_units='xy', scale=25, linewidth = 1, color = mycolors[i])
            i += 1

        hh = 100
        #plt.xlim([-hh, hh])
        #plt.ylim([-hh, hh])
        mytitle = str(freq[0]/1e3) + 'kHz rays at ' + str(ray_datenum.month) + '-' + str(ray_datenum.day) + '-' + str(ray_datenum.hour) + ':' + str(minutes)
        plt.title(mytitle)
        ax.legend(fontsize =  'x-small')

        savename = datadir + str(freq[0]/1e3) + 'kHz_' + str(ray_datenum.month) + str(ray_datenum.day) + str(ray_datenum.year) + str(ray_datenum.hour) + str(minutes) + '_all' + '.png'
        plt.savefig(savename, format='png')
        #plt.show()
        plt.close()
    """

dates = [dt.datetime(2020,9,14,22,53)]

fs = [3e3]
bs = ['kvecs2' for i in range(len(dates))]
plotraydir(dates, fs, bs)