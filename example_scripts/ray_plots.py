# import needed packages
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.patches import Circle
from matplotlib.offsetbox import (TextArea, DrawingArea, OffsetImage,
                                  AnnotationBbox)
from matplotlib.cbook import get_sample_data

import cartopy.crs as ccrs
import datetime as dt
from spacepy import coordinates as coord
from spacepy.time import Ticktock

from constants_settings import *
from coordinates import create_spc, convert_spc
from raytracer_utils import read_rayfile,readdump, get_yearmiliday 
from bfield import Bfieldinfo, trace_fieldline_ODE


#from plotrefractivesurface import getLshell, stix_parameters

#---------------------------------------------------------------------------
def rotateplane(plane_long, rc, tvec_datetime, crs, carsph, units):
    # rotate long to be in same plane
    rot_crs = convert_spc(rc, tvec_datetime, 'GEO', 'sph', units=['Re','deg','deg'])
    rot_lon = [plane_long for new_lon in rot_crs.long]
    rot_crs_new = create_spc(list(zip(rot_crs.radi, rot_crs.lati, rot_lon)), tvec_datetime, 'GEO', 'sph', units=['Re','deg','deg'])
    final_crs = convert_spc(rot_crs_new, tvec_datetime, crs, carsph, units)
    return final_crs
#---------------------------------------------------------------------------

#---------------------------------------------------------------------------
def get_Lshells(lshells, tvec_datetime, crs, carsph, units):
    Lshell_fline = []
    for lsh in lshells:

        L = Bfieldinfo()
        L.time = tvec_datetime[0]
        stpos = create_spc(cor_array=[lsh,0,0], dt_array=L.time, crs='GEO', carsph='car', units=['Re','Re','Re'])
        L.pos = stpos

        # trace a fieldline to the north + south hemisphere
        trace_fieldline_ODE(L, 'north', crs, carsph, units)
        Lshell_fline.append(L.fieldline)
        trace_fieldline_ODE(L, 'south', crs, carsph, units)
        Lshell_fline.append(L.fieldline)

    return Lshell_fline
#---------------------------------------------------------------------------

def get_plasmasphere():
    # need to figure out how to change this coordinate system! TODO!
    path2plasma = '/modeldumps/'
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

#---------------------------------------------------------------------------
def plotray2D(ray_datenum, raylist, ray_out_dir, crs, carsph, units, plot_kvec=True, show_plot=True):

    # convert to desired coordinate system into vector list rays
    ray_coords = []
    k_coords = []
    for r in raylist:
        w = r['w']

        tmp_coords = coord.Coords(list(zip(r['pos'].x, r['pos'].y, r['pos'].z)), 'SM', 'car', units=['m', 'm', 'm'])
        tmp_kcoords = coord.Coords(list(zip((w/C) * r['n'].x, (w/C) * r['n'].y, (w/C) * r['n'].z)), 'SM', 'car', units=['m', 'm', 'm'])
        print(tmp_kcoords)
        # convert to a unit vector first
        unitk = [(float(tmp_kcoords[s].x), float(tmp_kcoords[s].y), float(tmp_kcoords[s].z)) / np.sqrt(tmp_kcoords[s].x**2 + tmp_kcoords[s].y**2 + tmp_kcoords[s].z**2) for s in range(len(tmp_kcoords))]

        tvec_datetime = [ray_datenum + dt.timedelta(seconds=s) for s in r['time']]
        tmp_kcoords = create_spc(unitk, tvec_datetime, 'SM', 'car', units=['m','m','m'])

        tmp_coords.ticks = Ticktock(tvec_datetime, 'UTC')  # add ticks
        tmp_kcoords.ticks = Ticktock(tvec_datetime, 'UTC') # add ticks

        new_kcoords = convert_spc(tmp_kcoords, tvec_datetime, crs, carsph, units)

        tmp_coords.sim_time = r['time']
        new_coords = convert_spc(tmp_coords, tvec_datetime, crs, carsph, units)

        # save it
        ray_coords.append(new_coords)
        k_coords.append(new_kcoords)
        
    # -------------------------------- PLOTTING --------------------------------
    
    # make this into a an image (sneaky)
    ray1 = ray_coords[0][0]
    ray1_sph = convert_spc(ray1, tvec_datetime, 'GEO', 'sph', units=['m','deg','deg'])

    ax = plt.axes(projection=ccrs.Orthographic(central_longitude=float(ray1_sph.long)-90, central_latitude=0.0))
    ax.coastlines()
    plt.savefig(ray_out_dir+'/ccrs_proj.png')
    plt.close()

    img = plt.imread(ray_out_dir+'/ccrs_proj.png', format='png')
    fig, ax = plt.subplots()
    ax.imshow(img, extent=[-1.62,1.62,-1.3,1.3]) # not the most scientific

    # plot the rays and kvec
    for rc, kc in zip(ray_coords, k_coords):
        rotated_rcoords = rotateplane(float(ray1_sph.long), rc, tvec_datetime, crs, carsph, units)
    if plot_kvec:
        int_plt = len(rotated_rcoords)//10
        ax.quiver(rotated_rcoords.x[::int_plt], rotated_rcoords.z[::int_plt], kc.x[::int_plt], kc.z[::int_plt], color='Black', zorder=104)
        #for ii,i in enumerate(range(0,len(rotated_rcoords),int_plt)): # not the prettiest but it will work, just trying to annoate
            #ax.text(rotated_rcoords.x[ii], rotated_rcoords.z[ii], str(ii))
    else:
        ax.scatter(rotated_rcoords.x, rotated_rcoords.z, c = 'Black', s = 1, zorder = 103)

    # plot field lines (from IGRF13 model)
    L_shells = [2, 3, 4]  # Field lines to draw
    Lshell_flines = get_Lshells(L_shells, tvec_datetime, crs, carsph, units)
    
    for lfline in Lshell_flines:
        ax.plot(lfline.x, lfline.z, color='Black', linewidth=1, linestyle='dashed')

    plt.xlabel('L (R$_E$)')
    plt.ylabel('L (R$_E$)')
    plt.xlim([0, max(L_shells)])
    plt.ylim([-2, 2])

    if show_plot:
        plt.show()

    return fig
# --------------------------------------------------------------------------------------

# ---------------------------------------- STIX PARAM --------------------------------------------
def stix_parameters(ray, t, w, tvec_datetime, rspot,crs,carsph):

    rb = Bfieldinfo()
    rb.time = tvec_datetime[t]
    rb.pos = rspot
    rb.Bfield_direction('north',crs,carsph)
    newb = rb.bvec
    Bmag = np.sqrt(newb.x**2+newb.y**2+newb.z**2) * 1e-9

    B   =  ray['B0'].iloc[t]
    Bmag2 = np.linalg.norm(B)

    print(Bmag, Bmag2)
    Q    = np.abs(np.array(ray['qs'].iloc[t,:]))
    M    = np.array(ray['ms'].iloc[t,:])
    Ns   = np.array(ray['Ns'].iloc[t,:])

    Wcs   = Q*Bmag/M
    Wps2  = Ns*pow(Q,2)/EPS0/M

    R = 1.0 - np.sum(Wps2/(w*(w + Wcs)))
    L = 1.0 - np.sum(Wps2/(w*(w - Wcs)))
    P = 1.0 - np.sum(Wps2/(w*w))
    S = (R+L)/2.0
    D = (R-L)/2.0

    return R, L, P, S, D
# ---------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------
def plotrefractivesurface(ray_datenum, ray,crs, carsph, units):
    # set up phi vec
    phi_vec = np.linspace(0,360,int(1e5))*D2R
    w = ray['w']

    tmp_kcoords = coord.Coords(list(zip((w/C) * ray['n'].x, (w/C) * ray['n'].y, (w/C) * ray['n'].z)), 'SM', 'car', units=['m', 'm', 'm'])
        
    # convert to a unit vector first
    unitk = [(float(tmp_kcoords[s].x), float(tmp_kcoords[s].y), float(tmp_kcoords[s].z)) / np.sqrt(tmp_kcoords[s].x**2 + tmp_kcoords[s].y**2 + tmp_kcoords[s].z**2) for s in range(len(tmp_kcoords))]

    tvec_datetime = [ray_datenum + dt.timedelta(seconds=s) for s in ray['time']]
    kcoords = create_spc(unitk, tvec_datetime, 'SM', 'car', units=['m','m','m'])
    kcoords.ticks = Ticktock(tvec_datetime)
    kcoords_rot = convert_spc(kcoords, tvec_datetime, crs, carsph, units)

    tmp_coords = coord.Coords(list(zip(ray['pos'].x, ray['pos'].y, ray['pos'].z)), 'SM', 'car', units=['m', 'm', 'm'])

    tvec_datetime = [ray_datenum + dt.timedelta(seconds=s) for s in ray['time']]

    tmp_coords.ticks = Ticktock(tvec_datetime, 'UTC')  # add ticks
    tmp_coords.sim_time = ray['time']
    r_coords = convert_spc(tmp_coords, tvec_datetime, crs, carsph, units)
    
    int_plt = len(kcoords)//10

    for ti, t in enumerate(ray['time']):
        if ti % int_plt == 0:
            print('plotting a refractive surface')
            # get stix param
            R, L, P, S, D = stix_parameters(ray, ti, w, tvec_datetime, r_coords[ti],crs,carsph)
            root = -1 # whistler solution

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
                
                k_vec[phi_ind] = w*n2/C
                eta_vec[phi_ind] = n2

            # plot it --------------------------
            fig, ax = plt.subplots(1,1)

            # plot the surface
            ax.plot(eta_vec*np.sin(phi_vec), eta_vec*np.cos(phi_vec), 'gray', LineWidth = 1, label = 'e + ions')

            B   =  ray['B0'].iloc[ti]
            Bmag = np.linalg.norm(B)

            rb = Bfieldinfo()
            rb.time = tvec_datetime[ti]
            rb.pos = r_coords[ti]
            rb.Bfield_direction('north',crs,carsph)
            newb = rb.unit_vec
            
            kmag = np.sqrt(kcoords_rot.x[ti]**2 + kcoords_rot.z[ti]**2)
            kunit = [kcoords_rot.x[ti]/kmag, kcoords_rot.z[ti]/kmag]
            print(kcoords_rot.y[ti])
            
            bmag = np.sqrt(newb.x**2 + newb.z**2)
            bunit = [newb.x/bmag, newb.z/bmag]

            alpha = np.arccos(bunit[0]*kunit[0]+bunit[1]*kunit[1])
            alphaedg = float(alpha)*R2D

            ax.quiver(0,0, 50*np.sin(float(alpha)), 50*np.cos(float(alpha)))
            #ax.quiver(0,0,bunit[0], bunit[1])


            xlim1 = -100
            xlim2 = -xlim1
            ylim1 = xlim1
            ylim2 = -ylim1

            ax.set_xlim([xlim1, xlim2])
            ax.set_ylim([ylim1, ylim2])
            ax.set_xlabel('Transverse Refractive Component')
            plt.legend(loc='upper right')
            resonanceangle = float(R2D*resangle)
            ax.text(25, -75,'resonance cone < '+ str(round(resonanceangle,2)))
            ax.text(25, 70, 'theta = '+ str(round(alphaedg,2)))

            
            rayfile_directory = '/home/rileyannereid/workspace/SR-output/rayfiles'
            ray_out_dir = rayfile_directory + '/'+ dt.datetime.strftime(ray_datenum, '%Y-%m-%d %H:%M:%S')
            plt.savefig(ray_out_dir+'/refractive_surface'+str(ti)+'.png')

        # note to self -- left off correctly plotting the kvector (yay! ) but still  not getting the line crossing
        # correct - maybe try a different method or use an intersection method? 
        # next, get the wavenormal figured out, but this is looking good

        # find normal at that point
        #f1 = eta_vec[etaind-1]*np.cos(phi_vec[etaind-1]+bang)
        #f2 = eta_vec[etaind+1]*np.cos(phi_vec[etaind+1]+bang)
        #x1 = eta_vec[etaind-1]*np.sin(phi_vec[etaind-1]+bang)
        #x2 = eta_vec[etaind+1]*np.sin(phi_vec[etaind+1]+bang)

        # this was all Sam
        #taneta = (f2-f1)/(x2-x1)
        #normeta = -1/taneta
        #ntheta = np.arctan2(x1-x2, f2-f1)
        #if ntheta < 0: 
        #    ntheta = ntheta + np.pi

        #pp = ntheta
        # bottom surface
        # for JB thesis, uncomment this and comment out the next condition
        #if ang > np.pi/2:
        #    ntheta = ntheta + np.pi

        #if ntheta > np.pi/2 and ang < np.pi/2:
        #    ntheta = ntheta + np.pi

        # make a line
        # y = mx + b
        #intercept = (etaang * np.cos(ang + bang)) - (normeta * etaang * np.sin(ang + bang))
        #ax.scatter(x1, f1, label='1')
        #ax.scatter(x2, f2, label='2')
        #ax.plot([etaang * np.sin(ang + bang), etaang * np.cos(ang + bang), 1, normeta * 1 + intercept) 
        #ax.quiver(etaang*np.sin(ang + bang), etaang*np.cos(ang + bang), np.cos(ntheta),  np.sin(ntheta))
        #ax.quiver(etaang*np.sin(ang + bang), etaang*np.cos(ang + bang)

    return

# ------------------------------------------- END --------------------------------------------