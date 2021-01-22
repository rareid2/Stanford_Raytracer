import numpy as np
import spacepy.irbempy as irbem
import spacepy.coordinates as coord
from spacepy.time import Ticktock
import datetime as dt
from scipy.integrate import ode

from coordinates import create_spc, convert_spc
from constants_settings import *

# -------------------------------- Bfieldinfo CLASS -----------------------------------------
# just a nice way to mangage getting bfield info

# uses irbempy from spacepy
# set time and position 
#
# see https://spacepy.github.io/irbempy.html for opts description:
# extfield:
#   ‘0’ = No external field model
#   ‘MEAD’ = Mead and Fairfield
#   ‘T87SHORT’ = Tsyganenko 1987 short (inner magnetosphere)
#   ‘T87LONG’ = Tsyganenko 1987 long (valid in extended tail region)
#   ‘T89’ = Tsyganenko 1989
#   ‘OPQUIET’ = Olsen-Pfitzer static model for quiet conditions
#   ‘OPDYN’ = Olsen-Pfitzer static model for active conditions
#   ‘T96’ = Tsyganenko 1996
#   ‘OSTA’ = Ostapenko and Maltsev
#   ‘T01QUIET’ = Tsyganenko 2001 model for quiet conditions
#   ‘T01STORM’ = Tsyganenko 2001 model for active conditions
#   ‘T05’ = Tsyganenko and Sitnov 2005 model
#   ‘ALEX’ = Alexeev model
#   ‘TS07’ = Tsyganenko and Sitnov 2007 model
#
# set opts for irbempy bmodels:
#            - 0 = IGRF13
#            - 1 = Eccentric tilted dipole
#            - 2 = Jensen&Cain 1960
#            - 3 = GSFC 12/66 updated to 1970
#            - 4 = User-defined model (Default: Centred dipole + uniform [Dungey open model] )
#            - 5 = Centred dipole
# see end of script for example call
# -------------------------------------------------------------------------------------

def B_dir(t, x, ray_datenum, extfield, opts, dir_val):
    pos = coord.Coords([x[0], x[1], x[2]], 'GEO', 'car')
    tv = Ticktock(ray_datenum, 'UTC')
    B = irbem.get_Bfield(tv, pos, extMag=extfield, options=opts, omnivals=None)
    Bmags = dir_val * B['Bvec'] / B['Blocal']
    return [Bmags[0][0], Bmags[0][1], Bmags[0][2]]

def trace_fieldline_ODE(bfieldinfo_obj, hemis, crs, carsph, units):

    # set direction
    if hemis == 'north':
        dir_val = 1
    else:
        dir_val = -1

    if bfieldinfo_obj.pos.dtype != 'GEO' or bfieldinfo_obj.pos.units != ['Re','Re','Re'] or bfieldinfo_obj.pos.carsph == 'sph':
        p0 = convert_spc(bfieldinfo_obj.pos, bfieldinfo_obj.time, 'GEO', 'car', ['Re','Re','Re'])
    else:
        p0 = bfieldinfo_obj.pos
    
    x = []
    y = []
    z = []
    dt = 5e-2
    r = ode(B_dir)
    r.set_integrator('vode')

    r.set_initial_value([float(p0.x),float(p0.y),float(p0.z)], 0)
    r.set_f_params(bfieldinfo_obj.time, bfieldinfo_obj.extfield, bfieldinfo_obj.opts, dir_val)
    counts = 0
    while r.successful():
        r.integrate(r.t + dt)
        x.append(r.y[0])
        y.append(r.y[1])
        z.append(r.y[2])

        counts += 1
        if np.linalg.norm(r.y) < 1:
            # hit the earth
            final_coords = [[xx,yy,zz] for xx, yy, zz in zip(x,y,z)]
            dt_array = [bfieldinfo_obj.time for i in range(len(x))]
            fieldlines_spc = create_spc(final_coords, dt_array, 'GEO', 'car', ['Re', 'Re', 'Re'])
            fieldline_con = convert_spc(fieldlines_spc, dt_array, crs, carsph, units)
            bfieldinfo_obj.fieldline = fieldline_con

            break

        if counts > 500:
            print('max count - error')
            break
    return x, y, z


class Bfieldinfo:
    def __init__(self, time=None, pos=None, extfield=None, opts=None, unit_vec=None, bvec=None, fieldline=None, footpoint=None):
        self.time = None
        self.pos = None

        # set based on raytracer settings
        if use_tsyg == 1:
            self.extfield = '0'
        else:
            self.extfield = '0'
        if use_IGRF == 1:
            self.opts = [0,0,0,0,0]
        else: # dipole (assume centered?)
            self.opts = [0,0,0,0,5]

        self.unit_vec = None
        self.bvec = None
        self.fieldline = None 
        self.footpoint = None 

    def Bfield_direction(self, hemis, crs, carsph, units=None):
        if self.time == None:
            print('need time info')
            return
        elif self.pos == None:
            print('need pos info')
            return

        tv = Ticktock(self.time)

        # make sure correct input coordinates
        if self.pos.dtype != 'GEO' or self.pos.units != ['Re','Re','Re'] or self.pos.carsph == 'sph':
            bpos = convert_spc(self.pos, self.time, 'GEO', 'car', ['Re','Re','Re'])
            B = irbem.get_Bfield(tv, bpos, extMag=self.extfield, options=self.opts, omnivals=None)
        else:
            B = irbem.get_Bfield(tv, self.pos, extMag=self.extfield, options=self.opts, omnivals=None)

        if hemis == 'north':
            dir_val = 1
        else:
            dir_val = -1

        B_dir = dir_val * B['Bvec'] / B['Blocal'] # in geo car units are in nT
        B_dir_spc = create_spc(B_dir, self.time, 'GEO', 'car', units=None)
        B_dir_con = convert_spc(B_dir_spc, self.time, crs, carsph, units=None)

        B_vec_spc = create_spc(B['Bvec'], self.time, 'GEO', 'car', units=None)
        B_vec_con = convert_spc(B_vec_spc, self.time, crs, carsph, units=None)

        self.unit_vec = B_dir_con
        self.bvec = B_vec_con

    def findFootpoints(self, hemis, stopalt, crs, carsph, units):

        tv = Ticktock(self.time,'UTC')
        if self.pos.dtype != 'GEO' or self.pos.units != ['Re','Re','Re'] or self.pos.carsph == 'sph':
            bpos = convert_spc(self.pos, self.time, 'GEO', 'car', ['Re','Re','Re'])
            footpoint = irbem.find_footpoint(tv, bpos, extMag=self.extfield, options=self.opts, hemi=hemis, alt=stopalt, omnivals=None)
        else:
            footpoint = irbem.find_footpoint(tv, self.pos, extMag=self.extfield, options=self.opts, hemi=hemis, alt=stopalt, omnivals=None)
        
        # already spc
        gdz = footpoint['loci']
        gdz.ticks = Ticktock(self.time, 'UTC')

        footpoint_coords = convert_spc(gdz, self.time, crs, carsph, units)

        self.footpoint = footpoint_coords
# ---------------------------------------------------------------------

# EXAMPLE CALL

# create a bfield object - how about Lshell3 at the equator?

# let's use today's date
#L3 = Bfieldinfo()
#L3.time = dt.datetime.now().replace(tzinfo=dt.timezone.utc) 

# set the position -- first create a spacepy coord obj
#mypos = create_spc(cor_array=[3,0,0], dt_array=L3.time, crs='GEO', carsph='car', units=['Re','Re','Re'])
#L3.pos = mypos

# find the bfield unit vector here towards the north hemisphere, in GEO car coords
#L3.Bfield_direction(hemis='north', crs='GEO', carsph='car')
#print(L3.unit_vec)

# trace a fieldline to the northern hemisphere
#trace_fieldline_ODE(bfieldinfo_obj=L3, hemis='north', crs='GEO', carsph='car', units=['km','km','km'])
#print(L3.fieldline)

# find the footpoint of the fieldline in the southern hemisphere at 400km
#L3.findFootpoints(hemis='south', stopalt=400, crs='MAG', carsph='sph', units=['km','deg','deg'])
#print(L3.footpoint)