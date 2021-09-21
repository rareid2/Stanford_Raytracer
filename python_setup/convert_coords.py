import numpy as np
import datetime as dt 
import os

from constants import *

# point to coordinate transform library
from libxformd import xflib
cwd = os.getcwd()
xf = xflib.xflib(lib_path=cwd+'/libxformd/libxformd.so')

def convert2(pos_array, dt_array, crs1, carsph1, units1, crs2, carsph2, units2):
    # supports SM, GEI, and GEO cartesian/spherical

    if carsph1 == 'sph': # came in as sph
        units1 = [units1[0],units1[0],units1[0]] # redefine units in car
        converted_coords = []
        for pi, pos in enumerate(pos_array):
            new_coords = xf.s2c(pos)
            converted_coords.append(new_coords)
        pos_array = converted_coords # now in cartesian
    # now all in cartesian, proceed
    converted_coords = []
    if crs1 == 'GEI':
        if crs2 == 'SM':
            for pi, pos in enumerate(pos_array):
                new_coords = xf.gei2sm(pos, dt_array[pi])
                converted_coords.append(new_coords)
        if crs2 == 'GEO':
            for pi, pos in enumerate(pos_array):
                new_coords = xf.gei2geo(pos, dt_array[pi])
                converted_coords.append(new_coords)
    elif crs1 == 'GEO':
        if crs2 == 'SM':
            for pi, pos in enumerate(pos_array):
                new_coords = xf.geo2sm(pos, dt_array[pi])
                converted_coords.append(new_coords)
        if crs2 == 'GEI':
            for pi, pos in enumerate(pos_array):
                new_coords = xf.geo2gei(pos, dt_array[pi])
                converted_coords.append(new_coords)
    elif crs1 == 'SM':
        if crs2 == 'GEO':
            for pi, pos in enumerate(pos_array):
                new_coords = xf.sm2geo(pos, dt_array[pi])
                converted_coords.append(new_coords)
        if crs2 == 'GEI':
            for pi, pos in enumerate(pos_array):
                new_coords = xf.sm2gei(pos, dt_array[pi])
                converted_coords.append(new_coords)
    else:
        pass

    if not converted_coords:
        converted_coords = pos_array # not converting between frames

    if carsph2 == 'sph':
        units2 = [units2[0],units2[0],units2[0]] # redefine units in car

    if units1 != units2:
        if units1 == ['m','m','m']:
            if units2 == ['km','km','km']:
                converted_coords2 = [[cr[0]/1e3,cr[1]/1e3,cr[2]/1e3] for cr in converted_coords]
        if units1 == ['m','m','m']:
            if units2 == ['Re','Re','Re']:
                converted_coords2 = [[cr[0]/R_E,cr[1]/R_E,cr[2]/R_E] for cr in converted_coords]
        if units1 == ['km','km','km']:
            if units2 == ['m','m','m']:
                converted_coords2 = [[cr[0]*1e3,cr[1]*1e3,cr[2]*1e3] for cr in converted_coords]
        if units1 == ['km','km','km']:
            if units2 == ['Re','Re','Re']:
                converted_coords2 = [[cr[0]*1e3/R_E,cr[1]*1e3/R_E,cr[2]*1e3/R_E] for cr in converted_coords]
        if units1 == ['Re','Re','Re']:
            if units2 == ['m','m','m']:
                converted_coords2 = [[cr[0]*R_E,cr[1]*R_E,cr[2]*R_E] for cr in converted_coords]
        if units1 == ['Re','Re','Re']:
            if units2 == ['km','km','km']:
                converted_coords2 = [[cr[0]*R_E/1e3,cr[1]*R_E/1e3,cr[2]*R_E/1e3] for cr in converted_coords]
    else:
        converted_coords2 = converted_coords
    
    if carsph2 == 'sph':
        converted_coords = []
        for pi, pos in enumerate(converted_coords2):
                new_coords = xf.c2s(pos)
                converted_coords.append(new_coords)
        converted_coords2 = converted_coords


    return converted_coords2