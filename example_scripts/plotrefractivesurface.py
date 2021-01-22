"""
plot refractive surface @ ray starting point!
"""

import numpy as np
import datetime as dt
import os
import sys
import tempfile

from spacepy import coordinates as coord
from spacepy.time import Ticktock
from spacepy import irbempy

import matplotlib.pyplot as plt
import matplotlib
import matplotlib.patches as Patch
import seaborn as sns
sns.set(style="whitegrid")

from raytracer_utils import readdump, read_rayfile, read_rayfiles
from run_rays import run_rays
from constants_settings import *
from IGRF_funcs import B_dir, trace_fieldline_ODE, findFootprints, B_direasy
#from TLE_funcs import TLE2pos
from run_model_dump import modeldump 

# --------------- CONSTANTS --------------------------
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

# ---------------------------------------- STIX PARAM --------------------------------------------
def stix_parameters(ray, t, w):

    B   =  ray['B0'].iloc[t]
    Bmag = np.linalg.norm(B)
    Q    = np.abs(np.array(ray['qs'].iloc[t,:]))
    M    = np.array(ray['ms'].iloc[t,:])
    Ns   = np.array(ray['Ns'].iloc[t,:])

    Wcs   = Q*Bmag/M
    Wps2  = Ns*pow(Q,2)/eo/M

    R = 1.0 - np.sum(Wps2/(w*(w + Wcs)))
    L = 1.0 - np.sum(Wps2/(w*(w - Wcs)))
    P = 1.0 - np.sum(Wps2/(w*w))
    S = (R+L)/2.0
    D = (R-L)/2.0

    return R, L, P, S, D
# ---------------------------------------------------------------------------------------------

# -------------------------------------- FIND L SHELL -----------------------------------------
def getLshell(ray, t, ray_datenum):

    pos = coord.Coords([ray['pos'].iloc[t,:].x, ray['pos'].iloc[t,:].y, ray['pos'].iloc[t,:].z], 'SM', 'car', units = ['m', 'm', 'm']) 
    pos.ticks = Ticktock(ray_datenum, 'UTC') # add ticks
    MAG_pos = pos.convert('MAG', 'sph')
    Lshell = (float(MAG_pos.radi) / R_E) / (pow(np.cos(float(MAG_pos.lati)),2))
    
    return Lshell

# ---------------------------------------------------------------------------------------------