# import packages needed to define settings
import numpy as np
import os
import sys
import json

# Constants (matching the Fortran raytracer)
EPS0 = 8.854187817e-12        # C^2/Nm^2
MU0 = np.pi * 4e-7
C = np.sqrt(1.0/(EPS0*MU0))   # m/s
R_E = 6371.2e3                # m
D2R = np.pi / 180.0           # rad
R2D = 180.0 / np.pi           # deg
H_IONO = 1000e3               # m
Hz2Rad = 2.*np.pi             # rad
Rad2Hz = 1./Hz2Rad            # Hz
Q_EL = 1.602e-19              # C
M_EL = 9.1e-31                # kg
M_P = 1.67e-27                # kg

file_location = os.path.dirname(__file__)
raytracer_dir = os.path.join(file_location, '../bin')

# Opening JSON file
fpath = sys.argv[1]
# returns JSON object as a dictionary
with open(fpath) as f:
    settings = json.load(f)


minalt = settings["minalt"]+R_E
Kp = settings['Kp']
AE = settings['AE']
Pdyn = settings['Pdyn']
Dst = settings['Dst']
ByIMF = settings['ByIMF']
BzIMF = settings['BzIMF']
W = settings['W']
use_tsyg = settings['use_tsyg']
use_IGRF = settings['use_IGRF']

t_max = settings['t_max']
dt0 = settings['dt0']
dtmax = settings['dtmax']
root = settings['root']
fixedstep = settings['fixedstep']
maxerr = settings['maxerr']
maxsteps = settings['maxsteps']
outputper = settings['outputper']
MLT = settings['MLT']

# the name of the json config file
# runName = test_rays if file is test_rays.json
runName = __file__.split('/')[-1].split('.')[0]
