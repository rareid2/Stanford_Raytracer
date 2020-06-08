"""
here is a script to define ALL settings for the ray tracer and models
to keep everything
central
"""

# import packages needed to define settings
import numpy as np               # for math
import os                        # for running commands
import datetime as dt            # for coordinate transforms
import tempfile
import random, string

# Constants
D2R = (np.pi / 180.0)
R2D = 1.0 / D2R
R_E = 6371e3  # m
H_IONO = 1000e3

#  --------------------------- CHANGE ENV SETTINGS HERE  --------------------------
# Environmental parameters
Kp = 2
AE = 1.6
Pdyn = 4
Dst = 1.0
ByIMF = 0.0
BzIMF = -5
# Tsyganenko correction parameters
W = [0.132, 0.303, 0.083, 0.070, 0.211, 0.308]  # Doesn't matter if we're not using Tsyg

# Simulation parameters
t_max = 30       # Maximum duration in seconds
dt0 = 1e-3       # Initial timestep in seconds
dtmax = 0.1      # Maximum allowable timestep in seconds
root = 2         # Which root of the Appleton-Hartree equation
                 # (1 = negative, 2 = positive)
                 # (2=whistler in magnetosphere)
fixedstep = 0    # Don't use fixed step sizes, that's a bad idea.
maxerr = 5.0e-4  # Error bound for adaptive timestepping
maxsteps = 2e3   # Max number of timesteps (abort if reached)
use_IGRF = 1     # Magnetic field model (1 for IGRF, 0 for dipole)
use_tsyg = 0     # Use the Tsyganenko magnetic field model corrections
minalt = R_E + 400e3   # cutoff altitude in meters
# TODO: make a max cutoff alt or check if inside the plasmasphere

# Which plasmasphere models should we run?
#   1 - Legacy (Ngo) model
#   2 - GCPM (Accurate, but * s l o w * )
#   3 - Uniformly interpolated precomputed model
#   4 - Randomly interpolated precomputed model
#   5 - (not real)
#   6 - Simplified GCPM from Austin Sousa's thesis

modes_to_do = [1]

# Should we include a geometric focusing term in the damping?
include_geom_factor = 1  # 1 for yes

# interpolation parameters for mode 4
scattered_interp_window_scale = 1.2
scattered_interp_order = 2
scattered_interp_exact = 0  # Try 0 if there's weirdness at discontinuities
scattered_interp_local_window_scale = 5

# change this someday
configfile = 'ngoconfig.in'

#  --------------------------- END CHANGE ENV SETTINGS HERE  --------------------------