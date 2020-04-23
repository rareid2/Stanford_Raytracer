"""

here is a script to define ALL settings for the ray tracer and models
to keep everything
central

"""

# import packages needed to define settings
import numpy as np               # for math
import os                        # for running commands
import datetime as dt            # for coordinate transforms

# Constants
D2R = (np.pi / 180.0)
R2D = 1.0 / D2R
R_E = 6371e3  # m
H_IONO = 1000e3

#  --------------------------- CHANGE ENV SETTINGS HERE  --------------------------
# Environmental parameters

# change time information here - use UTC -
year = 2020
month = 4
day = 6
hours = 22
minutes = 4
seconds = 30

# convert for raytracer settings
days_in_the_year = (dt.date(year, month, day) - dt.date(year,1,1)).days + 1
days_in_the_year = format(days_in_the_year, '03d')

# yearday and miliseconds day are used by raytracer
yearday = str(year)+ str(days_in_the_year)   # YYYYDDD
milliseconds_day = hours*3.6e6 + minutes*6e4 + seconds*1e3

# used for plotting the ray
ray_datenum = dt.datetime(int(yearday[0:4]), 1, 1) + dt.timedelta(int(yearday[4:]) - 1)

# space weather settings
Kp = 2
AE = 1.6
Pdyn = 4
Dst = 1.0
ByIMF = 0.0
BzIMF = -5
# Tsyganenko correction parameters
W = [0.132, 0.303, 0.083, 0.070, 0.211, 0.308]  # Doesn't matter if we're not using Tsyg

# Simulation parameters
t_max = 10       # Maximum duration in seconds
dt0 = 1e-3       # Initial timestep in seconds
dtmax = 0.1      # Maximum allowable timestep in seconds
root = 2         # Which root of the Appleton-Hartree equation
                 # (1 = negative, 2 = positive)
                 # (2=whistler in magnetosphere)
fixedstep = 0    # Don't use fixed step sizes, that's a bad idea.
maxerr = 5.0e-4  # Error bound for adaptive timestepping
maxsteps = 2e5   # Max number of timesteps (abort if reached)
use_IGRF = 1     # Magnetic field model (1 for IGRF, 0 for dipole)
use_tsyg = 1     # Use the Tsyganenko magnetic field model corrections
minalt = R_E     # cutoff altitude in meters
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

#  --------------------------- END CHANGE ENV SETTINGS HERE  --------------------------


#  ---------------------------- SET UP FILE SETTINGS  ------------------------------
project_root = os.getcwd()  # grabs current full path

# Set input file path
ray_inpfile = os.path.join(project_root, "ray_inpfile.txt")

# Set config file for Ngo Model (Mode 1)
configfile = os.path.join(project_root, "newray_default.in")

# Set config file for mode 3
mode3_interpfile = os.path.join(project_root, 'precomputed_grids',
                                            'gcpm_kp4_2001001_L10_80x80x80_noderiv.txt')
# Set config file for mode 4
mode4_modelfile = os.path.join(project_root,
                               'precomputed_grids', 'precomputed_model_gcpm_2010001_0_kp2_L10_random.dat')

# Set up the output directory
ray_out_dir = os.path.join(project_root, "test_outputs")

# Create directory for outputs if doesn't already exist
# print("output directory:", ray_out_dir)
if not os.path.exists(ray_out_dir):
    os.mkdir(ray_out_dir)

#  ---------------------------- END SET UP FILE SETTINGS  ------------------------------

