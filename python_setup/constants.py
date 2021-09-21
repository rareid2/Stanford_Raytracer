# import packages needed to define settings
import numpy as np

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