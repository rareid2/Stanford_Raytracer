# python ray tracer driver so that I can read this code

# import statements are a bit diff
import math

# first we need constants -- fortran specifies this as
# double precision, p sure that means Float64, which 
# is python's default for floats

EPS0 = 8.854187817e-12
PI = 3.141592653589793238462643
MU0 = PI * 4e-7
C = math.sqrt(1.0/EPS0/MU0)
R_E = 6371.2e3
D2R = PI / 180.0
R2D = 180.0 / PI
VERSION = 1.17
REkm = R_E*1e-3

#def raytracer_driver

# set a LOT of variables

