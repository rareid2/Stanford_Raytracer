"""
plot refractive surface @ ray starting point!
"""

import numpy as np
import datetime as dt
from spacepy import coordinates as coord
from spacepy.time import Ticktock
import matplotlib.pyplot as plt
import matplotlib.patches as Patch
import seaborn as sns
sns.set(style="whitegrid")
from raytracer_utils import readdump, read_rayfile, read_rayfiles

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

# change time information here - use UTC -
year = 2020
month = 5
day = 20
hours = 0
minutes = 0
seconds = 0

ray_datenum = dt.datetime(year, month, day, hours, minutes, seconds)
# --------------- CONSTANTS --------------------------

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
    S = 1.0/(2.0*(R+L))
    D = 1.0/(2.0*(R-L))

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

# Read in a rayfile -- get the plasma density parameters from within
rf = read_rayfile('/var/folders/51/h992wgvj4kld4w4yhw1vx5600000gn/T/tmppv9fxqqu/example_ray_mode1.ray')

# get an entire ray lets just try one for now
ray = rf[0]
w = ray['w']

f = w/(2*np.pi)
print('frequency of ray is:', w/(2*np.pi), ' Hz')

# set up phi vec
phi_vec = np.linspace(0,360,int(1e4))*D2R

# grab only t = 0
t = 0

# for earlier use
Lshell = getLshell(ray, t, ray_datenum)

# get stix param
R, L, P, S, D = stix_parameters(ray, t, w)

root = 1 # why ??

k_vec = np.zeros_like(phi_vec)
eta_vec=np.zeros_like(phi_vec)

for phi_ind, phi  in enumerate(phi_vec):

    # Solve the cold plasma dispersion relation
    cos2phi = pow(np.cos(phi),2)
    sin2phi = pow(np.sin(phi),2)

    A = S*sin2phi + P*cos2phi
    B = R*L*sin2phi + P*S*(1.0+cos2phi)

    discriminant = B*B - 4.0*A*R*L*P

    n1sq = (B + np.sqrt(discriminant))/(2.0*A)
    n2sq = (B - np.sqrt(discriminant))/(2.0*A)

    n1 = np.sqrt(n1sq)
    n2 = np.sqrt(n2sq)
    print(n1, n2)

    # Order the roots
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

# grab for ONLY electrons
Ns = float(ray['Ns'].iloc[t,0])
Q = float(ray['qs'].iloc[t,0])
M = float(ray['ms'].iloc[t,0])
B   =  ray['B0'].iloc[t]
Bmag = np.linalg.norm(B)

# makes it easier to square here
w2 = w*w
wp2 = Ns*pow(Q,2)/eo/M
wh = Q*Bmag/M
wh2 = wh*wh

root = 1

# solve appleton-hartree eq
numerator = wp2/w2
denom1 = (wh2*pow(np.sin(phi_vec),2))/(2*(w2 - wp2))
denom2 = np.sqrt(pow(denom1,2) + wh2*pow(np.cos(phi_vec), 2)/w2)
eta2_AH   = 1 - (numerator/(1 - denom1 + root*denom2))
eta_AH = np.sqrt(-eta2_AH)

# plot it 
fig, ax = plt.subplots(1,1)

# ax.plot(eta_AH*np.sin(phi_vec), eta_AH*np.cos(phi_vec), LineWidth = 1, label = 'e only')
ax.plot(eta_vec*np.sin(phi_vec), eta_vec*np.cos(phi_vec), LineWidth = 1, label = 'e + ions')

findcone = eta_vec*np.sin(phi_vec)
for cone in findcone:
    if cone > 100:
        wherecone = np.where(findcone == cone)
        conetheta = phi_vec[wherecone]
        conetheta = conetheta * R2D
        break

archeight = float(eta_vec[wherecone])

# formatting
xlim1 = -100
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

resonanceangle = float(90 - conetheta)

pac = Patch.Arc([0, 0], archeight, archeight, angle=0, theta1=0, theta2=float(resonanceangle), edgecolor = 'r')
#ax.add_patch(pac)

ax.annotate('${\Theta}$ < ' + str(round(float(conetheta), 2)) + 'deg', xy=(xlim2 - (scale*200), scale*50))

ax.set_title(str(round(f/1e3, 2)) + ' kHz Refractive Surface at ' + str(ray_datenum))
plt.legend()
plt.savefig('/Users/rileyannereid/Desktop/refractivesurface' + str(round(f/1e3, 2)) + 'kHz.png', format='png')
plt.show()
