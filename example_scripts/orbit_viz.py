"""
============================================================
Getting orbital data and saving
============================================================
"""

# don't need this script anymore - consider moving to scratches

import numpy as np
import datetime as dt
from get_TLE import get_pos
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation

R_E = 6371e3  # m

# change time information here - use UTC -
# start date
year = 2020
month = 4
day = 29
hours = 0
minutes = 0
seconds = 0
datenum = dt.datetime(year, month, day, hours, minutes, seconds)

# initialize
x_DSX = []
y_DSX = []
z_DSX = []
x_VPM = []
y_VPM = []
z_VPM = []

# let's get for the next 2weeks
secperweek = 604800
secperday = 86400
minperweek = 10080
hrperweek = 168
time_period = np.linspace(0, int(secperweek*2), int(secperweek*2)+1)
datenum_o = [datenum + dt.timedelta(seconds=s) for s in time_period]

s = 0
for daten in datenum_o:
    # DSX TLE:
    line1 = '1 44344U 19036F   20117.92375283 -.00000009 +00000-0 +00000-0 0  9998'
    line2 = '2 44344 042.2535 091.4116 1974932 131.9492 246.6867 04.54371554013931'
    DSX_pos, DSX_t = get_pos(line1, line2, 'DSX', daten)
    x_DSX.append(DSX_pos[0])
    y_DSX.append(DSX_pos[1])
    z_DSX.append(DSX_pos[2])

    # VPM TLE:
    line1 = '1 45120U 19071K   20119.07726718  .00004309  00000-0  14125-3 0  9993'
    line2 = '2 45120  51.6431 246.5927 0011612 229.4528 130.5444 15.33735479013284'
    VPM_pos, VPM_t = get_pos(line1, line2, 'VPM', daten)
    x_VPM.append(VPM_pos[0])
    y_VPM.append(VPM_pos[1])
    z_VPM.append(VPM_pos[2])

    s += 1
    print(s)

print('finished generating data')

# save data to a txtfile
orbitdata = list(zip(x_DSX, y_DSX, z_DSX, x_VPM, y_VPM, z_VPM))
afile = open('orbit_pos.txt', 'w')
np.savetxt(afile, orbitdata)
afile.close()

# save time to a txtfile
afile = open('orbit_time.txt', 'w')
np.savetxt(afile, datenum_o)
afile.close()


"""
============================================================
3D animation
============================================================
"""
"""

# get data
data = np.genfromtxt('orbit_data.txt')

# Attaching 3D axis to the figure
fig = plt.figure()
ax = p3.Axes3D(fig)

# Setting the fig properties
ax.set_xlim3d([-3, 3])
ax.set_xlabel('X')
ax.set_ylim3d([-3, 3])
ax.set_ylabel('Y')
ax.set_zlim3d([-3, 3])
ax.set_zlabel('Z')
figname = 'orbits_' + str(datenum)
ax.set_title(figname)

line, = ax.plot([], [], [], '-bo', lw=2, zorder = 101)
line2, = ax.plot([], [], [], '-ro', lw=2, zorder = 102)

# add the earth
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
earthx = np.outer(np.cos(u), np.sin(v))
earthy = np.outer(np.sin(u), np.sin(v))
earthz = np.outer(np.ones(np.size(u)), np.cos(v))
ax.plot_surface(earthx, earthy, earthz, color='g', alpha = 0.5, zorder = 100)

# add in orbit trajectories
plt.plot(data[:,0] / R_E, data[:,1] / R_E, data[:,2] / R_E, 'y')
plt.plot(data[:,3] / R_E, data[:,4] / R_E, data[:,5] / R_E, 'y')

# not sure if this is needed
def init():
    line.set_data([], [])
    line2.set_data([], [])
    return line, line2,

def animate(i):
    x = data[i,0] / R_E
    y = data[i,1] / R_E
    z = data[i,2] / R_E
    line.set_data(x, y)
    line.set_3d_properties(z)

    x2 = data[i,3] / R_E
    y2 = data[i,4] / R_E
    z2 = data[i,5] / R_E

    line2.set_data(x2, y2)
    line2.set_3d_properties(z2)

    print(i)
    return line, line2,

#animate it
anim = animation.FuncAnimation(fig, animate, init_func=init, frames=int(len(data)), interval=2, blit=True)

#name and save it
framerate = 300 # fps
savename = 'orbits_' + str(datenum) + '.mp4'
FFwriter=animation.FFMpegWriter(fps=int(framerate), extra_args=['-vcodec', 'libx264'])
anim.save(savename, writer=FFwriter)

print('we done')

"""
"""
############################################################
# figure out when at apogee and perigee over the next 2 weeks
# get data
data = np.genfromtxt('orbit_data.txt')
dsxdata = data[:,0:3]
apa = []

for i in range(len(dsxdata)):
    dat = dsxdata[i]
    x = dat[0] / R_E
    if x<-1 or x>1:
        apa.append(i)

new_orbitdata = [data[ap] for ap in apa]

# save data to a txtfile
afile = open('new_orbit_data.txt', 'w')
np.savetxt(afile, new_orbitdata)
afile.close()

# save apa to a txtfile
afile = open('apa.txt', 'w')
np.savetxt(afile, apa)
afile.close()
"""


"""
#distances = [np.sqrt(data[i,0]**2 + data[i,1]**2 + data[i,2]**2) for i in range(len(data))]
#for dist in distances:
#    if dist < 12350540.7 + 500e3 or dist > 18444281.95 - 500e3:
#        apa.append(distances.index(dist))

new_data = [data[i,0:3] for i in apa]
xray = [dat[0]/R_E for dat in new_data]
yray = [dat[1]/R_E for dat in new_data]
zray = [dat[2]/R_E for dat in new_data]

#apa_new = []
xray = []
yray = []
zray = []

for i in range(len(new_data)):
    dat = new_data[i]
    x = dat[0] / R_E
    if x<-1 or x>1:
        xray.append(dat[0])
        yray.append(dat[1])
        zray.append(dat[2])
        apa_new.append(i)



plt.xlim([-3,3])
plt.ylim([-3,3])
plt.scatter(xray, zray)

savename = 'plots/rayUGHS.png'
plt.savefig(savename, format='png')
"""
