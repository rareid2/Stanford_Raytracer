"""
============================================================
Getting orbital data and saving
============================================================
"""

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

secperweek = 604800
minperweek = 10080
hrperweek = 168
time_period = np.linspace(0, 2*hrperweek, (2*hrperweek)+1)

for s in time_period:

    #iterate by min
    datenum_o = datenum + dt.timedelta(minutes=s)

    # DSX TLE:
    line1 = '1 44344U 19036F   20116.38366941 -.00000011 +00000-0 +00000+0 0  9990'
    line2 = '2 44344 042.2529 091.9758 1974961 131.2888 247.5136 04.54371596013862'
    DSX_pos, DSX_t = get_pos(line1, line2, 'DSX', datenum_o)
    x_DSX.append(DSX_pos[0])
    y_DSX.append(DSX_pos[1])
    z_DSX.append(DSX_pos[2])

    # VPM TLE:
    line1 = '1 45120U 19071K   20116.53609353  .00004260  00000-0  13987-3 0  9997'
    line2 = '2 45120  51.6427 258.8751 0011734 219.6893 140.3231 15.33710671012897'
    VPM_pos, VPM_t = get_pos(line1, line2, 'VPM', datenum_o)
    x_VPM.append(VPM_pos[0])
    y_VPM.append(VPM_pos[1])
    z_VPM.append(VPM_pos[2])

    print(s)

print('finished generating data')

# save data to a txtfile
orbitdata = list(zip(x_DSX, y_DSX, z_DSX, x_VPM, y_VPM, z_VPM))
afile = open('orbit_data.txt', 'w')
np.savetxt(afile, orbitdata)
afile.close()



"""
============================================================
3D animation
============================================================
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