import numpy as np
import datetime as datetime
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from dateutil import parser
from example_scripts.raytracer_settings import *

# read in orbital data
# comes in as GCRS car coordinate system in m == ECI == GEI car??????
orbitdata = np.genfromtxt('orbit_pos.txt')
dsxdata = orbitdata[:,0:3]

# convert back to time data - time is in UTC
f = open('orbit_time.txt', 'r')
times = f.readlines()
timestr = [t.strip() for t in times]
timedata = [parser.parse(t) for t in timestr]

print('got data')

# get when in apo/peri
apa = []
for i in range(len(dsxdata)):
    dat = dsxdata[i]
    x = dat[0] / R_E
    if x<-1 or x>1:
        apa.append(i)

dsxpositions = []
vpmpositions = []
raytime = []

for ap in apa:
    dsxpositions.append(orbitdata[ap][0:3])
    vpmpositions.append(orbitdata[ap][3:6])
    raytime.append(timedata[ap])

# only do for 3 days (let's do next week) - 150,000 time points
# do every 5 seconds - 10,000 rays

dsxpositions = dsxpositions[300000:450000]
vpmpositions = vpmpositions[300000:450000]
raytime = raytime[300000:450000]

print('got all positions')

looplen = len(dsxpositions)
n = 1500


# get data needed here:
vpmx = np.genfromtxt('vpmx.txt')
vpmy = np.genfromtxt('vpmy.txt')
vpmz = np.genfromtxt('vpmz.txt')

vpm = [np.hstack([x, y, z]) for x,y,z in zip(vpmx, vpmy, vpmz)]

xray = np.genfromtxt('xray.txt')
yray = np.genfromtxt('yray.txt')
zray = np.genfromtxt('zray.txt')

ray = [np.hstack([x, y, z]) for x,y,z in zip(xray, yray, zray)]

footx = np.genfromtxt('footx.txt')
footy = np.genfromtxt('footy.txt')
footz = np.genfromtxt('footz.txt')

foot = [np.hstack([x, y, z]) for x,y,z in zip(footx, footy, footz)]

rayt = [raytime[int(numcount)] for numcount in np.linspace(0,looplen-1,n)]

dist = []
distf = []
# find distances:
for r,v,f in zip(ray, vpm, foot):
    d = np.sqrt((v[0] - r[0]) ** 2 + (v[1] - r[1]) ** 2 + (v[2] - r[2]) ** 2)
    df = np.sqrt((v[0] - f[0]) ** 2 + (v[1] - f[1]) ** 2 + (v[2] - f[2]) ** 2)
    d = d/1e3
    df = df/1e3
    dist.append(d)
    distf.append(df)

# easy
fig, ax = plt.subplots(1, 1)
ax.plot(rayt,dist, label='field-aligned ray')
ax.plot(rayt,distf, label='fieldline footprint')

# formatting
for label in ax.get_xticklabels():
    label.set_rotation(40)
    label.set_horizontalalignment('right')

#ax.set_ylim([0,6000])
ax.set_xlabel('UTC Time')
ax.set_ylabel('Distance in km')
plt.legend()
ax.set_title('Distance from Launched Rays on DSX to VPM  \n for 8.2kHz field-aligned ray')
plt.savefig('distvstime.png')
plt.show()
