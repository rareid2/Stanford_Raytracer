"""
TLE to position functions to plot nicely with ray tracer
get TLE's from http://celestrak.com/satcat/search.php
DSX NORAD ID: 44344
VPM NORAD ID: 45120
"""

# import packages
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
from spacepy.coordinates import Coords
from spacepy.time import Ticktock
from sgp4.api import Satrec, SatrecArray, jday
import julian

"""

FUNCTION TLE2pos

INPUTS:  line 1 and 2 = list of TLE lines strings
         satnames = list of satnames strings
         plen = length to propagate the orbit in SECONDS AS AN INT
         FYI - one orbit is about 5 hours for DSX and about 1.5 hours for VPM
OUTPUTS: orbital position in ECI cart km and time vector corresponding to orbit (UTC) 


"""

def TLE2pos(lines1, lines2, satnames, plen, ray_datenum):

    # fast function uses Julian time (weird)
    jd, fr = jday(ray_datenum.year, ray_datenum.month, ray_datenum.day, ray_datenum.hour, ray_datenum.minute,
                ray_datenum.second)
    np.set_printoptions(precision=2)

    frac = 1/86400

    # create time arrays in desired format
    frarray = np.arange(fr, fr + float(plen*frac), frac)
    jdarray = jd * np.ones(len(frarray))

    # set up dict
    sats = {}
    i = 0

    # get satellite objects from TLEs
    for name in satnames:
        line1 = lines1[i]
        line2 = lines2[i]
        sats[name] = Satrec.twoline2rv(line1, line2)
        i +=1
    
    # get positions over time
    satlist = list(sats.values())
    a = SatrecArray(satlist)
    e, r, v = a.sgp4(jdarray, frarray)

    tvec = []
    # convert back to datetime
    for f, j in zip(frarray, jdarray):
        myjd = j + f
        dti = julian.from_jd(myjd, fmt='jd')
        tvec.append(dti)

    return r, tvec # returned in km from Earth center in ECI coordinates

# --------------------------------- END FUNCTION -------------------------------------


"""
# example call

# change time information here - use UTC -
year = 2020
month = 6
day = 11
hours = 21
minutes = 58
seconds = 0

ray_datenum = dt.datetime(year, month, day, hours, minutes, seconds)

# last updated 6/21

# DSX TLE
l11 = '1 44344U 19036F   20173.14565688 -.00000031  00000-0  00000-0 0  9999'
l21 = '2 44344  42.2760  71.1855 1973524 155.6114 215.1832  4.54371095 16448'
# VPM TLE
l12 = '1 45120U 19071K   20173.93473231  .00003239  00000-0  10800-3 0  9994'
l22 = '2 45120  51.6437 341.3758 0012446  71.4995 288.7339 15.34053724 21707'


lines1 = [l11, l12]
lines2 = [l21, l22]
satnames = ['DSX', 'VPM']

r, tvec = TLE2pos(lines1, lines2, satnames, 6*3600, ray_datenum)

# convert to meters
dsx = [rpos*1e3 for rpos in r[0]]
vpm = [rpos*1e3 for rpos in r[1]]

rx_d = [r[0][i][0] for i in range(int(len(r[0])))]
rz_d = [r[0][i][2] for i in range(int(len(r[0])))]

rx_v = [r[1][i][0] for i in range(int(len(r[1])))]
rz_v = [r[1][i][2] for i in range(int(len(r[1])))]

plt.plot(rx_d, rz_d)
plt.plot(rx_v, rz_v)
plt.show()
#print(tvec[0])
"""