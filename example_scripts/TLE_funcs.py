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
month = 5
day = 17
hours = 3
minutes = 0
seconds = 0

ray_datenum = dt.datetime(year, month, day, hours, minutes, seconds)

# DSX TLE
l11 = '1 44344U 19036F   20130.24435661 -.00000027 +00000-0 +00000+0 0  9994'
l21 = '2 44344 042.2583 086.8979 1974641 137.2296 239.9665 04.54371389014496'

# VPM TLE
l12 = '1 45120U 19071K   20132.49935632  .00001453  00000-0  55129-4 0  9998'
l22 = '2 45120  51.6416 181.7127 0011592 280.5137  79.4539 15.33820525015342'

lines1 = [l11, l12]
lines2 = [l21, l22]
satnames = ['DSX', 'VPM']
r, tvec = TLE2posfast(lines1, lines2, satnames, 3*3600, ray_datenum)

rx = [r[0][i][0] for i in range(int(len(r[0])))]
rz = [r[0][i][2] for i in range(int(len(r[0])))]
plt.plot(rx, rz)
plt.show()
print(tvec[0])
"""