import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import requests
import configparser
import json
import time
from sgp4.earth_gravity import wgs84
from sgp4.io import twoline2rv

from coordinates import create_spc, convert_spc

# -------------------------------- SAT CLASS -----------------------------------------
# just a nice way to mangage the satellites (VPM, DSX...)

# provide catnmbr (44344 or 45120 for DSX and VPM respectively)
# provide time in UTC

# return TLE, orbit (provide desired length to propagate the orbit and direction of propagation)
# and get velocity interesting for doppler

# see end of script for example call
# -------------------------------------------------------------------------------------

class sat:
    def __init__(self, TLE=None, time=None, pos=None, vel=None):
        self.catnmbr = None # catnmbr is NOARD ID
        self.TLE = None
        self.time = None
        self.pos = None
        self.vel = None

    def getTLE_ephem(self):
        if self.catnmbr == None:
            print('need catnmbr first')
            return
        if self.time == None:
            print('need time first')
            return

        class MyError(Exception):
            def __init___(self,args):
                Exception.__init__(self,"my exception was raised with arguments {0}".format(args))
                self.args = args

        uriBase                = "https://www.space-track.org"
        requestLogin           = "/ajaxauth/login"
        requestCmdAction       = "/basicspacedata/query" 
        requestFindSat         = "/class/tle/NORAD_CAT_ID/"+str(self.catnmbr)+"/"

        # Use configparser package to pull in the ini file (pip install configparser)
        config = configparser.ConfigParser()
        config.read("./SLTrack.ini")
        configUsr = config.get("configuration","username")
        configPwd = config.get("configuration","password")
        configOut = config.get("configuration","output")
        siteCred = {'identity': configUsr, 'password': configPwd}

        # use requests package to drive the RESTful session with space-track.org
        with requests.Session() as session:
            # run the session in a with block to force session to close if we exit

            # need to log in first. note that we get a 200 to say the web site got the data, not that we are logged in
            resp = session.post(uriBase + requestLogin, data = siteCred)
            if resp.status_code != 200:
                raise MyError(resp, "POST fail on login")

            # this query picks up all Starlink satellites from the catalog. Note - a 401 failure shows you have bad credentials 
            resp = session.get(uriBase + requestCmdAction + requestFindSat)
            if resp.status_code != 200:
                raise MyError(resp, "GET fail on request for satellite")

            # use the json package to break the json formatted response text into a Python structure (a list of dictionaries)
            retData = json.loads(resp.text)
            epoch_dt = [] # save epoch tiemdelta
            for n in range(len(retData)):
                tle_time = dt.datetime.strptime(retData[n]['EPOCH'],'%Y-%m-%d %H:%M:%S').replace(tzinfo=dt.timezone.utc)
                tle_us = dt.timedelta(microseconds=float(retData[n]['EPOCH_MICROSECONDS']))
                tle_fulltime = tle_time + tle_us
                epoch_dt.append(np.abs((self.time - tle_fulltime).total_seconds()))
            
            epoch_n = epoch_dt.index(min(epoch_dt))
            # get the data!
            self.TLE = [retData[epoch_n]['TLE_LINE1'], retData[n]['TLE_LINE2']]

        session.close()
        print("retrieved TLE from SpaceTrak for " + retData[epoch_n]['EPOCH'])
    
    def propagatefromTLE(self, sec, orbit_dir, crs, carsph, units):
        if self.TLE == None:
            print('need TLE')
            return

        if orbit_dir == 'future':
            dt_array = [self.time + dt.timedelta(seconds=s) for s in range(0,sec+1)]
        elif orbit_dir == 'past':
            dt_array = [self.time + dt.timedelta(seconds=s) for s in range(-sec-1,0)]
        else:
            dt_array = [self.time + dt.timedelta(seconds=s) for s in range(-sec-1,sec+1)] 

        # store pos and vel data
        pos_array = np.zeros((len(dt_array),3))
        vel_array = np.zeros((len(dt_array),3))

        for ti, tt in enumerate(dt_array):
            satellite = twoline2rv(self.TLE[0], self.TLE[1], wgs84)
            position, velocity = satellite.propagate(tt.year,tt.month,tt.day,tt.hour,tt.minute,tt.second)
            pos_array[ti,:] = position
            vel_array[ti,:] = velocity

        # finally, convert!
        pos_spc = create_spc(pos_array, dt_array, 'GEI', 'car', ['km','km','km'])
        #vel_spc = create_spc(vel_array, dt_array, 'GEI', 'car', ['km','km','km'])
    
        self.pos = convert_spc(pos_spc, dt_array, crs, carsph, units)
        #self.vel = convert_spc(vel_spc, dt_array, crs, carsph, units)

# ---------------------------------------------------------------------

# Example Call! 
# define an obj called dsx
#dsx = sat()

# set the NORAD ID for dsx
#dsx.catnmbr = 44344

# set the time we want orbit info
#dsx.time = dt.datetime.now().replace(tzinfo=dt.timezone.utc)

# get TLE
#dsx.getTLE_ephem()

# propagate for 10 seconds
#psec = 60.0

# propagate into the future
# alternatively, use 'past' for previous
# or 'both' for propagating in both directions in time
#pdir = 'past'

# select desired coordinate system
#crs = 'GEO'
#carsph = 'sph'
#units = ['m','m','m']

# propagate
#dsx.propagatefromTLE(psec, pdir, crs, carsph, units)

# look @ updated position and velocity vectors
#print(dsx.pos)
#dsx.vel