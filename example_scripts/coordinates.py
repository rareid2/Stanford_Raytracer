from spacepy import coordinates as coord
from spacepy.time import Ticktock
import datetime as dt

from constants_settings import *

# coordinate system transformations
# options: 

# crs:
# GDZ (Geodetic; WGS84),
# GEO (Geographic Coordinate System),
# GSM (Geocentric Solar Magnetospheric),
# GSE (Geocentric Solar Ecliptic),
# SM (Solar Magnetic),
# GEI (Geocentric Equatorial Inertial; True-of-Date),
# MAG (Geomagnetic Coordinate System),
# SPH (Spherical Coordinate System),
# RLL (Radius, Latitude, Longitude; Geodetic)

# carsph:
# car or sph 

# units
# Re, km, m 
# (spherical defaults to degrees)

# create a spacepy coordinate object
def create_spc(cor_array, dt_array, crs, carsph, units):
    cvals = coord.Coords(cor_array, crs, carsph, units=units)
    cvals.ticks = Ticktock(dt_array, 'UTC')
    return cvals

# convert between objects!
def convert_spc(cvals, dt_array, crs, carsph, units):

    # convert coords here
    newcoord = cvals.convert(crs, carsph)
    if units == None:
        return newcoord

    # there's a bug in spacepy coordinates
    # spacepy defaulted to Re despite me asking it not to
    if carsph == 'car':
        if newcoord.units[0] != 'Re' and np.abs(newcoord.x[0]) < 10:
            newcoord.units = ['Re', 'Re', 'Re']
    else:
        if newcoord.units[0] != 'Re' and np.abs(newcoord.radi[0]) < 10:
            newcoord.units = ['Re', 'deg', 'deg']

    if newcoord.units[0] == units[0]:
        conversionf = 1.0
    elif newcoord.units[0] == 'km' and units[0] == 'm':
        conversionf = 1.0e3
    elif newcoord.units[0] == 'm' and units[0] == 'km':
        conversionf = 1.0e-3
    elif newcoord.units[0] == 'Re' and units[0] == 'km':
        conversionf = R_E * 1.0e-3
    elif newcoord.units[0] == 'Re' and units[0] == 'm':
        conversionf = R_E 
    elif newcoord.units[0] == 'm' and units[0] == 'Re':
        conversionf = 1.0 / R_E
    elif newcoord.units[0] == 'km' and units[0] == 'Re':
        conversionf = 1.0 / (R_E * 1.0e-3)
    else:
        print('unit conversion not yet implemented')

    # convert units here
    if conversionf != 1:
        if carsph == 'car':
            newx = newcoord.x * conversionf
            newy = newcoord.y * conversionf
            newz = newcoord.z * conversionf
            cor_array = [[nx,ny,nz] for nx,ny,nz in zip(newx,newy,newz)]
            converted_coords = create_spc(cor_array, dt_array, crs, carsph, units)
        
        if carsph == 'sph':
            newr = newcoord.radi *conversionf
            cor_array = [[nr,ny,nz] for nr,ny,nz in zip(newr,newcoord.lati,newcoord.long)]
            converted_coords = create_spc(cor_array, dt_array, crs, carsph, units)
    else:
        converted_coords = newcoord
    
    return converted_coords