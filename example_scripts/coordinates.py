from spacepy import coordinates as coord
from spacepy.time import Ticktock
import datetime as dt

# coordinate systems

def create_spc(cor_array, dt_array, all_units, crs, carsph):
    cvals = coord.Coords(cor_array, crs, carsph, units=all_units)
    cvals.ticks = Ticktock(dt_array, 'UTC')
    return cvals


def convert_spc(cvals, crs, carsph):
    newcoord = cvals.convert(crs, carsph)
    return newcoord