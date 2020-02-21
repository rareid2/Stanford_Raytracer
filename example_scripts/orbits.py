def orbits(n):
    # n is the array size of time from 0 to orbit period

    # where are our lil spacecraft friends
    from PyAstronomy import pyasl
    import numpy as np
    from spacepy import coordinates as coord
    from spacepy.time import Ticktock
    from datetime import datetime

    R_E = 6371e3         # [m]
    mu = 3.986004418e14  # [m^3/s^2] gravitational parameter for Earth
    dsx_a = 9e6 + R_E    # [m] semi major axis from center of Earth
    vpm_a = 5e5 + R_E    # [m] semi major axis from center of Earth

    dsx_period = 2 * np.pi * np.sqrt(dsx_a ** 3 / mu)
    vpm_period = 2 * np.pi * np.sqrt(vpm_a ** 3 / mu)

    dsx = pyasl.KeplerEllipse(dsx_a, dsx_period, e=((12e6 + R_E) / dsx_a) - 1, i=120)
    vpm = pyasl.KeplerEllipse(vpm_a, vpm_period, e=0, i=51.6)

    t = np.linspace(0, dsx_period, n)  # Get a time vector

    dsx_pos = dsx.xyzPos(t)
    vpm_pos = vpm.xyzPos(t)

    dsx_x = dsx_pos[::, 0]
    dsx_y = dsx_pos[::, 1]
    dsx_z = dsx_pos[::, 2]
    vpm_x = vpm_pos[::, 0]
    vpm_y = vpm_pos[::, 1]
    vpm_z = vpm_pos[::, 2]

    dtime = datetime(2010, 1, 1, 0, 0, 0)  # date for 2010

    dsxorbit = np.zeros((t.size, 3))
    vpmorbit = np.zeros((t.size, 3))

    for i in range(0, t.size):
        dsxpos = coord.Coords([[dsx_x[i], dsx_y[i], dsx_z[i]]], 'GEO', 'car')
        dsxpos.ticks = Ticktock([dtime])  # add ticks
        dsxpos = dsxpos.convert('SM', 'car')
        dsxorbit[i, 0] = dsxpos.x
        dsxorbit[i, 1] = dsxpos.y
        dsxorbit[i, 2] = dsxpos.z

        vpmpos = coord.Coords([[vpm_x[i], vpm_y[i], vpm_z[i]]], 'GEO', 'car')
        vpmpos.ticks = Ticktock([dtime])  # add ticks
        vpmpos = vpmpos.convert('SM', 'car')
        vpmorbit[i, 0] = vpmpos.x
        vpmorbit[i, 1] = vpmpos.y
        vpmorbit[i, 2] = vpmpos.z

    return dsxorbit, vpmorbit
