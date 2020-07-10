# This example plots "gappy" data with and without explicit fill values.

import matplotlib.pyplot as plt
import numpy
import spacepy.datamanager
x = numpy.append(numpy.arange(0, 6, 0.1), numpy.arange(12, 18, 0.1))
y = numpy.sin(x)
xf, yf = spacepy.datamanager.insert_fill(x, y)
fig = plt.figure()
ax0 = fig.add_subplot(211)
ax0.plot(x, y)
ax1 = fig.add_subplot(212)
ax1.plot(xf, yf)
plt.show()
