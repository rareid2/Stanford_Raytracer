import spacepy.plot.utils
import matplotlib.pyplot as plt
import datetime
times = [datetime.datetime(2010, 1, 1) + datetime.timedelta(hours=i)
    for i in range(0, 48, 3)]
plt.plot(times, range(16))
# [<matplotlib.lines.Line2D object at 0x0000000>]
spacepy.plot.utils.annotate_xaxis(' UT') #mark that times are UT
# <matplotlib.text.Text object at 0x0000000>
