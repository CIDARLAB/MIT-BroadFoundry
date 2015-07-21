
#Practice plotting
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111)
x_points = [1,2,3,4,5,6,7,8,7]
y_points = [1,2,3,45,6,7,8,3,7]
p = ax.plot(x_points, y_points, 'ro')
plt.xlim(0,5)
plt.ylim(0,8)
ax.set_xlabel('x-points')
ax.set_ylabel('y-points')
ax.set_title('Simple XY point plot')
fig.show()

##from pylab import *
##t = linspace(0.0, 2.0, 100) # makes a list
##y = cos(2*pi*t) # makes a list with same length as t
##figure() # opens a new figure
##plot(t, y, ls=':', color='r', linewidth=5.0)
##xlabel('time (s)')
##ylabel('voltage (mV)')
##title('A Simple Plot')
##grid()
##show()
