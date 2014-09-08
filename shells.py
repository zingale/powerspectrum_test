#!/bin/env python

import numpy as np
import matplotlib.pyplot as plt


# draw a grid
N = 9

xmin = ymin = 0.0
xmax = ymax = 1.0

x = np.arange(N)*(xmax-xmin)/(N-1) + xmin
y = np.arange(N)*(ymax-ymin)/(N-1) + ymin

xc = (x[0:N-1] + x[1:N])/2
yc = (y[0:N-1] + y[1:N])/2

for n in range(N):
    plt.plot([xmin, xmax], [y[n], y[n]], color="k")
    plt.plot([x[n], x[n]], [ymin, ymax], color="k")



# now draw shells
print x

dr = (2.0*(x[1]-x[0]))

R = np.arange(0.0, 1.0+dr, dr)

npts = 200
theta = np.arange(npts)*(np.pi/2)/npts

for r in R:
    if r == 0: continue
    plt.plot(r*np.cos(theta), r*np.sin(theta), color="b")


ax = plt.gca()
ax.set_aspect("equal")    


plt.savefig("shells.png")
