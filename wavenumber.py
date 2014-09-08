#!/bin/env python

import numpy as np
import matplotlib.pyplot as plt

plt.subplot(211)

x = np.linspace(0.0, 1.0, 200)

k = 1.0

plt.plot(x, np.sin(2.0*np.pi*k*x))

xlabel = [0, 0.5, 1]
xlabelnames = [r"$0$", r"$L/2$", r"$L$"]

ylabel = [-1, 1]
ylabelnames = [r"$-1$", r"$1$"]


ax = plt.gca()
ax.spines['left'].set_position('zero')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_position('zero')
ax.spines['top'].set_color('none')
ax.spines['left'].set_smart_bounds(True)
ax.spines['bottom'].set_smart_bounds(True)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

plt.xticks(xlabel, xlabelnames)
plt.yticks(ylabel, ylabelnames)

plt.ylim(-1.1,1.1)



plt.subplot(212)

x = np.linspace(0.0, 1.0, 200)

k = 4.0

plt.plot(x, np.sin(2.0*np.pi*k*x))

ax = plt.gca()
ax.spines['left'].set_position('zero')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_position('zero')
ax.spines['top'].set_color('none')
ax.spines['left'].set_smart_bounds(True)
ax.spines['bottom'].set_smart_bounds(True)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

plt.xticks(xlabel, xlabelnames)
plt.yticks(ylabel, ylabelnames)

plt.ylim(-1.1,1.1)

plt.tight_layout()

plt.savefig("wavenumber.eps")


