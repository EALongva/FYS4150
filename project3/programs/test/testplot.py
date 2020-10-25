# reading file from c++

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

csv = ".csv"
n = 10
T = 5.0
dt = T/np.floor(np.logspace(2,5,num=10))


""" 3d plot
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
ax.plot(earth[0,:], earth[1,:], earth[2,:], 'g')
plt.title("Test plot of earth orbit around the sun (heliocentric)")
plt.xlabel("x [AU]")
plt.ylabel("y [AU]")
"""
plt.figure()

plt.subplot(321)
plt.title("Euler's Forward Method")
for i in range(n):
    folder = "test/eulerdat/distanceN"
    filename = folder + str(i) + csv
    r = np.genfromtxt(filename, delimiter=',')
    plt.plot(np.linspace(0,1,np.size(r)),r)

plt.xlabel("time [yr]")
plt.ylabel("distance to sun [AU]")

plt.subplot(323)
for i in range(n):
    folder = "test/eulerdat/orbitN"
    filename = folder + str(i) + csv
    orbit = np.genfromtxt(filename, delimiter=',')
    plt.plot(orbit[0],orbit[1])

plt.xlabel("x [AU]")
plt.ylabel("y [AU]")

plt.subplot(322)
plt.title("Velocity Verlet method")
for i in range(n):
    folder = "test/vvdat/distanceN"
    filename = folder + str(i) + csv
    r = np.genfromtxt(filename, delimiter=',')
    plt.plot(np.linspace(0,1,np.size(r)),r)

plt.xlabel("time [yr]")
plt.ylabel("distance to sun [AU]")

plt.subplot(324)
for i in range(n):
    folder = "test/vvdat/orbitN"
    filename = folder + str(i) + csv
    orbit = np.genfromtxt(filename, delimiter=',')
    plt.plot(orbit[0],orbit[1])

plt.xlabel("x [AU]")
plt.ylabel("y [AU]")

plt.subplot(325) #euler

dr = np.zeros(n)

for i in range(n):
    folder = "test/eulerdat/distanceN"
    filename = folder + str(i) + csv
    r = np.genfromtxt(filename, delimiter=',')
    dr[i] = np.mean(r)


plt.plot(dt,dr,'bo')
plt.xscale("log")
plt.xlabel("timestep [yr]")
plt.ylabel("average distance to sun [AU]")

plt.subplot(326) #verlet

drvv = np.zeros(n)

for i in range(n):
    folder = "test/vvdat/distanceN"
    filename = folder + str(i) + csv
    r = np.genfromtxt(filename, delimiter=',')
    drvv[i] = np.mean(r)


plt.plot(dt,drvv,'bo')
plt.xscale("log")
plt.xlabel("timestep [yr]")
plt.ylabel("average distance to sun [AU]")






plt.show()
