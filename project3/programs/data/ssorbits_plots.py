# reading file from c++

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

csv = ".csv"
names = np.array(["Sun","Earth","Mercury", "Venus", "Mars", "Jupiter", \
                  "Saturn", "Uranus", "Neptune", "Pluto"])

"""# 3d plot
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
for i in range(10):
    folder = "data/solarsystem/planet_"
    filename = folder + names[i] + csv
    orbit = np.genfromtxt(filename, delimiter=',')
    ax.plot(orbit[0],orbit[1],orbit[2])
plt.title("Test plot of earth orbit around the sun (heliocentric)")
plt.xlabel("x [AU]")
plt.ylabel("y [AU]")
"""
# 2d plot
for i in range(10):
    folder = "data/solarsystem/planet_"
    filename = folder + names[i] + csv
    orbit = np.genfromtxt(filename, delimiter=',')
    plt.plot(orbit[0],orbit[1])

plt.title("Solar System Orbits, Velocity Verlet, T=250years")
plt.legend(names)
plt.xlabel("x [AU]")
plt.ylabel("y [AU]")



plt.show()
