# reading file from c++

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

earth = np.genfromtxt("data/EarthSun/earth.csv", delimiter=',')
#x = np.linspace(0, 1, X.shape[0]) # generating the x-arrays

fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
ax.plot(earth[0,:], earth[1,:], earth[2,:], 'g')
plt.title("Test plot of earth orbit around the sun (heliocentric)")
plt.xlabel("x [AU]")
plt.ylabel("y [AU]")

plt.show()
