# reading file from c++

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#importing data
SunOrbit = np.genfromtxt("data/ESJ/SunOrbit.csv", delimiter=',')
SunOrbitM1 = np.genfromtxt("data/ESJ/SunOrbitM1.csv", delimiter=',')
SunOrbitM10 = np.genfromtxt("data/ESJ/SunOrbitM10.csv", delimiter=',')
SunOrbitM1000 = np.genfromtxt("data/ESJ/SunOrbitM1000.csv", delimiter=',')

EarthOrbit = np.genfromtxt("data/ESJ/EarthOrbit.csv", delimiter=',')
EarthOrbitM1 = np.genfromtxt("data/ESJ/EarthOrbitM1.csv", delimiter=',')
EarthOrbitM10 = np.genfromtxt("data/ESJ/EarthOrbitM10.csv", delimiter=',')
EarthOrbitM1000 = np.genfromtxt("data/ESJ/EarthOrbitM1000.csv", delimiter=',')

JupOrbitM1 = np.genfromtxt("data/ESJ/JupOrbitM1.csv", delimiter=',')
JupOrbitM10 = np.genfromtxt("data/ESJ/JupOrbitM10.csv", delimiter=',')
JupOrbitM1000 = np.genfromtxt("data/ESJ/JupOrbitM1000.csv", delimiter=',')

#calculating earth - sun distance
r0 = np.linalg.norm(EarthOrbit - SunOrbit, axis=0)
r1 = np.linalg.norm(EarthOrbitM1 - SunOrbitM1, axis=0)
r2 = np.linalg.norm(EarthOrbitM10 - SunOrbitM10, axis=0)
r3 = np.linalg.norm(EarthOrbitM1000 - SunOrbitM1000, axis=0)

t = np.linspace(0,3,len(r0)) # have to update each time T is changed

#SORRY FOR NOT MAKING LOOPS D:<
plt.figure()
plt.subplot(211)
plt.title("Orbital Trajectories (Sun, Earth and Jupiter), Varying Mj (Relative Jupiter mass), N=10^4")
plt.plot(SunOrbit[0],SunOrbit[1],'y')
plt.plot(SunOrbitM1[0],SunOrbitM1[1],'y',alpha=0.9)
plt.plot(SunOrbitM10[0],SunOrbitM10[1],'y',alpha=0.7)
plt.plot(SunOrbitM1000[0],SunOrbitM1000[1],'y',alpha=0.5)
plt.plot(EarthOrbit[0],EarthOrbit[1],'b')
plt.plot(EarthOrbitM1[0],EarthOrbitM1[1],'b',alpha=0.9)
plt.plot(EarthOrbitM10[0],EarthOrbitM10[1],'b',alpha=0.7)
plt.plot(EarthOrbitM1000[0],EarthOrbitM1000[1],'b',alpha=0.5)
plt.plot(JupOrbitM1[0],JupOrbitM1[1],'r',alpha=0.9)
plt.plot(JupOrbitM10[0],JupOrbitM10[1],'r',alpha=0.7)
plt.plot(JupOrbitM1000[0],JupOrbitM1000[1],'r',alpha=0.5)

plt.xlabel("y position [AU]")
plt.xlabel("x position [AU]")

plt.legend(["Sun Mj=0", "Sun Mj=1", "Sun Mj=10", "Sun Mj=1000", \
            "Earth Mj=0", "Earth Mj=1", "Earth Mj=10", "Earth Mj=1000",\
            "Jupiter Mj=1", "Jupiter Mj=10", "Jupiter Mj=1000"])

plt.subplot(212)

plt.title("Earth Distance From Sun, Varying Jupiter Mass")

plt.plot(t,r0)
plt.plot(t,r1)
plt.plot(t,r2)
plt.plot(t,r3)

plt.legend(["Earth - Sun system", "Earth - Sun - Jupiter system", "ESJ, Mj=10",\
            "ESJ, Mj=1000"])

plt.xlabel("time [yr]")
plt.ylabel("distance from sun [AU]")

plt.show()



"""
csv = ".csv"
names = np.array(["Sun","Earth","Mercury", "Venus", "Mars", "Jupiter", \
                  "Saturn", "Uranus", "Neptune", "Pluto"])

# 3d plot
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
"""
