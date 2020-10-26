# reading file from c++

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

csv = ".csv"
n = 10
T = 5.0
N = np.floor(np.logspace(2,4,num=2))
G = 4*np.pi**2
M = (6.*10**24)/(2.*10**30)

""" 3d plot
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
ax.plot(earth[0,:], earth[1,:], earth[2,:], 'g')
plt.title("Test plot of earth orbit around the sun (heliocentric)")
plt.xlabel("x [AU]")
plt.ylabel("y [AU]")
"""
plt.figure()
plt.rc("font",size=14)

plt.subplot(211)
plt.title("Conservation of Energy, Forward Euler and Velocity Verlet, N=10^2 and N=10^4")

folder = "nrg/eulPosN"
filename = folder + str(0) + csv
r = np.genfromtxt(filename, delimiter=',')
r = np.linalg.norm(r,axis=0)
P = -G*M/r

folder = "nrg/eulVelN"
filename = folder + str(0) + csv
v = np.genfromtxt(filename, delimiter=',')
v = np.linalg.norm(v,axis=0)
K = 0.5*M*v**2

folder = "nrg/vvPosN"
filename = folder + str(0) + csv
vv_r = np.genfromtxt(filename, delimiter=',')
vv_r = np.linalg.norm(vv_r,axis=0)
vvP = -G*M/vv_r

folder = "nrg/vvVelN"
filename = folder + str(0) + csv
vvv = np.genfromtxt(filename, delimiter=',')
vvv = np.linalg.norm(vvv,axis=0)
vvK = 0.5*M*vvv**2

plt.plot(np.linspace(0,T,np.size(r)),P)
plt.plot(np.linspace(0,T,np.size(r)),K)
plt.plot(np.linspace(0,T,np.size(r)),vvP)
plt.plot(np.linspace(0,T,np.size(r)),vvK)
plt.plot(np.linspace(0,T,np.size(r)),P+K,"y--")
plt.plot(np.linspace(0,T,np.size(r)),vvP+vvK, "k--")
plt.legend(["P_euler","K_euler","P_vv", "K_vv", "Etot_euler","Etot_vv"])
plt.grid()
plt.xlabel("time [yr]")
plt.ylabel("Energy [M_sun][AU^2/yr^2]")


"""plot2"""

plt.subplot(212)
folder = "nrg/eulPosN"
filename = folder + str(1) + csv
r = np.genfromtxt(filename, delimiter=',')
r = np.linalg.norm(r,axis=0)
P = -G*M/r

folder = "nrg/eulVelN"
filename = folder + str(1) + csv
v = np.genfromtxt(filename, delimiter=',')
v = np.linalg.norm(v,axis=0)
K = 0.5*M*v**2

folder = "nrg/vvPosN"
filename = folder + str(1) + csv
vv_r = np.genfromtxt(filename, delimiter=',')
vv_r = np.linalg.norm(vv_r,axis=0)
vvP = -G*M/vv_r

folder = "nrg/vvVelN"
filename = folder + str(1) + csv
vvv = np.genfromtxt(filename, delimiter=',')
vvv = np.linalg.norm(vvv,axis=0)
vvK = 0.5*M*vvv**2

plt.plot(np.linspace(0,T,np.size(r)),P)
plt.plot(np.linspace(0,T,np.size(r)),K)
plt.plot(np.linspace(0,T,np.size(r)),vvP)
plt.plot(np.linspace(0,T,np.size(r)),vvK)
plt.plot(np.linspace(0,T,np.size(r)),P+K,"y--")
plt.plot(np.linspace(0,T,np.size(r)),vvP+vvK,"k--")
plt.legend(["P_euler","K_euler","P_vv", "K_vv", "Etot_euler","Etot_vv"])
plt.grid()
plt.xlabel("time [yr]")
plt.ylabel("Energy [M_sun][AU^2/yr^2]")


plt.show()
