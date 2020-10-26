# reading file from c++

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

dt = np.genfromtxt("test/test_time_dt.csv", delimiter=',')
eult = np.genfromtxt("test/test_time_eult.csv", delimiter=',')
vvt = np.genfromtxt("test/test_time_vvt.csv", delimiter=',')

plt.figure()

plt.subplot(221)
plt.title("CPU Time Usage, Forward Euler, 10 Runs Average")
plt.plot(dt,eult,'bo')
plt.plot(dt,eult,'k--',alpha=0.8)
plt.xlabel("timestep dt [yr]")
plt.ylabel("CPU time [s]")
plt.xscale("log")
plt.grid()

plt.subplot(222)
plt.title("CPU Time Usage, Velocity Verlet, 10 Runs Average")
plt.plot(dt,vvt,'bo')
plt.plot(dt,vvt,'k--',alpha=0.8)
plt.xlabel("timestep dt [yr]")
plt.ylabel("CPU time [s]")
plt.xscale("log")
plt.grid()

plt.subplot(223)
#plt.title("CPU Time Usage, Forward Euler, 10 Runs Average")
plt.plot(dt,eult,'bo')
plt.plot(dt,eult,'k--',alpha=0.8)
plt.xlabel("timestep dt [yr]")
plt.ylabel("CPU time [s]")
plt.yscale("log")
plt.xscale("log")
plt.grid()

plt.subplot(224)
#plt.title("CPU Time Usage, Velocity Verlet, 10 Runs Average")
plt.plot(dt,vvt,'bo')
plt.plot(dt,vvt,'k--',alpha=0.8)
plt.xlabel("timestep dt [yr]")
plt.ylabel("CPU time [s]")
plt.yscale("log")
plt.xscale("log")
plt.grid()

plt.show()
