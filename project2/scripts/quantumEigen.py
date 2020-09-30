# reading file from c++

import numpy as np
import matplotlib.pyplot as plt

X = np.genfromtxt("dat/QMeigen_Stddev.csv", delimiter=',',skip_header=1)
N = np.linspace(X[0,0], X[-1,0], X.shape[0]) # generating the x-arrays
#x = np.linspace(0.,40.,X[-1,0])

plt.figure()

plt.errorbar(N,X[:,2],yerr=X[:,1],capsize=7,fmt="o")
plt.xscale("log")
plt.yscale("log")
plt.title("Mean deviation from analytical eigenvalues, Harmonic Osc, rhoN = 40")
plt.ylabel("Mean deviation")
plt.xlabel("N")

plt.show()
