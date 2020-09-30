# reading file from c++

import numpy as np
import matplotlib.pyplot as plt

X = np.genfromtxt("dat/time.csv", delimiter=',',skip_header=1)
#x = np.linspace(0, 1, X.shape[0]) # generating the x-arrays

plt.figure()

plt.plot(X[:,0], X[:,1], 'b', X[:,0], X[:,2], 'y--')
plt.legend(["Jacobi", "Armadillo"])
plt.title("Buckling beam, CPU time comparison")
plt.xlabel("N (dimensionality of matrix)")
plt.ylabel("CPU time")

plt.show()
