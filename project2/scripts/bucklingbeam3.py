# reading file from c++

import numpy as np
import matplotlib.pyplot as plt

X = np.genfromtxt("dat/analyticalEigComparisonBB.csv", delimiter=',',skip_header=1)
x = np.linspace(1, X.shape[0], X.shape[0]) # generating the x-arrays
rho = np.linspace(0,1,X.shape[0])

plt.figure()

plt.subplot(211)
plt.plot(x, X[:,0], 'g', x, X[:,1], 'r--')
plt.legend(["Analytic", "Armadillo"])
plt.title("Buckling beam, eigenvalues, rhoN=1")
plt.xlabel("i")
plt.ylabel("eigenvalue i")


plt.subplot(212)
plt.plot(rho, X[:,2], 'g', rho, X[:,3], 'r--')
plt.legend(["Analytic", "Armadillo"])
plt.title("Eigenvector of smallest eigenvalue")
plt.xlabel("rho [m/L]")
plt.ylabel("u(rho)")

plt.show()
