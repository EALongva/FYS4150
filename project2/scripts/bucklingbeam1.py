# reading file from c++

import numpy as np
import matplotlib.pyplot as plt

X = np.genfromtxt("dat/SimTransCount.csv", delimiter=',',skip_header=1)
#x = np.linspace(0, 1, X.shape[0]) # generating the x-arrays

plt.figure()

plt.subplot(211)
plt.plot(X[1:,0], X[1:,1], 'g', X[1:,0], X[1:,0]**2, 'r--')
plt.xscale("log")
plt.yscale("log")
plt.title("Buckling beam, similarity transformations")
plt.xlabel("N (dimensionality of matrix)")
plt.ylabel("Similarity transformations")

plt.subplot(212)
plt.plot(X[1:,0], X[1:,2], 'b')
plt.title("CPU time as function of dimensionality")
plt.xlabel("N (dimensionality of matrix)")
plt.ylabel("CPU time")

plt.show()
