# reading file from c++

import numpy as np
import matplotlib.pyplot as plt

X = np.genfromtxt("dat/eigenvalues.csv", delimiter=',',skip_header=1)
x = np.linspace(1, X.shape[0], X.shape[0]) # generating the x-arrays

plt.figure()

plt.subplot(211)
plt.plot(x, X[:,0], 'g', x, X[:,1], 'r--')

plt.subplot(212)
plt.plot(x, X[:,2], 'g', x, X[:,3], 'r--')

plt.show()
