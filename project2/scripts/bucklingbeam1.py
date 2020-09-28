# reading file from c++

import numpy as np
import matplotlib.pyplot as plt

X = np.genfromtxt("dat/SimTransCount.csv", delimiter=',',skip_header=1)
#x = np.linspace(0, 1, X.shape[0]) # generating the x-arrays

plt.figure()

plt.subplot(211)
plt.plot(X[1:,0], X[1:,1], 'g', X[1:,0], X[1:,0]**2.1, 'r--')

plt.subplot(212)
plt.plot(X[1:,0], X[1:,2], 'g')

plt.show()
