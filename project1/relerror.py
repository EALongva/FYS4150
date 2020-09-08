# reading file from c++, relative error and creating plot

import numpy as np
import matplotlib.pyplot as plt

X = np.genfromtxt("relError.csv", delimiter=',',skip_header=1)
#x = np.linspace(0, 1, X.shape[0]) # generating the x-arrays

# comparing numerical and exact solutions
plt.plot(X[:,0], X[:,1], '-g')
plt.ylabel('log10(max relative error)')
plt.xlabel('n, number of steps [10^n]')
plt.title('maximum relative error (specialized algo) for decreasing step size')
plt.show()
