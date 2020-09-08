# reading file from c++

import numpy as np
import matplotlib.pyplot as plt

# importing data for n=10
X1 = np.genfromtxt("genAlgo1.csv", delimiter=',',skip_header=3)
t1 = np.genfromtxt("genAlgo1.csv", delimiter=',',skip_header=1,max_rows=1)
x1 = np.linspace(0, 1, X1.shape[0]) # generating the x-arrays

# importing data for n=10
X2 = np.genfromtxt("genAlgo2.csv", delimiter=',',skip_header=3)
t2 = np.genfromtxt("genAlgo2.csv", delimiter=',',skip_header=1,max_rows=1)
x2 = np.linspace(0, 1, X2.shape[0]) # generating the x-arrays

# importing data for n=10
X3 = np.genfromtxt("genAlgo3.csv", delimiter=',',skip_header=3)
t3 = np.genfromtxt("genAlgo3.csv", delimiter=',',skip_header=1,max_rows=1)
x3 = np.linspace(0, 1, X3.shape[0]) # generating the x-arrays

# comparing numerical and exact solutions
plt.figure()

plt.subplot(311)
plt.plot(x1, X1[:,0], 'g', label="numerical")
plt.plot(x1, X1[:,1], label="exact")
plt.title('Solution for n=10, CPU time=' + str(t1[1]))
plt.grid(True)

plt.subplot(312)
plt.plot(x2, X2[:,0], label="numerical")
plt.plot(x2, X2[:,1], label="exact")
plt.title('Solution for n=100, CPU time=' + str(t2[1]))
plt.grid(True)

plt.subplot(313)
plt.plot(x3, X3[:,0], label="numerical")
plt.plot(x3, X3[:,1], label="exact")
plt.title('Solution for n=1000, CPU time=' + str(t3[1]))
plt.grid(True)

plt.show()
