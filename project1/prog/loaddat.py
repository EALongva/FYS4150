# reading file from c++

import numpy as np
import matplotlib.pyplot as plt

'''
t = np.linspace(0.0,4.0,)

def y(t):
    return t**2 - 2*t

y = y(t)

plt.plot(t,y)
plt.show()

def loadcsv(filename):
    data = np.loadtxt(filename, delimiter=',',skiprows=1)
    return data
'''

X = np.genfromtxt("data.csv", delimiter=',',skip_header=3)
x = np.linspace(0, 1, X.shape[0])
plt.plot(x, X[:,0], label="num")
plt.plot(x, X[:,1], label="true")
plt.show()

plt.plot(x, X[:,2])
plt.show()
