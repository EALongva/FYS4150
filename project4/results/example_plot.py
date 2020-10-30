# example plot taking data from results/data and saving a figure to
# results/figures

import numpy as np
import matplotlib.pyplot as plt

X = np.genfromtxt("results/data/example.csv", delimiter=',')
x = np.linspace(0,1,np.size(X))

plt.title("example plot")
plt.plot(x,X)
plt.savefig("results/figures/example.png",dpi=300)
