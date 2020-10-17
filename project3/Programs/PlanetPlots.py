import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

dim = 2 # Spatial dimension
Nbodies = 2 # Number of celestial bodies

filename = ["out_positions_beta_2.000000", "out_positions_beta_2.250000", "out_positions_beta_2.500000", "out_positions_beta_2.750000", "out_positions_beta_3.000000"]
beta = [2, 2.25, 2.5, 2.75, 3]

plt.figure()
for b in range(0,5):
    raw = pd.read_csv(filename[b],header=None)
    n_steps = int(len(raw)/(dim*Nbodies+1))

    #t = np.zeros(n_steps)
    x_positions = np.zeros((n_steps, Nbodies,5))
    y_positions = np.zeros((n_steps, Nbodies,5))
    #z_positions = np.zeros((n_steps, Nbodies))
    planet_KE = np.zeros((n_steps, Nbodies))
    planet_PE = np.zeros((n_steps, Nbodies))

    for i in range(0, n_steps):
        int_t = i*(dim*Nbodies+1)
        #t[i]= raw.iloc[int_t,0]
        for n in range(0,Nbodies):
            x_positions[i,n,b] = raw.iloc[int_t+1+n*dim, 0]
            y_positions[i,n,b] = raw.iloc[int_t+2+n*dim, 0]


    plt.plot(x_positions[:,0,b], y_positions[:,0,b])#, label="Sun, $beta=%g$" %beta[b])
    plt.plot(x_positions[:,1,b], y_positions[:,1,b], label=r"Earth, $\beta=%g$" %beta[b])

plt.legend(loc='center')
plt.xlabel('x')
plt.ylabel('y')
plt.title(r"The Orbit of the Earth with varying $\beta$-parameter")
plt.show()
