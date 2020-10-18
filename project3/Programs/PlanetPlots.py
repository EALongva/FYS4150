import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

dim = 2 # Spatial dimension
Nbodies = 2 # Number of celestial bodies
Msun = 1.989E+30 #kg

raw = pd.read_csv("out_positions.dat",header=None)
n_steps = int(len(raw)/(1+(dim+3)*Nbodies))

t = np.zeros(n_steps)
x_positions = np.zeros((n_steps, Nbodies))
y_positions = np.zeros((n_steps, Nbodies))
#z_positions = np.zeros((n_steps, Nbodies))
planet_KE = np.zeros((n_steps, Nbodies))
planet_PE = np.zeros((n_steps, Nbodies))
planet_L = np.zeros((n_steps, Nbodies))

for i in range(0, n_steps):
    int_t = i*(1+(dim+3)*Nbodies)
    t[i]= raw.iloc[int_t,0]
    for n in range(0,Nbodies):
        x_positions[i,n] = raw.iloc[int_t+1+n*dim, 0]
        y_positions[i,n] = raw.iloc[int_t+2+n*dim, 0]
        planet_KE[i,n] = raw.iloc[int_t+dim*Nbodies+2+n*dim, 0]
        planet_PE[i,n] = raw.iloc[int_t+dim*Nbodies+3+n*dim, 0]
        planet_L[i,n] = raw.iloc[int_t+dim*Nbodies+4+n*dim, 0]

# planet_L = planet_L * Msun

plt.figure()
plt.plot(x_positions[0:,0], y_positions[0:,0], label="Sun")
plt.plot(x_positions[0:,1], y_positions[0:,1], label="Earth")
plt.xlabel('x')
plt.ylabel('y')
plt.title('The Orbit of the Earth around the Sun')
ax = plt.gca()
ax.set_facecolor('midnightblue')
plt.legend()
plt.show()

plt.figure()
plt.plot(t, planet_KE[:,1], label="Kinetic Energy")
plt.plot(t, planet_PE[:,1], label="Potential Enery")
plt.plot(t, planet_KE[:,1]+planet_PE[:,1], label="Total Energy")
plt.xlabel('Time [years]')
plt.ylabel('Energy')
plt.title('The Mechanical Energy of the Earth')
#ax = plt.gca()
#ax.set_facecolor('midnightblue')
plt.legend()
plt.show()

plt.figure()
plt.plot(t, planet_L[:,1], label="Angular Momentum")
plt.xlabel('Time [years]')
plt.ylabel('kg AU²/year')
#plt.ylabel(r'M$_{sun}$AU²/year')
plt.title('The Angular Momentum of the Earth')
#ax = plt.gca()
#ax.set_facecolor('midnightblue')
plt.legend()
plt.show()
"""
# Varying beta

filename = ["out_positions_beta_2.000000", "out_positions_beta_2.250000", "out_positions_beta_2.500000", "out_positions_beta_2.750000", "out_positions_beta_3.000000"]
beta = [2, 2.25, 2.5, 2.75, 3]

plt.figure()
for b in range(0,5):
    raw = pd.read_csv(filename[b],header=None)
    n_steps = int(len(raw)/(dim*Nbodies+1))

    t = np.zeros(n_steps)
    x_positions = np.zeros((n_steps, Nbodies,5))
    y_positions = np.zeros((n_steps, Nbodies,5))
    #z_positions = np.zeros((n_steps, Nbodies))
    planet_KE = np.zeros((n_steps, Nbodies))
    planet_PE = np.zeros((n_steps, Nbodies))

    for i in range(0, n_steps):
        int_t = i*(dim*Nbodies+1)
        t[i]= raw.iloc[int_t,0]
        for n in range(0,Nbodies):
            x_positions[i,n,b] = raw.iloc[int_t+1+n*dim, 0]
            y_positions[i,n,b] = raw.iloc[int_t+2+n*dim, 0]


    plt.plot(x_positions[:,0,b], y_positions[:,0,b])#, label="Sun, $beta=%g$" %beta[b])
    plt.plot(x_positions[:,1,b], y_positions[:,1,b])#, label="Earth, $beta=%g$" %beta[b])
    #distance = np.sqrt((x_positions[:,1,b] - x_positions[:,0,b])**2 + (y_positions[:,1,b] - y_positions[:,0,b])**2)
    #plt.plot(t, distance, label=r"$\beta=%g$" %beta[b])

plt.legend(loc='center')
#plt.legend(loc='lower left')
plt.xlabel('x')
plt.ylabel('y')
#plt.xlabel('Time [years]')
#plt.ylabel('Distance from Sun to Earth [AU]')
plt.title(r'The Orbit of the Earth with varying $\beta$-parameter')
#plt.title(r"The Earth-Sun distance with varying $\beta$-parameter")
plt.show()
"""
