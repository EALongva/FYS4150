import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14})
inputpath = "../results/output/"
figpath = "../results/figures/"

dim = 2 # Spatial dimension
Nbodies = 2 # Number of celestial bodies

beta = [2, 2.25, 2.5, 2.75, 3]
name = ["2", "2.25", "2.5", "2.75", "3"]
# Initial velocity in data
v_fname = "7" # filename format
v = "7" # latex format

raw = pd.read_csv(inputpath+"ES_beta2_v"+v_fname"_t2N10e-5.dat", header=None)
n_steps = int(len(raw)/(1+(dim+3)*Nbodies))

t = np.zeros(n_steps)
x_positions = np.zeros((n_steps, Nbodies,len(beta)))
y_positions = np.zeros((n_steps, Nbodies,len(beta)))
planet_KE = np.zeros((n_steps, Nbodies, len(beta)))
planet_PE = np.zeros((n_steps, Nbodies, len(beta)))
planet_L = np.zeros((n_steps, Nbodies, len(beta)))

for b in range(0,len(beta)):
    raw = pd.read_csv(inputpath+"ES_beta"+name[b]+"_v"+v_fname"_t2N10e-5.dat",header=None)

    for i in range(0, n_steps):
        int_t = i*(1+(dim+3)*Nbodies)
        t[i]= raw.iloc[int_t,0]
        for n in range(0,Nbodies):
            x_positions[i,n,b] = raw.iloc[int_t+1+n*dim, 0]
            y_positions[i,n,b] = raw.iloc[int_t+2+n*dim, 0]
            planet_KE[i,n,b] = raw.iloc[int_t+dim*Nbodies+2+n*dim, 0]
            planet_PE[i,n,b] = raw.iloc[int_t+dim*Nbodies+3+n*dim, 0]
            planet_L[i,n,b] = raw.iloc[int_t+dim*Nbodies+4+n*dim, 0]

# Plotting orbits

plt.figure()
for b in range(0,len(beta)):
    #plt.plot(x_positions[:,0,b], y_positions[:,0,b], label="Sun, $beta=%g$" %beta[b])
    plt.plot(x_positions[:,1,b], y_positions[:,1,b], label=r"$\beta=%g$" %beta[b])
ax = plt.gca()
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
# Put a legend below current axis
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), fancybox=True, shadow=True, ncol=5)
plt.xlabel('x [AU]')
plt.ylabel('y [AU]')
plt.title(r'The Orbit of the Earth, $v_i=%s$ AU/year' %v)
plt.savefig(figpath+"ES_v"+v_fname+"_orbit_varyingbeta.png",bbox_inches = 'tight')

plt.figure()
for b in range(0,len(beta)):
    distance = np.sqrt((x_positions[:,1,b] - x_positions[:,0,b])**2 + (y_positions[:,1,b] - y_positions[:,0,b])**2)
    plt.plot(t, distance, label=r"$\beta=%g$" %beta[b])
plt.legend(loc='best')
plt.xlabel('Time [years]')
plt.ylabel('Distance from Sun to Earth [AU]')
plt.title(r"Earth-Sun distance, $v_i=%s$ AU/year" %v)
plt.savefig(figpath+"ES_v"+v_fname+"_distance_varyingbeta.png",bbox_inches = 'tight')

# Plotting mechanical energy

plt.figure()
plt.plot(t, planet_KE[:,1,0], label=r"Kinetic, $\beta=2$")
plt.plot(t, planet_PE[:,1,0], label=r"Potential, $\beta=2$")
plt.plot(t, planet_KE[:,1,0]+planet_PE[:,1,0], label=r"Total, $\beta=2$")
plt.xlabel('Time [years]')
plt.ylabel(r'Energy [$M_{\odot}$AU$^2$/year$^2$]')
plt.title(r"Earth's Mechanical Energy, $v_i=%s$ AU/year" %v)
plt.legend()
plt.savefig(figpath+"ES_v"+v_fname+"_energy.png",bbox_inches = 'tight')

plt.figure()
plt.plot(t, planet_KE[:,1,4], label=r"Kinetic, $\beta=3$")
plt.plot(t, planet_PE[:,1,4], label=r"Potential, $\beta=3$")
plt.plot(t, planet_KE[:,1,4]+planet_PE[:,1,4], label=r"Total Energy, $\beta=3$")
plt.xlabel('Time [years]')
plt.ylabel(r'Energy [$M_{\odot}$AU$^2$/year$^2$]')
plt.title(r"Earth's Mechanical Energy, $v_i=%s$ AU/year" %v)
plt.legend()
plt.savefig(figpath+"ES_v"+v_fname+"_beta3_energy.png",bbox_inches = 'tight')

# Plotting angular momentum

plt.figure()
plt.plot(t, planet_L[:,1,0], label=r"$\beta=2$")
plt.xlabel('Time [years]')
plt.ylabel('Specific Angular Momentum [AU²/year]')
plt.title(r"Earth's Angular Momentum, $v_i=%s$ AU/year" %v)
plt.legend()
plt.savefig(figpath+"ES_v"+v_fname+"_angmom.png",bbox_inches = 'tight')

plt.figure()
plt.plot(t, planet_L[:,1,4], label=r"$\beta=3$")
plt.xlabel('Time [years]')
plt.ylabel('Specific Angular Momentum [AU²/year]')
plt.title(r"Earth's Angular Momentum, $v_i=%s$ AU/year" %v)
plt.legend()
plt.savefig(figpath+"ES_v"+v_fname+"_beta3_angmom.png",bbox_inches = 'tight')

plt.show()
