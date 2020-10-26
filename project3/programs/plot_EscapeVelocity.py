import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14})
inputpath = "../results/output/"
figpath = "../results/figures/"

dim = 2 # Spatial dimension
Nbodies = 2 # Number of celestial bodies

escv = [0.9, 0.99, 1.0, 1.01, 1.1]
name = ["0.9", "0.99", "1.0", "1.01", "1.1"]

raw = pd.read_csv(inputpath+"ES_escv0.9_t10N10e-5.dat", header=None)
n_steps = int(len(raw)/(1+(dim+3)*Nbodies))

t = np.zeros(n_steps)
x_positions = np.zeros((n_steps, Nbodies,len(escv)))
y_positions = np.zeros((n_steps, Nbodies,len(escv)))
planet_KE = np.zeros((n_steps, Nbodies, len(escv)))
planet_PE = np.zeros((n_steps, Nbodies, len(escv)))
planet_L = np.zeros((n_steps, Nbodies, len(escv)))

for v in range(0,len(escv)):
    raw = pd.read_csv(inputpath+"ES_escv"+name[v]+"_t10N10e-5.dat",header=None)

    for i in range(0, n_steps):
        int_t = i*(1+(dim+3)*Nbodies)
        t[i]= raw.iloc[int_t,0]
        for n in range(0,Nbodies):
            x_positions[i,n,v] = raw.iloc[int_t+1+n*dim, 0]
            y_positions[i,n,v] = raw.iloc[int_t+2+n*dim, 0]
            planet_KE[i,n,v] = raw.iloc[int_t+dim*Nbodies+2+n*dim, 0]
            planet_PE[i,n,v] = raw.iloc[int_t+dim*Nbodies+3+n*dim, 0]
            planet_L[i,n,v] = raw.iloc[int_t+dim*Nbodies+4+n*dim, 0]

# Plotting orbits

plt.figure()
for v in range(0,len(escv)):
    #plt.plot(x_positions[:,0,b], y_positions[:,0,b], label="Sun, $beta=%g$" %beta[b])
    plt.plot(x_positions[:,1,v], y_positions[:,1,v], label=r"$v_i=%gv_e$" %escv[v])
#ax = plt.gca()
#box = ax.get_position()
#ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
# Put a legend below current axis
#ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), fancybox=True, shadow=True, ncol=5)
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlabel('x [AU]')
plt.ylabel('y [AU]')
plt.title(r'The Orbit of the Earth')
plt.savefig(figpath+"ES_escv_orbit.png",bbox_inches = 'tight')

plt.figure()
for v in range(0,len(escv)):
    distance = np.sqrt((x_positions[:,1,v] - x_positions[:,0,v])**2 + (y_positions[:,1,v] - y_positions[:,0,v])**2)
    plt.plot(t, distance, label=r"$v_i=%gv_e$" %escv[v])
plt.legend(loc='best')
plt.xlabel('Time [years]')
plt.ylabel('Distance from Sun to Earth [AU]')
plt.title(r"Earth-Sun distance")
plt.savefig(figpath+"ES_escv_distance.png",bbox_inches = 'tight')

plt.show()
"""
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
"""
"""
# Plotting orbits

plt.figure()
for b in range(0,len(beta)):
    #plt.plot(x_positions[:,0,b], y_positions[:,0,b], label="Sun, $beta=%g$" %beta[b])
    plt.plot(x_positions[:13200,1,b], y_positions[:13200,1,b], label=r"$\beta=%g$" %beta[b])
plt.legend(loc='best')
plt.xlabel('x')
plt.ylabel('y')
plt.title(r'The Orbit of the Earth')
plt.savefig(figpath+"earthsun_orbit_varyingbeta_zoom.png",bbox_inches = 'tight')

plt.figure()
for b in range(0,len(beta)):
    distance = np.sqrt((x_positions[:,1,b] - x_positions[:,0,b])**2 + (y_positions[:,1,b] - y_positions[:,0,b])**2)
    plt.plot(t[:13200], distance[:13200], label=r"$\beta=%g$" %beta[b])
plt.legend(loc='best')
plt.xlabel('Time [years]')
plt.ylabel('Distance from Sun to Earth [AU]')
plt.title(r"Earth-Sun distance")
plt.savefig(figpath+"earthsun_distance_varyingbeta_zoom.png",bbox_inches = 'tight')

# Plotting mechanical energy

plt.figure()
plt.plot(t[:13200], planet_KE[:13200,1,0], label=r"Kinetic, $\beta=2$")
plt.plot(t[:13200], planet_PE[:13200,1,0], label=r"Potential, $\beta=2$")
plt.plot(t[:13200], planet_KE[:13200,1,0]+planet_PE[:13200,1,0], label=r"Total, $\beta=2$")
plt.xlabel('Time [years]')
plt.ylabel(r'Energy [$M_{\odot}$AU$^2$/year$^2$]')
plt.title(r'The Mechanical Energy of the Earth')
plt.legend()
plt.savefig(figpath+"earthsun_energy_zoom.png",bbox_inches = 'tight')

plt.figure()
plt.plot(t[:13200], planet_KE[:13200,1,4], label=r"Kinetic, $\beta=3$")
plt.plot(t[:13200], planet_PE[:13200,1,4], label=r"Potential, $\beta=3$")
plt.plot(t[:13200], planet_KE[:13200,1,4]+planet_PE[:13200,1,4], label=r"Total Energy, $\beta=3$")
plt.xlabel('Time [years]')
plt.ylabel(r'Energy [$M_{\odot}$AU$^2$/year$^2$]')
plt.title(r'The Mechanical Energy of the Earth')
plt.legend()
plt.savefig(figpath+"earthsun_beta3_energy_zoom.png",bbox_inches = 'tight')

# Plotting angular momentum

plt.figure()
plt.plot(t[:13200], planet_L[:13200,1,0], label=r"$\beta=2$")
plt.xlabel('Time [years]')
plt.ylabel('Specific Angular Momentum [AU²/year]')
plt.title(r'The Specific Angular Momentum of the Earth')
plt.legend()
plt.savefig(figpath+"earthsun_angmom_zoom.png",bbox_inches = 'tight')

plt.figure()
plt.plot(t[:13200], planet_L[:13200,1,4], label=r"$\beta=3$")
plt.xlabel('Time [years]')
plt.ylabel('Specific Angular Momentum [AU²/year]')
plt.title(r'The Specific Angular Momentum of the Earth')
plt.legend()
plt.savefig(figpath+"earthsun_angmom_beta3_zoom.png",bbox_inches = 'tight')

plt.show()
"""
