import numpy as np
import matplotlib.pyplot as plt

N = 100000 # Integration points
FinalTime = 10.0 # Years
GM_sun = 4*np.pi**2 # [AU³/yr²] Gravitational Constant * Mass of the Sun
dt = FinalTime/N
t = np.linspace(0,FinalTime,N)

M_sun = 2e+30 # kg
# Masses relative to solar mass (https://en.wikipedia.org/wiki/Planetary_mass#cite_note-maia.usno.navy.mil-15)
m_earth = 3.04e-06
m_jupiter = 954.79e-6

# Distance to Sun in AU
x_earth = 1.
x_jupiter = 5.2

# Initial positions
r_sun_0 = np.array([0,0])
r_earth_0 = np.array([x_earth,0])
r_jupiter_0 = np.array([x_jupiter,0])

# Initial velocities
v_sun_0 = np.array([0,0])
v_earth_0 = np.array([0,2*np.pi])
v_jupiter_0 = np.array([0,2*np.pi*x_jupiter/12])

# Setting up vectors for storing position and velocity
r_sun = np.zeros((N,2)); r_sun[0,:] = r_sun_0
r_earth = np.zeros((N,2)); r_earth[0,:] = r_earth_0
r_jupiter = np.zeros((N,2)); r_jupiter[0,:] = r_jupiter_0
v_sun = np.zeros((N,2)); v_sun[0,:] = v_sun_0
v_earth = np.zeros((N,2)); v_earth[0,:] = v_earth_0
v_jupiter = np.zeros((N,2)); v_jupiter[0,:] = v_jupiter_0

"""
# Euler Forward
for i in range(1,N):
    r_E_J = r_earth[i-1,:] - r_jupiter[i-1,:]
    distance_E_J = np.linalg.norm(r_E_J)
    a_earth = - GM_sun/x_earth**3*r_earth[i-1,:] + GM_sun*m_jupiter/distance_E_J**3*r_E_J
    a_jupiter = - GM_sun/x_jupiter**3*r_jupiter[i-1,:] - GM_sun*m_earth/distance_E_J**3*r_E_J
    r_earth[i,:] = r_earth[i-1,:] + dt*v_earth[i-1,:]
    v_earth[i,:] = v_earth[i-1,:] + dt*a_earth
    r_jupiter[i,:] = r_jupiter[i-1,:] + dt*v_jupiter[i-1,:]
    v_jupiter[i,:] = v_jupiter[i-1,:] + dt*a_jupiter

plt.plot(r_earth[:,0],r_earth[:,1],label='Earth')
plt.plot(r_jupiter[:,0], r_jupiter[:,1],label='Jupiter')
plt.scatter(r_sun[:,0],r_sun[:,0],label='Sun')
plt.legend()
plt.show()
"""

# Velocity Verlet
r_E_J = r_jupiter[0,:] - r_earth[0,:]
distance_E_J = np.linalg.norm(r_E_J)
a_earth_new = - GM_sun/x_earth**3*r_earth[0,:] - GM_sun*m_jupiter/distance_E_J**3*r_E_J
a_jupiter_new = - GM_sun/x_jupiter**3*r_jupiter[0,:] - GM_sun*m_earth/distance_E_J**3*(-r_E_J)

for i in range(1,N):
    a_earth = a_earth_new
    a_jupiter = a_jupiter_new
    r_earth[i,:] = r_earth[i-1,:] + dt*v_earth[i-1,:] + 0.5*dt**2*a_earth
    r_jupiter[i,:] = r_jupiter[i-1,:] + dt*v_jupiter[i-1,:] + 0.5*dt**2*a_jupiter
    r_E_J = r_earth[i,:] - r_jupiter[i,:]
    distance_E_J = np.linalg.norm(r_E_J)
    a_earth_new = - GM_sun/np.linalg.norm(r_earth[i,:])**3*r_earth[i,:] - GM_sun*m_jupiter/distance_E_J**3*r_E_J
    a_jupiter_new = - GM_sun/np.linalg.norm(r_jupiter[i,:])**3*r_jupiter[i,:] - GM_sun*m_earth/distance_E_J**3*(-r_E_J)
    v_earth[i,:] = v_earth[i-1,:] + 0.5*dt*(a_earth_new + a_earth)
    v_jupiter[i,:] = v_jupiter[i-1,:] + 0.5*dt*(a_jupiter_new + a_jupiter)

plt.plot(r_earth[:,0],r_earth[:,1],label='Earth')
plt.plot(r_jupiter[:,0], r_jupiter[:,1],label='Jupiter')
plt.scatter(r_sun[:,0],r_sun[:,0],label='Sun')
plt.legend()
plt.show()
