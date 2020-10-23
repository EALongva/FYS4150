import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14})
inputpath = "../results/output/"
figpath = "../results/figures/"

dim = 2 # Spatial dimension
Nbodies = 2 # Number of celestial bodies

#raw = pd.read_csv(inputpath+"mercurysun_t100N10e-6.dat",header=None)
raw_r = pd.read_csv(inputpath+"mercurysun_relcorr_t100N10e-6.dat",header=None)
n_steps = int(len(raw_r)/(1+(dim+3)*Nbodies))

t = np.zeros(n_steps)
# x_positions = np.zeros((n_steps, Nbodies))
# y_positions = np.zeros((n_steps, Nbodies))
# planet_KE = np.zeros((n_steps, Nbodies))
# planet_PE = np.zeros((n_steps, Nbodies))
# planet_L = np.zeros((n_steps, Nbodies))
x_positions_r = np.zeros((n_steps, Nbodies))
y_positions_r = np.zeros((n_steps, Nbodies))
planet_KE_r = np.zeros((n_steps, Nbodies))
planet_PE_r = np.zeros((n_steps, Nbodies))
planet_L_r = np.zeros((n_steps, Nbodies))

for i in range(0, n_steps):
    int_t = i*(1+(dim+3)*Nbodies)
    t[i]= raw.iloc[int_t,0]
    for n in range(0,Nbodies):
        # x_positions[i,n] = raw.iloc[int_t+1+n*dim, 0]
        # y_positions[i,n] = raw.iloc[int_t+2+n*dim, 0]
        # planet_KE[i,n] = raw.iloc[int_t+dim*Nbodies+2+n*dim, 0]
        # planet_PE[i,n] = raw.iloc[int_t+dim*Nbodies+3+n*dim, 0]
        # planet_L[i,n] = raw.iloc[int_t+dim*Nbodies+4+n*dim, 0]
        x_positions_r[i,n] = raw_r.iloc[int_t+1+n*dim, 0]
        y_positions_r[i,n] = raw_r.iloc[int_t+2+n*dim, 0]
        # planet_KE_r[i,n] = raw_r.iloc[int_t+dim*Nbodies+2+n*dim, 0]
        # planet_PE_r[i,n] = raw_r.iloc[int_t+dim*Nbodies+3+n*dim, 0]
        # planet_L_r[i,n] = raw_r.iloc[int_t+dim*Nbodies+4+n*dim, 0]

mercury_year = 88/365 # Length of a Mercury year in years (Earth years)
n_mercury_years = int(round(100/mercury_year)) # Number of Mercury years in 100 years
n_steps_year = int(round(n_steps / 100 * mercury_year)) # Number of steps each Mercury year
rad_to_arcsec = 206264.806 # Converting radians to arcsec


#perihelion_angle = np.zeros(n_mercury_years)
perihelion_angle_r = np.zeros(n_mercury_years)
for i in range(0,n_mercury_years):
    #theta = np.arctan(y_positions[i*n_steps_year,1]/x_positions[i*n_steps_year,1])
    #perihelion_angle[i] = theta#*rad_to_arcsec
    theta_r = np.arctan(y_positions_r[i*n_steps_year,1]/x_positions_r[i*n_steps_year,1])
    perihelion_angle_r[i] = theta_r#*rad_to_arcsec

t_vec = np.linspace(0,100,n_mercury_years)

#plt.plot(t_vec, perihelion_angle, label="Pure Newtonian")
plt.plot(t_vec, perihelion_angle_r, label="Relativistic")
plt.legend()
plt.xlabel('Time [years]')
#plt.ylabel('Angle [arc seconds]')
plt.ylabel('Angle [radians]')
plt.title('Perihelion Angle of Mercury')
plt.savefig(figpath+"perihelion_mercury.png", bbox_inches = 'tight')
plt.show()
