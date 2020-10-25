import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14})
inputpath = "../results/output/"
figpath = "../results/figures/"

dim = 2 # Spatial dimension
rad_to_arcsec = 206264.806 # Converting radians to arcsec

#raw = pd.read_csv(inputpath+"mercurysun_t100N10e-6.dat",header=None)
#raw_r = pd.read_csv(inputpath+"mercurysun_relcorr_t100N10e-6.dat",header=None)
# raw = pd.read_csv(inputpath+"mercurysun_t10N10e-8.dat",header=None)
# raw_r = pd.read_csv(inputpath+"mercurysun_relcorr_t10N10e-8.dat",header=None)
raw = pd.read_csv(inputpath+"mercurysun_t2N10e-8.dat",header=None)
raw_r = pd.read_csv(inputpath+"mercurysun_relcorr_t2N10e-8.dat",header=None)
n_steps = int(len(raw_r)/(1+dim))

t = np.zeros(n_steps)
x_p = np.zeros(n_steps)
y_p = np.zeros(n_steps)
theta_p = np.zeros(n_steps)
x_p_r = np.zeros(n_steps)
y_p_r = np.zeros(n_steps)
theta_p_r = np.zeros(n_steps)

for i in range(0, n_steps):
    int_t = i*(1+dim)
    t[i]= raw_r.iloc[int_t,0]
    x_p[i] = raw.iloc[int_t+1,0]
    y_p[i] = raw.iloc[int_t+2,0]
    theta_p[i] = np.arctan(y_p[i]/x_p[i])
    x_p_r[i] = raw_r.iloc[int_t+1,0]
    y_p_r[i] = raw_r.iloc[int_t+2,0]
    theta_p_r[i] = np.arctan(y_p_r[i]/x_p_r[i])


# mercury_year = 88/365 # Length of a Mercury year in years (Earth years)
# n_mercury_years = int(round(100/mercury_year)) # Number of Mercury years in 100 years
# n_steps_year = int(round(n_steps / 100 * mercury_year)) # Number of steps each Mercury year
# rad_to_arcsec = 206264.806 # Converting radians to arcsec


plt.plot(t, theta_p*rad_to_arcsec, label="Pure Newtonian")
plt.plot(t, theta_p_r*rad_to_arcsec, label="Relativistic")
plt.legend()
plt.xlabel('Time [years]')
plt.ylabel('Angle [arc seconds]')
#plt.ylabel('Angle [radians]')
plt.title('Perihelion Angle of Mercury')
plt.savefig(figpath+"perihelion_mercury.png", bbox_inches = 'tight')
plt.show()
