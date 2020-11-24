import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
import pandas as pd
plt.rcParams.update({'font.size': 14})

lattice40=pd.read_csv("data/lattice40", index_col=False, sep="\s+",names=["T","E","Cv","M","Chi","|M|"])
lattice60=pd.read_csv("data/lattice60", index_col=False, sep="\s+",names=["T","E","Cv","M","Chi","|M|"])
lattice80=pd.read_csv("data/lattice80", index_col=False, sep="\s+",names=["T","E","Cv","M","Chi","|M|"])
lattice100=pd.read_csv("data/lattice100", index_col=False, sep="\s+",names=["T","E","Cv","M","Chi","|M|"])

f_chi40 = CubicSpline(lattice40['T'], lattice40['Chi'])
f_chi60 = CubicSpline(lattice60['T'], lattice60['Chi'])
f_chi80 = CubicSpline(lattice80['T'], lattice80['Chi'])
f_chi100 = CubicSpline(lattice100['T'], lattice100['Chi'])

Tc40 = lattice40['T'].iloc[np.argmax(f_chi40(lattice40['T']))]
Tc60 = lattice60['T'].iloc[np.argmax(f_chi60(lattice60['T']))]
Tc80 = lattice80['T'].iloc[np.argmax(f_chi80(lattice80['T']))]
Tc100 = lattice100['T'].iloc[np.argmax(f_chi100(lattice100['T']))]

invL = np.linspace(0, 1./30, 100)
Tc = np.array([Tc100, Tc80, Tc60, Tc40])
invLc = np.array([1./100, 1./80, 1./60, 1./40])

# Make a linear fit using the critical temperatures
poly = np.polyfit(invLc, Tc, deg=1, full=True)

plt.figure()
plt.plot(lattice40['T'], f_chi40(lattice40['T']))
plt.plot(lattice60['T'], f_chi60(lattice60['T']))
plt.plot(lattice80['T'], f_chi80(lattice80['T']))
plt.plot(lattice100['T'], f_chi100(lattice100['T']))
plt.scatter(Tc40, max(f_chi40(lattice40['T'])), label="L = 40")
plt.scatter(Tc60, max(f_chi60(lattice60['T'])), label="L = 60")
plt.scatter(Tc80, max(f_chi80(lattice80['T'])), label="L = 80")
plt.scatter(Tc100, max(f_chi100(lattice100['T'])), label="L = 100")
plt.xlabel(r'Temperature [$J/k_B$]')
plt.ylabel(r'Susceptibility $\chi$ [$1/J$]')
plt.title('Maxima of Interpolated Susceptibility')
plt.legend()
plt.grid()
plt.savefig("figures/maximachi.png", bbox_inches = 'tight')

plt.figure()
plt.title("Critical temperature estimation")
plt.scatter(invLc, Tc, label='Simulations')
plt.scatter(0,2.269, label="Exact")
plt.plot(invL, poly[0][1] + poly[0][0]*invL, color='r', linestyle='--', label='Linear fit')
plt.legend()
plt.grid()
plt.ylabel(r'Critical temperature $T_C$  []$J/k_B$]')
plt.xlabel(r"Inverse lattice size $L^{-1}$")
plt.savefig("figures/criticalT.png", bbox_inches = 'tight')

plt.show()
