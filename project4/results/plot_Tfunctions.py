import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
plt.rcParams.update({'font.size': 14})

lattice20=pd.read_csv("data/lattice20", index_col=False, sep="\s+",names=["T","E","Cv","M","Chi","|M|"])
lattice40=pd.read_csv("data/lattice40", index_col=False, sep="\s+",names=["T","E","Cv","M","Chi","|M|"])
lattice60=pd.read_csv("data/lattice60", index_col=False, sep="\s+",names=["T","E","Cv","M","Chi","|M|"])
lattice80=pd.read_csv("data/lattice80", index_col=False, sep="\s+",names=["T","E","Cv","M","Chi","|M|"])
lattice100=pd.read_csv("data/lattice100", index_col=False, sep="\s+",names=["T","E","Cv","M","Chi","|M|"])


plt.figure()
plt.scatter(lattice20['T'], lattice20['Chi'], label=r'$L=20$')
plt.scatter(lattice40['T'], lattice40['Chi'], marker="x", label=r'$L=40$')
plt.scatter(lattice60['T'], lattice60['Chi'], marker="+", label=r'$L=60$')
plt.scatter(lattice80['T'], lattice80['Chi'], marker="h", label=r'$L=80$')
plt.scatter(lattice100['T'], lattice100['Chi'], marker="d", label=r'$L=100$')
plt.xlabel(r'Temperature [$J/k_B$]')
plt.ylabel(r'Susceptibility $\chi$ [$1/J$]')
plt.grid()
plt.legend(loc="best")
plt.savefig("figures/susceptibility.png",bbox_inches = 'tight')

plt.figure()
plt.scatter(lattice20['T'], lattice20['E'], label=r'$L=20$')
plt.scatter(lattice40['T'], lattice40['E'], marker="x", label=r'$L=40$')
plt.scatter(lattice60['T'], lattice60['E'], marker="+", label=r'$L=60$')
plt.scatter(lattice80['T'], lattice80['E'], marker="h", label=r'$L=80$')
plt.scatter(lattice100['T'], lattice100['E'], marker="d", label=r'$L=100$')
plt.xlabel(r'Temperature [$J/k_B$]')
plt.ylabel(r'Energy [$J$]')
plt.grid()
plt.legend(loc="best")
plt.savefig("figures/meanenergy.png",bbox_inches = 'tight')


plt.figure()
plt.scatter(lattice20['T'], lattice20['Cv'], label=r'$L=20$')
plt.scatter(lattice40['T'], lattice40['Cv'], marker="x", label=r'$L=40$')
plt.scatter(lattice60['T'], lattice60['Cv'], marker="+", label=r'$L=60$')
plt.scatter(lattice80['T'], lattice80['Cv'], marker="h", label=r'$L=80$')
plt.scatter(lattice100['T'], lattice100['Cv'], marker="d", label=r'$L=100$')
plt.xlabel(r'Temperature [$J/k_B$]')
plt.ylabel(r'Heat Capacity [$k_B$]')
plt.grid()
plt.legend(loc="best")
plt.savefig("figures/heatcapacity.png",bbox_inches = 'tight')


plt.figure()
plt.scatter(lattice20['T'], lattice20['|M|'], label=r'$L=20$')
plt.scatter(lattice40['T'], lattice40['|M|'], marker="x", label=r'$L=40$')
plt.scatter(lattice60['T'], lattice60['|M|'], marker="+", label=r'$L=60$')
plt.scatter(lattice80['T'], lattice80['|M|'], marker="h", label=r'$L=80$')
plt.scatter(lattice100['T'], lattice100['|M|'], marker="d", label=r'$L=100$')
plt.xlabel(r'Temperature [$J/k_B$]')
plt.ylabel(r'Absolute magnetic moment')
plt.grid()
plt.legend(loc="best")
plt.savefig("figures/meanmagnetisation.png",bbox_inches = 'tight')

plt.show()
