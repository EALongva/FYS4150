import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

lattice20=pd.read_csv("data/lattice20", index_col=False, sep="\s+",names=["T","E","Cv","M","Chi","|M|"])
lattice40=pd.read_csv("data/lattice40", index_col=False, sep="\s+",names=["T","E","Cv","M","Chi","|M|"])
lattice60=pd.read_csv("data/lattice60", index_col=False, sep="\s+",names=["T","E","Cv","M","Chi","|M|"])
lattice80=pd.read_csv("data/lattice80", index_col=False, sep="\s+",names=["T","E","Cv","M","Chi","|M|"])
lattice100=pd.read_csv("data/lattice100", index_col=False, sep="\s+",names=["T","E","Cv","M","Chi","|M|"])


plt.figure()
plt.scatter(lattice20['T'], lattice20['Chi']/(20*20), label=r'$L=20$')
plt.scatter(lattice40['T'], lattice40['Chi']/(40*40), marker="x", label=r'$L=40$')
plt.scatter(lattice60['T'], lattice60['Chi']/(60*60), marker="+", label=r'$L=60$')
plt.scatter(lattice80['T'], lattice80['Chi']/(80*80), marker="h", label=r'$L=80$')
plt.scatter(lattice100['T'], lattice100['Chi']/(100*100), marker="d", label=r'$L=100$')
plt.xlabel(r'Temperature [$k_BT/J$]')
plt.ylabel(r'Susceptibility $\chi$')
plt.grid()
plt.legend(loc="best")
plt.savefig("figures/susceptibility.png",bbox_inches = 'tight')
plt.show()
