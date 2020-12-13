import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
plt.rcParams.update({'font.size': 14})

datadir = "results/data/"
figdir = "results/figures/"

wave = "bounded_sine"

data_centered = pd.read_csv(datadir+"psi_"+wave+"_centered.csv", sep=',', header=None)

x = np.linspace(0, 1, len(data_centered.iloc[:,0]))
t = np.linspace(0, 150, len(data_centered.iloc[0, :]))

levels = np.linspace(-1.1, 1.1, 12)
fig, ax = plt.subplots()
c = ax.contourf(t,x, np.array(data_centered))#, levels=levels, vmin=-1.1, vmax=1.1)
ax.set_xlabel("Time")
ax.set_ylabel("Spatial extent")
fig.colorbar(c, label="Amplitude")
plt.savefig(figdir + "hovmuller_"+wave+".png",)
plt.show()
