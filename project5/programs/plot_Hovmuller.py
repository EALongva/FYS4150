import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
plt.rcParams.update({'font.size': 14})

datadir = "../results/data/"
figdir = "../results/figures/"

wave = "bounded_sine"

data_centered = pd.read_csv(datadir+"psi_"+wave+"_centered.csv", sep=',', header=None)

x = np.linspace(0, 1, len(data_centered.iloc[:,0]))
t = np.linspace(0, 150, len(data_centered.iloc[0, :]))
data_flipped = np.flip(np.transpose(np.array(data_centered)),0)

fig, ax = plt.subplots()
c = ax.contourf(x, t, data_flipped)
ax.set_xlabel("East-west spatial extent")
ax.set_ylabel("Time")
fig.colorbar(c, label="Amplitude")
#plt.savefig(figdir + "hovmuller_"+wave+".png",bbox_inches = 'tight')
plt.show()
