import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
plt.rcParams.update({'font.size': 14})

datadir = "results/data/"
figdir = "results/figures/"

wave = "periodic_sine"

data_centered = pd.read_csv(datadir+"psi_"+wave+"_centered.csv", sep=',', header=None)
data_forward = pd.read_csv(datadir+"psi_"+wave+"_forward.csv", sep=',', header=None)

x = np.linspace(0, 1, len(data_centered.iloc[:,0]))
t = np.linspace(0, 150, len(data_centered.iloc[0, :]))
dt = t[-1]/len(data_centered.iloc[0, :])

fig, ax = plt.subplots(2, 1, sharex=True)
fig.add_subplot(111, frameon=False)

plt.tick_params(labelcolor="none", top=False, bottom=False, left=False, right=False)
plt.grid()
plt.xlabel("Spatial extent")
plt.ylabel("Amplitude")

times = [50, 100, 150]
labels = ["t = 0"]
ax[0].plot(x, data_forward.iloc[:,0])
ax[1].plot(x, data_centered.iloc[:,0])

for time in times:
    labels.append("t = {}".format(time))
    timeIndex = int(time/dt) - 1
    ax[0].plot(x, data_forward.iloc[:,timeIndex])
    ax[1].plot(x, data_centered.iloc[:,timeIndex])
fig.legend(labels, frameon=False, ncol=4, bbox_to_anchor=(1.0, 0.96), fontsize=14)
plt.savefig(figdir + "compare_forward_centered.png",bbox_inches = 'tight')
plt.show()
