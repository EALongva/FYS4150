import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, FormatStrFormatter
import pandas as pd
plt.rcParams.update({'font.size': 14})

datadir = "../results/data/"
figdir = "../results/figures/"

wave = "periodic_sine"
times = [50, 100, 150]

data1 = pd.read_csv(datadir+"psi_2d_"+wave+"_forward0.csv", sep=',', header=None)
data2 = pd.read_csv(datadir+"psi_2d_"+wave+"_forward50.csv", sep=',', header=None)
data3 = pd.read_csv(datadir+"psi_2d_"+wave+"_forward100.csv", sep=',', header=None)
data4 = pd.read_csv(datadir+"psi_2d_"+wave+"_forward150.csv", sep=',', header=None)

x = np.linspace(0, 1, len(data1.iloc[:,0]))
y = np.linspace(0, 1, len(data1.iloc[0, :]))
data_flipped1 = np.transpose(np.array(data1))
data_flipped2 = np.transpose(np.array(data2))
data_flipped3 = np.transpose(np.array(data3))
data_flipped4 = np.transpose(np.array(data4))

fig, axs = plt.subplots(2, 2, sharex='all', sharey='all')
vmax = 1.6; vmin = -vmax

c1 = axs[0,0].contourf(x, y, data_flipped1, vmin=vmin, vmax=vmax)
ymin, ymax = axs[0,0].get_ylim()
axs[0,0].set_yticks(np.round(np.linspace(ymin, ymax, 3), 2)) # Reducing number of ticks
axs[0,0].set_title("t = 0",fontsize=14)

c2 = axs[0,1].contourf(x, y, data_flipped2, vmin=vmin, vmax=vmax)
axs[0,1].set_title("t = %d" %times[0],fontsize=14)

c3 = axs[1,0].contourf(x, y, data_flipped3, vmin=vmin, vmax=vmax)
axs[1,0].set_title("t = %d" %times[1],fontsize=14)

c4 = axs[1,1].contourf(x, y, data_flipped4, vmin=vmin, vmax=vmax)
axs[1,1].set_title("t = %d" %times[2],fontsize=14)

# Adding common xlabels and ylabels
fig.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.xlabel("East-west spatial extent")
plt.ylabel("North-south spatial extent")

# Adding common colorbar
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(c4, cax=cbar_ax, label="Amplitude")

plt.savefig(figdir + "hovmuller_2d_"+wave+".png",bbox_inches = 'tight')
plt.show()
