import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

temp1rand=pd.read_csv("data/temp1random", sep="\s+",names=["MC","E","|M|"])
temp2rand=pd.read_csv("data/temp2random", sep="\s+",names=["MC","E","|M|"])
temp1ord=pd.read_csv("data/temp1ordered", sep="\s+",names=["MC","E","|M|"])
temp2ord=pd.read_csv("data/temp2ordered", sep="\s+",names=["MC","E","|M|"])

plt.figure()
plt.plot(temp1rand["MC"], temp1rand['E'], label=r"$T=1$ (random)")
plt.plot(temp2rand["MC"], temp2rand['E'], label=r"$T=2.4$ (random)")
plt.plot(temp1ord["MC"], temp1ord['E'], label=r"$T=1$ (ordered)")
plt.plot(temp2ord["MC"], temp2ord['E'], label=r"$T=2.4$ (ordered)")
plt.grid()
plt.title("Expected energy")
plt.xlabel("Monte Carlo Cycles")
plt.legend()
plt.xscale("log")
plt.savefig("figures/expnrg_10e6.png", dpi=300)

plt.figure()
plt.plot(temp1rand["MC"], temp1rand['|M|'], label=r"$T=1$ (random)")
plt.plot(temp2rand["MC"], temp2rand['|M|'], label=r"$T=2.4$ (random)")
plt.plot(temp1ord["MC"], temp1ord['|M|'], label=r"$T=2.4$ (ordered)")
plt.plot(temp2ord["MC"], temp2ord['|M|'], label=r"$T=1$ (ordered)")
plt.grid()
plt.xlabel("Monte Carlo Cycles")
plt.xscale("log")
plt.legend()
plt.title("Expected magnetisation")
plt.savefig("figures/expmag_10e6.png", dpi=300)

plt.show()


"""AverageEati = np.zeros(N)
AverageMati = np.zeros(N)
sumE = np.zeros(N)
sumM = np.zeros(N)
sumE[0] = temp1rand["E"][0]
sumM[0] = temp1rand["|M|"][0]
AverageEati[0] = temp1rand["E"][0]
AverageMati[0] = temp1rand["|M|"][0]
for i in range(N-1):
    sumE[i+1] = sumE[i] + temp1rand["E"][i+1]
    sumM[i+1] = sumM[i] + temp1rand["|M|"][i+1]
    AverageEati[i+1] = sumE[i+1]/(i+1)
    AverageMati[i+1] = sumM[i+1]/(i+1)
"""
