import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
plt.rcParams.update({'font.size': 14})


temp1rand=pd.read_csv("data/temp1random", sep="\s+",names=["MC","expE","|M|","E"])
temp2rand=pd.read_csv("data/temp2random", sep="\s+",names=["MC","expE","|M|","E"])
temp1ord=pd.read_csv("data/temp1ordered", sep="\s+",names=["MC","expE","|M|","E"])
temp2ord=pd.read_csv("data/temp2ordered", sep="\s+",names=["MC","expE","|M|","E"])

plt.figure()
plt.plot(temp1rand["MC"], temp1rand['expE'], label=r"$T=1 J/k_B$ (random)")
plt.plot(temp2rand["MC"], temp2rand['expE'], label=r"$T=2.4 J/k_B$ (random)")
plt.plot(temp1ord["MC"], temp1ord['expE'], label=r"$T=1 J/k_B$ (ordered)")
plt.plot(temp2ord["MC"], temp2ord['expE'], label=r"$T=2.4 J/k_B$ (ordered)")
plt.grid()
plt.title("Expected Energy")
plt.ylabel(r"Energy [$J$]")
plt.xlabel("Monte Carlo Cycles")
plt.legend()
plt.xscale("log")
plt.savefig("figures/burninenergy.png",bbox_inches = 'tight')


plt.figure()
plt.plot(temp1rand["MC"], temp1rand['|M|'], label=r"$T=1 J/k_B$ (random)")
plt.plot(temp2rand["MC"], temp2rand['|M|'], label=r"$T=2.4 J/k_B$ (random)")
plt.plot(temp1ord["MC"], temp1ord['|M|'], label=r"$T=2.4 J/k_B$ (ordered)")
plt.plot(temp2ord["MC"], temp2ord['|M|'], label=r"$T=1 J/k_B$ (ordered)")
plt.grid()
plt.xlabel("Monte Carlo Cycles")
plt.xscale("log")
plt.ylabel("Absolute magnetisation")
plt.legend()
plt.title("Expected Magnetisation")
plt.savefig("figures/burninmagnet.png",bbox_inches = 'tight')

plt.show()
