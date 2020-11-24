import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
plt.rcParams.update({'font.size': 14})

temp1rand=pd.read_csv("data/temp1random", sep="\s+",names=["MC","expE","|M|","E"])
temp2rand=pd.read_csv("data/temp2random", sep="\s+",names=["MC","expE","|M|","E"])


steadystateE1 = temp1rand['E'].iloc[100000:]
steadystateE2 = temp2rand['E'].iloc[100000:]

var1 = np.var(steadystateE1)
var2 = np.var(steadystateE2)

print("The variance for T = 1.0: %g" %var1)
print("The variance for T = 2.4: %g" %var2)

plt.figure()
plt.hist(steadystateE1,bins=20,color='orange')
plt.title(r'P(E) for $T=1.0 \quad J/k_B$')
plt.xlabel('Energy [$J$]')
plt.ylabel('Frequency')
plt.grid()
plt.savefig('figures/probdist_T1.png',bbox_inches = 'tight')

plt.figure()
plt.title(r'P(E) for $T=2.4 \quad J/k_B$')
plt.hist(steadystateE2,bins=20, color='orange')
plt.xlabel('Energy [$J$]')
plt.ylabel('Frequency')
plt.grid()
plt.savefig('figures/probdist_T24.png',bbox_inches = 'tight')

plt.show()
