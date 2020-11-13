import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

lattice2=pd.read_csv("data/lattice2", index_col=False, sep="\s+",names=["T","E","Cv","M","Chi","|M|"])

plt.figure()
plt.plot(lattice2['T'], lattice2['Chi'])
plt.show()

lattice20=pd.read_csv("data/lattice20", index_col=False, sep="\s+",names=["T","E","Cv","M","Chi","|M|"])

plt.figure()
plt.plot(lattice20['T'], lattice20['Chi'])
plt.show()
