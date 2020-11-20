import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14})

MC = np.array([100, 1000, 10000, 100000, 1000000])
T1 = np.array([904, 2444, 4788, 32605, 291390])
T15 = np.array([4924, 24993, 53575, 530220, 5364380])
T2 = np.array([8014, 33109, 326008, 3211813, 31678952])
T25 = np.array([13219, 129370, 1271547, 12758606, 127559658])

plt.plot(MC, T1, '.-',label=r"$T=1.0 J/k_B$")
plt.plot(MC, T15,'.-', label=r"$T=1.5 J/k_B$")
plt.plot(MC, T2, '.-', label=r"$T=2.0 J/k_B$")
plt.plot(MC, T25, '.-', label=r"$T=2.5 J/k_B$")
plt.xscale("log")
plt.yscale("log")
plt.grid()
plt.legend()
plt.ylabel('# Accepted configurations')
plt.xlabel('Monte Carlo cycles')
plt.savefig('figures/acceptedconfigs.png', bbox_inches = 'tight')

plt.show()
