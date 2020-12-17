# test plot rossby

import numpy as np
import matplotlib.pyplot as plt

analytical_Psi = np.genfromtxt("../results/data/test_analytical_periodic.csv", delimiter=",")
numerical_Psi = np.genfromtxt("../results/data/test_numerical_periodic.csv", delimiter=",")
relerror = np.genfromtxt("../results/data/test_relerror_periodic.csv", delimiter=",")

analytical_Psi0 = np.genfromtxt("../results/data/test_analytical_periodic0.csv", delimiter=",")
numerical_Psi0 = np.genfromtxt("../results/data/test_numerical_periodic0.csv", delimiter=",")
analytical_PsiHalf = np.genfromtxt("../results/data/test_analytical_periodicHalf.csv", delimiter=",")
numerical_PsiHalf = np.genfromtxt("../results/data/test_numerical_periodicHalf.csv", delimiter=",")

x = np.linspace(0,1,len(analytical_Psi))

plt.figure()

plt.subplot(411)
plt.title("Initial time")
plt.plot(x,analytical_Psi0,"g--")
plt.plot(x,numerical_Psi0, "r")
plt.legend(["Analytic solution", "Numerical result"])
plt.grid()

plt.subplot(412)
plt.title("t = 8$\pi^2$")
plt.plot(x,analytical_PsiHalf,"g--")
plt.plot(x,numerical_PsiHalf, "r")
plt.grid()

plt.subplot(413)
plt.title("Final time, t = 16$\pi^2$")
plt.plot(x,analytical_Psi,"g--")
plt.plot(x,numerical_Psi, "r")
plt.grid()

plt.subplot(414)
plt.title("Absolute error (final time)")
plt.ylim((0,1))
plt.plot(x,relerror,"b--")
plt.grid()

plt.tight_layout()
plt.savefig("../results/figures/test_fig.png", dpi=400)
plt.show()
