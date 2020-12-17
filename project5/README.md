# Project 5

In this project we study the behavior and evolution of Rossby waves in one and two spatial dimensions. The Rossby waves are described by a partial differential equation, which we solve using both forward and centered finite difference methods and Jacobi's method. We found the centered-time-centered-space (CTCS) method to be stable for $\Delta t \leq \Delta x$. The Rossby wave equation was solved using both sine waves and centered Gaussians as initial stream functions, and the simulations were done in both a periodic and a bounded domain. We reproduced the expected behavior for the wave propagation in the westward direction in all domains and in both one and two dimensions. The phase speed was found to be $|c| \simeq 0.00625$ with a relative error of 0.0128 for an initial sine wave in a periodic domain (1D). 

The programs used for generating the results are placed in the programs folder. These include
* rossby_test.cpp: Simulates an initial sine wave in one dimensions in a periodic domain. The resulting evolution of the wave is then compared with analytical closed form solutions. testplot.py then plots the results and compares the analytical and numerical results as well as the absolute error. 
* rossby_periodic.cpp: Simulates in one dimension in a periodic domain. Input arguments: 'gridStep' 'timeStep' 'endTime' 'sineWave=True/False' 'forwardStep=True/False'
* rossby_bounded.cpp: Simulates in one dimension in a bounded domain. Input arguments: 'gridStep' 'timeStep' 'endTime' 'sineWave=True/False' 'forwardStep=True/False'
* rossby_periodic_2d.cpp: Simulates in two dimensions in a periodic domain. Input arguments: 'gridStep' 'timeStep' 'endTime' 'sineWave=True/False' 'forwardStep=True/False'
* rossby_bounded_2d.cpp: Simulates in two dimensions in a bounded domain. Input arguments: 'gridStep' 'timeStep' 'endTime' 'sineWave=True/False' 'forwardStep=True/False'

The .hpp-file for the rossby classes is placed in the folder "include", while the .cpp-file is placed in "src".

The generated output data from the programs is placed in "results/data". The generated plots are placed in "results/figures". The "results"-folder also include Python files used to make figures of the data.
