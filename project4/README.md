
# Project 4

In this project we have performed numerical studies of a two-dimensional Ising model using the Metropolis algorithm. 
Our report can be found in the folder "report". 

The programs used for generating the results are placed in the main folder. These include
* IsingBurnIn.cpp: this program used for calculating the expected energy and absolute magnetisation of the Ising model given a specific temperature and initial orientation. It writes out these properties, in addition to the total energy energy of the system and the number of accepted spins flips, for each Monte Carlo cycle is performs. It takes 5 arguments: the output file name, the lattice dimension, the number of Monte Carlo cycles, the temperature and the orientation (given by either 'random' or 'ordered'). 
* IsingPara.cpp: this program is used for calculating the average expectated energy, expected absolute magnetisation, heat capacity and magnetic susceptibility for a range of different temperatures. It writes these average thermodynamical properties to file for each temperature. It takes 7 arguments: the output file name, the lattice dimension, the number of burn-in Monte Carlo cycles to be run before calculating average values, the number of Monte Carlo cycles used in the averages, and the temperature range start, end and step. 

The .hpp-file for the Ising class is placed in the folder 'include', while the .cpp-file is placed in 'src'.

The generated output data from the programs is placed in 'results/data'. The generated plots are placed in 'results/figures. The 'results'-folder also include Python files used to make figures of the data. 
