## FYS4150 - PROJECT 2
Computational Physics

Most of the code needed to run the scripts for this project is located as
functions in src/utils.cpp, but the bisection method has to be run indepen-
dently using "bisection.cpp".

To obtain all the plots and results use "main.py" and run the different tasks.
The tasks are described at the bottom of the file but the same description is
included here:

# cmd line argument:  Description:
- test
- task1             buckling beam, iterations for N=5,15,25,...,125")
- task2             buckling beam, eigvals and vec using jacobimethod")
- task3             solves the QM problem and plots eig- vals & vecs")
- task4             QM, Arguments required: Nstart, Nstop, Nstep, rhoN")
- task5             BB, comparing analytical eigvals/vecs with armadillo")
- task6             BB, comparing CPU time, jacobi vs armadillo")

example: "python3 task1 -compile
note: -compile flag is only required on the first run :)

Corresponding scripts:
- test    -> test.cpp
- task1   -> scripts/bucklingbeam1.cpp
- task2   -> scripts/bucklingbeam2.cpp
- task3   -> scripts/quantumdots.cpp
- task4   -> scripts/quantumEigen.cpp
- task5   -> scripts/bucklingbeam3.cpp
- task6   -> scripts/bucklingbeam4.cpp
