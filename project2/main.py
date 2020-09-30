"""
Assuming you have scripts
    taskn.cpp (takes i arguments); i = 0,1,2,..i

Here I'm using "-compile" as an optional argument for
compiling whatever program we're running.
"""
import sys
import subprocess

# check task name has been provided
if len(sys.argv) < 1:
    print("Not enough arguments")
    sys.exit(1)

# setup
program, args = sys.argv[1], sys.argv[2:]
COMPILE = True if "-compile" in args else False
if COMPILE: args.remove("-compile")

# execute task1
if program == "test":

    if COMPILE:
        res = subprocess.run(["g++", "-std=c++17", "-c", "test.cpp"])
        res = subprocess.run(["g++", "-std=c++17", "-c", "src/utils.cpp"])
        res = subprocess.run(["g++", "-std=c++17", "test.o", "utils.o" \
        , "-o", "test.x", "-l", "armadillo"])
        res = subprocess.run(["rm", "test.o", "utils.o"])

    res = subprocess.run(["./test.x"] + args)

elif program == "task1":

    if COMPILE:
        res = subprocess.run(["g++", "-std=c++17", "-c", "scripts/bucklingbeam1.cpp"])
        res = subprocess.run(["g++", "-std=c++17", "-c", "src/utils.cpp"])
        res = subprocess.run(["g++", "-std=c++17", "bucklingbeam1.o", "utils.o" \
        , "-o", "task1.x", "-l", "armadillo"])
        res = subprocess.run(["rm", "bucklingbeam1.o", "utils.o"])

    res = subprocess.run(["./task1.x"] + args)
    res = subprocess.run(["python", "scripts/bucklingbeam1.py"])

elif program == "task2":

    if COMPILE:
        res = subprocess.run(["g++", "-std=c++17", "-c", "scripts/bucklingbeam2.cpp"])
        res = subprocess.run(["g++", "-std=c++17", "-c", "src/utils.cpp"])
        res = subprocess.run(["g++", "-std=c++17", "bucklingbeam2.o", "utils.o" \
        , "-o", "task2.x", "-l", "armadillo"])
        res = subprocess.run(["rm", "bucklingbeam2.o", "utils.o"])

    res = subprocess.run(["./task2.x"] + args)
    res = subprocess.run(["python", "scripts/bucklingbeam2.py"])

elif program == "task3":

    if COMPILE:
        res = subprocess.run(["g++", "-std=c++17", "-c", "scripts/quantumdots.cpp",\
        "-Ofast"])
        res = subprocess.run(["g++", "-std=c++17", "-c", "src/utils.cpp", \
        "-Ofast"])
        res = subprocess.run(["g++", "-std=c++17", "quantumdots.o", "utils.o" \
        , "-o", "task3.x", "-l", "armadillo", "-Ofast"])
        res = subprocess.run(["rm", "quantumdots.o", "utils.o"])

    res = subprocess.run(["./task3.x"] + args)
    res = subprocess.run(["python", "scripts/quantumdots.py"])

elif program == "task4":

    if COMPILE:
        res = subprocess.run(["g++", "-std=c++17", "-c", "scripts/quantumEigen.cpp",\
         "-Ofast", "-l", "armadillo"])
        res = subprocess.run(["g++", "-std=c++17", "-c", "src/utils.cpp", \
        "-Ofast", "-l", "armadillo"])
        res = subprocess.run(["g++", "-std=c++17", "quantumEigen.o", "utils.o" \
        , "-o", "task4.x", "-l", "armadillo", "-Ofast"])
        res = subprocess.run(["rm", "quantumEigen.o", "utils.o"])

    res = subprocess.run(["./task4.x"] + args)
    res = subprocess.run(["python", "scripts/quantumEigen.py"])

# catch failure to provide a valid program
else:
    print("There's no program called " + program + ". Your options are:")
    print("    'test'")
    print("    'task1'    - buckling beam, iterations for N=5,15,25,...,125")
    print("    'task2'    - buckling beam, eigvals and vec using jacobimethod")
    print("    'task3'    - solves the QM problem and plots eig- vals & vecs")
    print("    'task4'    - QM, Arguments required: Nstart, Nstop, Nstep, rhoN")

"""

# execute task2
elif program == "task2":
    if COMPILE:
        #compile task2.cpp
    res = subprocess.run(["./task2.x"] + args)

# execute task3
elif program == "task3":
    if COMPILE:
        #compile task3.cpp
    res = subprocess.run(["./task3.x"] + args)

"""
