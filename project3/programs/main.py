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
if program == "test_methods":

    if COMPILE:
        res = subprocess.run(["g++", "-std=c++17", "-c", "test/test.cpp", "-O3"])
        res = subprocess.run(["g++", "-std=c++17", "-c", "src/planet.cpp", "-O3"])
        res = subprocess.run(["g++", "-std=c++17", "-c", "src/solarSystem.cpp", "-O3"])
        res = subprocess.run(["g++", "-std=c++17", "test.o", "planet.o", "solarSystem.o" \
        , "-o", "test.x", "-l", "armadillo", "-O3"])
        res = subprocess.run(["rm", "test.o", "planet.o", "solarSystem.o"])

    res = subprocess.run(["./test.x"] + args)
    res = subprocess.run(["python", "test/testplot.py"])

elif program == "sstest":

    if COMPILE:
        res = subprocess.run(["g++", "-std=c++17", "-c", "scripts/sstest.cpp"])
        res = subprocess.run(["g++", "-std=c++17", "-c", "src/ss.cpp"])
        res = subprocess.run(["g++", "-std=c++17", "sstest.o", "ss.o" \
        , "-o", "sstest.x", "-l", "armadillo"])
        res = subprocess.run(["rm", "sstest.o", "ss.o"])

    res = subprocess.run(["./sstest.x"] + args)
    #res = subprocess.run(["python", "scripts/bucklingbeam1.py"])

elif program == "crp":

    if COMPILE:
        res = subprocess.run(["g++", "-std=c++17", "-c", "crp.cpp", "-O3"])
        res = subprocess.run(["g++", "-std=c++17", "-c", "src/planet.cpp", "-O3"])
        res = subprocess.run(["g++", "-std=c++17", "-c", "src/ss.cpp", "-O3"])
        res = subprocess.run(["g++", "-std=c++17", "-c", "src/utils.cpp", "-O3"])
        res = subprocess.run(["g++", "-std=c++17", "crp.o", \
        "planet.o", "ss.o", "utils.o", "-o", "crp.x", "-l", "armadillo", "-O3"])
        res = subprocess.run(["rm", "crp.o", "ss.o", "planet.o"])

    res = subprocess.run(["./crp.x"] + args)




# catch failure to provide a valid program
else:
    print("There's no program called " + program + ". Your options are:")
    print("    'test'")
    print("    'task1'    - buckling beam, iterations for N=5,15,25,...,125")

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
