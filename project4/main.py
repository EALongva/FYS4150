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
if program == "crp":

    if COMPILE:
        res = subprocess.run(["mpic++", "-std=c++17", "-c", "crp.cpp", "-O3"])
        res = subprocess.run(["mpic++", "-std=c++17", "-c", "src/ising.cpp", "-O3"])
        res = subprocess.run(["mpic++", "-std=c++17", "-c", "src/utils.cpp", "-O3"])
        res = subprocess.run(["mpic++", "-std=c++17", "crp.o", "ising.o", "utils.o" \
        , "-o", "crp.x", "-l", "armadillo", "-O3"])
        res = subprocess.run(["rm", "crp.o", "ising.o", "utils.o"])

    res = subprocess.run(["./crp.x"] + args)


elif program == "test":

    if COMPILE:
        res = subprocess.run(["mpic++", "-std=c++17", "-c", "test.cpp", "-O3"])
        res = subprocess.run(["mpic++", "-std=c++17", "-c", "src/ising.cpp", "-O3"])
        res = subprocess.run(["mpic++", "-std=c++17", "-c", "src/utils.cpp", "-O3"])
        res = subprocess.run(["mpic++", "-std=c++17", "test.o", "ising.o", "utils.o" \
        , "-o", "test.x", "-l", "armadillo", "-O3"])
        res = subprocess.run(["rm", "test.o", "ising.o", "utils.o"])

    res = subprocess.run(["./test.x"] + args)
    res = subprocess.run(["python", "results/example_plot.py"])

elif program == "IsingPara":

    if COMPILE:
        res = subprocess.run(["mpic++", "-std=c++17", "-c", "IsingPara.cpp", "-O3"])
        res = subprocess.run(["mpic++", "-std=c++17", "-c", "src/ising.cpp", "-O3"])
        res = subprocess.run(["mpic++", "-std=c++17", "-c", "src/utils.cpp", "-O3"])
        res = subprocess.run(["mpic++", "-std=c++17", "IsingPara.o", "ising.o", "utils.o" \
        , "-o", "IsingPara.x", "-l", "armadillo", "-O3"])
        res = subprocess.run(["rm", "IsingPara.o", "ising.o", "utils.o"])

    res = subprocess.run(["mpirun", "-np", "4", "./IsingPara.x"] + args)
    #res = subprocess.run(["python", "results/example_plot.py"])




# catch failure to provide a valid program
else:
    print("There's no program called " + program + ". Your options are:")
    print("    'crp'    ")
    print("    'test'      ")
    print("    'IsingPara'    - please provide, time in years T, power of timesteps N")
