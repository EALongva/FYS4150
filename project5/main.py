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
if program == "example":

    if COMPILE:
        res = subprocess.run(["g++", "-std=c++17", "-c", "example.cpp"])
        res = subprocess.run(["g++", "-std=c++17", "-c", "src/exampleClass.cpp"])
        res = subprocess.run(["g++", "-std=c++17", "-c", "src/utils.cpp"])
        res = subprocess.run(["g++", "-std=c++17", "example.o", "exampleClass.o", "utils.o" \
        , "-o", "example.x", "-l", "armadillo"])
        res = subprocess.run(["rm", "example.o", "exampleClass.o", "utils.o"])

    res = subprocess.run(["./example.x"] + args)
    res = subprocess.run(["python", "results/example_plot.py"])


elif program == "test_time":

    if COMPILE:
        res = subprocess.run(["g++", "-std=c++17", "-c", "test/test_time.cpp"])
        res = subprocess.run(["g++", "-std=c++17", "-c", "src/planet.cpp"])
        res = subprocess.run(["g++", "-std=c++17", "-c", "src/solarSystem.cpp"])
        res = subprocess.run(["g++", "-std=c++17", "test_time.o", "planet.o", "solarSystem.o" \
        , "-o", "test_time.x", "-l", "armadillo"])
        res = subprocess.run(["rm", "test_time.o", "planet.o", "solarSystem.o"])

    res = subprocess.run(["./test_time.x"] + args)
    res = subprocess.run(["python", "test/timeplot.py"])

elif program == "nrg":

    if COMPILE:
        res = subprocess.run(["g++", "-std=c++17", "-c", "nrg/nrg.cpp"])
        res = subprocess.run(["g++", "-std=c++17", "-c", "src/planet.cpp"])
        res = subprocess.run(["g++", "-std=c++17", "-c", "src/solarSystem.cpp"])
        res = subprocess.run(["g++", "-std=c++17", "nrg.o", "planet.o", "solarSystem.o" \
        , "-o", "nrg.x", "-l", "armadillo"])
        res = subprocess.run(["rm", "nrg.o", "planet.o", "solarSystem.o"])

    res = subprocess.run(["./nrg.x"] + args)
    res = subprocess.run(["python", "nrg/nrg.py"])




# catch failure to provide a valid program
else:
    print("There's no program called " + program + ". Your options are:")
    print("    'test_methods'    ")
    print("    'test_time'      ")
    print("    'nrg'            ")
    print("    'solarsys_orbits'    - please provide, time in years T, power of timesteps N")
    print("    'SunEarthJupiter'    - please provide, time in years T, power of timesteps N")
    print("    'varying_force'      - please provide, filename + init_velocity (single float) \
    + init_distance (single float) + endtime + exponent (1/r^exponent)")
    print("    'MercuryPerihelion'      - please provide, filename, endtime (years), power of timesteps N")
    print("    'escape_velocity'      - please provide, filename, endtime (years), power of timesteps N")
