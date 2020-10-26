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
        res = subprocess.run(["g++", "-std=c++17", "-c", "test/test.cpp"])
        res = subprocess.run(["g++", "-std=c++17", "-c", "src/planet.cpp"])
        res = subprocess.run(["g++", "-std=c++17", "-c", "src/solarSystem.cpp"])
        res = subprocess.run(["g++", "-std=c++17", "test.o", "planet.o", "solarSystem.o" \
        , "-o", "test.x", "-l", "armadillo"])
        res = subprocess.run(["rm", "test.o", "planet.o", "solarSystem.o"])

    res = subprocess.run(["./test.x"] + args)
    res = subprocess.run(["python", "test/testplot.py"])


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

elif program == "solarsys_orbits":

    if COMPILE:
        res = subprocess.run(["g++", "-std=c++17", "-c", "ssorbits.cpp"])
        res = subprocess.run(["g++", "-std=c++17", "-c", "src/planet.cpp"])
        res = subprocess.run(["g++", "-std=c++17", "-c", "src/solarSystem.cpp"])
        res = subprocess.run(["g++", "-std=c++17", "-c", "src/utils.cpp"])
        res = subprocess.run(["g++", "-std=c++17", "ssorbits.o", "planet.o", "solarSystem.o" \
        ,"utils.o", "-o", "ssorbits.x", "-l", "armadillo"])
        res = subprocess.run(["rm", "ssorbits.o", "planet.o", "solarSystem.o", "utils.o"])

    res = subprocess.run(["./ssorbits.x"] + args)
    res = subprocess.run(["python", "data/ssorbits_plots.py"])

elif program == "SunEarthJupiter":

    if COMPILE:
        res = subprocess.run(["g++", "-std=c++17", "-c", "SunEarthJupiter.cpp"])
        res = subprocess.run(["g++", "-std=c++17", "-c", "src/planet.cpp"])
        res = subprocess.run(["g++", "-std=c++17", "-c", "src/solarSystem.cpp"])
        res = subprocess.run(["g++", "-std=c++17", "-c", "src/utils.cpp"])
        res = subprocess.run(["g++", "-std=c++17", "SunEarthJupiter.o", "planet.o", "solarSystem.o" \
        ,"utils.o", "-o", "SunEarthJupiter.x", "-l", "armadillo"])
        res = subprocess.run(["rm", "SunEarthJupiter.o", "planet.o", "solarSystem.o", "utils.o"])

    res = subprocess.run(["./SunEarthJupiter.x"] + args)
    res = subprocess.run(["python", "data/SunEarthJupiter_plots.py"])

elif program == "varying_force":

    if COMPILE:
        res = subprocess.run(["g++", "-std=c++17", "-c", "varying_forceform.cpp"])
        res = subprocess.run(["g++", "-std=c++17", "-c", "src/planet.cpp"])
        res = subprocess.run(["g++", "-std=c++17", "-c", "src/solarSystem.cpp"])
        res = subprocess.run(["g++", "-std=c++17", "varying_forceform.o", "planet.o", "solarSystem.o" \
    , "-o", "varying_forceform.x", "-l", "armadillo"])
        res = subprocess.run(["rm", "varying_forceform.o", "planet.o", "solarSystem.o"])

    res = subprocess.run(["./varying_forceform.x"] + args)

elif program == "MercuryPerihelion":

    if COMPILE:
        res = subprocess.run(["g++", "-std=c++17", "-c", "perihelion_mercury.cpp"])
        res = subprocess.run(["g++", "-std=c++17", "-c", "src/planet.cpp"])
        res = subprocess.run(["g++", "-std=c++17", "-c", "src/solarSystem.cpp"])
        res = subprocess.run(["g++", "-std=c++17", "perihelion_mercury.o", "planet.o", "solarSystem.o" \
    , "-o", "perihelion_mercury.x", "-l", "armadillo"])
        res = subprocess.run(["rm", "perihelion_mercury.o", "planet.o", "solarSystem.o"])

    res = subprocess.run(["./perihelion_mercury.x"] + args)

elif program == "escape_velocity":

    if COMPILE:
        res = subprocess.run(["g++", "-std=c++17", "-c", "escape_velocity.cpp"])
        res = subprocess.run(["g++", "-std=c++17", "-c", "src/planet.cpp"])
        res = subprocess.run(["g++", "-std=c++17", "-c", "src/solarSystem.cpp"])
        res = subprocess.run(["g++", "-std=c++17", "escape_velocity.o", "planet.o", "solarSystem.o" \
    , "-o", "escape_velocity.x", "-l", "armadillo"])
        res = subprocess.run(["rm", "escape_velocity.o", "planet.o", "solarSystem.o"])

    res = subprocess.run(["./escape_velocity.x"] + args)




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
