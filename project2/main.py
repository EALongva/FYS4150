"""
Assuming you have scripts
    task1.cpp (takes 1 argument)
    task2.cpp (takes 2 arguments)
    task3.cpp (takes no arguments)

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
if program == "task1":

    if COMPILE:
        res = subprocess.run(["g++", "-std=c++17", "-c", "scripts/bucklingbeam1.cpp"])
        res = subprocess.run(["g++", "-std=c++17", "-c", "src/utils.cpp"])
        res = subprocess.run(["g++", "-std=c++17", "bucklingbeam1.o", "utils.o" \
        , "-o", "task1.x", "-l", "armadillo"])
        res = subprocess.run(["rm", "bucklingbeam1.o", "utils.o"])

    res = subprocess.run(["./task1.x"] + args)
    res = subprocess.run(["python", "scripts/bucklingbeam1.py"])

# catch failure to provide a valid program
else:
    print("There's no program called " + program + ". Your options are:")
    print("    task1")

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
