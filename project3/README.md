
# Project 3

Please run the programs from the /programs folder.
The main.py master file can be used to run most of the programs more easily. A small description of this file can be found below.

Please note that the main components of this project is the planet and solarsystem classes. The .cpp files are as always located in the source folder /programs/src and the .hpp headers are found in /programs/include. 

Running programs using main.py:
>>python main.py 'task' 'args' -compile

The compile flag is only required if running for the first time.
The 'task's that can be run through this program are:
 - 'test_methods'
 - 'test_time'
 - 'nrg'
 - 'solarsys_orbits'   args: 'Time' + 'N_timestepPower'
 - 'SunEarthJupiter'   args: 'Time' + 'N_timestepPower'
 -'varying_force'      args: 'filename' + 'init_velocity' (single float) + 'init_distance' (single float) + 'endtime' + 'exponent' (1/r^exponent)")
 - 'MercuryPerihelion' args: 'filename' + 'endtime' (years) + 'N_timestepPower'
 - 'escape_velocity'   args: 'filename' + 'endtime' (years) + 'N_timestepPower'
