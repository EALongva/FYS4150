#include <iostream>
#include <fstream>
#include <armadillo>
#include <cmath>
#include <vector>
#include "time.h"
#include <new>
#include "include/solarSystem.hpp"
#include "include/planet.hpp"

using namespace arma;
using namespace std;

int main(){
  /*
  double *x, *y, *z, *vx, *vy, *vz; //To store initial conditions for each particle.
  double *mass; //Store mass of particles.
  int Nparticles = 10;
  x = new double[Nparticles];
  y = new double[Nparticles];
  z = new double[Nparticles];
  vx = new double[Nparticles];
  vy = new double[Nparticles];
  vz = new double[Nparticles];
  mass = new double[Nparticles];
  char* filename_pos_and_vel = "positions_and_velocities.txt";   //Each line of file gives initial condition for a particle on the form: x y z vx vy vz
  char* filename_mass = "masses.txt"; //Each line of the file contains a mass for a given particle.

  //Open files
  FILE *fp_init = fopen(filename_pos_and_vel, "r"); //Open file to read, specified by "r".
  FILE *fp_mass = fopen(filename_mass, "r"); //Open file to read.

  //Loop over each particle and extract its mass and initial conditions:
  for (int i = 0; i < Nparticles; i++){
  	fscanf(fp_init, "%lf %lf %lf %lf %lf %lf", &x[i], &y[i], &z[i], &vx[i], &vy[i], &vz[i]); // One %lf (lf=long float or double) for each floating point number on each line of the file.
  	fscanf(fp_mass, "%lf", &mass[i]); //Extract mass for particle i.
  }

  fclose(fp_init); //Close file with initial conditions
  fclose(fp_mass); //Close file with masses.
  */

  vec v0Sun(2,fill::zeros);
  vec p0Sun("0,0");
  vec v0Earth("0,6.2832");
  vec p0Earth("1.0,0");

  vec beta_list = {2, 2.25, 2.5, 2.75, 3};

  // planet Earth("Earth", 0.0003, p0Earth,v0Earth);
  // planet Sun("Sun", 1.0, p0Sun, v0Sun);
  //
  // solarSystem SolarSystem(2, 4*9.87, 100);
  // SolarSystem.add_planet(Sun);
  // SolarSystem.add_planet(Earth);
  // SolarSystem.velocityVerlet(5.6,10000,"out_positions.dat");

  for (double beta:beta_list){
    planet Earth("Earth", 0.0003, p0Earth,v0Earth, beta);
    planet Sun("Sun", 1.0, p0Sun, v0Sun, beta);

    solarSystem SolarSystem(2, 4*9.87, 100);
    SolarSystem.add_planet(Sun);
    SolarSystem.add_planet(Earth);
    string filename = "out_positions_beta_";
    filename.append(to_string(beta));
    SolarSystem.velocityVerlet(5.6,10000,filename);

    cout << SolarSystem.allPlanets[1].position << endl;
  }

}
