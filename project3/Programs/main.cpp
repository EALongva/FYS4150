#include <iostream>
#include <fstream>
#include <armadillo>
#include <cmath>
#include <vector>
#include "time.h"
#include <new>

int main(int argc,char* argv[]){

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
  FILE *fp_init = fopen(filename_pos_vel, "r"); //Open file to read, specified by "r".
  FILE *fp_mass = fopen(filename_mass, "r") //Open file to read.

  //Loop over each particle and extract its mass and initial conditions:
  for (int i = 0; i < Nparticles; i++){
  	fscanf(fp_init, "%lf %lf %lf %lf %lf %lf", &x[i], &y[i], &z[i], &vx[i], &vy[i], &vz[i]); // One %lf (lf=long float or double) for each floating point number on each line of the file.
  	fscanf(fp_mass, "%lf", &mass[i]); //Extract mass for particle i.
  }

  fclose(fp_init); //Close file with initial conditions
  fclose(fp_mass); //Close file with masses.
}
