#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <armadillo>
#include <cmath>
#include <vector>
#include <math.h>
#include "time.h"
#include <new>
#include "include/solarSystem.hpp"
#include "include/planet.hpp"

using namespace arma;
using namespace std;

int main(int argc, char *argv[]){

  int exponent;
  double v0E;
  string fname;
  double endtime;

  if(argc == 5){
    fname = argv[1];
    v0E = atof(argv[2]);
    endtime = atof(argv[3]);
    exponent = atoi(argv[4]);
  }
  else{
    fname = "v2pi_t2N10e-5";
    v0E = 2*M_PI;
    endtime = 2;
    exponent = 5;
  }

  string outputpath = "../results/output/";

  vec v0Sun = {0, 0};
  vec p0Sun = {0, 0};
  vec v0Earth = {0, v0E};
  vec p0Earth = {1, 0};

  vec beta_list = {2, 2.25, 2.5, 2.75, 3};
  string name_list[5] = {"2", "2.25", "2.5", "2.75", "3"};
  for (int i = 0; i < beta_list.n_elem; i++){
    planet Earth("Earth", 3.00348959632E-6, p0Earth,v0Earth, beta_list(i));
    planet Sun("Sun", 1.0, p0Sun, v0Sun, beta_list(i));

    solarSystem SolarSystem(2, 4*M_PI*M_PI, 100);
    SolarSystem.add_planet(Sun);
    SolarSystem.add_planet(Earth);
    string filename = outputpath + "ES_beta" + name_list[i] + "_"+fname+".dat";
    SolarSystem.fixOriginCentreOfMass();
    SolarSystem.velocityVerlet(endtime,pow(10,exponent),filename);
  }
}
