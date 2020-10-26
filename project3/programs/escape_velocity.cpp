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
  string fname;
  double endtime;

  if(argc == 4){
    fname = argv[1];
    endtime = atof(argv[2]);
    exponent = atoi(argv[3]);
  }
  else{
    fname = "t2N10e-5";
    endtime = 2;
    exponent = 5;
  }

  string outputpath = "../results/output/";

  double v_e = 2*M_PI*sqrt(2);
  vec v0Sun = {0, 0};
  vec p0Sun = {0, 0};
  vec p0Earth = {1, 0};
  planet Sun("Sun", 1.0, p0Sun, v0Sun);

  vec v_list = {0.9*v_e, 0.99*v_e, v_e, 1.01*v_e, 1.1*v_e};
  string name_list[5] = {"0.9", "0.99", "1.0", "1.01", "1.1"};
  for (int i = 0; i < v_list.n_elem; i++){
    vec v0Earth = {0, v_list(i)};
    planet Earth("Earth", 3.00348959632E-6, p0Earth,v0Earth);

    solarSystem SolarSystem(2, 4*M_PI*M_PI, 100);
    SolarSystem.add_planet(Sun);
    SolarSystem.add_planet(Earth);
    string filename = outputpath + "ES_escv" + name_list[i] + "_"+fname+".dat";
    SolarSystem.fixOriginCentreOfMass();
    SolarSystem.velocityVerlet(endtime,pow(10,exponent),filename);
  }
}
