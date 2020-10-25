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

int main(int argc, char *argv[]){

  int exponent;
  string outputfilename;
  double endtime;

  if(argc <= 2){
    cout << "Two few arguments in " << argv[0] <<
    ". Read in name of output file and the exponent of the maximum number of gridpoints" << endl;
    exit(1);
  }
  else{
    outputfilename = argv[1];
    endtime = atoi(argv[2]);
    exponent = atoi(argv[3]);
  }

  string outputpath = "../results/output/";

  // Initial positions and velocities
  vec v0Sun("0, 0");
  vec p0Sun("0, 0");
  vec v0Merc("0, 12.44");
  vec p0Merc("0.3075,0");

  // Contruct celestial bodies
  planet Mercury("Mercury",0.16601E-6, p0Merc,v0Merc);
  planet Sun("Sun", 1.0, p0Sun, v0Sun);

  // Construct solar system and add planets
  solarSystem SolarSystem(2, 4*9.869604401, 100);
  SolarSystem.add_planet(Sun);
  SolarSystem.add_planet(Mercury);

  // Fix center of mass
  SolarSystem.fixOriginCentreOfMass();

  SolarSystem.perihelionAngle(endtime, pow(10,exponent), outputpath+"MS_"+outputfilename);

  // Reset planet positions and velocities
  planet Mercury_r("Mercury",0.16601E-6, p0Merc,v0Merc);
  planet Sun_r("Sun", 1.0, p0Sun, v0Sun);
  solarSystem SolarSystem_r(2, 4*9.869604401, 100);
  SolarSystem_r.add_planet(Sun_r);
  SolarSystem_r.add_planet(Mercury_r);
  SolarSystem_r.fixOriginCentreOfMass();

  SolarSystem_r.perihelionAngle_relcorr(endtime, pow(10,exponent), outputpath+"MS_relcorr_"+outputfilename);

}
