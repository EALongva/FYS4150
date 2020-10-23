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

  // Fix cent
  SolarSystem.fixOriginCentreOfMass();

  SolarSystem.velocityVerlet_relcorr(100,1000000, outputpath+"mercurysun_relcorr_t100N10e-6.dat");

}
