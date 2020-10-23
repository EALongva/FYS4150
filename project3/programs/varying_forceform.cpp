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

  vec v0Sun("0, 0");
  vec p0Sun("0, 0");
  vec v0Earth("0, 5.0")
  //vec v0Earth("0, 6.28318");
  vec p0Earth("1.0, 0");

  planet Earth("Earth", 3.00348959632E-6, p0Earth,v0Earth);
  planet Sun("Sun", 1.0, p0Sun, v0Sun);
  solarSystem SolarSystem(2, 4*9.869604401, 100);
  SolarSystem.add_planet(Sun);
  SolarSystem.add_planet(Earth);
  SolarSystem.fixOriginCentreOfMass();
  SolarSystem.velocityVerlet(2.0,50000,outputpath+"earthsun_t2N10e-5.dat");

  vec beta_list = {2, 2.25, 2.5, 2.75, 3};
  string name_list[5] = {"2", "2.25", "2.5", "2.75", "3"};
  for (int i = 0; i < beta_list.n_elem; i++){
    planet Earth("Earth", 3.00348959632E-6, p0Earth,v0Earth, beta_list(i));
    planet Sun("Sun", 1.0, p0Sun, v0Sun, beta_list(i));

    solarSystem SolarSystem(2, 4*9.869604401, 100);
    SolarSystem.add_planet(Sun);
    SolarSystem.add_planet(Earth);
    string filename = outputpath + "earthsun_beta" + name_list[i] + "_t2N10e-5.dat";
    SolarSystem.fixOriginCentreOfMass();
    SolarSystem.velocityVerlet(2.0,100000,filename);
  }

  planet Earth2("Earth2", 3.00348959632E-6, p0Earth,v0Earth,2);
  planet Sun2("Sun2", 1.0, p0Sun, v0Sun,3);
  solarSystem SolarSystem2(2, 4*9.869604401, 100);
  SolarSystem2.add_planet(Sun2);
  SolarSystem2.add_planet(Earth2);
  SolarSystem2.fixOriginCentreOfMass();
  SolarSystem2.velocityVerlet(2.0,100000,outputpath+"earthsun_beta3_t2N10e-5.dat");

}
