// UNIT TESTS PROJECT 3

#include <iostream>
#include <armadillo>
#include <vector>
#include <cstdlib>
#include <cmath>
#include "time.h"
#include "include/solarSystem.hpp"
#include "include/planet.hpp"
#include "include/utils.hpp"

int main(int argc,char* argv[]){

  if (argc < 3){
    std::cout << "please provide, time in years T, power of timesteps N" << std::endl;
  }

  double T = std::atof(argv[1]);

  int N_in;
  N_in = std::atoi(argv[2]);
  int N = std::pow(10,N_in);

  double G = 4*M_PI*M_PI;
  double C = 365.0;

  // importing planets from nasa data
  solarSystem ss(3,G,0);

  // SUN
  std::string name_sun;
  name_sun = "Sun";
  double M_sun;
  double M_sun0; // the mass for the sun that we use!!
  M_sun = 2*std::pow(10,30);
  M_sun0 = 1.0;

  arma::vec sun_pos;
  arma::vec sun_vel;
  sun_pos = {0.0,0.0,0.0};
  sun_vel = {0.0,0.0,0.0};
  planet sun(name_sun, M_sun0, sun_pos, sun_vel);

  // EARTH
  std::string name_earth;
  name_earth = "Earth";
  double M_earth;
  M_earth = 6*std::pow(10,24)/M_sun;

  arma::vec earth_pos;
  arma::vec earth_vel;
  earth_pos = {-1.663457613546699*std::pow(10,-1),  9.691203921746886*std::pow(10,-1), -4.125583172010008*std::pow(10,-5)};
  earth_vel = {-1.723919408870981*std::pow(10,-2)*C, -2.981520896064708*std::pow(10,-3)*C,  4.254600200473125*std::pow(10,-7)*C}; //<- multiply by C=365
  planet earth(name_earth, M_earth, earth_pos, earth_vel);

  // JUPITER
  std::string name_jupiter;
  name_jupiter = "Jupiter";
  double M_jupiter;
  M_jupiter = 1.9*std::pow(10,27)/M_sun;

  arma::vec jupiter_pos;
  arma::vec jupiter_vel;
  jupiter_pos = {5.261470562232079*std::pow(10,-1), -5.201022508399864,  9.830503793253315*std::pow(10,-3)};
  jupiter_vel = {7.423674451944973*std::pow(10,-3)*C,  1.116865602956755*std::pow(10,-3)*C, -1.707572410053752*std::pow(10,-4)*C};
  planet jupiter(name_jupiter, M_jupiter, jupiter_pos, jupiter_vel);

  // first we look at earths orbital trajectory around the sun (without jupiter)
  ss.add_planet(sun);
  ss.add_planet(earth);
  ss.vv_(T,N);
  ss.pos_evo.slice(0).save("data/ESJ/SunOrbit.csv", arma::csv_ascii);
  ss.pos_evo.slice(1).save("data/ESJ/EarthOrbit.csv", arma::csv_ascii);
  ss.reset_evolution();

  // now we include jupiter in our system as well
  ss.add_planet(jupiter);

  // jupiter mass normal
  ss.vv_(T,N);
  ss.pos_evo.slice(0).save("data/ESJ/SunOrbitM1.csv", arma::csv_ascii);
  ss.pos_evo.slice(1).save("data/ESJ/EarthOrbitM1.csv", arma::csv_ascii);
  ss.pos_evo.slice(2).save("data/ESJ/JupOrbitM1.csv", arma::csv_ascii);
  ss.reset_evolution();

  // jupiter mass normal x 10
  ss.allPlanets[2].mass *= 10;
  ss.vv_(T,N);
  ss.pos_evo.slice(0).save("data/ESJ/SunOrbitM10.csv", arma::csv_ascii);
  ss.pos_evo.slice(1).save("data/ESJ/EarthOrbitM10.csv", arma::csv_ascii);
  ss.pos_evo.slice(2).save("data/ESJ/JupOrbitM10.csv", arma::csv_ascii);
  ss.reset_evolution();

  // jupiter mass normal x 1000
  ss.allPlanets[2].mass *= 100;
  ss.vv_(T,N);
  ss.pos_evo.slice(0).save("data/ESJ/SunOrbitM1000.csv", arma::csv_ascii);
  ss.pos_evo.slice(1).save("data/ESJ/EarthOrbitM1000.csv", arma::csv_ascii);
  ss.pos_evo.slice(2).save("data/ESJ/JupOrbitM1000.csv", arma::csv_ascii);
  ss.reset_evolution();

}
