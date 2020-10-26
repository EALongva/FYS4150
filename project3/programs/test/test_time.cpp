// UNIT TESTS PROJECT 3

#include <iostream>
#include <armadillo>
#include <vector>
#include <cstdlib>
#include <cmath>
#include "time.h"
#include "../include/solarSystem.hpp"
#include "../include/planet.hpp"
//#include "../include/utils.hpp"

int main(int argc,char* argv[]){

  double T = 5;
  int N;
  int range = 11; // range of timesteps
  int runs = 10; // number of timing runs per N

  //arma::vec timesteps = arma::linspace(std::pow(10,2),std::pow(10,6),range);

  arma::vec timesteps = arma::linspace(2,6,range);
  timesteps = arma::exp10(timesteps);
  timesteps = arma::floor(timesteps);

  double M_sun = 2*std::pow(10,30);
  double G = 4*M_PI*M_PI;


  // defining the sun
  arma::vec sunpos = {0.0,0.0};
  arma::vec sunvel = {0.0,0.0};
  double sunmass = 1.0;
  std::string sunname = "sun";
  planet sun(sunname, sunmass, sunpos, sunvel);

  // definin earth for circular orbit
  arma::vec earthpos = {-1.0,0.0};
  arma::vec earthvel = {0.0,2*M_PI};
  double earthmass = 6*std::pow(10,24)/M_sun;
  std::string earthname = "earth";
  planet earth(earthname, earthmass, earthpos, earthvel);

  solarSystem ss(2, G, 0);
  ss.add_planet(sun);
  ss.add_planet(earth);

  std::cout << "now we cookin" << std::endl;

  arma::vec dt = T/timesteps;
  arma::vec eult(range, arma::fill::zeros);
  arma::vec vvt(range, arma::fill::zeros);

  clock_t start;
  clock_t stop;

  // simulating the orbit using euler and vv, saving the position data and
  // distance from the sun

  for (int n=0; n < range; n++){

    N = timesteps(n);

    // timing
    start = clock();
    for (int i=0; i<runs; i++){}

      ss.euler_(T,N);
      ss.reset_evolution();

    stop = clock();
    eult(n) = ((double) (stop-start)/CLOCKS_PER_SEC)/((double)runs);

  }


  for (int n=0; n < range; n++){

    N = timesteps(n);

    // timing
    start = clock();
    for (int i=0; i<runs; i++){}

      ss.vv_(T,N);
      ss.reset_evolution();

    stop = clock();
    vvt(n) = ((double) (stop-start)/CLOCKS_PER_SEC)/((double)runs);

  }

dt.save("test/test_time_dt.csv", arma::csv_ascii);
eult.save("test/test_time_eult.csv", arma::csv_ascii);
vvt.save("test/test_time_vvt.csv", arma::csv_ascii);

}
