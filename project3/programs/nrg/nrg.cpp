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
  /*
  if (argc < 3){
    std::cout << "please provide, time in years T, max power in timesteps N, range of timesteps" << std::endl;
  }

  double T = std::atof(argv[1]);

  int N_in;
  N_in = std::atoi(argv[2]);
  int N = std::pow(10,N_in);

  int range = std::atoi(argv[3]);
  */

  double T = 5;
  int N;
  int range = 2;
  arma::vec timesteps = arma::linspace(2,4,range);
  timesteps = arma::exp10(timesteps);
  timesteps = arma::floor(timesteps);


  double M_sun = 2*std::pow(10,30);
  double G = 4*M_PI*M_PI;

  /*
  arma::mat g;
  arma::vec r;
  g = {{1,2,2,4,5},
        {1,2,2,3,2}};
  r = arma::trans(arma::sqrt(arma::square(g.row(0)) + arma::square(g.row(1))));
  r.print();
  */



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

  arma::vec dr;
  std::string filename;
  std::string nstr;
  std::string csv = ".csv";

  // simulating the orbit using euler and vv, saving the position data and
  // distance from the sun

  for (int n=0; n < range; n++){

    N = timesteps(n);

    ss.euler_(T,N);

    // setting up filename
    nstr = std::to_string(n);
    filename = "nrg/eulPosN" + nstr + csv;
    // saving orbit positions
    ss.pos_evo.slice(1).save(filename, arma::csv_ascii);
    filename = "nrg/eulVelN" + nstr + csv;
    ss.vel_evo.slice(1).save(filename, arma::csv_ascii);

    ss.reset_evolution();

  }


  for (int n=0; n < range; n++){

    N = timesteps(n);

    ss.vv_(T,N);

    // setting up filename
    nstr = std::to_string(n);
    filename = "nrg/vvPosN" + nstr + csv;
    // saving orbit positions
    ss.pos_evo.slice(1).save(filename, arma::csv_ascii);
    filename = "nrg/vvVelN" + nstr + csv;
    ss.vel_evo.slice(1).save(filename, arma::csv_ascii);

    ss.reset_evolution();

  }


  /*
  ss.euler_(T,N);
  ss.vv_(T,N*10);
  ss.pos_evo.slice(1).save("test/test.csv", arma::csv_ascii);
  */

}
