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

  double M_sun = 2*std::pow(10,30);
  double G = 4*M_PI*M_PI;

  // importing planets using the utils function import_solarsys
  solarSystem ss(3,G,0);
  ss = import_solarsys();

  std::cout << "now we cookin" << std::endl;

  ss.fixOriginCentreOfMass();
  ss.vv_(T,N);

  arma::vec dr;
  std::string filename;
  std::string nstr;
  std::string csv = ".csv";

  for (int i=0; i<10; i++){
    // setting up filename
    nstr = ss.allPlanets[i].name;
    filename = "data/solarsystem/planet_" + nstr + csv;
    // saving orbit positions
    ss.pos_evo.slice(i).save(filename, arma::csv_ascii);
  }

  /*
  // calculating the average distance to the sun
  dr = arma::vec(N, arma::fill::zeros);

  for (int i=0; i<N; i++){
  dr(i) = arma::norm(ss.pos_evo.slice(1).col(i)-ss.pos_evo.slice(0).col(i));
  }

  filename = "test/vvdat/distanceN" + nstr + csv;
  dr.save(filename, arma::csv_ascii);
  */


}
