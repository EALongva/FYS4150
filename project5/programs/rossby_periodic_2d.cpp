#include <iostream>
#include <armadillo>
#include <vector>
#include <cstdlib>
#include <cmath>
#include "time.h"
#include "include/rossby_2d.hpp"

using namespace arma;
using namespace std;

int main(int argc, char *argv[]){
  string psiname = "../results/data/psi_2d_periodic";

  // Takes size of grid space, size of time step and length of time period as input
  double deltapos = atof(argv[1]);
  double deltatime = atof(argv[2]);
  double endtime = atof(argv[3]);
  double sigma = 0.1; double x0 = 0.5; double y0 = 0.5;

  // Takes two additional input arguments that determines initial wave and time step method
  bool sineWave;
  if(atof(argv[4])==0){
    sineWave = true;
    psiname += "_sine";
  }
  else{
    sineWave = false;
    psiname += "_gaussian";
  }

  bool forwardStep;
  if(atof(argv[5])==0){
    forwardStep = true;
    psiname += "_forward";
  }
  else{
    forwardStep = false;
    psiname += "_centered";
  }

  vec times = {0, 50/deltatime-1, 100/deltatime-1, 150/deltatime-1};
  string times_list[4] = {"0", "50", "100", "150"};
  string psiname_temp;

  // Rossby instance with periodic domain
  rossby ross(deltapos, deltatime, endtime);
  ross.initialize_wave(sineWave, sigma, x0, y0);
  ross.evolve_periodic(forwardStep);

  // Stores lattice for four different times
  for (int i = 0; i < 4; i++){
    psiname_temp = psiname;
    psiname_temp = psiname_temp+times_list[i]+".csv";
    ross.Psi.slice((int) times(i)).save(psiname_temp, arma::csv_ascii);
  }

  return 0;
}
