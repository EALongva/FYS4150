#include <iostream>
#include <armadillo>
#include <vector>
#include <cstdlib>
#include <cmath>
#include "time.h"
#include "include/rossby.hpp"

using namespace arma;
using namespace std;

int main(int argc, char *argv[]){
  string psiname = "../results/data/psi_bounded";

  // Takes size of grid space, size of time step and length of time period as input
  double deltapos = atof(argv[1]);
  double deltatime = atof(argv[2]);
  double endtime = atof(argv[3]);
  double sigma = 0.1; double x0 = 0.5;

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
  psiname += ".csv";

  // Rossby instance with bounded domain
  rossby ross(deltapos, deltatime, endtime);
  ross.initialize_wave(sineWave, sigma, x0);
  ross.evolve_bounded(forwardStep);
  ross.Psi.save(psiname, arma::csv_ascii);

  return 0;
}
