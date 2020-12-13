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
  // string zetaname = "results/data/zeta_bounded";
  // string psiname = "results/data/psi_bounded";
  //
  // // funksjonen tar tre cmd argument, dt, dx og slutt tid
  // double deltapos = atof(argv[1]);
  // double deltatime = atof(argv[2]);
  // double endtime = atof(argv[3]);
  // double sigma = 0.1; double x0 = 0.5;
  //
  // bool sineWave;
  // if(atof(argv[4])==0){
  //   sineWave = true;
  //   zetaname += "_sine";
  //   psiname += "_sine";
  // }
  // else{
  //   sineWave = false;
  //   zetaname += "_gaussian";
  //   psiname += "_gaussian";
  //   if (argc > 5){
  //     sigma = atof(argv[6]);
  //     x0 = atof(argv[7]);
  //   }
  // }
  //
  // bool forwardStep;
  // if(atof(argv[5])==0){
  //   forwardStep = true;
  //   zetaname += "_forward";
  //   psiname += "_forward";
  // }
  // else{
  //   forwardStep = false;
  //   zetaname += "_centered";
  //   psiname += "_centered";
  // }
  // zetaname += ".csv";
  // psiname += ".csv";

  // // testing rossby class
  // rossby ross(deltapos, deltatime, endtime);
  // ross.initialize_wave(sineWave, sigma, x0);
  // ross.evolve_bounded(forwardStep);
  // ross.Psi.save(psiname, arma::csv_ascii);

  rossby ross(0.1, 0.1, 10);
  cout << ross.Psi.col(0) << endl;

  return 0;
}
