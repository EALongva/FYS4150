/* script solving project 2 b) Finding the similarity transformations needed
to diagonalize the buckling beam matrix */

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <cmath>
#include <vector>
#include "time.h"
#include "../include/utils.hpp"

void lambda(int N, arma::vec& analyticEig){

    analyticEig = arma::zeros(N);

  for (int i=0;i<N;i++){
    analyticEig(i) = 3 + i*4;

  }
}

int main(int argc,char* argv[]){

  if (argc < 5){
    std::cout << "Arguments required: Nstart, Nstop, Nstep, rhoN" << std::endl;
    exit(1);

  }

  double   Nstart = std::atoi(argv[1]);
  double   Nstop  = std::atof(argv[2]);
  double   Nstep  = std::atoi(argv[3]);
  double  rhoN = std::atof(argv[4]);
  double  rho0 = 0;
  int      eps = -5; // we are looking for 4 decimal precision

  arma::vec N;
  arma::vec u;
  arma::vec relerr;
  arma::vec analyticEig;

  //N = arma::linspace(Nstart, Nstop, Nstep);
  //N = arma::ceil(N);

  N = arma::linspace(Nstart, Nstop, Nstep);
  N = arma::exp10(N);
  N = arma::ceil(N);

  arma::mat A;
  arma::mat R;
  arma::mat B;

  int iterations; // not necessary for this problem
  arma::vec dev(Nstep, arma::fill::zeros); // std deviation vector
  arma::vec mean(Nstep, arma::fill::zeros); // mean deviation vector

  for (int i=0; i<Nstep; i++){

    A = HOmatrix(rho0, rhoN, N(i)); //initial qm harmonic oscillator matrix
    B = jacobimethod(A, R, N(i), eps, iterations);

    u = B.diag();
    arma::vec eigval = arma::sort(u,"ascend");

    lambda(N(i), analyticEig);
    relerr = arma::abs( (analyticEig - eigval)/analyticEig );
    dev(i) = arma::stddev(relerr);
    mean(i)= arma::mean(relerr);

  }

  arma::mat M(3,Nstep,arma::fill::zeros);
  std::vector<std::string> v = {"N steps", "stddev", "mean"};
  std::string filename = "dat/QMeigen_Stddev.csv";

  M.row(0) = N.t();
  M.row(1) = dev.t();
  M.row(2) = mean.t();

  ToFile(M, v, filename);

  return 0;
}
