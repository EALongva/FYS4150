// testing ising class

#include <iostream>
#include <armadillo>
#include <vector>
#include <cstdlib>
#include <cmath>
#include "time.h"
#include "include/ising.hpp"

double E_variance(double T, double J);
double M_variance(double T, double J);

int main(){

  // testing ising class
  int N       = 2;
  int seed    = 1337;
  int ordered = 1;
  double T    = 1.0;
  double J    = 1.0;
  int burnin  = 2;

  // monte carlo cycles to check for
  int MC_start = 1;
  int MC_end   = 7;
  int MC_steps = 20;
  double dMC   = (MC_end - MC_start)/((double) MC_steps);
  //std::cout << dMC << std::endl;

  // analytical values

  double Evar = E_variance(T, J);
  double Mvar = M_variance(T, J);

  // checking mag
  double beta     = 1/T; // note: Boltzmann constant = 1
  double fac      = 8 * J * beta;
  double Z        = 12 + 4 * std::cosh(fac);

  double M        = (16 + 8 * std::exp(fac))/Z; // absolute value of magnetisation |M|
  double M2       = (32 + 32 * std::exp(fac))/Z;

  arma::mat relError(2, MC_steps, arma::fill::zeros);

  Ising Isi(N, T, J, seed, ordered);
  //Isi.burnin(burnin);

  // running the monte carlo simulation for MC_steps number of MC cycles
  for(int i = 0; i <= MC_steps; i++){
    double cyc = MC_start + i*dMC;
    Isi.MC(cyc);

    double numEvar = Isi.expvals[1] - Isi.expvals[0]*Isi.expvals[0];
    double numMvar = Isi.expvals[3] - Isi.expvals[2]*Isi.expvals[2];

    //std::cout << numEvar << " --- " << Evar << std::endl;
    //std::cout << numMvar << " --- " << Mvar << std::endl;
    //std::cout << Isi.expvals[2] << " --- " << M << std::endl;
    //std::cout << Isi.expvals[3] << " --- " << M2 << std::endl;

    relError[0,i] = std::fabs((numEvar - Evar)/Evar);
    relError[1,i] = std::fabs((numMvar - Mvar)/Mvar);

    Isi.reset2init();

  }

  // comparing numerical and analytical values
  relError.print();
  relError.save("results/data/relative_error.csv", arma::csv_ascii);

  return 0;
}


double E_variance(double T, double J){

  double beta     = 1/T; // note: Boltzmann constant = 1
  double fac      = 8 * J * beta;
  double Z        = 12 + 4 * std::cosh(fac);

  double E        = -(32 * J * std::sinh(fac))/Z;
  double E2       = (256 * J*J * std::cosh(fac))/Z;

  double Evar     = E2 - E*E;
  return Evar;
}


double M_variance(double T, double J){

  double beta     = 1/T; // note: Boltzmann constant = 1
  double fac      = 8 * J * beta;
  double Z        = 12 + 4 * std::cosh(fac);

  double M        = (16 + 8 * std::exp(fac))/Z; // absolute value of magnetisation |M|
  double M2       = (32 + 32 * std::exp(fac))/Z;

  double Mvar     = M2 - M*M;
  return Mvar;
}
