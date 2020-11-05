#ifndef METRO_H
#define METRO_H

#include <iostream>
#include <cmath>
#include <armadillo>

class metro
{
public:
  // Properties
  int N; // dimensions of lattice (NxN)
  double J; // coupling constant
  double T; // temperature, same as energy with K_b=1
  int seed; // seed for randomgenerators
  double E;
  arma::mat state; // state of system

  // Initializer
  metro(int N_in);
  metro(int N_in, double T_in, double J_in, int seed );

  // Functions
  void genState();
  void energy();

};

#endif
