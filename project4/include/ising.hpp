#ifndef ISING_H
#define ISING_H

#include <iostream>
#include <cmath>
#include <armadillo>
#include <stdlib.h>

class Ising
{
public:
  // Properties
  int N; // dimensions of lattice (LxL)
  double J; // coupling constant
  double T; // temperature, same as energy with K_b=1
  arma::mat state; // state spin configurations of system

  // Initializer
  Ising(int N);
  Ising(int N_in, double T_in, double J_in, int seed );

  // Functions
  void genState();
  void energy();
  void Metropolis();
};

#endif
