#ifndef METRO_H
#define METRO_H

#include <iostream>
#include <cmath>
#include <armadillo>
#include <stdlib.h>

class Ising
{
public:
  // Properties
  int L; // dimensions of lattice (LxL)
  double J; // coupling constant
  double T; // temperature, same as energy with K_b=1
  arma::mat spin_state; // state spin configurations of system

  // Initializer
  Ising(int L);
  Ising(int N_in, double T_in, double J_in, int seed );

  // Functions
  void genState();
  void energy();
  void metropolis();
};

#endif
