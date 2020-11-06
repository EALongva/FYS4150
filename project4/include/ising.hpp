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
  double T; // temperature, same as energy with k_b=1
  arma::mat state; // state spin configurations of system
  arma::vec w;
  int E;
  int M;

  // Initializer
  Ising(int N_in, double T_in);
  Ising(int N_in, double T_in, double J_in);

  // Functions
  void genState();
  void energy();
  void Metropolis();
};

#endif
