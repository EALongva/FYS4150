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
  int seed; // seeding our rngs
  int E; // energy of model
  int M; // magnetization

  double J; // coupling constant
  double T; // temperature, same as energy with K_b=1

  arma::mat state; // state spin configurations of system
  arma::vec w; // probability array (boltzmann<3)

  // Initializer
  Ising(int N, double T_in);
  Ising(int N_in, double T_in, double J_in, int seed_in );

  // Functions
  void init();
  void genState();
  void energy();
  void Metropolis();

};

#endif
