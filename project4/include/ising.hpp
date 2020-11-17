#ifndef ISING_H
#define ISING_H

//#include "mpi.h"
#include <iostream>
#include <cmath>
#include <armadillo>
#include <stdlib.h>

class Ising
{
public:
  // Properties

  int N; // dimensions of lattice (LxL)
  int Nspins; // total spins
  int seed; // seeding our rngs
  int ordered; // controls the initial matrix 0 = ordered, 1 = random

  double E; // energy of model
  double M; // magnetization
  double J; // coupling constant
  double T; // temperature, same as energy with K_b=1

  std::mt19937_64 gen; // RNG object
  std::uniform_real_distribution<double> URD; // uniform real dist [0,1]

  arma::mat state; // state spin configurations of system
  arma::vec w; // probability array (boltzmann<3)

  arma::vec avgE;
  arma::vec avgM;
  arma::vec expE; // energy expectation value
  arma::vec expM; // magnetisation expectation value

  // Initializer
  Ising(int N, double T_in);
  Ising(int N_in, double T_in, double J_in, int seed_in, int ordered_in );

  // Functions
  void init();
  void genState();
  void energy();
  void Metropolis();
  void MC(int cycles_in, int burnin);
  void paraMC(int cycles_in);
  void burnin(int cycles_in);
  void save(std::string filename);

};

#endif
