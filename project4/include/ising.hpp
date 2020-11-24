#ifndef ISING_H
#define ISING_H

#include "mpi.h"
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

  // system after burnin
  arma::mat burnin_state; // state after burn in
  double burnin_E; // energy after burnin
  double burnin_M; // magnetisation after burnin

  // system after init
  arma::mat init_state; // state after init
  double init_E; // energy after init
  double init_M; // magnetisation after init

  arma::vec avgE; // averages used to find number of burnin cycles
  arma::vec avgM;

  arma::vec expvals; // array of expectation values E, E^2, |M|, M^2
  arma::mat totexpvals;

  // used to benchmark time consumption of runs
  double time;     // total time
  arma::vec ptime; // time of each process

  // Initializer
  Ising(int N_in, int seed_in);
  Ising(int N_in, double T_in, double J_in, int seed_in, int ordered_in );

  // Functions
  void init();
  void genState();
  void energy();
  void Metropolis();
  void burninMC(int cycles_in);
  void MC(int cycles_in);
  void paraMC(int cycles_in);
  void burnin(int cycles_in);
  void observables(int cycles_in, int burnin_cycles_in, double T_start, double T_end, int T_N);
  void reset2burnin();
  void reset2init();
  void save(std::string filename);
  void hello(int &argc, char** &argv);
  void parallel(int &argc, char** &argv);

};

#endif
