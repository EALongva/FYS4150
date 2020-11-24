// example file project 4

#include <iostream>
#include <armadillo>
#include <vector>
#include <cstdlib>
#include <cmath>
#include "time.h"
//#include "include/exampleClass.hpp"
//#include "include/utils.hpp"
//#include "include/metro.hpp"
#include "include/ising.hpp"

int main(int argc, char* argv[]){

  // declaring variables
  /*
  int N, burnin, MC;
  double T_init, T_fin, dT;
  std::string filename;
  std::string outpath = "results/data/";

  if(argc < 8){

    // standard values (aka lazy values)
    N       = 20;
    burnin  = 5;
    MC      = 5;

    T_init     = 2.0;
    T_fin      = 2.4;
    dT         = 0.05;

    filename = "Simulation_N" + std::to_string(N);

    std::cout << "Bad usage: " << argv[0] << std::endl;
    std::cout << " read output file, Number of spins, burn-in MC cycles, MC cycles,"
    << "initial and final temperature and tempurate step" << std::endl;
    std::cout << "-------------------------------------------------" << std::endl;
    std::cout << "Implementing standard values: " << filename << ", " << N <<
    ", " << burnin << ", " << MC << ", " << T_init << ", " << T_fin << ", "<<
    dT << std::endl;
  }
  else{

    // command line args must be passed to class method (so MPI_Init can remove
    // the unecessary argvs)

  }

  std::cout << N << std::endl;
  */

  Ising M(1, 1337);
  M.parallel(argc, argv);



  return 0;
}
